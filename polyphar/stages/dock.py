"""Stage 07: dock the ligand library (e-Drug3D) into every selected target with smina.

Inputs: 06_targets/final_targets.csv (Protein, Pocket, Receptor_Path, Pocket_Path; paths
relative to the run folder) and the library SDF (`--library`, default
inputs/e-Drug3D_2083.sdf; a library outside <run>/inputs is copied there first).

Receptors are protonated once by pdbfixer at pH 7.0 and never re-protonated:
  - a PDB chain's Receptor_Path is the fix-stage output (03_pdbfixer/..., pdbfixer pH 7.0);
  - an AlphaFold model (Receptor_Path in 02_structures, not pdbfixer'd) gets hydrogens only,
    pdbfixer addMissingHydrogens(7.0) (seeded, Reference platform), into
    07_docking/<label>/receptor_h.pdb.
The PDBQT (07_docking/<label>/receptor.pdbqt) is written by Open Babel keeping those
hydrogens (no CorrectForPH, no addh; non-polar H merged as AutoDock does), Gasteiger
charges, rigid. This replaces the legacy CorrectForPH(7.4) + addh, which re-placed the
hydrogens of ionisable groups at random on every conversion (scores moved by up to
0.9 kcal/mol between two conversions of one receptor). Conversion is byte-deterministic;
each receptor PDBQT's sha256 is recorded in the manifest.

Kept legacy parameters (smina_screen_P_2_2023_08_15.py): box = the pocket's atom extent,
unpadded (centre (max+min)/2, size max-min); `--exhaustiveness 12 --num_modes 12`; the
ligand records as supplied. New: a fixed recorded `--seed` (default 42), and ligand
identity from the library `ID` field: every record's title line is set to EDRUG3D_<ID>
(the library's titles are blank or placeholders, and 35 different drugs are named
"Molecule-1"), and poses are mapped back to the library by it.

Work units are (target x chunk of --chunk-size ligands, in library order), docked with
--cpu --threads-per-job (default 1) on nprocs // threads_per_job spawn processes. A
ligand's result depends neither on --cpu nor on its chunk or position in it (pinned by
test_cpu_and_chunk_invariance), so the output is the same for any chunking and nprocs.
A chunk is done when 07_docking/<label>/chunks/<chunk>.done holds its job key (receptor
and pocket sha256, chunk sha256, seed, exhaustiveness, num_modes, smina version, receptor
preparation); outputs are atomic, so a killed run redoes only unfinished chunks.

Outputs per target, label = f"{Protein}_{Pocket}": 07_docking/<label>/poses.sdf and
07_docking/df_results_<label>_fixed.csv (Pose, SMINA_Score, Ligand, ID, Pocket; Ligand is
the library `name`, ID the library ID, SMINA_Score smina's minimizedAffinity). A target
with any ligand ID missing from its poses (a failed chunk, or a ligand smina skipped) is
incomplete: failed.tsv lists the missing IDs and NO df_results/poses.sdf is written for it.
manifest.json is written last; a rerun on the same inputs changes nothing.
"""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from polyphar.stages.fix import H_SEED, reference_platform, seeded_placement

STAGE = "dock"
DEFAULT_SEED = 42
DEFAULT_CHUNK = 50
DEFAULT_TPJ = 1
EXHAUSTIVENESS = 12
NUM_MODES = 12
AF_H_SEED = H_SEED   # seeds OpenMM's random hydrogen placement (deterministic output)
AF_PH = 7.0   # AlphaFold receptors: pdbfixer addMissingHydrogens, as the fix stage
TITLE_PREFIX = "EDRUG3D_"
RESULT_COLUMNS = ["Pose", "SMINA_Score", "Ligand", "ID", "Pocket"]
LIBRARY_COLUMNS = ["ID", "name", "cas", "status", "record_index"]
END = "$$$$"


def smina_bin() -> str:
    return os.environ.get("SMINA_BIN") or "smina"


# --- SDF text handling ---------------------------------------------------------------

_FIELD = re.compile(r"^>.*<([^>]+)>")


@dataclass
class Record:
    lines: list[str]                     # verbatim, line endings kept, incl. the $$$$ line
    fields: dict = field(default_factory=dict)


def _records(text: str) -> list[Record]:
    """Split SDF text into records, keeping every byte; parse the `> <field>` data items
    (first value line, stripped)."""
    out, cur = [], []
    for line in text.splitlines(keepends=True):
        cur.append(line)
        if line.rstrip("\r\n") == END:
            out.append(Record(cur))
            cur = []
    if any(line.strip() for line in cur):
        out.append(Record(cur))   # last record without a $$$$ line
    for rec in out:
        lines = rec.lines
        for i, line in enumerate(lines):
            m = _FIELD.match(line)
            if m and i + 1 < len(lines):
                rec.fields.setdefault(m.group(1), lines[i + 1].strip())
    return out


def _titled(rec: Record) -> str:
    """The record with its title line set to EDRUG3D_<ID> (the rest verbatim)."""
    first = rec.lines[0]
    eol = first[len(first.rstrip("\r\n")):] or "\n"
    return f"{TITLE_PREFIX}{rec.fields['ID']}{eol}" + "".join(rec.lines[1:])


def _read_text(path) -> str:
    with open(path, encoding="utf-8", errors="surrogateescape", newline="") as fh:
        return fh.read()


def _atomic_text(path: Path, text: str) -> None:
    path = Path(path)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", errors="surrogateescape", newline="") as fh:
            fh.write(text)
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _library_records(sdf) -> list[Record]:
    recs = _records(_read_text(sdf))
    missing = [i for i, r in enumerate(recs) if not r.fields.get("ID")]
    if missing:
        raise ValueError(f"{sdf}: {len(missing)} record(s) without an ID field "
                         f"(record index {', '.join(map(str, missing[:10]))})")
    seen: dict[str, int] = {}
    dups = []
    for i, r in enumerate(recs):
        if r.fields["ID"] in seen:
            dups.append(r.fields["ID"])
        seen.setdefault(r.fields["ID"], i)
    if dups:
        raise ValueError(f"{sdf}: duplicate ID(s): {', '.join(sorted(set(dups))[:10])}")
    return recs


def read_library(sdf):
    """DataFrame[ID, name, cas, status, record_index] (ID as text). Raises on a missing or
    duplicate ID."""
    import pandas as pd

    rows = [[r.fields["ID"], r.fields.get("name", ""), r.fields.get("cas", ""),
             r.fields.get("status", ""), i] for i, r in enumerate(_library_records(sdf))]
    return pd.DataFrame(rows, columns=LIBRARY_COLUMNS)


def write_titled_library(sdf_in, sdf_out) -> None:
    """Copy the records verbatim, each title line set to EDRUG3D_<ID>."""
    _atomic_text(Path(sdf_out), "".join(_titled(r) for r in _library_records(sdf_in)))


# --- receptor and box ----------------------------------------------------------------

def _strip_date_remark(text: str) -> str:
    """Drop OpenMM's dated header line, so the same input gives the same bytes any day."""
    return "".join(l for l in text.splitlines(keepends=True)
                   if not l.startswith("REMARK   1 CREATED WITH OPENMM"))


def protonate_alphafold(pdb_in, pdb_out, ph: float = AF_PH) -> None:
    """Add hydrogens to an AlphaFold model with pdbfixer (addMissingHydrogens(ph) only:
    no missing residues/atoms, no heterogen removal). Written atomically, without the
    dated OpenMM header line.

    OpenMM places new hydrogens at random (Python `random`) and then minimises them; the
    RNG is seeded (AF_H_SEED, state restored afterwards) and the Reference platform is
    used, so the same model always gives the same bytes."""
    import io

    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(pdb_in), platform=reference_platform())
    with seeded_placement(AF_H_SEED):
        fixer.addMissingHydrogens(ph)
    buf = io.StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, buf, keepIds=True)
    _atomic_text(Path(pdb_out), _strip_date_remark(buf.getvalue()))


def receptor_to_pdbqt(pdb, pdbqt) -> None:
    """Rigid receptor PDBQT with Open Babel, keeping the input's hydrogens: no
    CorrectForPH, no addh (the receptor was protonated once by pdbfixer at pH 7.0).
    Gasteiger partial charges as in the legacy script; non-polar hydrogens are merged
    (AutoDock convention, as the legacy PDBQT), polar ones written as HD at their input
    coordinates. The `REMARK Name =` line carries the file's base name only, so the bytes
    do not depend on where the run folder is mounted. Written atomically."""
    from openbabel import pybel

    pybel.ob.obErrorLog.SetOutputLevel(0)   # no Open Babel chatter on stderr
    pdb, pdbqt = Path(pdb), Path(pdbqt)
    mol = next(iter(pybel.readfile("pdb", str(pdb))), None)
    if mol is None or not mol.atoms:
        raise ValueError(f"{pdb}: Open Babel read no molecule")
    for atom in mol.atoms:
        atom.OBAtom.GetPartialCharge()   # Gasteiger (Open Babel default)
    text = mol.write("pdbqt", opt={"r": None})
    if not any(l.startswith(("ATOM", "HETATM")) for l in text.splitlines()):
        raise ValueError(f"{pdb}: Open Babel wrote no atoms")
    lines = [f"REMARK  Name = {pdb.name}" if l.startswith("REMARK  Name =") else l
             for l in text.splitlines()]
    _atomic_text(pdbqt, "\n".join(lines) + "\n")


def pocket_box(pocket_pdb) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """(centre, size) of the pocket's ATOM/HETATM coordinates, unpadded."""
    xyz = []
    with open(pocket_pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not xyz:
        raise ValueError(f"{pocket_pdb}: no atoms")
    lo = [min(c[k] for c in xyz) for k in range(3)]
    hi = [max(c[k] for c in xyz) for k in range(3)]
    return (tuple((hi[k] + lo[k]) / 2 for k in range(3)),
            tuple(hi[k] - lo[k] for k in range(3)))


# --- smina ---------------------------------------------------------------------------

class SminaError(RuntimeError):
    pass


def smina_version() -> str:
    """First line of `smina --version`, or "unavailable"."""
    try:
        p = subprocess.run([smina_bin(), "--version"], capture_output=True, text=True,
                           timeout=60)
    except (OSError, subprocess.SubprocessError):
        return "unavailable"
    lines = (p.stdout or p.stderr).strip().splitlines()
    return lines[0].strip() if p.returncode == 0 and lines else "unavailable"


def run_smina(receptor, ligands, center, size, out_sdf, log, seed: int, cpu: int,
              exhaustiveness: int = EXHAUSTIVENESS, num_modes: int = NUM_MODES) -> None:
    """Dock `ligands` into the box; out_sdf and log are written atomically. Raises
    SminaError on a non-zero exit."""
    out_sdf, log = Path(out_sdf), Path(log)
    tok = f"{os.getpid()}.{os.urandom(4).hex()}"
    # smina picks the output format from the extension
    tmp_out = out_sdf.with_name(f".{out_sdf.name}.{tok}.tmp.sdf")
    tmp_log = log.with_name(f".{log.name}.{tok}.tmp")
    cmd = [smina_bin(), "-r", str(receptor), "-l", str(ligands)]
    for axis, c in zip("xyz", center):
        cmd += [f"--center_{axis}", repr(float(c))]
    for axis, s in zip("xyz", size):
        cmd += [f"--size_{axis}", repr(float(s))]
    cmd += ["--out", str(tmp_out), "--log", str(tmp_log), "--num_modes", str(num_modes),
            "--exhaustiveness", str(exhaustiveness), "--cpu", str(cpu), "--seed", str(seed)]
    try:
        p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                           text=True, errors="replace")
        if p.returncode != 0:
            tail = " | ".join((p.stderr or "").strip().splitlines()[-5:])
            raise SminaError(f"smina exited with {p.returncode}: {tail}")
        if not tmp_out.exists():
            raise SminaError("smina wrote no output file")
        if tmp_log.exists():
            os.replace(tmp_log, log)
        os.replace(tmp_out, out_sdf)
    finally:
        tmp_out.unlink(missing_ok=True)
        tmp_log.unlink(missing_ok=True)


def parse_poses(out_sdf, library, label: str):
    """DataFrame[Pose, SMINA_Score, Ligand, ID, Pocket] from a smina output SDF whose
    titles are EDRUG3D_<ID>: Pose counts 1.. within each ligand's consecutive records,
    SMINA_Score is minimizedAffinity, Ligand the library name."""
    import pandas as pd

    names = dict(zip(library["ID"].astype(str), library["name"]))
    rows, prev, pose = [], None, 0
    for rec in _records(_read_text(out_sdf)):
        title = rec.lines[0].strip()
        if not title.startswith(TITLE_PREFIX) or title[len(TITLE_PREFIX):] not in names:
            raise ValueError(f"{out_sdf}: pose title {title!r} is not a library ID")
        lid = title[len(TITLE_PREFIX):]
        pose = pose + 1 if lid == prev else 1
        prev = lid
        if "minimizedAffinity" not in rec.fields:
            raise ValueError(f"{out_sdf}: pose {pose} of {title} has no minimizedAffinity")
        rows.append([pose, float(rec.fields["minimizedAffinity"]), names[lid], lid, label])
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


# --- stage -------------------------------------------------------------------------

DEFAULT_LIBRARY = "inputs/e-Drug3D_2083.sdf"
PREP_PDB = "pdbfixer pH 7.0 (fix stage); Open Babel PDBQT, H kept, Gasteiger, rigid"
PREP_AF = (f"pdbfixer addMissingHydrogens({AF_PH}) seed {AF_H_SEED} Reference; "
           "Open Babel PDBQT, H kept, Gasteiger, rigid")
PREP_VERSION = 1   # bump when the receptor preparation changes


def _positive_int(value: str) -> int:
    import argparse

    try:
        n = int(value)
    except (TypeError, ValueError):
        n = 0
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be an integer >= 1 (got {value!r})")
    return n


def add_arguments(parser) -> None:
    parser.add_argument("--library", default=DEFAULT_LIBRARY,
                        help=f"ligand SDF, relative to the run folder (default {DEFAULT_LIBRARY})")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"smina --seed, recorded (default {DEFAULT_SEED})")
    parser.add_argument("--chunk-size", type=_positive_int, default=DEFAULT_CHUNK,
                        help=f"ligands per job (default {DEFAULT_CHUNK})")
    parser.add_argument("--threads-per-job", type=_positive_int, default=DEFAULT_TPJ,
                        help=f"smina --cpu per job; nprocs // this jobs run at once "
                             f"(default {DEFAULT_TPJ})")


@dataclass
class Target:
    label: str
    receptor: Path        # absolute
    pocket: Path          # absolute
    receptor_rel: str
    pocket_rel: str
    alphafold: bool

    @property
    def prep(self) -> str:
        return PREP_AF if self.alphafold else PREP_PDB


def _rel(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path)


def read_targets(layout) -> list[Target]:
    """One Target per final_targets.csv row; raises on a missing file or column, or on
    two rows with the same label."""
    import pandas as pd

    path = layout.final_targets_csv
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; run `polyphar select` first")
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    need = ["Protein", "Pocket", "Receptor_Path", "Pocket_Path"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
    out, seen = [], set()
    for row in df.itertuples(index=False):
        label = f"{row.Protein}_{row.Pocket}"
        if label in seen:
            raise ValueError(f"{path}: two targets with the label {label!r}")
        seen.add(label)
        rel = Path(row.Receptor_Path)
        out.append(Target(label, layout.run_dir / rel, layout.run_dir / row.Pocket_Path,
                          row.Receptor_Path, row.Pocket_Path,
                          alphafold=rel.parts[:1] == (layout.structures.name,)))
    return out


def _library_path(layout, args, copy: bool = True) -> Path:
    """The library inside <run>/inputs; one given elsewhere is copied there (never
    modified; with copy=False, as for the ETA, it is only located). Raises when it does
    not exist."""
    import shutil

    lib = Path(getattr(args, "library", DEFAULT_LIBRARY))
    if not lib.is_absolute():
        lib = layout.run_dir / lib
    if not lib.is_file():
        raise FileNotFoundError(f"ligand library {lib} not found")
    inputs = layout.inputs.resolve()
    if lib.resolve().parent == inputs or not copy:
        return lib
    dst = layout.inputs / lib.name
    if not (dst.is_file() and _md5(dst) == _md5(lib)):
        if dst.exists():
            raise FileExistsError(f"{dst} exists and differs from {lib}; not overwritten")
        layout.inputs.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_name(f".{dst.name}.{os.getpid()}.tmp")
        shutil.copyfile(lib, tmp)
        os.replace(tmp, dst)
    return dst


def _md5(path) -> str:
    import hashlib

    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _sha_text(text: str) -> str:
    import hashlib

    return hashlib.sha256(text.encode("utf-8", "surrogateescape")).hexdigest()


@dataclass
class Chunk:
    name: str
    ids: list[str]
    text: str             # titled records
    sha256: str


def make_chunks(records: list[Record], size: int) -> list[Chunk]:
    out = []
    for k in range(0, len(records), size):
        part = records[k:k + size]
        text = "".join(_titled(r) for r in part)
        out.append(Chunk(f"c{size}_{k // size:05d}", [r.fields["ID"] for r in part], text,
                         _sha_text(text)))
    return out


def _target_dir(layout, t: Target) -> Path:
    return layout.docking / t.label


def _df_path(layout, t: Target) -> Path:
    return layout.docking / f"df_results_{t.label}_fixed.csv"


def _src_key(t: Target) -> str:
    """What the receptor PDBQT depends on (it is byte-deterministic given these)."""
    from polyphar.common import sha256_file

    return _sha_text(_json({"receptor_sha256": sha256_file(t.receptor),
                                 "prep": t.prep, "prep_version": PREP_VERSION}))


def _json(obj) -> str:
    import json

    return json.dumps(obj, sort_keys=True)


def _job_key(src_key: str, pocket_sha: str, chunk: Chunk, opts: dict) -> str:
    return _sha_text(_json({"receptor": src_key, "pocket_sha256": pocket_sha,
                                 "chunk_sha256": chunk.sha256, **opts}))


def _options(args, smina: str) -> dict:
    return {"seed": int(getattr(args, "seed", DEFAULT_SEED)),
            "exhaustiveness": EXHAUSTIVENESS, "num_modes": NUM_MODES, "smina": smina}


def _chunk_paths(layout, t: Target, chunk: Chunk) -> tuple[Path, Path, Path]:
    d = _target_dir(layout, t) / "chunks"
    return d / f"{chunk.name}.sdf", d / f"{chunk.name}.log", d / f"{chunk.name}.done"


def _chunk_done(layout, t: Target, chunk: Chunk, key: str) -> bool:
    sdf, _, done = _chunk_paths(layout, t, chunk)
    try:
        return done.read_text().strip() == key and sdf.is_file()
    except OSError:
        return False


@dataclass
class _Plan:
    targets: list[Target]
    chunks: list[Chunk]
    keys: dict            # (label, chunk name) -> job key
    usable: list[Target]  # receptor and pocket present
    missing: list[tuple[str, str]]   # (label, reason) for unusable targets

    def todo(self, layout) -> list[tuple[Target, Chunk]]:
        return [(t, c) for t in self.usable for c in self.chunks
                if not _chunk_done(layout, t, c, self.keys[(t.label, c.name)])]


def _plan(layout, args, lib: Path, smina: str) -> _Plan:
    from polyphar.common import sha256_file

    targets = read_targets(layout)
    chunks = make_chunks(_library_records(lib), int(getattr(args, "chunk_size", DEFAULT_CHUNK)))
    opts = _options(args, smina)
    keys, usable, missing = {}, [], []
    for t in targets:
        absent = [p for p, what in ((t.receptor, "receptor"), (t.pocket, "pocket"))
                  if not p.is_file()]
        if absent:
            missing.append((t.label, "missing file: " + ", ".join(
                _rel(p, layout.run_dir) for p in absent)))
            continue
        usable.append(t)
        src, pocket = _src_key(t), sha256_file(t.pocket)
        for c in chunks:
            keys[(t.label, c.name)] = _job_key(src, pocket, c, opts)
    return _Plan(targets, chunks, keys, usable, missing)


def _prepare_receptor(layout, t: Target) -> Path:
    """07_docking/<label>/receptor.pdbqt, reused while its source key is unchanged."""
    d = _target_dir(layout, t)
    d.mkdir(parents=True, exist_ok=True)
    pdbqt, stamp = d / "receptor.pdbqt", d / "receptor.key"
    key = _src_key(t)
    try:
        if pdbqt.is_file() and stamp.read_text().strip() == key:
            return pdbqt
    except OSError:
        pass
    stamp.unlink(missing_ok=True)
    src = t.receptor
    if t.alphafold:
        src = d / "receptor_h.pdb"
        protonate_alphafold(t.receptor, src)
    receptor_to_pdbqt(src, pdbqt)
    _atomic_text(stamp, key + "\n")
    return pdbqt


# --- workers (spawn processes) -------------------------------------------------------

def _dock_chunk(receptor, ligands, center, size, sdf, log, done, seed, cpu, key) -> int:
    run_smina(receptor, ligands, center, size, sdf, log, seed=seed, cpu=cpu)
    _atomic_text(Path(done), key + "\n")
    return 0


def _chunk_worker(*args) -> int:
    from polyphar.stages.similarity import _task

    return _task(lambda: _dock_chunk(*args))


# --- merge -----------------------------------------------------------------------------

def _merge(layout, t: Target, chunks: list[Chunk], library, failed_chunks: set) -> list[str]:
    """Merge a target's done chunks. Returns the library IDs without poses ([] when the
    target is complete; then poses.sdf and df_results are written atomically). An
    incomplete target gets neither, and stale ones are removed."""
    d = _target_dir(layout, t)
    poses, df_path = d / "poses.sdf", _df_path(layout, t)
    texts = []
    for c in chunks:
        if c.name in failed_chunks:
            continue
        sdf, _, _ = _chunk_paths(layout, t, c)
        if sdf.is_file():
            texts.append(_read_text(sdf))
    text = "".join(texts)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".poses.sdf.", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp)
    try:
        tmp.write_text(text, encoding="utf-8", errors="surrogateescape")
        df = parse_poses(tmp, library, t.label)
        have = set(df["ID"])
        missing = [i for i in library["ID"] if i not in have]
        if missing:
            poses.unlink(missing_ok=True)
            df_path.unlink(missing_ok=True)
            return missing
        os.chmod(tmp, 0o644)
        os.replace(tmp, poses)
    finally:
        tmp.unlink(missing_ok=True)
    fd, tmpc = tempfile.mkstemp(dir=layout.docking, prefix=f".{df_path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            df.to_csv(fh, index=False)
        os.chmod(tmpc, 0o644)
        os.replace(tmpc, df_path)
    except BaseException:
        Path(tmpc).unlink(missing_ok=True)
        raise
    return []


# --- fingerprint / skip --------------------------------------------------------------

def _stage_fingerprint(layout, plan: _Plan, lib_md5: str, opts: dict) -> str | None:
    """What the finished outputs depend on (not nprocs, chunk size or threads): the
    final_targets.csv, each target's receptor, pocket and receptor PDBQT sha256, the
    library md5, seed, exhaustiveness, num_modes and smina version. None when a receptor
    PDBQT is not prepared yet."""
    from polyphar.common import sha256_file

    rows = []
    for t in plan.targets:
        pdbqt = _target_dir(layout, t) / "receptor.pdbqt"
        if t in plan.usable:
            if not pdbqt.is_file():
                return None
            rows.append([t.label, sha256_file(t.receptor), sha256_file(t.pocket), t.prep,
                         sha256_file(pdbqt)])
        else:
            rows.append([t.label, None, None, None, None])
    return _sha_text(_json({"final_targets_sha256": sha256_file(layout.final_targets_csv),
                                 "targets": rows, "library_md5": lib_md5, **opts}))


def _outputs(layout) -> dict:
    d = layout.docking
    files = sorted(d.glob("df_results_*_fixed.csv")) + sorted(d.glob("*/poses.sdf"))
    return {p.relative_to(d).as_posix(): p.stat().st_size for p in files}


def _finished(layout, plan: _Plan, lib_md5: str, opts: dict) -> dict | None:
    import json

    try:
        manifest = json.loads((layout.docking / "manifest.json").read_text())
    except (OSError, ValueError):
        return None
    fp = _stage_fingerprint(layout, plan, lib_md5, opts)
    if fp is None or manifest.get("fingerprint") != fp or manifest.get("outputs") != _outputs(layout):
        return None
    if plan.todo(layout):   # a chunk lost its .done marker: redo it
        return None
    return dict(manifest["counts"])


def units(layout, args):
    """Ligands still to dock: the ligands of every (target, chunk) without a valid .done
    marker (0 when the stage finished on exactly these inputs)."""
    from polyphar.common import WorkUnits

    lib = _library_path(layout, args, copy=False)
    smina = smina_version()
    plan = _plan(layout, args, lib, smina)
    if _finished(layout, plan, _md5(lib), _options(args, smina)) is not None:
        return WorkUnits(0, "ligand")
    return WorkUnits(sum(len(c.ids) for _, c in plan.todo(layout)), "ligand")


def run(layout, args, progress) -> dict:
    from polyphar.common import append_failure, sha256_file, write_manifest
    from polyphar.stages.similarity import _run_pool

    lib = _library_path(layout, args)
    smina = smina_version()
    if smina == "unavailable":
        raise FileNotFoundError(f"smina ({smina_bin()}) not runnable")
    opts = _options(args, smina)
    lib_md5 = _md5(lib)
    plan = _plan(layout, args, lib, smina)
    done = _finished(layout, plan, lib_md5, opts)
    if done is not None:   # finished on exactly these inputs: touch nothing
        return done

    d = layout.docking
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    library = read_library(lib)
    rows: list[tuple[str, str]] = list(plan.missing)

    # receptors (once per target) and boxes, in the parent
    prepared: dict[str, tuple[Path, tuple, tuple]] = {}
    for t in plan.usable:
        try:
            pdbqt = _prepare_receptor(layout, t)
            prepared[t.label] = (pdbqt, *pocket_box(t.pocket))
        except Exception as exc:
            rows.append((t.label, f"receptor/pocket preparation failed: "
                                  f"{type(exc).__name__}: {exc}"))
    ligdir = d / "ligands"
    ligdir.mkdir(exist_ok=True)
    for c in plan.chunks:
        f = ligdir / f"{c.name}.sdf"
        if not (f.is_file() and sha256_file(f) == c.sha256):
            _atomic_text(f, c.text)

    tpj = int(getattr(args, "threads_per_job", DEFAULT_TPJ))
    slots = max(1, args.nprocs // tpj)
    todo = [(t, c) for t, c in plan.todo(layout) if t.label in prepared]
    n_of = {c.name: len(c.ids) for c in plan.chunks}
    tasks = []
    for t, c in todo:
        pdbqt, center, size = prepared[t.label]
        sdf, log, marker = _chunk_paths(layout, t, c)
        sdf.parent.mkdir(parents=True, exist_ok=True)
        tasks.append(((t.label, c.name), _chunk_worker,
                      (str(pdbqt), str(ligdir / f"{c.name}.sdf"), center, size, str(sdf),
                       str(log), str(marker), opts["seed"], tpj,
                       plan.keys[(t.label, c.name)])))
    failed_chunks: dict[str, set] = {}
    try:
        outcome = _run_pool(tasks, slots,
                            on_result=lambda item, _r: progress.advance(n_of[item[1]]))
        for (label, chunk), reason in outcome.failures:
            failed_chunks.setdefault(label, set()).add(chunk)
            rows.append((f"{label}/{chunk}", reason))
            progress.advance(n_of[chunk])

        complete = incomplete = 0
        for t in plan.usable:
            if t.label not in prepared:
                incomplete += 1
                _df_path(layout, t).unlink(missing_ok=True)
                continue
            missing = _merge(layout, t, plan.chunks, library, failed_chunks.get(t.label, set()))
            if missing:
                incomplete += 1
                rows.append((t.label, f"incomplete: {len(missing)} ligand ID(s) without "
                                      f"poses: {','.join(missing)}"))
            else:
                complete += 1
    finally:
        for item, reason in rows:   # appended from the parent only
            append_failure(d, item, STAGE, reason)

    incomplete += len(plan.missing)
    counts = {"targets": len(plan.targets), "complete": complete, "incomplete": incomplete,
              "ligands": len(library), "chunks": len(plan.chunks) * len(plan.usable),
              "chunks_docked": len(tasks) - sum(map(len, failed_chunks.values())),
              "skipped_existing": len(plan.chunks) * len(prepared) - len(tasks),
              "failed": len(rows)}
    receptors = {t.label: {"receptor": t.receptor_rel, "pocket": t.pocket_rel,
                           "preparation": t.prep,
                           "pdbqt_sha256": sha256_file(prepared[t.label][0])}
                 for t in plan.usable if t.label in prepared}
    fingerprint = _stage_fingerprint(layout, plan, lib_md5, opts)
    write_manifest(d, STAGE,
                   args={"nprocs": args.nprocs, "seed": opts["seed"],
                         "chunk_size": int(getattr(args, "chunk_size", DEFAULT_CHUNK)),
                         "threads_per_job": tpj, "exhaustiveness": EXHAUSTIVENESS,
                         "num_modes": NUM_MODES, "library": _rel(lib, layout.run_dir)},
                   inputs=[layout.final_targets_csv, lib], versions=_versions(smina),
                   counts=counts, timing={},
                   extra={"fingerprint": fingerprint, "outputs": _outputs(layout),
                          "library_md5": lib_md5, "seed": opts["seed"], "smina": smina,
                          "receptors": receptors})
    return counts


def _versions(smina: str) -> dict:
    import importlib.metadata

    out = {"smina": smina}
    for pkg in ("openbabel", "pdbfixer", "openmm"):
        try:
            out[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            out[pkg] = "unknown"
    return out
