"""Stage 07: dock the ligand library (e-Drug3D) into every selected target.

The engine is selectable (`--engine`): **smina** on CPU, the single-host / development
path, or **Uni-Dock >= 1.2.0** on a CUDA GPU through `unidock_tools`
(polyphar/stages/dock_unidock.py). The engine name and its search settings are part of
every unit key, so changing engine supersedes the recorded draws rather than mixing two
populations under one target. Uni-Dock needs a GPU and says so before it does any work.

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
ligand records as supplied, with smina's ligand hydrogen addition pinned explicitly to its
default `--addH 1` (recorded in the manifest and the unit key). New: a recorded per-replicate
`--seed` (base `--seed`, default 42, plus the replicate index), and ligand
identity from the library `ID` field: every record's title line is set to EDRUG3D_<ID>
(the library's titles are blank or placeholders, and 35 different drugs are named
"Molecule-1"), and poses are mapped back to the library by it.

Work unit = (target, replicate): ONE docking engine invocation over the WHOLE prepared
library, repeated --replicates times (default 3). Replicate i uses seed base_seed + i
(--seed is the base, default 42), recorded per unit.

**The replicate index is part of a unit's identity, and resume must skip, never redo.**
The engine this stage is being moved to (Uni-Dock) is not deterministic: identical
same-seed repeats move 95-98% of ligand scores. So a unit's key --

    (receptor sha256, pocket sha256, library sha256, replicate i, seed_i,
     search parameters, engine version)

-- does NOT assert "re-running would give these numbers". It asserts "this (target,
replicate) was run under this identity and its result file is the record of it". The
recorded seed is provenance, not a reproducibility claim (the engine's seeding is broken
upstream). Replicate 2 is a distinct observation, not a repeat of replicate 1: re-running
a finished unit would not reproduce it, it would draw a different sample. Hence:

* `unit_state` == "done" (marker holds the key, the result file is intact) -> SKIPPED. The
  stage never re-runs it, for any reason short of the key changing.
* "corrupt" (marker holds the key but the result file is gone, empty or not terminated) ->
  the draw is unrecoverable, so it is DISCARDED and re-drawn under the SAME replicate
  index. This is the only legitimate repeat. It is appended to 07_docking/discarded.tsv
  (append-only, across runs) because the discarded draw was a real sample and dropping it
  is a small selection effect. A unit whose engine output is systematically unusable is
  therefore re-drawn, and logged, on every invocation -- discarded.tsv is where that shows.
* "missing" (no marker, or a marker for another key) -> run it. When a result file from a
  superseded key is about to be overwritten that, too, is logged as a discard.

**Discarding is DEFERRED to the moment a replacement is actually about to be written** --
inside the worker, immediately before the engine call that replaces the record, and inside
`import_bundle` immediately before the bundle's copy lands. Planning what would be thrown
away (`plan_discards`) throws nothing away, so a pass that is filtered (`--max-units`, a
target whose receptor preparation failed), aborted, or not executing at all
(`dock-merge`, which plans with execute=False) costs no draws.

**And a discarded record is MOVED, never deleted**, to
07_docking/discarded/<timestamp>/<label>/r<i>.sdf, the policy orphan labels already follow;
discarded.tsv names the archived file in its `archived` column.

Per unit, beside the record, the host that DREW it writes r<i>.prov.json: engine build,
unidock_tools version, CUDA runtime, GPU model, search settings, replicate and seed, the
library / receptor / pocket sha256s, the host path, and the RESULT's sha256 -- which is
what lets `dock-merge` refuse a record edited after the engine wrote it, and what lets it
warn when one target's replicates span GPU models.

Markers live in 07_docking/<label>/replicates/r<i>.done, beside the record they certify --
inside the directory a copy, a tar or an `rm -rf` would treat as the unit, so the two
cannot be separated. Results are 07_docking/<label>/replicates/r<i>.sdf with r<i>.log
beside them. An engine
call running longer than --unit-timeout (default 12 h) is killed and the unit failed; its
log is named in failed.tsv. counts["failed"] is the number of failed units,
counts["discarded"] the number of discarded draws; incomplete targets are counted in
counts["incomplete_targets"].

Outputs per target, label = f"{Protein}_{Pocket}":

* 07_docking/df_results_<label>.csv -- the SUMMARY over every replicate (SUMMARY_COLUMNS:
  Ligand, ID, Pocket, n, score_mean, score_sd, score_se, score_min, score_max). One row per
  ligand; the observation summarised is that ligand's BEST score in each replicate.
* 07_docking/<label>/poses.sdf -- replicate 0's poses, under the default `--poses first`
  (`all` adds poses_r<i>.sdf, `none` writes none). Replicate 0's poses are what the
  control stage's geometric readouts read.

The raw per-replicate records (<label>/replicates/r<i>.sdf) are kept: they cost little,
they let any other summary be recomputed without re-docking, and they are the only way to
revisit the jitter question later.

A target with a replicate not done, or with any ligand ID missing from ANY replicate, is
incomplete: failed.tsv says why and no df_results/poses file is written for it -- a
replicate that fails validation can never enter the summary.
Outputs of labels no longer in final_targets.csv (df_results_<label>.csv, the retired
df_results_<label>_fixed.csv under any label, <label>/)
are moved, never deleted, to 07_docking/obsolete/<timestamp>/ and listed in the manifest
(`obsoleted`); the manifest's `outputs` lists the current targets' files only.
manifest.json is written last; a rerun on the same inputs changes nothing.

`dock` is the single-host path: it plans, runs the outstanding units in place and merges.
When the engine is on a different machine from the run folder, the same work is split
into `dock-export` (bundle the outstanding units), `dock-worker` (run them, no run folder
needed) and `dock-merge` (fold the results back in). See `bundle_*` below and the
polyphar.stages.dock_export / dock_worker / dock_merge modules.
"""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from polyphar.stages import dock_unidock as unidock
from polyphar.stages import ligands
from polyphar.stages.fix import H_SEED, reference_platform, seeded_placement, strip_date_remark

STAGE = "dock"
DEFAULT_SEED = 42         # BASE seed: replicate i runs with base_seed + i
DEFAULT_REPLICATES = 3
DEFAULT_TPJ = 1
EXHAUSTIVENESS = 12
NUM_MODES = 12
LIGAND_ADDH = True   # smina --addH: add ligand hydrogens (smina's default, pinned explicitly)
DEFAULT_UNIT_TIMEOUT = 43200   # seconds per engine call (12 h)
AF_H_SEED = H_SEED   # seeds OpenMM's random hydrogen placement (deterministic output)
AF_PH = 7.0   # AlphaFold receptors: pdbfixer addMissingHydrogens, as the fix stage
TITLE_PREFIX = "EDRUG3D_"
# Per-replicate raw poses. `RESULT_COLUMNS = ["Pose", "SMINA_Score", ...]` is RETIRED: the
# engine-specific column name does not survive a selectable engine, and one score column
# cannot express a distribution over replicates (see SUMMARY_COLUMNS).
POSE_COLUMNS = ["Pose", "Score", "Ligand", "ID", "Pocket"]
# df_results_<label>.csv: the summary over replicates. `n` is how many replicates
# contributed, `score_sd` the SAMPLE standard deviation (ddof=1) of the per-replicate best
# scores and `score_se` = score_sd / sqrt(n). Both are empty at n == 1, because a single
# draw has no spread -- not 0, which would claim it had none.
SUMMARY_COLUMNS = ["Ligand", "ID", "Pocket", "n", "score_mean", "score_sd", "score_se",
                   "score_min", "score_max"]
LIBRARY_COLUMNS = ["ID", "name", "cas", "status", "record_index"]
END = "$$$$"

# score field written by each engine (smina's own, and the one unidock_tools adds). Uni-Dock
# results are keyed on <ID>, never on the title: RDKit may rewrite the title line.
#
# This is NOT a first-match preference list, and must not become one again. `unidock_tools`
# clears only five of its own SD properties and passes every other INPUT field through, so a
# library whose records already carry <minimizedAffinity> (any SDF that has been through
# smina, `obabel --minimize`, or an earlier docking round) would hand every Uni-Dock record a
# stale score next to the real <docking_score>. Taking the first present field replaced every
# engine score with the library's, invisibly to the absurd-score filter. The score field is
# therefore chosen by ENGINE (SCORE_FIELD_OF), a record carrying both is refused rather than
# guessed, and `dock_unidock.split_library` refuses such a library before any GPU time.
SCORE_FIELDS = ("minimizedAffinity", "docking_score")
SCORE_FIELD_OF = {"smina": "minimizedAffinity", "unidock": "docking_score"}

ENGINE_SMINA = "smina"
ENGINE_UNIDOCK = "unidock"
ENGINES = (ENGINE_SMINA, ENGINE_UNIDOCK)
# Uni-Dock on a GPU is the production route: 0.176 s/ligand measured on a Colab L4 and
# flat across a 5x box range, against ~200 core-seconds on 64 CPU cores for smina. smina
# stays available with --engine smina as the single-host / development path.
#
# `control` inherits this default, and must: its library-percentile readout compares the
# co-crystal ligand against the PRODUCTION library at the same target, so a control run on
# a different engine than `dock` compares two score scales (Uni-Dock reads +0.326 kcal/mol
# less negative than smina) and biases the percentile.
DEFAULT_ENGINE = ENGINE_UNIDOCK

POSES_NONE, POSES_FIRST, POSES_ALL = "none", "first", "all"
POSE_POLICIES = (POSES_NONE, POSES_FIRST, POSES_ALL)
# Replicate 0's poses feed the control stage's redocking RMSD (the gate) and its
# pose-in-pocket diagnostic -- whether the best-scoring pose is the one the
# crystallographer saw, and whether it sits in the selected cavity at all, the assumption
# the pocket-similarity premise rests on. Three replicates x 482 targets x 1,960 ligands
# of poses is far more than useful, so
# only r0's are kept by default.
DEFAULT_POSES = POSES_FIRST


def smina_bin() -> str:
    return os.environ.get("SMINA_BIN") or "smina"


def _engine_bin(args) -> str:
    return (unidock.tools_bin() if engine_name(args) == ENGINE_UNIDOCK else smina_bin())


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
    _atomic_text(Path(pdb_out), strip_date_remark(buf.getvalue()))


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


def _engine_reason(reason: str) -> str:
    """A failure message without the engine's exception class in front of it: failed.tsv
    reads the same whichever engine produced the failure."""
    for prefix in ("SminaError: ", "UniDockError: "):
        reason = reason.removeprefix(prefix)
    return reason


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
              exhaustiveness: int = EXHAUSTIVENESS, num_modes: int = NUM_MODES,
              add_h: bool = LIGAND_ADDH, timeout: float | None = None) -> None:
    """Dock `ligands` into the box; out_sdf and log are written atomically. Raises
    SminaError on a non-zero exit or after `timeout` seconds (smina is killed). On
    failure no out_sdf is written and the log holds smina's log so far plus its stderr."""
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
            "--exhaustiveness", str(exhaustiveness), "--cpu", str(cpu), "--seed", str(seed),
            "--addH", "1" if add_h else "0"]
    stderr = ""
    try:
        try:
            p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                               text=True, errors="replace", timeout=timeout)
        except subprocess.TimeoutExpired as exc:   # subprocess.run killed smina
            err = exc.stderr or b""
            stderr = err.decode(errors="replace") if isinstance(err, bytes) else err
            raise SminaError(f"timeout after {timeout:g} s") from None
        stderr = p.stderr or ""
        if p.returncode != 0:
            tail = " | ".join(stderr.strip().splitlines()[-5:])
            raise SminaError(f"smina exited with {p.returncode}: {tail}")
        if not tmp_out.exists():
            raise SminaError("smina wrote no output file")
        if tmp_log.exists():
            os.replace(tmp_log, log)
        os.replace(tmp_out, out_sdf)
    except SminaError as exc:
        so_far = tmp_log.read_text(errors="replace") if tmp_log.exists() else ""
        _atomic_text(log, f"{so_far}\n--- polyphar: {exc} ---\n--- smina stderr ---\n{stderr}")
        raise
    finally:
        tmp_out.unlink(missing_ok=True)
        tmp_log.unlink(missing_ok=True)


def _pose_id(rec: Record, names: dict, out_sdf) -> str:
    """The library ID of one pose record.

    The `<ID>` SD field comes FIRST: `unidock_tools` runs the ligand through RDKit and MAY
    rewrite the title line, so the EDRUG3D_<ID> title that identifies a smina pose cannot be
    relied on. The input SD fields are preserved. The title is the fallback, for engines
    (smina) that pass the record through verbatim and write no `<ID>`."""
    lid = rec.fields.get("ID")
    if lid and str(lid) in names:
        return str(lid)
    title = rec.lines[0].strip()
    if title.startswith(TITLE_PREFIX) and title[len(TITLE_PREFIX):] in names:
        return title[len(TITLE_PREFIX):]
    raise ValueError(f"{out_sdf}: pose {title!r} has no <ID> field and no library-ID title")


def _pose_score(rec: Record, out_sdf, what: str, score_field: str | None = None) -> float:
    """The pose's score. `score_field` pins it to the engine that wrote the file; without
    one, exactly one of SCORE_FIELDS must be present -- a record carrying both is ambiguous
    and is refused, never silently resolved in favour of either."""
    if score_field is not None:
        if score_field not in rec.fields:
            raise ValueError(f"{out_sdf}: {what} has no <{score_field}> field")
        return float(rec.fields[score_field])
    present = [f for f in SCORE_FIELDS if f in rec.fields]
    if len(present) == 1:
        return float(rec.fields[present[0]])
    if not present:
        raise ValueError(f"{out_sdf}: {what} has none of the score fields "
                         f"{', '.join(SCORE_FIELDS)}")
    raise ValueError(f"{out_sdf}: {what} carries more than one score field "
                     f"({', '.join(present)}); which one the engine wrote cannot be told "
                     f"apart from a stale value the input library carried in")


def parse_poses(out_sdf, library, label: str, max_poses: int = NUM_MODES,
                score_field: str | None = None):
    """DataFrame[Pose, Score, Ligand, ID, Pocket] from an engine's output SDF.

    Engine-independent: a ligand is identified by its `<ID>` field (Uni-Dock) or, failing
    that, by an EDRUG3D_<ID> title (smina), and the score is `minimizedAffinity` (smina) or
    `docking_score` (unidock_tools) -- NOT the engine's own `ENERGY=` remark, which
    `unidock_tools` does not write. Pose counts 1.. within each ligand's consecutive
    records; each ligand's poses must be one contiguous run of at most max_poses records,
    otherwise ValueError."""
    import pandas as pd

    names = dict(zip(library["ID"].astype(str), library["name"]))
    rows, prev, pose, seen = [], None, 0, set()
    for rec in _records(_read_text(out_sdf)):
        lid = _pose_id(rec, names, out_sdf)
        if lid != prev and lid in seen:
            raise ValueError(f"{out_sdf}: the poses of ligand ID {lid} are not contiguous")
        seen.add(lid)
        pose = pose + 1 if lid == prev else 1
        prev = lid
        if pose > max_poses:
            raise ValueError(f"{out_sdf}: ligand ID {lid} has more than {max_poses} poses")
        rows.append([pose, _pose_score(rec, out_sdf, f"pose {pose} of ligand ID {lid}",
                                       score_field),
                     names[lid], lid, label])
    return pd.DataFrame(rows, columns=POSE_COLUMNS)


def summarize(frames: dict, library, label: str):
    """DataFrame[SUMMARY_COLUMNS] over the replicates in `frames` ({replicate: poses}).

    One observation per (ligand, replicate): that replicate's BEST (lowest) score for the
    ligand. The summary is over those draws, so `n` is the number of replicates that
    contributed, never the number of poses.

    `score_sd` is the sample standard deviation (ddof=1) and `score_se` = sd / sqrt(n). A
    previous attempt used mean|delta| / sqrt(N), which is NOT a standard error: for a normal,
    E|X - Y| = 1.128 sigma, so that figure is biased by ~13 % and is not the sd of anything.
    The raw per-replicate draws stay on disk precisely so any other summary can be
    recomputed without re-running the engine."""
    import math

    import pandas as pd

    names = dict(zip(library["ID"].astype(str), library["name"]))
    best: dict[str, list[float]] = {}
    for _i, df in sorted(frames.items()):
        for lid, score in df.groupby("ID")["Score"].min().items():
            best.setdefault(str(lid), []).append(float(score))
    rows = []
    for lid in library["ID"].astype(str):
        draws = best.get(lid, [])
        if not draws:
            continue
        n = len(draws)
        mean = sum(draws) / n
        if n > 1:
            sd = math.sqrt(sum((x - mean) ** 2 for x in draws) / (n - 1))
            se = sd / math.sqrt(n)
        else:
            sd = se = float("nan")
        rows.append([names[lid], lid, label, n, mean, sd, se, min(draws), max(draws)])
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


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


def add_engine_options(parser) -> None:
    """The engine's own options, with no ligand set attached to them.

    Split out of `add_arguments` so `polyphar control` -- which runs the SAME engine over
    its OWN ligand set (the co-crystal ligands) and its own replicate count -- can offer
    them with byte-identical argparse signatures. `polyphar run` merges the stages'
    options and REFUSES two stages that define one option differently, so sharing the
    definition rather than copying it is what keeps `run` working."""
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"BASE seed; replicate i runs with base + i and records it "
                             f"(provenance, not reproducibility) (default {DEFAULT_SEED})")
    parser.add_argument("--threads-per-job", type=_positive_int, default=DEFAULT_TPJ,
                        help=f"engine threads per unit; nprocs // this units run at once "
                             f"(default {DEFAULT_TPJ})")
    parser.add_argument("--unit-timeout", type=_positive_int, default=DEFAULT_UNIT_TIMEOUT,
                        help=f"seconds before one unit's engine call is killed and the unit "
                             f"failed (default {DEFAULT_UNIT_TIMEOUT})")
    parser.add_argument("--engine", choices=ENGINES, default=DEFAULT_ENGINE,
                        help=f"docking engine (default {DEFAULT_ENGINE}); 'unidock' is "
                             f"Uni-Dock >= 1.2.0 on a CUDA GPU through unidock_tools, "
                             f"'smina' the CPU single-host path. The engine is part of "
                             f"every unit key, so changing it supersedes recorded draws")
    parser.add_argument("--search-mode", default=unidock.DEFAULT_SEARCH_MODE,
                        choices=unidock.SEARCH_MODES,
                        help=f"Uni-Dock search mode (default "
                             f"{unidock.DEFAULT_SEARCH_MODE} = exhaustiveness 384); "
                             f"ignored by --engine smina")
    parser.add_argument("--batch-size", type=_positive_int,
                        default=unidock.DEFAULT_BATCH_SIZE,
                        help=f"Uni-Dock ligands per GPU batch (default "
                             f"{unidock.DEFAULT_BATCH_SIZE}); ignored by --engine smina")


def add_arguments(parser) -> None:
    """Options shared by `dock`, `dock-export`, `dock-worker` and `dock-merge`: they all
    have to agree on what a unit is, or the keys they compute would not match."""
    add_engine_options(parser)
    parser.add_argument("--library", default=DEFAULT_LIBRARY,
                        help=f"ligand SDF, relative to the run folder (default {DEFAULT_LIBRARY})")
    parser.add_argument("--replicates", type=_positive_int, default=DEFAULT_REPLICATES,
                        help=f"independent draws per target (default {DEFAULT_REPLICATES}); "
                             f"raising it adds units, it never re-runs existing ones")
    parser.add_argument("--poses", choices=POSE_POLICIES, default=DEFAULT_POSES,
                        help=f"which replicates' poses to write out (default "
                             f"{DEFAULT_POSES}: replicate 0 only, which is what the "
                             f"pose-in-pocket validation needs)")


def engine_name(args) -> str:
    """Which engine this invocation uses ('smina' / 'unidock')."""
    return str(getattr(args, "engine", None) or DEFAULT_ENGINE)


def engine_settings(args) -> dict:
    """The engine's own search settings -- part of what a draw IS, so part of the unit key.

    smina's exhaustiveness/num_modes are a CPU figure and are deliberately NOT mirrored
    onto Uni-Dock, whose `--search_mode balance` is exhaustiveness 384."""
    if engine_name(args) == ENGINE_UNIDOCK:
        return {"search_mode": str(getattr(args, "search_mode", None)
                                   or unidock.DEFAULT_SEARCH_MODE),
                "batch_size": int(getattr(args, "batch_size", None)
                                  or unidock.DEFAULT_BATCH_SIZE)}
    return {"exhaustiveness": EXHAUSTIVENESS, "num_modes": NUM_MODES}


def calibration_tool(args) -> str:
    """Which calibration entry this invocation's ETA may use.

    'smina' for the CPU path; 'unidock:<device class>' for Uni-Dock, because the rate is
    per device -- an L4 figure is not an A100's. The existing single seconds-per-ligand
    smina calibration is therefore never borrowed by a GPU run (and the Uni-Dock rate never
    by a CPU one): a mismatch prints "uncalibrated for ..." instead of a wrong ETA."""
    if engine_name(args) != ENGINE_UNIDOCK:
        return ENGINE_SMINA
    from polyphar.benchmark import unidock_rate, unidock_tool_prefix

    return unidock_tool_prefix() + unidock_rate(unidock.gpu_model())[1]


def require_engine_host(args) -> None:
    """Refuse to start when this host cannot run the chosen engine.

    `--engine unidock` on a host with no GPU fails HERE, with the guidance for the split
    export/worker/merge path, rather than after the stage has planned, prepared receptors
    and launched a call that crashes. POSEIDON -- where the run folders are -- is exactly
    such a host, so this is the common case."""
    if engine_name(args) != ENGINE_UNIDOCK:
        return
    unidock.require_gpu()
    unidock.check_version()


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


def _library_path(layout, args, copy: bool = True, progress=None) -> Path:
    """The library to dock: the prepared one from the `ligands` stage when that stage has
    finished on this source library, else the source library itself (see
    _raw_library_path). The prepared file's md5 then goes into the stage fingerprint like
    any other library, so a re-preparation re-docks.

    Falling back to the raw library is never silent once a `ligands` manifest exists:
    `progress` (passed by `run`) gets a stderr warning and a priority-1 notification. The
    fallback changes library_md5, so every target is re-docked and its existing results
    replaced — by scores from an unprepared library, which are not comparable with the
    prepared ones."""
    from polyphar.stages.ligands import prepared_state

    raw = _raw_library_path(layout, args, copy=copy)
    prepared, why = prepared_state(layout, raw)
    if prepared is not None:
        return prepared
    if why:
        _warn_raw_library(progress, why, raw)
    return raw


def _warn_raw_library(progress, why: str, raw: Path) -> None:
    import sys

    message = (f"a `ligands` stage has run here but its prepared library cannot be used: "
               f"{why}. Docking the RAW library {raw.name} instead: every target will be "
               f"re-docked and its existing results REPLACED, and the new scores will NOT "
               f"be comparable with scores from the prepared library. Rerun "
               f"`polyphar ligands` (or restore 07_docking/library/) first if that matters.")
    print(f"{STAGE}: WARNING: {message}", file=sys.stderr, flush=True)
    notifier = getattr(progress, "notifier", None)
    if notifier is None:
        return
    try:
        notifier(f"{STAGE}: docking the UNPREPARED library", message, priority=1)
    except Exception:
        pass


def _library_is_prepared(layout, args, lib: Path) -> bool:
    """Whether `lib` is the `ligands` stage's prepared library (for dock's manifest, so
    Methods can say which library produced the scores without path archaeology).

    False -- not an error -- when the source library named by --library is not there at
    all: a caller that supplied its own library (`run(..., library=...)`, `polyphar
    control`) need not have e-Drug3D in the run folder, and the honest answer for a
    hand-supplied library is "no, the `ligands` stage did not prepare this"."""
    from polyphar.stages.ligands import prepared_library

    try:
        source = _raw_library_path(layout, args, copy=False)
    except FileNotFoundError:
        return False
    prepared = prepared_library(layout, source)
    return prepared is not None and Path(prepared) == Path(lib)


def _raw_library_path(layout, args, copy: bool = True) -> Path:
    """The source library inside <run>/inputs; one given elsewhere is copied there (never
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


# --- units: (target, replicate) --------------------------------------------------------

STATE_DONE = "done"           # ran under this key, the record of it is intact: SKIP
STATE_MISSING = "missing"     # never ran under this key: run it
STATE_CORRUPT = "corrupt"     # ran under this key, the record is unusable: discard + re-draw

LIB_DIRNAME = "ligands"       # 07_docking/ligands/ (this stage's engine input, see RESERVED_DIRS)
TITLED_LIBRARY = "library_titled.sdf"
DISCARD_LOG = "discarded.tsv"
DISCARD_DIRNAME = "discarded"   # 07_docking/discarded/<ts>/<label>/r<i>.sdf (archive)
EXPORT_DIRNAME = "export"       # 07_docking/export/<ts>/ (dock_export's default bundle dir)


@dataclass(frozen=True)
class Unit:
    """One (target, replicate): one engine invocation over the whole library.

    `seed` is base_seed + replicate, recorded for provenance. It is NOT a claim that the
    unit can be reproduced from it: the GPU engine this stage targets scores the same
    ligand differently across identical same-seed repeats."""

    label: str
    replicate: int
    seed: int
    key: str

    @property
    def name(self) -> str:
        return f"r{self.replicate}"


def seed_for(base_seed: int, replicate: int) -> int:
    return int(base_seed) + int(replicate)


def _target_dir(layout, t) -> Path:
    label = t if isinstance(t, str) else t.label
    return layout.docking / label


def unit_paths(layout, u: Unit) -> tuple[Path, Path, Path]:
    """(result SDF, log, done marker), all three in <label>/replicates/.

    The marker lives beside the record it certifies, not in a sibling directory: the spec's
    return step is "pack only new results", and packing <label>/replicates/ alone must not
    turn a drawn unit into an undrawn one. A copy, a tar or an `rm -rf` that treats the
    directory as the unit now moves the record and its certificate together."""
    d = _target_dir(layout, u.label) / "replicates"
    return d / f"{u.name}.sdf", d / f"{u.name}.log", d / f"{u.name}.done"


def result_intact(sdf) -> bool:
    """Whether a unit's result file is a usable record of its draw: present, non-empty and
    terminated by the SDF end marker. Engine output is written atomically, so what this
    catches is a file deleted, truncated or half-copied out of a bundle."""
    sdf = Path(sdf)
    try:
        if not sdf.is_file() or sdf.stat().st_size == 0:
            return False
        with open(sdf, "rb") as fh:
            fh.seek(max(0, sdf.stat().st_size - 4096))
            tail = fh.read().decode("utf-8", "replace")
    except OSError:
        return False
    lines = [l for l in tail.splitlines() if l.strip()]
    return bool(lines) and lines[-1].strip() == END


def unit_state(layout, u: Unit) -> str:
    """STATE_DONE / STATE_MISSING / STATE_CORRUPT. Reads only; it never discards anything
    (so `units()` can call it for the ETA without side effects)."""
    sdf, _, done = unit_paths(layout, u)
    try:
        marked = done.read_text().strip() == u.key
    except OSError:
        marked = False
    if not marked:
        return STATE_MISSING
    return STATE_DONE if result_intact(sdf) else STATE_CORRUPT


def log_discard(stage_dir, u: Unit, reason: str, archived: str = "") -> None:
    """Append one discarded draw to <stage_dir>/discarded.tsv. Append-only and never
    truncated, across runs: a discarded sample was a real draw and dropping it is a (small)
    selection effect, so the record of it has to outlive the run that did it."""
    import time

    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    path = stage_dir / DISCARD_LOG
    clean = lambda x: str(x).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    _ensure_discard_log(stage_dir)
    # One O_APPEND write of one short line: the discards of a run now happen in the worker
    # processes (deferred to the moment each replacement is written), so several of them
    # append here at once. The header is created exclusively by _ensure_discard_log, or two
    # workers would each write one and the second header would count as a discarded draw.
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')}\t{clean(u.label)}\t{u.replicate}\t"
                 f"{u.seed}\t{u.key}\t{clean(reason)}\t{clean(archived)}\n")


def _ensure_discard_log(stage_dir) -> Path:
    """Create <stage_dir>/discarded.tsv with its header, exactly once, whoever gets there
    first (O_EXCL). Never truncates an existing log."""
    path = Path(stage_dir) / DISCARD_LOG
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return path
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(DISCARD_COLUMNS + "\n")
    return path


DISCARD_COLUMNS = "timestamp\tlabel\treplicate\tseed\tunit_key\treason\tarchived"


def _archive_record(layout, u: Unit) -> str:
    """Move a unit's result file out of the way instead of deleting it, to
    07_docking/discarded/<timestamp>/<label>/r<i>.sdf (its log goes with it). Returns the
    archived SDF's path relative to 07_docking, or "" when there was no file.

    A docking result is a draw from a non-deterministic engine: destroying one is
    irreversible. Orphan labels are already "moved, never deleted"; a per-unit draw is a
    smaller thing to lose, not a less valuable one."""
    import time

    sdf, log, _ = unit_paths(layout, u)
    if not sdf.exists():
        return ""
    root = Path(layout.docking) / DISCARD_DIRNAME
    stamp, k = time.strftime("%Y%m%dT%H%M%S"), 1
    dest = root / stamp / u.label
    while (dest / f"{u.name}.sdf").exists():
        k += 1
        dest = root / f"{stamp}_{k}" / u.label
    dest.mkdir(parents=True, exist_ok=True)
    os.replace(sdf, dest / f"{u.name}.sdf")
    if log.is_file():
        os.replace(log, dest / f"{u.name}.log")
    prov = prov_path(layout, u)
    if prov.is_file():   # the provenance travels with the draw it describes
        os.replace(prov, dest / (u.name + PROV_SUFFIX))
    return (dest / f"{u.name}.sdf").relative_to(Path(layout.docking)).as_posix()


def discard(layout, u: Unit, reason: str, archived: str | None = None) -> None:
    """Take a unit's record off the books so it can be re-drawn under the SAME replicate
    index, and log it. The only legitimate repeat of a completed replicate.

    The record itself is ARCHIVED (moved under 07_docking/discarded/), never deleted; only
    the marker goes, because a marker with no record is a false certificate. `archived`
    overrides where the draw can still be found (for a bundle's draw that was dropped: it
    stays in the bundle)."""
    _, _, done = unit_paths(layout, u)
    where = _archive_record(layout, u) if archived is None else archived
    done.unlink(missing_ok=True)
    log_discard(layout.docking, u, reason, where)


def _df_path(layout, t) -> Path:
    """The summary over replicates. The name drops the legacy `_fixed` suffix: the file is
    no longer replicate 0's poses but a distribution (SUMMARY_COLUMNS), and a downstream
    reader must not pick it up as if it were the old one. Old `_fixed` files are swept into
    obsolete/ by the orphan sweep (_DF_RE matches both spellings)."""
    label = t if isinstance(t, str) else t.label
    return layout.docking / f"df_results_{label}.csv"


def _openbabel_version() -> str:
    import importlib.metadata

    try:
        return importlib.metadata.version("openbabel")
    except Exception:
        return "unknown"


def _prep_stamp(t: Target) -> str:
    """What 07_docking/<label>/receptor.pdbqt was built from: the source key PLUS the Open
    Babel version that wrote it.

    Open Babel is deliberately NOT in `_src_key` and so not in any unit key: putting it
    there would invalidate every recorded draw on a Babel upgrade. It IS in this stamp, so
    an upgrade re-prepares the receptor instead of a resumed folder and a fresh one feeding
    the engine different PDBQTs under identical unit keys. A changed PDBQT sha256 is
    warned about against the previous manifest (see `_warn_receptor_drift`)."""
    return _sha_text(_json({"src": _src_key(t), "openbabel": _openbabel_version()}))


def _src_key(t: Target) -> str:
    """What the receptor PDBQT depends on (it is byte-deterministic given these)."""
    from polyphar.common import sha256_file

    return _sha_text(_json({"receptor_sha256": sha256_file(t.receptor),
                                 "prep": t.prep, "prep_version": PREP_VERSION}))


def _json(obj) -> str:
    import json

    return json.dumps(obj, sort_keys=True)


def _unit_key(src_key: str, pocket_sha: str, library_sha: str, replicate: int, seed: int,
              opts: dict) -> str:
    """The identity of one (target, replicate) draw: receptor, pocket, library, REPLICATE
    INDEX, that replicate's seed, the search parameters and the engine version.

    Two replicates of one target differ in `replicate` and `seed`, so they are different
    units with different keys and different result files. The key says which draw this is,
    not what it scored."""
    return _sha_text(_json({"receptor": src_key, "pocket_sha256": pocket_sha,
                            "library_sha256": library_sha, "replicate": int(replicate),
                            "seed": int(seed),
                            # how MANY draws the run wants is not part of what ONE draw is:
                            # raising --replicates must add units, never invalidate the
                            # recorded ones (it is still in the stage fingerprint).
                            **{k: v for k, v in opts.items() if k != "replicates"}}))


def _options(args, engine: str) -> dict:
    """What every unit of this run shares. `base_seed`, not `seed`: the seed is per unit.

    `engine_name` and the engine's own search settings are in here, so they are in every
    unit key: changing the engine, or its search mode, makes the recorded draws draws of
    something else, and they are correctly superseded rather than mixed in."""
    return {"base_seed": int(getattr(args, "seed", DEFAULT_SEED)),
            "replicates": int(getattr(args, "replicates", DEFAULT_REPLICATES)),
            "ligand_addH": LIGAND_ADDH, "engine": engine,
            "engine_name": engine_name(args), **engine_settings(args)}


@dataclass
class _Plan:
    targets: list[Target]
    library_sha: str
    n_ligands: int
    usable: list[Target]             # receptor and pocket present
    missing: list[tuple[str, str]]   # (label, reason) for unusable targets
    units: list[Unit]                # every (usable target, replicate)
    opts: dict

    def units_of(self, label: str) -> list[Unit]:
        return [u for u in self.units if u.label == label]

    def outstanding(self, layout) -> list[tuple[Unit, str]]:
        """(unit, state) for every unit that is not STATE_DONE. A done unit never appears
        here: the stage does not re-run a recorded draw."""
        states = ((u, unit_state(layout, u)) for u in self.units)
        return [(u, s) for u, s in states if s != STATE_DONE]


def _plan(layout, args, lib: Path, engine: str) -> _Plan:
    from polyphar.common import sha256_file

    targets = read_targets(layout)
    records = _library_records(lib)
    library_sha = _sha_text("".join(_titled(r) for r in records))
    opts = _options(args, engine)
    base = opts["base_seed"]
    usable, missing, units = [], [], []
    for t in targets:
        absent = [p for p in (t.receptor, t.pocket) if not p.is_file()]
        if absent:
            missing.append((t.label, "missing file: " + ", ".join(
                _rel(p, layout.run_dir) for p in absent)))
            continue
        usable.append(t)
        src, pocket = _src_key(t), sha256_file(t.pocket)
        for i in range(opts["replicates"]):
            seed = seed_for(base, i)
            units.append(Unit(t.label, i, seed,
                              _unit_key(src, pocket, library_sha, i, seed, opts)))
    return _Plan(targets, library_sha, len(records), usable, missing, units, opts)


def _prepare_receptor(layout, t: Target) -> Path:
    """07_docking/<label>/receptor.pdbqt, reused while its source key is unchanged."""
    d = _target_dir(layout, t)
    d.mkdir(parents=True, exist_ok=True)
    pdbqt, stamp = d / "receptor.pdbqt", d / "receptor.key"
    key = _prep_stamp(t)
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

def _dock_unit(job: dict) -> int:
    """Run one unit and, only once its result file AND its provenance are in place, write
    the marker. The marker is the claim "this draw was made and recorded", so it must never
    precede the record.

    When a record of a superseded or uncertified draw is standing in this unit's place, it
    is archived and logged HERE -- at the last moment before the engine call that replaces
    it. That is the whole of P1: a pass that plans this unit but never reaches this point
    (filtered, aborted, or merging rather than executing) costs no draws. The archive is
    taken before the engine runs rather than after, because `run_smina` renames its output
    over the record; a failed engine call therefore leaves the old draw in
    07_docking/discarded/ rather than destroyed."""
    unit, docking_dir = job["unit"], job["docking_dir"]
    if job.get("discard_reason") is not None:
        discard(_ArchiveLayout(docking_dir), unit, job["discard_reason"])
    sdf, log = Path(job["sdf"]), Path(job["log"])
    if job["engine_name"] == ENGINE_UNIDOCK:
        unidock.run_unidock(job["receptor"], job["ligands"], job["center"], job["size"],
                            sdf, log, seed=unit.seed,
                            search_mode=job["settings"]["search_mode"],
                            batch_size=job["settings"]["batch_size"],
                            timeout=job["timeout"])
    else:
        run_smina(job["receptor"], job["ligands"], job["center"], job["size"], sdf, log,
                  seed=unit.seed, cpu=job["cpu"], timeout=job["timeout"])
    _write_provenance(Path(job["prov"]), job, sdf)
    _atomic_text(Path(job["done"]), unit.key + "\n")
    return 0


# --- per-unit provenance ------------------------------------------------------------------
#
# Written beside the record, by whoever actually ran the engine, and carried with it into
# the run folder by `import_bundle`. It is the only place that can know the GPU model, the
# CUDA runtime and the engine build: the planning host (`dock-export`) and the merging host
# (`dock-merge`) need have none of them. It also closes L3 -- without `result_sha256` a
# result SDF edited to different scores but still `$$$$`-terminated was copied in unchecked.

PROV_SUFFIX = ".prov.json"


def prov_path(layout, u: Unit) -> Path:
    sdf = unit_paths(layout, u)[0]
    return sdf.with_name(u.name + PROV_SUFFIX)


def read_provenance(layout, u: Unit) -> dict:
    import json

    try:
        info = json.loads(prov_path(layout, u).read_text())
    except (OSError, ValueError):
        return {}
    return info if isinstance(info, dict) else {}


def _write_provenance(path: Path, job: dict, sdf: Path) -> None:
    """Everything the Methods need per unit, recorded where the draw was made."""
    import platform

    from polyphar.common import sha256_file

    unit = job["unit"]
    name = job["engine_name"]
    info = {"label": unit.label, "replicate": unit.replicate, "seed": unit.seed,
            "unit_key": unit.key, "engine_name": name,
            "engine_version": job.get("engine_version") or "unknown",
            "settings": job.get("settings") or {},
            "library_sha256": job.get("library_sha256") or "",
            "receptor_sha256": job.get("receptor_sha256") or "",
            "pocket_sha256": job.get("pocket_sha256") or "",
            "result_sha256": sha256_file(sdf),
            "host": job.get("host") or host_path(),
            "hostname": platform.node()}
    if name == ENGINE_UNIDOCK:
        info.update({"gpu_model": unidock.gpu_model() or "unknown",
                     "cuda_driver_max": unidock.cuda_driver_max(),
                     "unidock_tools": unidock.tools_version()})
    else:
        info.update({"gpu_model": "", "cuda_driver_max": "", "unidock_tools": ""})
    _atomic_text(path, _json_indent(info))


def host_path() -> str:
    """Where the draw was made: Colab / TRUBA / a container image digest, as the host
    declares it in $POLYPHAR_HOST_PATH, else the container digest or the hostname."""
    import platform

    return (os.environ.get("POLYPHAR_HOST_PATH")
            or os.environ.get("POLYPHAR_IMAGE_DIGEST") or platform.node() or "unknown")


def verify_result_digest(layout, u: Unit) -> str | None:
    """None when the record matches the digest its provenance recorded, else why not.

    No provenance (a draw made before this existed) is not a mismatch: it is only an
    absence of evidence, and refusing those would discard every legacy record."""
    from polyphar.common import sha256_file

    info = read_provenance(layout, u)
    want = info.get("result_sha256")
    if not want:
        return None
    sdf = unit_paths(layout, u)[0]
    if not sdf.is_file():
        return "the result file is missing"
    got = sha256_file(sdf)
    if got != want:
        return (f"the result file does not match the sha256 its provenance recorded "
                f"({got[:12]} vs {want[:12]}): it was changed after the engine wrote it")
    return None


class _ArchiveLayout:
    """The minimum `discard` needs (a 07_docking directory), picklable for a worker
    process: neither RunLayout nor BundleLayout has to cross the process boundary."""

    def __init__(self, docking):
        self.docking = Path(docking)
        self.run_dir = self.docking


def _unit_worker(job) -> int:
    from polyphar.pool import task

    return task(lambda: _dock_unit(job))


def _slots(args, tpj: int) -> int:
    """How many units run at once. ONE for Uni-Dock: it batches the whole library onto the
    GPU and already saturates it, so a second concurrent unit contends for the same device
    rather than adding throughput (and the authors' own guidance is one invocation per
    device with as many ligands as possible)."""
    if engine_name(args) == ENGINE_UNIDOCK:
        return 1
    return max(1, int(getattr(args, "nprocs", 1)) // max(1, tpj))


def slots(args, threads_per_job: int = 1) -> int:
    """The stage's parallel slots, for `cli.gate` (eta.stage_slots) as well as execution.

    THE SAME number that `_run_pool` runs at: it scales the ETA and the seconds-per-unit
    written back to calibration.json. Deriving it from nprocs alone (the engine-blind
    eta.PARALLEL_EFFICIENCY table) made `--engine unidock --nprocs 32` quote an ETA 32x too
    short and then store a 32x-inflated rate under the device's calibration key."""
    return _slots(args, int(threads_per_job or 1))


def _unit_dirs(layout, u: Unit) -> None:
    for p in unit_paths(layout, u):
        p.parent.mkdir(parents=True, exist_ok=True)


def _engine_input(layout, args, plan, titled_sdf: Path) -> Path:
    """What the engine is handed: the one titled SDF for smina, the split library's index
    file for Uni-Dock (`unidock_tools` writes one output per INPUT FILE, so it cannot take
    a single multi-record library)."""
    if engine_name(args) != ENGINE_UNIDOCK:
        return titled_sdf
    return unidock.split_library(titled_sdf, Path(titled_sdf).parent / "split",
                                 plan.library_sha)


def _unit_job(u: Unit, args, engine: str, receptor, engine_input, center, size, docking_dir,
              timeout, cpu: int, discard_reason, library_sha: str, pocket_sha: str) -> dict:
    """Everything one unit's engine call needs, in one picklable dict.

    One dict rather than a dozen positional arguments: the provenance the worker has to
    record (engine build, GPU model, settings, input digests) grew the payload past what a
    positional signature can carry readably, and every field here ends up in the unit's
    provenance file."""
    from polyphar.common import sha256_file

    return {"unit": u, "engine_name": engine_name(args), "engine_version": engine,
            "settings": engine_settings(args), "receptor": str(receptor),
            "receptor_sha256": sha256_file(receptor), "ligands": str(engine_input),
            "library_sha256": library_sha, "pocket_sha256": pocket_sha,
            "center": center, "size": size, "cpu": int(cpu), "timeout": timeout,
            "docking_dir": str(docking_dir), "discard_reason": discard_reason,
            "host": host_path(),
            "sdf": str(Path(docking_dir) / u.label / "replicates" / f"{u.name}.sdf"),
            "log": str(Path(docking_dir) / u.label / "replicates" / f"{u.name}.log"),
            "prov": str(Path(docking_dir) / u.label / "replicates"
                        / (u.name + PROV_SUFFIX)),
            "done": str(Path(docking_dir) / u.label / "replicates" / f"{u.name}.done")}


# --- merge -----------------------------------------------------------------------------

def _check_replicate(layout, t: Target, u: Unit, library):
    """(DataFrame, reason). Parse one replicate's record and check that EVERY library
    ligand has poses in it.

    `result_intact` only looks at the last line, so a record truncated on a clean record
    boundary passes it: 6 of 2,083 ligands with a final `$$$$` is "intact" and, before
    this check existed, entered the mean over replicates as a complete draw. Completeness
    is per replicate, not a property of r0."""
    sdf, log = unit_paths(layout, u)[:2]
    # The unit's provenance records which engine drew it, so the score field is chosen by
    # the engine rather than by which field happens to appear first in the record.
    score_field = SCORE_FIELD_OF.get(read_provenance(layout, u).get("engine_name"))
    try:
        df = parse_poses(sdf, library, t.label, score_field=score_field)
    except ValueError as exc:
        return None, (f"{u.name}: unusable result: {exc}; "
                      f"log: {_rel(log, layout.run_dir)}")
    missing = [i for i in library["ID"] if i not in set(df["ID"])]
    if missing:
        return None, (f"incomplete: {len(missing)} ligand ID(s) without poses in {u.name}: "
                      f"{','.join(missing)}; log: {_rel(log, layout.run_dir)}")
    return df, None


def _pose_files(layout, t) -> list[Path]:
    """Every poses file a target could have, under any --poses policy (so a policy change
    can clear the ones it no longer wants)."""
    d = _target_dir(layout, t)
    return sorted(set([d / "poses.sdf", *d.glob("poses_r*.sdf")]))


def _clear_outputs(layout, t) -> None:
    for p in [*_pose_files(layout, t), _df_path(layout, t)]:
        p.unlink(missing_ok=True)


def _write_poses(layout, t, units: list[Unit], policy: str) -> None:
    """Copy out the replicates' pose records according to --poses.

    `first` (the default) writes replicate 0's poses to <label>/poses.sdf, which is what
    the pose-in-pocket validation reads. `all` adds <label>/poses_r<i>.sdf for the others;
    `none` writes nothing. Files a previous policy left are removed first, so lowering the
    policy does not leave a stale r1 behind looking current."""
    for p in _pose_files(layout, t):
        p.unlink(missing_ok=True)
    if policy == POSES_NONE:
        return
    d = _target_dir(layout, t)
    for u in sorted(units, key=lambda u: u.replicate):
        if u.replicate == 0:
            dest = d / "poses.sdf"
        elif policy == POSES_ALL:
            dest = d / f"poses_{u.name}.sdf"
        else:
            continue
        _atomic_text(dest, _read_text(unit_paths(layout, u)[0]))


def _merge(layout, t: Target, units: list[Unit], library, policy: str = DEFAULT_POSES
           ) -> list[str]:
    """Build the target's poses file(s) and its summary over EVERY replicate. Returns the
    reasons the target is incomplete ([] when it is complete); an incomplete target gets no
    files at all, and stale ones are removed.

    Every replicate must be done AND complete, and every replicate's scores enter the
    summary. A target whose r1 holds 6 of 2,083 ligands is not finished however good r0
    looks: averaging a 0.3 %-complete draw in silently is the failure this guards. A
    replicate that fails `_check_replicate` can therefore never reach `summarize`."""
    ordered = sorted(units, key=lambda u: u.replicate)
    undrawn = [u.name for u in ordered if unit_state(layout, u) != STATE_DONE]
    if undrawn:
        _clear_outputs(layout, t)
        return [f"replicate(s) not recorded: {', '.join(undrawn)}"]
    if not any(u.replicate == 0 for u in units):
        _clear_outputs(layout, t)
        return ["replicate r0 has no usable result"]
    reasons, frames = [], {}
    for u in ordered:
        frame, why = _check_replicate(layout, t, u, library)
        if why:
            reasons.append(why)
        else:
            frames[u.replicate] = frame
    if reasons or 0 not in frames:
        _clear_outputs(layout, t)
        return reasons or ["replicate r0 has no usable result"]
    _write_poses(layout, t, ordered, policy)
    df_path = _df_path(layout, t)
    fd, tmpc = tempfile.mkstemp(dir=layout.docking, prefix=f".{df_path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            summarize(frames, library, t.label).to_csv(fh, index=False)
        os.chmod(tmpc, 0o644)
        os.replace(tmpc, df_path)
    except BaseException:
        Path(tmpc).unlink(missing_ok=True)
        raise
    return []


# --- fingerprint / skip --------------------------------------------------------------

def _stage_fingerprint(layout, plan: _Plan, lib_md5: str, opts: dict,
                       policy: str = DEFAULT_POSES) -> str | None:
    """What the finished outputs depend on (not nprocs or threads): the final_targets.csv,
    each target's receptor, pocket and receptor PDBQT sha256, the library md5, the base
    seed, the replicate count, exhaustiveness, num_modes, ligand addH and the engine
    version. None when a receptor PDBQT is not prepared yet.

    This decides whether the STAGE is finished, never whether a unit must be re-drawn:
    a changed fingerprint makes the stage re-plan and re-merge, and the units whose key is
    unchanged are still skipped."""
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
    # --poses is in here and NOT in the unit key: which replicates' poses are written out
    # is not part of what a draw is (changing it must not supersede a recorded draw), but
    # it does change the stage's outputs, so a changed policy has to re-merge.
    return _sha_text(_json({"final_targets_sha256": sha256_file(layout.final_targets_csv),
                            "targets": rows, "library_md5": lib_md5, "poses": policy,
                            **opts}))


def _outputs(layout, plan: _Plan) -> dict:
    """Sizes of the current targets' result files (never of an orphan label's)."""
    d = layout.docking
    files = [p for t in plan.targets
             for p in (_df_path(layout, t), *_pose_files(layout, t)) if p.is_file()]
    return {p.relative_to(d).as_posix(): p.stat().st_size for p in sorted(files)}


# Subdirectories of 07_docking that are NOT target labels and must never be swept into
# obsolete/. The two first names look swapped and are not: "ligands" is THIS stage's
# titled-library directory (07_docking/ligands/library_titled.sdf), while
# ligands.DIRNAME == "library" is the `ligands` STAGE's own directory
# (07_docking/library/library_prepared.sdf).
#
# "export" and "discarded" are here because `_obsolete` runs on every `dock` AND every
# `dock-merge`: without them a worked-but-unmerged bundle under the default
# 07_docking/export/<ts>/, or the archive of already-discarded draws, would be relocated
# into obsolete/ behind the user's back (and `dock-merge --bundle` would then fail to find
# its own bundle). Pinned by test_orphans_never_sweep_the_export_or_discard_directories.
RESERVED_DIRS = (LIB_DIRNAME, "obsolete", ligands.DIRNAME, EXPORT_DIRNAME, DISCARD_DIRNAME)
# Both spellings: `df_results_<label>.csv` is the current summary, `_fixed` the retired
# replicate-0 table. A leftover `_fixed` file is an orphan under every label, so the sweep
# moves it to obsolete/ instead of leaving a stale file with the old schema beside the new.
_DF_RE = re.compile(r"^df_results_(.+?)(_fixed)?\.csv$")


def _orphans(layout, plan: _Plan) -> list[Path]:
    """df_results_<label>_fixed.csv files and <label>/ dirs in 07_docking whose label is
    not a current target.

    It looks only at the TOP level of 07_docking and only at whole labels, so a current
    target's replicate records (<label>/replicates/r<i>.sdf and their r<i>.done markers)
    can never be swept: they are inside the label's directory, and the shape of what a
    label contains plays no part in the decision. Adding or removing replicates therefore
    does not turn anything into an orphan. Pinned by
    test_orphans_never_sweep_replicate_results.

    An EMPTY target list sweeps NOTHING. With no labels every df_results_*.csv and every
    <label>/ matches, so a header-only final_targets.csv -- a `select` run that selected
    nothing -- would relocate the entire docking output to obsolete/<ts>/ in one `polyphar
    dock`. It is recoverable and counted, but it is a loud event that arrived silently.
    Nothing to dock is not evidence that everything already docked is obsolete."""
    d = layout.docking
    if not d.is_dir():
        return []
    labels = {t.label for t in plan.targets}
    if not labels:
        if any(_DF_RE.match(p.name) or (p.is_dir() and p.name not in RESERVED_DIRS)
               for p in d.iterdir()):
            print(f"dock: {layout.final_targets_csv} lists no targets, so the orphan sweep "
                  f"is SKIPPED: with no labels every result would count as an orphan. "
                  f"Fix the target list, or move the old results yourself", flush=True)
        return []
    out = []
    for p in sorted(d.iterdir()):
        m = _DF_RE.match(p.name)
        if p.is_file() and m and (m.group(2) or m.group(1) not in labels):
            out.append(p)
        elif p.is_dir() and p.name not in labels and p.name not in RESERVED_DIRS:
            out.append(p)
    return out


def _obsolete(layout, orphans: list[Path]) -> list[dict]:
    """Move orphans (never delete) to 07_docking/obsolete/<timestamp>/."""
    import time

    if not orphans:
        return []
    root = layout.docking / "obsolete"
    stamp = time.strftime("%Y%m%dT%H%M%S")
    dest, k = root / stamp, 1
    while dest.exists():
        k += 1
        dest = root / f"{stamp}_{k}"
    dest.mkdir(parents=True)
    moved = []
    for p in orphans:
        os.replace(p, dest / p.name)
        moved.append({"from": p.relative_to(layout.docking).as_posix(),
                      "to": (dest / p.name).relative_to(layout.docking).as_posix()})
    return moved


def _finished(layout, plan: _Plan, lib_md5: str, opts: dict,
              policy: str = DEFAULT_POSES) -> dict | None:
    import json

    try:
        manifest = json.loads((layout.docking / "manifest.json").read_text())
    except (OSError, ValueError):
        return None
    fp = _stage_fingerprint(layout, plan, lib_md5, opts, policy)
    if (fp is None or manifest.get("fingerprint") != fp
            or manifest.get("outputs") != _outputs(layout, plan) or _orphans(layout, plan)):
        return None
    if plan.outstanding(layout):   # a unit lost its marker or its record: it is not done
        return None
    return dict(manifest["counts"])


def engine_version(layout, args) -> str:
    """The engine version that goes into the unit keys.

    `dock-merge` sets args.engine_version from the bundle's run_key.json, because the host
    that merges need not have the engine at all; without that, the local binary decides,
    and when it is missing the engine recorded in an existing manifest is used so that a
    merge-only host still reports its finished units as finished."""
    pinned = getattr(args, "engine_version", None)
    if pinned:
        return str(pinned)
    local = (unidock.unidock_version() if engine_name(args) == ENGINE_UNIDOCK
             else smina_version())
    if local != "unavailable":
        return local
    return _recorded_engine(layout) or local


def _recorded_engine(layout) -> str | None:
    import json

    try:
        manifest = json.loads((layout.docking / "manifest.json").read_text())
    except (OSError, ValueError):
        return None
    engine = manifest.get("engine") or manifest.get("smina")
    return str(engine) if engine and engine != "unavailable" else None


def units(layout, args, library: Path | None = None):
    """Ligands still to draw: n_ligands for every (target, replicate) that is not done
    (0 when the stage finished on exactly these inputs). A done replicate is never counted:
    it is not re-run.

    `library`, when given, is docked AS GIVEN: `_library_path` is not consulted, so the
    `ligands` stage's prepared e-Drug3D is never substituted for it and no "docking the
    UNPREPARED library" warning is raised. It is for a caller that supplies its own ligand
    set -- `polyphar control`, which redocks co-crystal ligands -- and writes elsewhere
    (see polyphar.stages.control.ControlLayout)."""
    from polyphar.common import WorkUnits

    lib = Path(library) if library is not None else _library_path(layout, args, copy=False)
    engine = engine_version(layout, args)
    plan = _plan(layout, args, lib, engine)
    policy = str(getattr(args, "poses", None) or DEFAULT_POSES)
    if _finished(layout, plan, _md5(lib), plan.opts, policy) is not None:
        return WorkUnits(0, "ligand")
    return WorkUnits(plan.n_ligands * len(plan.outstanding(layout)), "ligand")


def titled_library(layout, lib: Path, plan: _Plan) -> Path:
    """07_docking/ligands/library_titled.sdf: the whole library, titles set to
    EDRUG3D_<ID>, written once and reused while its sha matches the plan's."""
    from polyphar.common import sha256_file

    ligdir = layout.docking / LIB_DIRNAME
    ligdir.mkdir(parents=True, exist_ok=True)
    out = ligdir / TITLED_LIBRARY
    text = "".join(_titled(r) for r in _library_records(lib))
    if not (out.is_file() and sha256_file(out) == _sha_text(text)
            and _sha_text(text) == plan.library_sha):
        _atomic_text(out, text)
    return out


def prepare_receptors(layout, plan: _Plan, rows: list) -> dict:
    """{label: (receptor pdbqt, box centre, box size)} for every usable target; failures
    are appended to `rows` as (label, reason)."""
    out: dict[str, tuple[Path, tuple, tuple]] = {}
    for t in plan.usable:
        try:
            out[t.label] = (_prepare_receptor(layout, t), *pocket_box(t.pocket))
        except Exception as exc:
            rows.append((t.label, f"receptor/pocket preparation failed: "
                                  f"{type(exc).__name__}: {exc}"))
    return out


CORRUPT_REASON = ("result file missing, empty or not terminated; "
                  "re-drawn under the same replicate index")
SUPERSEDED_REASON = ("superseded: the unit key changed (receptor, pocket, library, seed, "
                     "parameters or engine); the recorded draw is replaced")
UNCERTIFIED_REASON = ("uncertified: a result file with no marker (the run was interrupted "
                      "between writing it and recording it); it is replaced")


def pending_discard(layout, u: Unit, state: str) -> str | None:
    """Why running `u` would cost a recorded draw, or None when it would cost nothing.

    PURE: it decides, it does not discard. The discard happens where the replacement is
    written (`_dock_unit`, `import_bundle`), so planning a unit that is then filtered out,
    aborted or never executed destroys nothing. This is the fix for the bug where
    `dock-merge` on a fully docked run folder deleted every draw before establishing that
    it had nothing to put in their place."""
    sdf, _, marker = unit_paths(layout, u)
    if state == STATE_CORRUPT:
        return CORRUPT_REASON
    if result_intact(sdf):
        return SUPERSEDED_REASON if marker.exists() else UNCERTIFIED_REASON
    # A result file that exists but is not terminated, with no marker, is overwritten
    # without a discard row: it is not a record of a draw (nothing certifies it and it is
    # not even a complete file), so there is no sample to account for. Note that
    # "unterminated" and "incomplete" are different things -- a record truncated on a
    # clean $$$$ boundary IS terminated, and _check_replicate is what catches that.
    return None


def plan_discards(layout, plan: _Plan) -> list[tuple[Unit, str | None]]:
    """[(unit to run, the discard its replacement would cost or None)] for every unit that
    is not STATE_DONE. STATE_DONE units are not here at all: they are skipped, never
    redone. Discards nothing; see `pending_discard`."""
    return [(u, pending_discard(layout, u, state))
            for u, state in plan.outstanding(layout)]


def run(layout, args, progress, execute: bool = True, discard_baseline: int | None = None,
        extra_counts: dict | None = None, extra_manifest: dict | None = None,
        library: Path | None = None) -> dict:
    """Plan, run the outstanding units (unless execute is False: `dock-merge` only folds
    in results produced elsewhere) and merge.

    `discard_baseline` is the discarded.tsv row count to measure counts["discarded"]
    against; `dock-merge` takes it BEFORE it imports (importing is what discards there),
    so the merge's own discards are not under-reported as 0. `extra_counts` /
    `extra_manifest` let `dock-merge` record that a merge happened -- which bundle, how
    many units, which engine actually drew them -- in the manifest this function writes.

    `library` overrides the library to dock; see `units`."""
    from polyphar.common import append_failure, sha256_file, write_manifest
    from polyphar.pool import run_pool as _run_pool

    lib = Path(library) if library is not None else _library_path(layout, args,
                                                                  progress=progress)
    if execute:
        # fail fast and with guidance, before any planning or receptor preparation: a
        # GPU-less host cannot run `--engine unidock` at all.
        require_engine_host(args)
    engine = engine_version(layout, args)
    if execute and engine == "unavailable":
        raise FileNotFoundError(f"{engine_name(args)} ({_engine_bin(args)}) not runnable")
    lib_md5 = _md5(lib)
    plan = _plan(layout, args, lib, engine)
    opts = plan.opts
    policy = str(getattr(args, "poses", None) or DEFAULT_POSES)
    done = _finished(layout, plan, lib_md5, opts, policy)
    if done is not None:   # finished on exactly these inputs: touch nothing
        done.update(extra_counts or {})   # ... but say what THIS invocation did
        return done

    d = layout.docking
    d.mkdir(parents=True, exist_ok=True)
    previous = _read_manifest(layout)   # read before it goes: receptor drift is vs this
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")   # discarded.tsv is NOT reset
    obsoleted = _obsolete(layout, _orphans(layout, plan))   # targets dropped by `select`
    library = read_library(lib)
    rows: list[tuple[str, str]] = list(plan.missing)

    prepared = prepare_receptors(layout, plan, rows)   # receptors and boxes, in the parent
    ligands_sdf = titled_library(layout, lib, plan)

    tpj = int(getattr(args, "threads_per_job", DEFAULT_TPJ))
    timeout = int(getattr(args, "unit_timeout", None) or DEFAULT_UNIT_TIMEOUT)
    slots = _slots(args, tpj)
    n_outstanding = len(plan.outstanding(layout))
    # BEFORE anything that discards. `dock-merge` discards in import_bundle, before it
    # calls this, and passes its own earlier baseline.
    discarded = _count_discards(d) if discard_baseline is None else int(discard_baseline)
    todo = [(u, why) for u, why in plan_discards(layout, plan)
            if execute and u.label in prepared]
    engine_input = _engine_input(layout, args, plan, ligands_sdf)
    pocket_sha = {t.label: sha256_file(t.pocket) for t in plan.usable}
    tasks = []
    for u, why in todo:
        pdbqt, center, size = prepared[u.label]
        job = _unit_job(u, args, engine, pdbqt, engine_input, center, size, layout.docking,
                        timeout, tpj, why, plan.library_sha, pocket_sha.get(u.label, ""))
        _unit_dirs(layout, u)
        tasks.append(((u.label, u.name), _unit_worker, (job,)))
    unit_by_item = {(u.label, u.name): u for u in plan.units}

    failed_units: set = set()
    try:
        outcome = _run_pool(tasks, slots,
                            on_result=lambda item, _r: progress.advance(plan.n_ligands))
        for item, reason in outcome.failures:
            failed_units.add(item)
            reason = _engine_reason(reason)
            log = _rel(unit_paths(layout, unit_by_item[item])[1], layout.run_dir)
            rows.append((f"{item[0]}/{item[1]}", f"{reason}; log: {log}"))
            progress.advance(plan.n_ligands)

        complete = incomplete = 0
        for t in plan.usable:
            if t.label not in prepared:
                incomplete += 1
                _clear_outputs(layout, t)
                continue
            reasons = _merge(layout, t, plan.units_of(t.label), library, policy)
            if reasons:
                incomplete += 1
                rows.extend((t.label, r) for r in reasons)
            else:
                complete += 1
    finally:
        for item, reason in rows:   # appended from the parent only
            append_failure(d, item, STAGE, reason)

    incomplete += len(plan.missing)
    counts = {"targets": len(plan.targets), "complete": complete,
              "incomplete_targets": incomplete,
              "ligands": plan.n_ligands, "replicates": opts["replicates"],
              "units": len(plan.units),
              "units_docked": len(tasks) - len(failed_units),
              "skipped_existing": len(plan.units) - n_outstanding,
              "discarded": _count_discards(d) - discarded,
              # records of replicate indices this run does not plan (e.g. after lowering
              # --replicates): they are KEPT, never swept, and this is where that shows.
              "units_present_outside_plan": _records_outside_plan(layout, plan),
              "obsoleted": len(obsoleted), "failed": len(failed_units)}
    counts.update(extra_counts or {})
    receptors = {t.label: {"receptor": t.receptor_rel, "pocket": t.pocket_rel,
                           "preparation": t.prep,
                           "pdbqt_sha256": sha256_file(prepared[t.label][0])}
                 for t in plan.usable if t.label in prepared}
    _warn_receptor_drift(progress, previous, receptors)
    fingerprint = _stage_fingerprint(layout, plan, lib_md5, opts, policy)
    write_manifest(d, STAGE,
                   args={"nprocs": args.nprocs, "seed": opts["base_seed"],
                         "replicates": opts["replicates"],
                         "threads_per_job": tpj, "ligand_addH": LIGAND_ADDH,
                         "engine": opts["engine_name"], "poses": policy,
                         **engine_settings(args),
                         "unit_timeout": timeout, "library": _rel(lib, layout.run_dir)},
                   inputs=[layout.final_targets_csv, lib],
                   versions=_versions(engine, opts["engine_name"]),
                   counts=counts, timing={},
                   extra={"fingerprint": fingerprint, "outputs": _outputs(layout, plan),
                          "library_md5": lib_md5, "library_sha256": plan.library_sha,
                          # True: the scores come from the `ligands` stage's prepared
                          # library; False: from the raw SDF (Methods needs to say which).
                          "library_prepared": _library_is_prepared(layout, args, lib),
                          "base_seed": opts["base_seed"], "engine": engine, "smina": engine,
                          # per unit: which draw it is, what seed it recorded, its key,
                          # and the provenance the host that DREW it recorded beside it
                          # (engine build, GPU model, CUDA runtime, settings, digests)
                          "units": [{"label": u.label, "replicate": u.replicate,
                                     "seed": u.seed, "key": u.key,
                                     "state": unit_state(layout, u),
                                     "provenance": read_provenance(layout, u)}
                                    for u in plan.units],
                          "receptors": receptors, "obsoleted": obsoleted,
                          **(extra_manifest or {})})
    return counts


def _read_manifest(layout) -> dict:
    import json

    try:
        return json.loads((Path(layout.docking) / "manifest.json").read_text())
    except (OSError, ValueError):
        return {}


_REPLICATE_RE = re.compile(r"^r(\d+)\.sdf$")


def _records_outside_plan(layout, plan: _Plan) -> int:
    """How many <label>/replicates/r<i>.sdf files belong to a replicate index this run
    does not plan. Lowering --replicates keeps r3, r4 ... (a recorded draw is never thrown
    away for being surplus); without this count nothing said so."""
    want = {(u.label, u.replicate) for u in plan.units}
    labels = {t.label for t in plan.targets}
    n = 0
    for label in labels:
        for f in sorted((Path(layout.docking) / label / "replicates").glob("r*.sdf")):
            m = _REPLICATE_RE.match(f.name)
            if m and (label, int(m.group(1))) not in want:
                n += 1
    return n


def _warn_receptor_drift(progress, previous: dict, receptors: dict) -> None:
    """Warn when a target's receptor PDBQT is not byte-identical to the one the previous
    manifest recorded.

    The PDBQT is what the engine actually reads, and it is deliberately NOT in the unit key
    (an Open Babel upgrade must not invalidate every recorded draw). The cost of leaving it
    out is that draws under one key could come from two different receptors, so the change
    is reported rather than silently absorbed."""
    import sys

    old = (previous.get("extra") or previous).get("receptors") or {}
    drift = sorted(label for label, info in receptors.items()
                   if label in old and old[label].get("pdbqt_sha256")
                   and old[label]["pdbqt_sha256"] != info["pdbqt_sha256"])
    if not drift:
        return
    message = (f"the prepared receptor PDBQT changed for {len(drift)} target(s) since the "
               f"last run ({', '.join(drift[:5])}"
               f"{', ...' if len(drift) > 5 else ''}). The PDBQT is what the engine reads "
               f"but is not part of a unit key, so draws recorded before and after this "
               f"change came from DIFFERENT receptors under the same keys. Check the "
               f"Open Babel version and the receptor inputs before using them together.")
    print(f"{STAGE}: WARNING: {message}", file=sys.stderr, flush=True)
    notifier = getattr(progress, "notifier", None)
    if notifier is not None:
        try:
            notifier(f"{STAGE}: receptor PDBQT changed", message, priority=1)
        except Exception:
            pass


def _count_discards(stage_dir) -> int:
    """Rows in discarded.tsv (0 without one). The difference across a run is how many
    recorded draws that run threw away."""
    path = Path(stage_dir) / DISCARD_LOG
    try:
        return max(0, sum(1 for _ in open(path, encoding="utf-8")) - 1)
    except OSError:
        return 0


def _versions(engine: str, name: str = ENGINE_SMINA) -> dict:
    import importlib.metadata

    out = {"engine": engine, "engine_name": name,
           # kept under its old key so a reader that knows only `smina` still finds the
           # engine string; it is the ENGINE's version, whichever engine that is.
           "smina": engine}
    if name == ENGINE_UNIDOCK:
        out.update({"unidock_tools": unidock.tools_version(),
                    "cuda_driver_max": unidock.cuda_driver_max(),
                    "gpu_model": unidock.gpu_model() or "none"})
    for pkg in ("openbabel", "pdbfixer", "openmm"):
        try:
            out[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            out[pkg] = "unknown"
    return out


# --- bundles: export / worker / merge ---------------------------------------------------
#
# `dock` needs a GPU (soon) and the run folder on one machine. When they are on different
# machines the same work is split in three, around a bundle directory:
#
#   dock-export   plan here, write the outstanding units and everything they need
#   dock-worker   run them there; no run folder, no plan, no key arithmetic
#   dock-merge    fold the results back in here
#
# The bundle's on-disk layout mirrors the run folder's (results/<label>/replicates/r<i>.sdf
# with its r<i>.done marker beside it), so the same unit_state / discard / resume code runs
# on both sides and merging is a copy. Packing it into a tarball with a verified run_key.json
# is a later step; what it needs from this one is that the unit KEY is computed once, at
# export, and carried through: the worker never recomputes it and the merge accepts only
# results whose key is the one the current plan expects.

BUNDLE_UNITS = "units.tsv"
BUNDLE_KEY = "run_key.json"
BUNDLE_LIBRARY = "library.sdf"
BUNDLE_UNIT_COLUMNS = ["label", "replicate", "seed", "key", "receptor",
                       "center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]


@dataclass(frozen=True)
class BundleUnit:
    unit: Unit
    receptor: str          # relative to the bundle
    center: tuple
    size: tuple


class BundleLayout:
    """Enough of a RunLayout for the unit helpers: a bundle's results/ directory stands in
    for 07_docking, so unit_paths / unit_state / discard work unchanged on both sides."""

    def __init__(self, bundle):
        self.bundle = Path(bundle)
        self.run_dir = self.bundle

    @property
    def docking(self) -> Path:
        return self.bundle / "results"


def write_bundle(layout, args, dest, plan: _Plan, prepared: dict, ligands_sdf: Path,
                 max_units: int | None = None) -> list[BundleUnit]:
    """Write the outstanding units and their inputs into `dest`. Returns the units written
    (at most max_units, in replicate-major order so a lost session costs one draw of many
    targets rather than every draw of a few).

    Exporting touches nothing in the run folder: a record that is going to be replaced is
    discarded by whoever actually replaces it (`dock-merge`, or `dock` itself), not here,
    so a session that is never worked costs no draws at all."""
    import shutil

    dest = Path(dest)
    (dest / "receptors").mkdir(parents=True, exist_ok=True)
    (dest / "results").mkdir(parents=True, exist_ok=True)
    todo = [u for u, _state in plan.outstanding(layout) if u.label in prepared]
    todo.sort(key=lambda u: (u.replicate, u.label))
    if max_units is not None:
        todo = todo[:max_units]
    out = []
    for u in todo:
        pdbqt, center, size = prepared[u.label]
        rel = f"receptors/{u.label}.pdbqt"
        if not (dest / rel).is_file():
            shutil.copyfile(pdbqt, dest / rel)
        out.append(BundleUnit(u, rel, center, size))
    shutil.copyfile(ligands_sdf, dest / BUNDLE_LIBRARY)
    lines = ["\t".join(BUNDLE_UNIT_COLUMNS)]
    for b in out:
        u = b.unit
        lines.append("\t".join([u.label, str(u.replicate), str(u.seed), u.key, b.receptor]
                               + [repr(float(c)) for c in b.center]
                               + [repr(float(s)) for s in b.size]))
    _atomic_text(dest / BUNDLE_UNITS, "\n".join(lines) + "\n")
    from polyphar.common import sha256_file

    _atomic_text(dest / BUNDLE_KEY, _json_indent({
        "run": Path(layout.run_dir).name, "library_sha256": plan.library_sha,
        # the panel this bundle was planned against: a bundle exported before `select` was
        # rerun describes targets that may no longer exist, and the per-unit key check
        # would say so one opaque rejection at a time
        "panel_sha256": sha256_file(layout.final_targets_csv),
        "exported_from": host_path(),
        "n_ligands": plan.n_ligands, "options": plan.opts,
        "unit_timeout": int(getattr(args, "unit_timeout", None) or DEFAULT_UNIT_TIMEOUT),
        "units": len(out), "outstanding_total": len(plan.outstanding(layout))}))
    return out


def _json_indent(obj) -> str:
    import json

    return json.dumps(obj, indent=2, sort_keys=True) + "\n"


def read_run_key(bundle) -> dict:
    import json

    path = Path(bundle) / BUNDLE_KEY
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found: {bundle} is not a dock bundle")
    return json.loads(path.read_text())


def read_bundle_units(bundle) -> list[BundleUnit]:
    """The bundle's units.tsv. The key is read, never recomputed: it was decided at export
    and is what the merge will check against the plan."""
    path = Path(bundle) / BUNDLE_UNITS
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found: {bundle} is not a dock bundle")
    rows = [l.split("\t") for l in path.read_text().splitlines() if l.strip()]
    if not rows or rows[0] != BUNDLE_UNIT_COLUMNS:
        raise ValueError(f"{path}: expected columns {', '.join(BUNDLE_UNIT_COLUMNS)}")
    out = []
    for r in rows[1:]:
        if len(r) != len(BUNDLE_UNIT_COLUMNS):
            raise ValueError(f"{path}: malformed row {r!r}")
        out.append(BundleUnit(Unit(r[0], int(r[1]), int(r[2]), r[3]), r[4],
                              tuple(float(x) for x in r[5:8]),
                              tuple(float(x) for x in r[8:11])))
    return out


def bundle_outstanding(bundle) -> list[tuple[BundleUnit, str | None]]:
    """[(unit still to draw, the discard re-drawing it would cost or None)]. Done units are
    skipped, exactly as in the run folder.

    The reason is only PLANNED here; it is acted on in `_dock_unit`, where the replacement
    is written. Both the corrupt case AND the uncertified one (an intact result whose
    marker was dropped on the way back from the engine host) are accounted for: before,
    only the first was, so a worker handed results without their markers re-docked over
    them with nothing in discarded.tsv."""
    bl = BundleLayout(bundle)
    out = []
    for b in read_bundle_units(bundle):
        state = unit_state(bl, b.unit)
        if state == STATE_DONE:
            continue
        out.append((b, pending_discard(bl, b.unit, state)))
    return out


RUN_KEY_SETTINGS = ("engine_name", "search_mode", "batch_size", "base_seed", "ligand_addH",
                    "exhaustiveness", "num_modes")


def verify_run_key(bundle, args):
    """Check the bundle's run_key.json against this host's settings BEFORE any work, and
    return the arguments the units must actually be run with.

    The bundle's settings win, because the unit keys were computed from them at export: a
    worker that docked with its own engine or search mode would produce draws that
    `dock-merge` then rejects on every key, after the GPU time was spent. A setting the
    CALLER gave explicitly and that contradicts the bundle is refused outright rather than
    silently overridden -- that is a user who believes they are running something else."""
    import copy

    opts = dict(read_run_key(bundle).get("options") or {})
    out = copy.copy(args)
    conflicts = []
    for field, attr in (("engine_name", "engine"), ("search_mode", "search_mode"),
                        ("batch_size", "batch_size")):
        if field not in opts:
            continue
        given = getattr(args, attr, None)
        if given is not None and str(given) != str(opts[field]):
            conflicts.append(f"--{attr.replace('_', '-')}: bundle {opts[field]!r}, "
                             f"this host {given!r}")
        setattr(out, attr, opts[field])
    if conflicts:
        raise ValueError(
            f"{bundle}: run_key.json disagrees with this host's settings, so its units "
            f"would be drawn under an identity the merge will reject -- "
            f"{'; '.join(conflicts)}. Drop the option, or export a new bundle with it")
    return out


def work_bundle(bundle, args, progress) -> dict:
    """Run the bundle's outstanding units in place. No run folder is needed and no key is
    computed: the keys come from units.tsv."""
    from polyphar.pool import run_pool as _run_pool

    bundle = Path(bundle)
    key = read_run_key(bundle)
    bundle_args = verify_run_key(bundle, args)
    bl = BundleLayout(bundle)
    todo = bundle_outstanding(bundle)
    name = engine_name(bundle_args)
    if todo:
        require_engine_host(bundle_args)
    engine = (unidock.unidock_version() if name == ENGINE_UNIDOCK else smina_version())
    if todo and engine == "unavailable":
        raise FileNotFoundError(f"{name} ({_engine_bin(bundle_args)}) not runnable")
    n_lig = int(key.get("n_ligands") or 0)
    tpj = int(getattr(args, "threads_per_job", DEFAULT_TPJ))
    timeout = int(getattr(args, "unit_timeout", None) or key.get("unit_timeout")
                  or DEFAULT_UNIT_TIMEOUT)
    engine_input = (bundle / BUNDLE_LIBRARY if name != ENGINE_UNIDOCK else
                    unidock.split_library(bundle / BUNDLE_LIBRARY, bundle / "split",
                                          str(key.get("library_sha256") or "")))
    tasks = []
    for b, why in todo:
        _unit_dirs(bl, b.unit)
        job = _unit_job(b.unit, bundle_args, engine, bundle / b.receptor, engine_input,
                        b.center, b.size, bl.docking, timeout, tpj, why,
                        str(key.get("library_sha256") or ""), "")
        tasks.append(((b.unit.label, b.unit.name), _unit_worker, (job,)))
    outcome = _run_pool(tasks, _slots(bundle_args, tpj),
                        on_result=lambda item, _r: progress.advance(n_lig))
    for item, reason in outcome.failures:
        progress.advance(n_lig)
        append_bundle_failure(bundle, f"{item[0]}/{item[1]}",
                              _engine_reason(reason))
    _atomic_text(bundle / "worker.json", _json_indent(
        {"engine": engine, "engine_name": name, "units": len(read_bundle_units(bundle)),
         "units_docked": len(tasks) - len(outcome.failures),
         "failed": len(outcome.failures), "versions": _versions(engine, name),
         "host": host_path(),
         "gpu_models": sorted({(read_provenance(bl, b.unit).get("gpu_model") or "")
                               for b in read_bundle_units(bundle)} - {"", None})}))
    return {"units": len(read_bundle_units(bundle)),
            "units_docked": len(tasks) - len(outcome.failures),
            "skipped_existing": len(read_bundle_units(bundle)) - len(todo),
            "failed": len(outcome.failures)}


def append_bundle_failure(bundle, item, reason) -> None:
    from polyphar.common import append_failure

    append_failure(Path(bundle), item, "dock-worker", reason)


def import_bundle(layout, plan: _Plan, bundle):
    """Copy the bundle's finished unit records into the run folder.
    Returns (taken, rejected, dropped).

    A result is taken only when the current plan has a unit with the SAME label, replicate
    AND key: a bundle made under other settings (another receptor, library, seed or engine)
    describes different draws, and mixing them silently is the thing this check exists to
    stop. Rejections are returned, not swallowed, and a bundle draw dropped because this
    unit is already on record is written to discarded.tsv -- the selection has to be
    auditable even when nothing failed."""
    import shutil

    bl = BundleLayout(bundle)
    want = {(u.label, u.replicate): u for u in plan.units}
    taken, rejected, dropped = [], [], []
    logged = _discard_log_text(layout.docking)   # re-merging a bundle must not re-log
    for b in read_bundle_units(bundle):
        u = b.unit
        mine = want.get((u.label, u.replicate))
        if mine is None:
            rejected.append((f"{u.label}/{u.name}", "not a unit of this run"))
            continue
        if mine.key != u.key:
            rejected.append((f"{u.label}/{u.name}",
                             "unit key differs from this run's: the bundle was made under "
                             "different settings and its draws are not these draws"))
            continue
        if unit_state(bl, u) != STATE_DONE:
            continue   # not drawn yet in the bundle; a later worker pass will
        # L3: the marker and the `$$$$` terminator say nothing about the CONTENT, so a
        # result edited to different scores and re-terminated used to be copied in
        # unchecked. The worker recorded the result's sha256 in its provenance; a record
        # that no longer matches it is not the draw that was made.
        tampered = verify_result_digest(bl, u)
        if tampered:
            rejected.append((f"{u.label}/{u.name}", tampered))
            continue
        src_sdf, src_log, src_done = unit_paths(bl, u)
        dst_sdf, dst_log, dst_done = unit_paths(layout, mine)
        if unit_state(layout, mine) == STATE_DONE:
            # Already recorded here; a recorded draw is never replaced. But the bundle's
            # draw was a real sample too, and it is now being dropped -- which happens
            # whenever two bundles were exported before either was merged, since export
            # only excludes what is already ON RECORD. Dropping it is a (small) selection
            # effect, so it is logged like any other discarded draw. Its bytes stay in the
            # bundle: `archived` says where.
            dropped.append((f"{u.label}/{u.name}", str(bundle)))
            if f"{u.key}\t" not in logged or str(src_sdf) not in logged:
                log_discard(layout.docking, u,
                            f"dropped on merge: a draw for this unit is already on record "
                            f"here, so the bundle's draw was not merged (bundle: {bundle})",
                            str(src_sdf))
            continue
        if dst_sdf.exists() or dst_done.exists():
            # Something is standing in this unit's place and is about to be replaced by
            # the bundle's copy: archive and log it here, at the moment of replacement.
            discard(layout, mine, "replaced by a draw merged from a bundle")
        dst_sdf.parent.mkdir(parents=True, exist_ok=True)
        dst_done.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src_sdf, dst_sdf)
        if src_log.is_file():
            shutil.copyfile(src_log, dst_log)
        src_prov = prov_path(bl, u)
        if src_prov.is_file():   # the provenance of the draw travels with the draw
            shutil.copyfile(src_prov, prov_path(layout, mine))
        _atomic_text(dst_done, mine.key + "\n")   # written last: the record is in place
        taken.append(mine)
    _carry_discards(bundle, layout)
    return taken, rejected, dropped


def _discard_log_text(stage_dir) -> str:
    try:
        return (Path(stage_dir) / DISCARD_LOG).read_text(encoding="utf-8")
    except OSError:
        return ""


def _carry_discards(bundle, layout) -> None:
    """Append the bundle's discarded draws to the run folder's log: a draw thrown away on
    the GPU host is still a draw thrown away.

    The bundle's own log is LEFT IN PLACE (it used to be deleted, so a bundle merged into a
    second run folder arrived with no discard history at all), and repeat merges are
    idempotent on CONTENT: a row already present here is not appended again."""
    src = BundleLayout(bundle).docking / DISCARD_LOG
    if not src.is_file():
        return
    rows = [l for l in src.read_text().splitlines()[1:] if l.strip()]
    if not rows:
        return
    dst = Path(layout.docking) / DISCARD_LOG
    dst.parent.mkdir(parents=True, exist_ok=True)
    have = set()
    if dst.is_file():
        have = {l for l in dst.read_text().splitlines()[1:] if l.strip()}
    rows = [r for r in rows if r not in have]
    if not rows:
        return
    _ensure_discard_log(layout.docking)
    with open(dst, "a", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")


def plan_for(layout, args, engine: str | None = None):
    """(plan, library path, engine version) for the bundle subcommands."""
    lib = _library_path(layout, args, copy=False)
    engine = engine or engine_version(layout, args)
    return _plan(layout, args, lib, engine), lib, engine
