"""Uni-Dock engine adapter for the dock stage (GPU; `unidock_tools` pipeline).

`smina` stays the single-host / development engine. This module is the other one, selected
with `--engine unidock`. What it has to get right, and what the rest of the stage must not
assume, is collected here:

* **One SDF per input ligand.** `unidock_tools` does not write one multi-record output; it
  writes `<savedir>/<ligand stem>.sdf` per input ligand, preserving the input SD fields.
  The engine input therefore has to be the library SPLIT into one file per ligand plus an
  index file, not the single titled SDF smina takes. `split_library` does that once per
  library and caches it.

* **Key on `<ID>`, never on the title.** The ligand goes through RDKit inside
  `unidock_tools`, which MAY rewrite the title line -- upstream `_postprocessing` does not
  touch `_Name` itself, and an `SDMolSupplier`/`SDWriter` round-trip preserves it, but the
  title is not something we control and nothing upstream promises to keep it. The input SD
  fields are explicitly preserved (`_postprocessing` clears only five of its own
  properties), so `<ID>` is the identifier. A previous investigation lost time to exactly
  this. The stub DOES rewrite the title, which pins the decision: the pipeline must work
  whether or not the title survives.

* **The score is `<docking_score>`, not `ENERGY=`.** `unidock_tools` adds its own
  `<docking_score>` field; the raw engine's `ENERGY=` remark is not what it writes. Looking
  for `ENERGY=` made working runs look as if they had scored nothing.

* **`--search_mode balance`** (exhaustiveness 384). smina's exhaustiveness 12 is a CPU
  figure and is deliberately NOT mirrored.

* **Uni-Dock < 1.2.0 is refused, and the gate guards the binary that docks.** Earlier
  versions silently corrupted about 20 % of results for boxes >= 30 A, which is most of
  this study's pockets. A wrong number that looks fine is worse than no number, so this is
  a refusal, not a warning. `unidock_tools` spawns a hardcoded `cmd = ["unidock"]` and runs
  no version check of its own, so the gate resolves `unidock` on $PATH and refuses a
  `$UNIDOCK_BIN` that is a different file -- otherwise it would certify one build while
  another ran. The banner must also say `Uni-Dock` (see BANNER) and name a released
  version, not a pre-release.

* **The library must not already carry a score field.** `_postprocessing` clears only five
  of its own properties and passes every other input SD field through, so a
  `<minimizedAffinity>` in the library would come back beside the real `<docking_score>`
  (see `_check_no_stale_score`).

* **Absurd scores are filtered and logged** (see ABSURD_SCORE).

The adapter normalises its output into the same single multi-record SDF the rest of the
stage reads, with each record carrying its input `<ID>` and a `<docking_score>` field, so
`parse_poses`, `result_intact` and the whole resume/merge machinery are engine-independent.

Nothing in this module has been executed against a real Uni-Dock: it is written against
the documented `unidock_tools` CLI and exercised through a stub (tests/fixtures/
fake_unidock). The argv is built in one place (`pipeline_cmd`) so a mismatch with the real
tool is a one-line fix.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

DEFAULT_SEARCH_MODE = "balance"      # = exhaustiveness 384
# Ligands per GPU batch. NOT the unidock_tools default, which is 18000 (its argparse
# default); 1200 is the figure its own README uses in the ligandprep example. Kept because
# this study's unit is 1,960 ligands: 18000 would make every unit a single batch whose peak
# GPU memory is set by the widest ligand in the whole library, and a batch that does not fit
# fails the whole unit rather than one batch of it.
DEFAULT_BATCH_SIZE = 1200
SEARCH_MODES = ("fast", "balance", "detail")
MIN_VERSION = (1, 2, 0)
SCORE_FIELD = "docking_score"
ID_FIELD = "ID"
INDEX_NAME = "ligand_index.txt"

# Scores beyond this (kcal/mol) are not observations, they are the known upstream failure
# mode: Uni-Dock issue #205 reports ~+300 on a Meeko-prepared ligand -- against 1.2.0
# ITSELF, not a version the 1.2.0 gate above excludes, and the figure carries no units
# upstream. That makes this filter necessary, not redundant with the version gate. 1.2.0
# also fixed out-of-bounds poses being reported with FLT_MAX (3.4e38) affinities, which
# this catches too.
#
# Why +-100 is unreachable by a genuine pose, correctly: a Vina-family score for a drug-like
# ligand in a real pocket lies between about -15 and +10. The repulsion term is NOT clamped
# -- it is d^2 for d < 0 -- but it is bounded by geometry: |d| <= R_i + R_j ~ 3.6 A, so one
# interpenetrating pair contributes ~11 kcal at most within the 8 A cutoff. Summed over
# hundreds of pairs that is NOT bounded "in the tens", so that is not the argument. What
# actually puts +100 out of reach is that Vina-family engines report MINIMISED poses: the
# optimiser drives a clashing pose apart before it is scored, and the result is divided by
# 1/(1 + 0.0585 * N_rot). A pose that still scores +100 after minimisation is not a pose.
# Non-finite values go the same way. Dropped ligands are NOT silently discarded: they are
# written to the unit's <unit>.filtered.tsv and named in its log, and the ligand then has no
# poses -- which makes the replicate incomplete, exactly as a ligand the engine skipped does.
ABSURD_SCORE = 100.0


class UniDockError(RuntimeError):
    pass


def unidock_bin() -> str:
    """The Uni-Dock binary that will ACTUALLY dock.

    `unidock_tools` hardcodes `cmd = ["unidock"]` (modules/docking/unidock.py): no
    shutil.which of its own, no env var, and no version check. So the engine that runs is
    whatever `unidock` resolves to on $PATH, and a gate on $UNIDOCK_BIN alone would be
    guarding a different file. $UNIDOCK_BIN is still honoured -- it is how you point at a
    specific build -- but `check_version` refuses when the two disagree instead of
    certifying one and running the other."""
    return os.environ.get("UNIDOCK_BIN") or path_unidock() or "unidock"


def path_unidock() -> str | None:
    """What `unidock` resolves to on $PATH: the binary `unidock_tools` will spawn."""
    import shutil

    return shutil.which("unidock")


def tools_bin() -> str:
    return os.environ.get("UNIDOCK_TOOLS_BIN") or "unidocktools"


# --- versions and host ------------------------------------------------------------------


# `suffix` captures a PRE-RELEASE marker only ('-rc1', 'rc1', '-beta.2'); a fourth numeric
# component ('1.2.0.1') and build metadata ('1.2.0+cuda12') are not pre-releases and are
# left alone.
_VERSION_RE = re.compile(
    r"(\d+)\.(\d+)\.(\d+)(?P<suffix>-[A-Za-z0-9][A-Za-z0-9.]*|[A-Za-z][A-Za-z0-9.]*)?")
# The banner must say what engine it is. A bare `\d+\.\d+\.\d+` match accepted
# 'AutoDock Vina v1.2.5' -- a different engine entirely, whose 1.2.5 has nothing to do
# with Uni-Dock's 1.2.0 grid fix.
BANNER = "uni-dock"


def _first_version(text: str) -> tuple[int, ...] | None:
    m = _VERSION_RE.search(text or "")
    return tuple(int(g) for g in m.groups()[:3]) if m else None


def _version_suffix(text: str) -> str:
    m = _VERSION_RE.search(text or "")
    return (m.group("suffix") or "") if m else ""


def _same_binary(a: str, b: str) -> bool:
    try:
        return os.path.samefile(a, b)
    except OSError:
        return os.path.realpath(a) == os.path.realpath(b)


def check_binary() -> str:
    """The binary that will actually dock, refusing when $UNIDOCK_BIN names another one.

    `unidock_tools` spawns a hardcoded `cmd = ["unidock"]`, so $PATH decides the engine and
    it runs no version check of its own. Gating $UNIDOCK_BIN while `unidock` on $PATH is a
    1.1.x build would certify one binary and run the other -- which is the exact silent
    corruption (~20 % of results for boxes >= 30 A) the gate exists to prevent."""
    on_path = path_unidock()
    override = os.environ.get("UNIDOCK_BIN")
    if not override:
        if on_path is None:
            raise UniDockError(
                "no `unidock` on $PATH. unidock_tools spawns a hardcoded `unidock`, so that "
                "is the binary that docks: install Uni-Dock >= 1.2.0 on $PATH (setting "
                "$UNIDOCK_BIN alone is not enough -- unidock_tools never reads it)")
        return on_path
    if not Path(override).is_file():
        raise UniDockError(
            f"$UNIDOCK_BIN={override} is not runnable (no such file). Set it to the "
            f"Uni-Dock >= {'.'.join(map(str, MIN_VERSION))} binary, which must also be the "
            f"`unidock` on $PATH -- unidock_tools spawns that one and never reads "
            f"$UNIDOCK_BIN")
    if on_path is None:
        raise UniDockError(
            f"$UNIDOCK_BIN={override} is set but there is no `unidock` on $PATH. "
            f"unidock_tools spawns a hardcoded `unidock` and never reads $UNIDOCK_BIN, so "
            f"nothing would run. Put that binary on $PATH")
    if not _same_binary(override, on_path):
        raise UniDockError(
            f"$UNIDOCK_BIN={override} is NOT the binary that will dock: unidock_tools "
            f"spawns a hardcoded `unidock`, which resolves to {on_path}. The version gate "
            f"would check one build and the run would use the other. Point $UNIDOCK_BIN at "
            f"{on_path}, or put {override} first on $PATH")
    return on_path


def unidock_version() -> str:
    """First line of `unidock --version`, or "unavailable"."""
    try:
        p = subprocess.run([unidock_bin(), "--version"], capture_output=True, text=True,
                           timeout=60)
    except (OSError, subprocess.SubprocessError):
        return "unavailable"
    lines = (p.stdout or p.stderr).strip().splitlines()
    return lines[0].strip() if p.returncode == 0 and lines else "unavailable"


def tools_version() -> str:
    import importlib.metadata

    try:
        return importlib.metadata.version("unidock_tools")
    except Exception:
        return "unknown"


def check_version(version: str | None = None) -> str:
    """Refuse Uni-Dock older than 1.2.0.

    Not a warning: versions before 1.2.0 silently corrupted about 20 % of results for boxes
    of 30 A and larger, which is most of this study's pockets. The failure is invisible in
    the output, so running anyway would produce numbers that cannot be told from good
    ones."""
    if version is None:
        check_binary()          # gate the binary that docks, not a different $UNIDOCK_BIN
        version = unidock_version()
    if version == "unavailable":
        raise UniDockError(
            f"Uni-Dock ({unidock_bin()}) not runnable (`{unidock_bin()} --version` failed). "
            f"Install Uni-Dock >= 1.2.0 and unidock_tools, or set $UNIDOCK_BIN")
    if BANNER not in version.lower():
        raise UniDockError(
            f"`{unidock_bin()} --version` says {version!r}, which does not identify itself "
            f"as Uni-Dock. A bare version number proves nothing -- 'AutoDock Vina v1.2.5' "
            f"would otherwise pass this gate as a 1.2.5. Point $UNIDOCK_BIN / $PATH at "
            f"Uni-Dock >= {'.'.join(map(str, MIN_VERSION))}")
    got = _first_version(version)
    if got is None:
        raise UniDockError(
            f"cannot read a version out of `{unidock_bin()} --version` ({version!r}); "
            f"Uni-Dock >= 1.2.0 is required and this cannot be confirmed to be one")
    suffix = _version_suffix(version)
    if suffix:
        raise UniDockError(
            f"Uni-Dock {'.'.join(map(str, got))}{suffix} is a pre-release build; a "
            f"{'.'.join(map(str, MIN_VERSION))}-rc may predate the grid fix that "
            f"{'.'.join(map(str, MIN_VERSION))} carries, and this stage will not take that "
            f"on trust. Install a released Uni-Dock >= {'.'.join(map(str, MIN_VERSION))}")
    if got < MIN_VERSION:
        raise UniDockError(
            f"Uni-Dock {'.'.join(map(str, got))} is too old: versions before "
            f"{'.'.join(map(str, MIN_VERSION))} silently corrupted ~20 % of results for "
            f"boxes >= 30 A, which is most of this study's pockets. Upgrade; this stage "
            f"will not run with it")
    return version


def gpu_model() -> str | None:
    """The first CUDA device's model name, or None when this host has no usable GPU.

    $POLYPHAR_GPU_MODEL overrides it (a container that cannot see nvidia-smi but has the
    device, and the tests)."""
    override = os.environ.get("POLYPHAR_GPU_MODEL")
    if override:
        return override.strip() or None
    try:
        p = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                           capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    if p.returncode != 0:
        return None
    names = [l.strip() for l in (p.stdout or "").splitlines() if l.strip()]
    return names[0] if names else None


def cuda_driver_max() -> str:
    """The HIGHEST CUDA version this driver supports, or "unknown".

    NOT the installed CUDA runtime: plain `nvidia-smi` prints the driver's maximum
    supported version ("CUDA Version: 12.4" next to the driver number), which a machine
    with CUDA 11.8 installed will also show. It goes into provenance as a property of the
    host, so it is named for what it is. $POLYPHAR_CUDA_DRIVER_MAX wins (and the older
    $POLYPHAR_CUDA_RUNTIME is still honoured)."""
    override = (os.environ.get("POLYPHAR_CUDA_DRIVER_MAX")
                or os.environ.get("POLYPHAR_CUDA_RUNTIME"))
    if override:
        return override.strip()
    try:
        p = subprocess.run(["nvidia-smi"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    m = re.search(r"CUDA Version:\s*([0-9.]+)", p.stdout or "")
    return m.group(1) if m else "unknown"


NO_GPU_MESSAGE = (
    "no CUDA GPU found (nvidia-smi reported none). `--engine unidock` needs one: Uni-Dock "
    "is GPU-only. Run the dock stage where the GPU is -- `polyphar dock-export` here, then "
    "`polyphar dock-worker --bundle <dir>` on the GPU host (Colab / TRUBA / a rented "
    "container) and `polyphar dock-merge --bundle <dir>` back here -- or use "
    "`--engine smina` for a CPU run. Set $POLYPHAR_GPU_MODEL if the device is present but "
    "nvidia-smi is not.")


def require_gpu() -> str:
    """The GPU model, or UniDockError with the guidance above. Fail fast: POSEIDON, where
    the run folders live, has no GPU, so this is the common case and not an edge one."""
    model = gpu_model()
    if model is None:
        raise UniDockError(NO_GPU_MESSAGE)
    return model


# --- library split ------------------------------------------------------------------------

SAFE_ID = re.compile(r"[A-Za-z0-9._-]+")
# `_postprocessing` names its output after the input file: `ligand.stem.rstrip('_out')`.
# rstrip takes a character SET, not a suffix, so every trailing o/u/t/_ is eaten: 'ABCt_out'
# and 'ABC_out' both become 'ABC', and an ID made only of those characters ('out', 'too')
# leaves an EMPTY stem. We choose the input filenames, so we refuse the IDs that collide
# instead of discovering it as a silently missing ligand. (Harmless for e-Drug3D's integer
# IDs; this is for the next library.)
MANGLED_BY_RSTRIP = re.compile(r"[_out]+$")
STALE_SCORE_FIELD = "minimizedAffinity"


def _check_id(lid: str, source) -> None:
    """Refuse a library ID that is not a safe single filename component.

    The ID becomes `dest/<ID>.sdf`. An ID of '../ligands' escaped `dest` entirely and
    overwrote 07_docking/ligands/ligands.sdf -- the prepared engine library -- and the
    escaped file then survived the next `shutil.rmtree(dest)`."""
    if not SAFE_ID.fullmatch(lid) or lid in (".", ".."):
        raise UniDockError(
            f"{source}: library ID {lid!r} is not usable as a filename (allowed: letters, "
            f"digits, '.', '_', '-'). Each ligand becomes one <ID>.sdf file for "
            f"unidock_tools, and an ID with a path separator writes outside the split "
            f"directory")
    if not MANGLED_BY_RSTRIP.sub("", lid):
        raise UniDockError(
            f"{source}: library ID {lid!r} is made only of the characters unidock_tools' "
            f"`stem.rstrip('_out')` strips, so its output file would have an empty stem "
            f"and the ligand would vanish. Rename it")


def _check_no_stale_score(fields: dict, lid: str, source) -> None:
    """Refuse a library record that already carries a score field.

    `unidock_tools._postprocessing` clears only its own five properties and passes every
    other INPUT SD field through, so a <minimizedAffinity> in the library would come back
    out beside the real <docking_score>. Caught here, before any GPU time, rather than as a
    wrong number in the results."""
    if STALE_SCORE_FIELD in fields:
        raise UniDockError(
            f"{source}: library ligand {lid} already carries a <{STALE_SCORE_FIELD}> field "
            f"({fields[STALE_SCORE_FIELD]!r}). unidock_tools passes input SD fields through, "
            f"so that stale score would come back beside the real <{SCORE_FIELD}> and could "
            f"be read as this run's result. Strip it from the library first")


def split_library(titled_sdf, dest, library_sha: str) -> Path:
    """One SDF per library record under `dest`, plus INDEX_NAME listing them in order.

    `unidock_tools` writes one output per INPUT FILE, so the engine input cannot be the
    single multi-record library. The split is cached: `dest/.library_sha` records which
    library it holds, and an unchanged library is not split again. An EMPTY `library_sha`
    never matches a stamp: it identifies nothing, so a caller that could not supply one
    (a bundle written without the field) must re-split rather than trust what is there."""
    from polyphar.stages.dock import _atomic_text, _library_records, _titled

    dest = Path(dest)
    stamp = dest / ".library_sha"
    index = dest / INDEX_NAME
    try:
        if library_sha and index.is_file() and stamp.read_text().strip() == library_sha:
            return index
    except OSError:
        pass
    import shutil

    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    paths = []
    for rec in _library_records(titled_sdf):
        lid = str(rec.fields[ID_FIELD])
        _check_id(lid, titled_sdf)
        _check_no_stale_score(rec.fields, lid, titled_sdf)
        path = dest / f"{lid}.sdf"
        _atomic_text(path, _titled(rec))
        paths.append(path)
    _atomic_text(index, "\n".join(str(p) for p in paths) + "\n")
    _atomic_text(stamp, library_sha + "\n")
    return index


# --- running one unit ----------------------------------------------------------------------

def pipeline_cmd(receptor, index, savedir, center, size, seed: int, search_mode: str,
                 batch_size: int) -> list[str]:
    """The whole `unidock_tools` argv, in one place: it is the only thing in this module
    that depends on the real CLI's spelling."""
    cmd = [tools_bin(), "unidock_pipeline", "--receptor", str(receptor),
           "--ligand_index", str(index), "--savedir", str(savedir)]
    for axis, c in zip("xyz", center):
        cmd += [f"--center_{axis}", repr(float(c))]
    for axis, s in zip("xyz", size):
        cmd += [f"--size_{axis}", repr(float(s))]
    cmd += ["--search_mode", str(search_mode), "--seed", str(int(seed)),
            "--batch_size", str(int(batch_size))]
    return cmd


def _score_of(fields: dict):
    try:
        return float(fields[SCORE_FIELD])
    except (KeyError, TypeError, ValueError):
        return None


def collect(savedir, out_sdf, filtered_tsv) -> tuple[int, list[tuple[str, str]]]:
    """Fold `unidock_tools`'s per-ligand SDFs into one multi-record SDF.

    Records are ordered by `<ID>` as they came, poses of one ligand kept contiguous (the
    whole stage's parsing contract), and every record keeps its input SD fields. Returns
    (records written, [(ID, reason) dropped]) and writes the dropped ones to
    `filtered_tsv`."""
    from polyphar.stages.dock import _atomic_text, _read_text, _records

    savedir = Path(savedir)
    out, dropped = [], []
    for path in sorted(savedir.glob("*.sdf"), key=lambda p: p.name):
        recs = _records(_read_text(path))
        if not recs:
            continue
        lid = recs[0].fields.get(ID_FIELD)
        if not lid:
            dropped.append((path.name, f"{path.name}: no <{ID_FIELD}> field; the ligand "
                                       f"cannot be mapped back to the library"))
            continue
        bad = None
        for rec in recs:
            if STALE_SCORE_FIELD in rec.fields:
                # split_library refuses such a library up front; this is the backstop for a
                # record that acquired the field some other way. Never resolved silently:
                # the stage cannot tell the engine's score from the library's leftover.
                bad = (f"carries a <{STALE_SCORE_FIELD}> field "
                       f"({rec.fields[STALE_SCORE_FIELD]!r}) beside <{SCORE_FIELD}>; "
                       f"unidock_tools passes input SD fields through, so which number is "
                       f"this run's result cannot be told")
                break
            score = _score_of(rec.fields)
            if score is None:
                bad = f"no <{SCORE_FIELD}> field"
                break
            if score != score or abs(score) > ABSURD_SCORE:   # NaN or absurd
                bad = (f"score {rec.fields.get(SCORE_FIELD)!r} is outside "
                       f"+-{ABSURD_SCORE:g} kcal/mol (upstream issue #205); the whole "
                       f"ligand is dropped")
                break
        if bad:
            dropped.append((str(lid), bad))
            continue
        out.extend("".join(r.lines) for r in recs)
    _atomic_text(Path(out_sdf), "".join(out))
    rows = ["id\treason"] + [f"{i}\t{r}" for i, r in dropped]
    _atomic_text(Path(filtered_tsv), "\n".join(rows) + "\n")
    return len(out), dropped


def run_unidock(receptor, index, center, size, out_sdf, log, seed: int,
                search_mode: str = DEFAULT_SEARCH_MODE,
                batch_size: int = DEFAULT_BATCH_SIZE,
                timeout: float | None = None) -> None:
    """Dock the indexed ligands into the box and write ONE multi-record SDF.

    Raises UniDockError on a non-zero exit, on a timeout, or when the run produced no
    usable record at all. The log holds the tool's output plus a summary of what was
    filtered."""
    import shutil
    import tempfile

    from polyphar.stages.dock import _atomic_text

    out_sdf, log = Path(out_sdf), Path(log)
    filtered = out_sdf.with_name(out_sdf.stem + ".filtered.tsv")
    # A re-draw that FAILS must not leave the previous draw's filtered.tsv beside a record
    # that no longer exists: stale provenance reading as this draw's. Removed before the
    # call; written again only when the draw succeeds.
    filtered.unlink(missing_ok=True)
    check_version()
    model = require_gpu()
    work = Path(tempfile.mkdtemp(dir=out_sdf.parent, prefix=f".{out_sdf.name}.unidock."))
    savedir = work / "poses"
    savedir.mkdir(parents=True, exist_ok=True)
    cmd = pipeline_cmd(receptor, index, savedir, center, size, seed, search_mode, batch_size)
    text = ""
    try:
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, errors="replace",
                               timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            text = _as_text(exc.stdout) + _as_text(exc.stderr)
            raise UniDockError(f"timeout after {timeout:g} s") from None
        text = (p.stdout or "") + (p.stderr or "")
        if p.returncode != 0:
            tail = " | ".join(text.strip().splitlines()[-5:])
            raise UniDockError(f"unidock_tools exited with {p.returncode}: {tail}")
        tmp_out = work / "merged.sdf"
        n, dropped = collect(savedir, tmp_out, work / "filtered.tsv")
        if n == 0:
            raise UniDockError("unidock_tools produced no usable record "
                               f"({len(dropped)} ligand(s) filtered)")
        _atomic_text(log, f"{text}\n--- polyphar: GPU {model}, search_mode {search_mode}, "
                          f"batch_size {batch_size}, seed {seed} ---\n"
                          f"--- polyphar: {len(dropped)} ligand(s) filtered ---\n"
                          + "".join(f"  dropped {i}: {r}\n" for i, r in dropped))
        shutil.copyfile(work / "filtered.tsv", filtered)
        os.replace(tmp_out, out_sdf)
    except UniDockError as exc:
        _atomic_text(log, f"{text}\n--- polyphar: {exc} ---\n")
        raise
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _as_text(x) -> str:
    if not x:
        return ""
    return x.decode(errors="replace") if isinstance(x, bytes) else str(x)
