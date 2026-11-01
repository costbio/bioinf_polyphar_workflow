"""Stage 03: repair the extracted PDB chains with pdbfixer (parallel, atomic).

Inputs are the non-AlphaFold rows of 02_structures/Protein_info.csv (column `File`,
relative to the run folder); AlphaFold models are not passed through pdbfixer and
nothing is globbed. Output: 03_pdbfixer/<UniProt>/<same file name>.

The legacy call sequence is kept exactly (it reproduces the 2023 receptors, although
removing heterogens after adding hydrogens is unusual):
  findMissingResidues, findMissingAtoms, addMissingAtoms, addMissingHydrogens(pH),
  removeHeterogens(keepWater=False)

Each output is written to a unique temp file `<out>.<random>.tmp` in the same directory and
renamed when complete, so a killed run never leaves a truncated file that looks done and
concurrent writers never share a temp name. A unit is done when its output exists and has
ATOM records. `*.tmp` files are ignored; those left by killed processes are harmless
leftovers (they are not deleted: they may belong to a live process). Inputs are never modified or deleted;
failures are listed in failed.tsv with the reason.
"""

from __future__ import annotations

import contextlib
import glob
import os
import random
import tempfile
from pathlib import Path

from polyphar.common import WorkUnits, append_failure, write_manifest

STAGE = "fix"
DEFAULT_PH = 7.0
H_SEED = 0   # seeds OpenMM's random hydrogen / missing-atom placement (deterministic output)


def reference_platform():
    """OpenMM's CPU-independent Reference platform: minimisation results do not depend on
    the host's CPU/GPU or thread count."""
    import openmm

    return openmm.Platform.getPlatformByName("Reference")


@contextlib.contextmanager
def seeded_placement(seed: int = H_SEED):
    """Make pdbfixer/OpenMM placement reproducible. OpenMM's Modeller.addHydrogens places new
    atoms with Python's `random` (unseeded); the energy minimiser that follows runs on the
    Reference platform (see reference_platform). The global `random` state is seeded for
    the block and restored afterwards, so callers' own random streams are untouched.
    pdbfixer and OpenMM use no numpy RNG here; addMissingAtoms takes `seed=` for its integrator."""
    state = random.getstate()
    try:
        random.seed(seed)
        yield
    finally:
        random.setstate(state)


def add_arguments(parser) -> None:
    parser.add_argument("--ph", type=float, default=DEFAULT_PH,
                        help=f"pH for adding hydrogens (default: {DEFAULT_PH})")


# --- inputs / outputs ----------------------------------------------------------------

def select_inputs(layout) -> list[Path]:
    """Extracted (non-AlphaFold) chain files listed in Protein_info.csv, in listed order."""
    import pandas as pd

    path = layout.protein_info
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; run `polyphar structures` first")
    info = pd.read_csv(path, dtype=str, keep_default_na=False)
    is_af = info["AlphaFold"].str.strip().str.lower() == "true"
    return [layout.run_dir / f for f in info.loc[~is_af, "File"]]


def _out_path(src: Path, structures_dir: Path, out_dir: Path) -> Path:
    return Path(out_dir) / Path(src).relative_to(structures_dir)


def _done(path: Path) -> bool:
    try:
        with open(path, "rb") as fh:
            return any(line.startswith(b"ATOM") for line in fh)
    except OSError:
        return False


def _size(path: Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def _pending(layout) -> tuple[list[Path], list[Path]]:
    inputs = select_inputs(layout)
    todo = [p for p in inputs
            if not _done(_out_path(p, layout.structures, layout.pdbfixer))]
    return inputs, todo


def units(layout, args) -> WorkUnits:
    _, todo = _pending(layout)
    return WorkUnits(len(todo), "bytes", [float(_size(p)) for p in todo])


# --- the work ------------------------------------------------------------------------

def fix_one(src, structures_dir, out_dir, ph: float = DEFAULT_PH, token: str = "") -> None:
    """Fix `src` and write it to out_dir/<path relative to structures_dir>, atomically.
    `token` (one per stage run) goes into the temp name so the stage can recognise its own
    leftovers. Raises on failure, leaving no output and the input untouched."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    src = Path(src)
    out = _out_path(src, Path(structures_dir), Path(out_dir))
    out.parent.mkdir(parents=True, exist_ok=True)
    # A unique temp file in the output directory: an orphaned worker, a rerun or a second
    # polyphar process on the same run folder never share a temp name.
    fd, tmp_name = tempfile.mkstemp(dir=out.parent, prefix=f"{out.name}.{token}." if token else out.name + ".", suffix=".tmp")
    tmp = Path(tmp_name)
    try:
        fixer = PDBFixer(filename=str(src), platform=reference_platform())
        with seeded_placement():
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms(seed=H_SEED)
            fixer.addMissingHydrogens(ph)
            fixer.removeHeterogens(keepWater=False)
        with os.fdopen(fd, "w") as fh:
            fd = None
            PDBFile.writeFile(fixer.topology, fixer.positions, fh)
        if not _done(tmp):
            raise ValueError("pdbfixer produced no ATOM records")
        os.chmod(tmp, 0o644)   # mkstemp creates 0600
        os.replace(tmp, out)
    except BaseException:
        if fd is not None:
            os.close(fd)
        tmp.unlink(missing_ok=True)
        raise


class FixError(RuntimeError):
    """Worker failure carried across the process boundary as a plain message, so it
    always pickles (an arbitrary exception that cannot be unpickled breaks the pool)."""


_CRASH_ENV = "POLYPHAR_FIX_TEST_CRASH_ON"   # test hook: worker dies hard on this file name


def _worker(src, structures_dir, out_dir, ph, token=""):
    if os.environ.get(_CRASH_ENV) == Path(src).name:
        os._exit(1)
    try:
        fix_one(src, structures_dir, out_dir, ph, token)
    except Exception as exc:   # fix_one already removed its .tmp
        raise FixError(f"{type(exc).__name__}: {exc}") from None


def _reason(exc) -> str:
    return str(exc) if isinstance(exc, FixError) else f"{type(exc).__name__}: {exc}"


REASON_DIED = "worker process died (BrokenProcessPool)"


def _new_token() -> str:
    import uuid

    return uuid.uuid4().hex[:8]


def _sweep_own_temps(units, structures_dir, out_dir, token) -> None:
    """Delete this run's temp files of `units`. Only called once the pool that wrote them
    is gone (a killed worker leaves its temp behind); other runs' temps (other token) are
    never touched."""
    for src in units:
        out = _out_path(src, structures_dir, out_dir)
        for tmp in out.parent.glob(f"{glob.escape(out.name)}.{token}.*.tmp"):
            tmp.unlink(missing_ok=True)


def _versions() -> dict:
    import importlib.metadata

    import openmm
    import pdbfixer

    pf = getattr(pdbfixer, "__version__", None) or importlib.metadata.version("pdbfixer")
    return {"pdbfixer": pf, "openmm": openmm.__version__, "platform": "Reference"}


def run(layout, args, progress) -> dict:
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool

    ph = getattr(args, "ph", DEFAULT_PH)
    d = layout.pdbfixer
    d.mkdir(parents=True, exist_ok=True)
    inputs, todo = _pending(layout)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    failed_tsv = d / "failed.tsv"
    failed_tsv.write_text("item\tstage\treason\n")

    fixed = 0
    failures: list[tuple[str, str]] = []
    token = _new_token()
    ctx = multiprocessing.get_context("spawn")   # workers start clean (no fork of a threaded parent)

    def item_of(src):
        return src.relative_to(layout.structures).as_posix()

    try:
        broken: list[Path] = []   # units whose pool died: cause unknown until run alone
        if todo:
            with ProcessPoolExecutor(max_workers=max(1, args.nprocs), mp_context=ctx) as ex:
                try:
                    futures = {}
                    for p in todo:
                        try:
                            futures[ex.submit(_worker, p, layout.structures, d, ph, token)] = p
                        except BrokenProcessPool:   # pool died while submitting: redone alone
                            broken.append(p)
                    for fut in as_completed(futures):
                        src = futures[fut]
                        try:
                            fut.result()
                        except BrokenProcessPool:
                            broken.append(src)   # a worker died hard; sort it out below
                            continue
                        except Exception as exc:   # one structure failing never stops the stage
                            failures.append((item_of(src), _reason(exc)))
                        else:
                            fixed += 1
                        progress.advance(1, size=_size(src))
                except BaseException:
                    ex.shutdown(wait=False, cancel_futures=True)
                    raise
        # The pool is dead for everything that was still in flight. Redo those one at a
        # time, each in a fresh single-worker pool, so a crashing structure is attributed
        # to itself and the rest are fixed.
        # (the pool and its workers are gone, so temps of those units are ours to delete)
        _sweep_own_temps(broken, layout.structures, d, token)
        for src in sorted(broken, key=todo.index):
            with ProcessPoolExecutor(max_workers=1, mp_context=ctx) as ex:
                try:
                    ex.submit(_worker, src, layout.structures, d, ph, token).result()
                except BrokenProcessPool:
                    _sweep_own_temps([src], layout.structures, d, token)
                    failures.append((item_of(src), REASON_DIED))
                except Exception as exc:
                    failures.append((item_of(src), _reason(exc)))
                except BaseException:
                    ex.shutdown(wait=False, cancel_futures=True)
                    raise
                else:
                    fixed += 1
            progress.advance(1, size=_size(src))
    finally:
        # appended from the parent only, in a fixed order, even when interrupted
        for item, reason in sorted(failures):
            append_failure(d, item, STAGE, reason)

    counts = {"inputs": len(inputs), "fixed": fixed,
              "skipped_existing": len(inputs) - len(todo), "failed": len(failures)}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE, args={"nprocs": args.nprocs, "ph": ph, "h_seed": H_SEED},
                   inputs=[layout.protein_info], versions=_versions(),
                   counts=counts, timing={})
    return counts
