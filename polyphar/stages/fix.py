"""Stage 03: repair the extracted PDB chains with pdbfixer (parallel, atomic).

Inputs are the non-AlphaFold rows of 02_structures/Protein_info.csv (column `File`,
relative to the run folder); AlphaFold models are not passed through pdbfixer and
nothing is globbed. Output: 03_pdbfixer/<UniProt>/<same file name>.

The legacy call sequence is kept exactly (it reproduces the 2023 receptors, although
removing heterogens after adding hydrogens is unusual):
  findMissingResidues, findMissingAtoms, addMissingAtoms, addMissingHydrogens(pH),
  removeHeterogens(keepWater=False)

Each output is written to `<out>.tmp` and renamed when complete, so a killed run never
leaves a truncated file that looks done. A unit is done when its output exists and has
ATOM records; a leftover `.tmp` is overwritten. Inputs are never modified or deleted;
failures are listed in failed.tsv with the reason.
"""

from __future__ import annotations

import os
from pathlib import Path

from polyphar.common import WorkUnits, append_failure, write_manifest

STAGE = "fix"
DEFAULT_PH = 7.0


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

def fix_one(src, structures_dir, out_dir, ph: float = DEFAULT_PH) -> None:
    """Fix `src` and write it to out_dir/<path relative to structures_dir>, atomically.
    Raises on failure, leaving no output and the input untouched."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    src = Path(src)
    out = _out_path(src, Path(structures_dir), Path(out_dir))
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    try:
        fixer = PDBFixer(filename=str(src))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        fixer.removeHeterogens(keepWater=False)
        with open(tmp, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh)
        if not _done(tmp):
            raise ValueError("pdbfixer produced no ATOM records")
        os.replace(tmp, out)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


class FixError(RuntimeError):
    """Worker failure carried across the process boundary as a plain message, so it
    always pickles (an arbitrary exception that cannot be unpickled breaks the pool)."""


_CRASH_ENV = "POLYPHAR_FIX_TEST_CRASH_ON"   # test hook: worker dies hard on this file name


def _worker(src, structures_dir, out_dir, ph):
    if os.environ.get(_CRASH_ENV) == Path(src).name:
        os._exit(1)
    try:
        fix_one(src, structures_dir, out_dir, ph)
    except Exception as exc:   # fix_one already removed its .tmp
        raise FixError(f"{type(exc).__name__}: {exc}") from None


def _reason(exc) -> str:
    return str(exc) if isinstance(exc, FixError) else f"{type(exc).__name__}: {exc}"


REASON_DIED = "worker process died (BrokenProcessPool)"


def _versions() -> dict:
    import importlib.metadata

    import openmm
    import pdbfixer

    pf = getattr(pdbfixer, "__version__", None) or importlib.metadata.version("pdbfixer")
    return {"pdbfixer": pf, "openmm": openmm.__version__}


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
    ctx = multiprocessing.get_context("spawn")   # workers start clean (no fork of a threaded parent)

    def item_of(src):
        return src.relative_to(layout.structures).as_posix()

    try:
        broken: list[Path] = []   # units whose pool died: cause unknown until run alone
        if todo:
            with ProcessPoolExecutor(max_workers=max(1, args.nprocs), mp_context=ctx) as ex:
                try:
                    futures = {ex.submit(_worker, p, layout.structures, d, ph): p for p in todo}
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
        for src in sorted(broken, key=todo.index):
            with ProcessPoolExecutor(max_workers=1, mp_context=ctx) as ex:
                try:
                    ex.submit(_worker, src, layout.structures, d, ph).result()
                except BrokenProcessPool:
                    out = _out_path(src, layout.structures, d)
                    out.with_name(out.name + ".tmp").unlink(missing_ok=True)
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
    write_manifest(d, STAGE, args={"nprocs": args.nprocs, "ph": ph},
                   inputs=[layout.protein_info], versions=_versions(),
                   counts=counts, timing={})
    return counts
