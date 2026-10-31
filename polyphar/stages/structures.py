"""Stage 02: PDB chains selected by SIFTS accession mapping, plus AlphaFold models.

For every (UniProt accession, PDB entry) pair listed in 01_uniprot/uniprot.csv
the chain is chosen from SIFTS (`SP_PRIMARY == accession`, strict: no ortholog
fallback); among matching chains the one with the largest UniProt coverage wins,
ties go to the smallest chain ID. Pairs without such a chain are reported as
unmapped and never fall back to another chain.

Chain IDs are AUTHOR chain IDs (mmCIF auth_asym_id), as in SIFTS. ProDy >= 2.6
matches `parseMMCIF(chain=...)` against label_asym_id, so extraction selects by
`segment=` (ProDy stores auth_asym_id as the segment name).

Outputs in 02_structures/:
  <UniProt>/{pdb}.cif            RCSB mmCIF (pdb id lowercase)
  <UniProt>/{pdb}_{chain}.pdb    extracted author chain, all models
  <UniProt>/AF-<acc>-F1.pdb      AlphaFold DB model (F1 only, no PAE json)
  sifts/                         pdb_chain_uniprot.tsv.gz + last_modified.txt
  Protein_info.csv               one row per extracted chain / AlphaFold model
  unmapped.tsv                   UniProt ID, PDB, reason
  alphafold.tsv                  accession, status (found|missing), version
  failed.tsv                     unmapped pairs and failed downloads/extractions
  manifest.json                  written last; SIFTS Last-Modified, AF versions
"""

from __future__ import annotations

import gzip
import io
import json
import os
from pathlib import Path

import requests

from polyphar.common import WorkUnits, append_failure, http_session, write_manifest

STAGE = "structures"
STAGE_DIR = "02_structures"
SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
SIFTS_NAME = "pdb_chain_uniprot.tsv.gz"
RCSB_URL = "https://files.rcsb.org/download/{}.cif"
AF_API = "https://alphafold.ebi.ac.uk/api/prediction/"
TIMEOUT = 120          # seconds per request (http_session has no default timeout)
SIFTS_TIMEOUT = 600
MAX_WORKERS = 8

REASON_NO_CHAIN = "no SIFTS chain for accession"
REASON_NO_PDB = "PDB absent from SIFTS"
INFO_COLUMNS = ["UniProt ID", "Gene Name", "PDB", "Chain ID", "AlphaFold", "Resolution",
                "Structure method", "SIFTS coverage", "File"]
PDB_RECORDS = ("HEADER", "OBSLTE", "TITLE ", "SPLIT ", "CAVEAT", "COMPND", "SOURCE",
               "KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE",
               "JRNL  ", "REMARK", "DBREF", "SEQADV", "SEQRES", "MODRES", "HET   ",
               "FORMUL", "HELIX ", "SHEET ", "CRYST1", "ORIGX", "SCALE", "MODEL ",
               "ATOM  ", "HETATM")


# --- small helpers ----------------------------------------------------------------

def _atomic_write(path: Path, data) -> None:
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    if isinstance(data, str):
        tmp.write_text(data, encoding="utf-8")
    else:
        tmp.write_bytes(data)
    os.replace(tmp, path)


def _head(path: Path, n: int) -> bytes:
    try:
        with open(path, "rb") as fh:
            return fh.read(n)
    except OSError:
        return b""


def _valid_cif(path: Path) -> bool:
    return _head(path, 5) == b"data_"


def _has_atoms(path: Path) -> bool:
    try:
        with open(path, "rb") as fh:
            return any(line.startswith(b"ATOM") for line in fh)
    except OSError:
        return False


def _valid_pdb(path: Path) -> bool:
    return _has_atoms(path)


def _is_pdb_text(data: bytes) -> bool:
    return data.decode("ascii", errors="replace").startswith(PDB_RECORDS) and (
        data.startswith(b"ATOM") or b"\nATOM" in data)


def _valid_af(path: Path) -> bool:
    return _head(path, 6).decode("ascii", errors="replace").startswith(PDB_RECORDS) \
        and _has_atoms(path)


def _prody():
    import prody
    prody.LOGGER.verbosity = "none"    # ProDy logs every parse to stderr
    return prody


# --- SIFTS ---------------------------------------------------------------------------

def download_sifts(dest_dir, session):
    """SIFTS pdb_chain_uniprot.tsv.gz into dest_dir. Returns (path, Last-Modified).

    An existing download (with its recorded Last-Modified) is reused, so a resumed
    run maps chains with the same SIFTS release. Any failure raises RuntimeError."""
    dest_dir = Path(dest_dir)
    path, lm_path = dest_dir / SIFTS_NAME, dest_dir / "last_modified.txt"
    if _head(path, 2) == b"\x1f\x8b" and lm_path.is_file():
        return path, lm_path.read_text(encoding="utf-8").strip()
    try:
        resp = session.get(SIFTS_URL, timeout=SIFTS_TIMEOUT)
        resp.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(f"SIFTS download failed ({SIFTS_URL}): {exc!r}") from exc
    try:
        head = gzip.decompress(resp.content)[:400].decode("utf-8", errors="replace")
    except (OSError, EOFError) as exc:
        raise RuntimeError(f"SIFTS download is not a valid gzip file: {exc}") from exc
    if "PDB\tCHAIN\tSP_PRIMARY" not in head:
        raise RuntimeError(f"SIFTS download has an unexpected header: {head[:200]!r}")
    last_modified = resp.headers.get("Last-Modified", "")
    dest_dir.mkdir(parents=True, exist_ok=True)
    _atomic_write(path, resp.content)
    _atomic_write(lm_path, last_modified + "\n")   # written last: marks the download done
    return path, last_modified


def load_sifts(path):
    """SIFTS pdb_chain_uniprot table (plain or gz); first line is a comment."""
    import pandas as pd
    return pd.read_csv(path, sep="\t", header=1, dtype=str, keep_default_na=False)


def select_chains(pairs, sifts):
    """Pick one author chain per (UniProt ID, PDB) pair.

    Returns (selected[UniProt ID, PDB, Chain ID, SP coverage],
             unmapped[UniProt ID, PDB, reason]). SP coverage is the number of
    UniProt residues the chain covers (sum of SP_END - SP_BEG + 1 over its
    segments); the largest wins, ties go to the smallest chain ID."""
    import pandas as pd

    accs = set(pairs["UniProt ID"])
    pdbs = {p.lower() for p in pairs["PDB"]}
    s = sifts[sifts["PDB"].isin(pdbs)].copy()
    present = set(s["PDB"])
    s = s[s["SP_PRIMARY"].isin(accs)]
    beg = pd.to_numeric(s["SP_BEG"], errors="coerce")
    end = pd.to_numeric(s["SP_END"], errors="coerce")
    s["cov"] = (end - beg + 1).fillna(0).clip(lower=0).astype(int)
    cov = s.groupby(["PDB", "SP_PRIMARY", "CHAIN"])["cov"].sum()
    chains = {}
    for (pdb, acc, chain), c in cov.items():
        chains.setdefault((pdb, acc), []).append((chain, int(c)))

    selected, unmapped = [], []
    for acc, pdb in zip(pairs["UniProt ID"], pairs["PDB"]):
        key = (pdb.lower(), acc)
        if key[0] not in present:
            unmapped.append((acc, pdb, REASON_NO_PDB))
        elif key not in chains:
            unmapped.append((acc, pdb, REASON_NO_CHAIN))
        else:
            chain, c = min(chains[key], key=lambda t: (-t[1], t[0]))
            selected.append((acc, pdb, chain, c))
    return (pd.DataFrame(selected, columns=["UniProt ID", "PDB", "Chain ID", "SP coverage"]),
            pd.DataFrame(unmapped, columns=["UniProt ID", "PDB", "reason"]))


# --- PDB entries -------------------------------------------------------------------------

def download_mmcif(pdb_id, dest_dir, session) -> Path:
    """RCSB mmCIF for pdb_id as <dest_dir>/<pdb lowercase>.cif; reuses a valid file."""
    dest_dir = Path(dest_dir)
    path = dest_dir / f"{pdb_id.lower()}.cif"
    if _valid_cif(path):
        return path
    resp = session.get(RCSB_URL.format(pdb_id.upper()), timeout=TIMEOUT)
    resp.raise_for_status()
    if not resp.content.startswith(b"data_"):
        raise ValueError(f"RCSB response for {pdb_id} is not mmCIF: {resp.content[:60]!r}")
    dest_dir.mkdir(parents=True, exist_ok=True)
    _atomic_write(path, resp.content)
    return path


def _header_text(cif) -> str:
    """The mmCIF text without the (large) atom_site loops."""
    out, skipping, lines = [], False, Path(cif).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
    for i, line in enumerate(lines):
        if skipping:
            if line.strip() == "#":
                skipping = False
                out.append(line)
            continue
        if line.strip() == "loop_" and i + 1 < len(lines) and \
                lines[i + 1].startswith(("_atom_site.", "_atom_site_anisotrop.")):
            skipping = True
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def cif_header(cif) -> dict:
    """Resolution (float, or "-") and Structure method (`_exptl.method`, ';'-joined)."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    d = MMCIF2Dict(io.StringIO(_header_text(cif)))
    resolution = "-"
    for key in ("_refine.ls_d_res_high", "_refine_hist.d_res_high",
                "_em_3d_reconstruction.resolution"):
        for value in d.get(key, []):
            try:
                resolution = float(value)
                break
            except (TypeError, ValueError):
                continue
        if resolution != "-":
            break
    methods = [m for m in d.get("_exptl.method", []) if m not in ("?", ".")]
    return {"Resolution": resolution, "Structure method": ";".join(methods) or "-"}


def extract_chain(cif, chain, out_pdb) -> dict:
    """Write author chain `chain` of `cif` (all models) to `out_pdb`.

    A multi-character chain ID is written as `A` (the PDB format has one column).
    Returns cif_header(cif). Raises ValueError when the chain has no atoms."""
    prody = _prody()
    out_pdb = Path(out_pdb)
    atoms = prody.parseMMCIF(str(cif), segment=chain)
    if atoms is None or atoms.numAtoms() == 0:
        raise ValueError(f"chain {chain} not found in {Path(cif).name}")
    atoms.setChids([chain if len(chain) == 1 else "A"] * atoms.numAtoms())
    tmp = out_pdb.with_name(out_pdb.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        prody.writePDBStream(fh, atoms)
    os.replace(tmp, out_pdb)
    return cif_header(cif)


# --- AlphaFold ---------------------------------------------------------------------------

def fetch_alphafold(acc, dest_dir, session):
    """AlphaFold DB F1 model of `acc` as <dest_dir>/AF-<acc>-F1.pdb.

    Returns (path, version), or (None, None) when AFDB has no F1 model (404 or no
    `AF-<acc>-F1` entry); nothing is written then. A body that is not a PDB file
    raises ValueError."""
    resp = session.get(AF_API + acc, timeout=TIMEOUT)
    if resp.status_code == 404:
        return None, None
    resp.raise_for_status()
    entry_id = f"AF-{acc}-F1"
    entry = next((e for e in resp.json() if e.get("entryId") == entry_id), None)
    if entry is None:
        return None, None
    version = int(entry["latestVersion"])
    model = session.get(entry["pdbUrl"], timeout=TIMEOUT)
    model.raise_for_status()
    if not _is_pdb_text(model.content):
        raise ValueError(f"AlphaFold model for {acc} is not a PDB file: {model.content[:60]!r}")
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / f"{entry_id}.pdb"
    _atomic_write(path, model.content)
    return path, version


# --- stage ---------------------------------------------------------------------------------

def _read_uniprot(layout):
    """(accessions, Entry Name by accession, pairs DataFrame[UniProt ID, PDB])."""
    import pandas as pd

    if not layout.uniprot_csv.is_file():
        raise FileNotFoundError(f"{layout.uniprot_csv} not found; run `polyphar uniprot` first")
    df = pd.read_csv(layout.uniprot_csv, dtype=str, keep_default_na=False)
    accessions, names, pairs = [], {}, []
    for acc, name, pdb in zip(df["Entry"], df["Entry Name"], df["PDB"]):
        if acc in names:
            continue
        accessions.append(acc)
        names[acc] = name
        ids = dict.fromkeys(p.strip().upper() for p in pdb.split(";") if p.strip())
        pairs.extend((acc, p) for p in ids)
    return accessions, names, pd.DataFrame(pairs, columns=["UniProt ID", "PDB"])


def _chain_pdb(d, acc, pdb, chain) -> Path:
    return d / acc / f"{pdb.lower()}_{chain}.pdb"


def _af_pdb(d, acc) -> Path:
    return d / acc / f"AF-{acc}-F1.pdb"


def _local_sifts(d):
    sifts_dir = d / "sifts"
    if _head(sifts_dir / SIFTS_NAME, 2) == b"\x1f\x8b" and (sifts_dir / "last_modified.txt").is_file():
        return sifts_dir / SIFTS_NAME
    return None


def _read_af_state(d) -> dict:
    """accession -> (status, version|None); the last line per accession wins."""
    path = d / "alphafold.tsv"
    state = {}
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 2:
                v = parts[2] if len(parts) > 2 else ""
                state[parts[0]] = (parts[1], int(v) if v.isdigit() else None)
    return state


def _af_done(d, acc, state) -> bool:
    return _valid_af(_af_pdb(d, acc)) or state.get(acc, ("",))[0] == "missing"


def units(layout, args) -> WorkUnits:
    accessions, _, pairs = _read_uniprot(layout)
    d = layout.structures
    sifts_path = _local_sifts(d)
    if sifts_path is not None:
        selected, _ = select_chains(pairs, load_sifts(sifts_path))
        n_pairs = sum(not _valid_pdb(_chain_pdb(d, a, p, c)) for a, p, c in
                      zip(selected["UniProt ID"], selected["PDB"], selected["Chain ID"]))
    else:   # not mapped yet: a pair is done if any extracted chain of it exists
        n_pairs = sum(not any(_valid_pdb(f) for f in (d / a).glob(f"{p.lower()}_*.pdb"))
                      for a, p in zip(pairs["UniProt ID"], pairs["PDB"]))
    state = _read_af_state(d)
    n_af = sum(not _af_done(d, a, state) for a in accessions)
    return WorkUnits(n_pairs + n_af, "download")


def _pdb_job(d, acc, pdb, chain, session) -> dict:
    cif = download_mmcif(pdb, d / acc, session)
    return extract_chain(cif, chain, _chain_pdb(d, acc, pdb, chain))


def _run_pool(jobs, nprocs, progress, on_result, on_error):
    """Run (key, fn, args) jobs in a thread pool; results are handled in the caller's
    thread, which is the only one that touches progress and failed.tsv."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    if not jobs:
        return
    ex = ThreadPoolExecutor(max_workers=max(1, min(nprocs, MAX_WORKERS)))
    try:
        futures = {ex.submit(fn, *a): key for key, fn, a in jobs}
        for fut in as_completed(futures):
            key = futures[fut]
            try:
                result = fut.result()
            except Exception as exc:   # one entry failing never stops the stage
                on_error(key, exc)
            else:
                on_result(key, result)
            progress.advance(1)
    except BaseException:
        ex.shutdown(wait=False, cancel_futures=True)
        raise
    ex.shutdown(wait=True)


def _reason(exc) -> str:
    return f"{type(exc).__name__}: {exc}"


def _uniprot_release(layout):
    path = layout.stage_dir("01_uniprot") / "manifest.json"
    try:
        return json.loads(path.read_text(encoding="utf-8"))["versions"]["uniprot"]
    except (OSError, ValueError, KeyError, TypeError):
        return None


def run(layout, args, progress) -> dict:
    import pandas as pd

    accessions, names, pairs = _read_uniprot(layout)
    d = layout.structures
    d.mkdir(parents=True, exist_ok=True)
    session = http_session()
    sifts_was_local = _local_sifts(d) is not None
    sifts_path, last_modified = download_sifts(d / "sifts", session)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    _atomic_write(d / "failed.tsv", "item\tstage\treason\n")

    selected, unmapped = select_chains(pairs, load_sifts(sifts_path))
    _atomic_write(d / "unmapped.tsv", "UniProt ID\tPDB\treason\n" + "".join(
        f"{a}\t{p}\t{r}\n" for a, p, r in unmapped.itertuples(index=False)))
    for a, p, r in unmapped.itertuples(index=False):
        append_failure(d, f"{a}/{p}", STAGE, r)
    if not sifts_was_local and len(unmapped):
        progress.advance(len(unmapped))   # units() counted them before SIFTS was known

    af_state = _read_af_state(d)
    af_log = d / "alphafold.tsv"
    if not af_log.is_file():
        _atomic_write(af_log, "accession\tstatus\tversion\n")
    headers, failed = {}, 0
    jobs = [(("pdb", a, p), _pdb_job, (d, a, p, c, session))
            for a, p, c, _ in selected.itertuples(index=False)
            if not _valid_pdb(_chain_pdb(d, a, p, c))]
    jobs += [(("af", a, None), fetch_alphafold, (a, d / a, session))
             for a in accessions if not _af_done(d, a, af_state)]

    def on_result(key, result):
        kind, acc, pdb = key
        if kind == "pdb":
            headers[(acc, pdb)] = result
        else:
            status, version = ("found", result[1]) if result[0] is not None else ("missing", None)
            af_state[acc] = (status, version)
            with open(af_log, "a", encoding="utf-8") as fh:
                fh.write(f"{acc}\t{status}\t{'' if version is None else version}\n")

    def on_error(key, exc):
        nonlocal failed
        kind, acc, pdb = key
        failed += 1
        append_failure(d, f"{acc}/{pdb if kind == 'pdb' else 'AlphaFold'}", STAGE, _reason(exc))

    _run_pool(jobs, args.nprocs, progress, on_result, on_error)

    rows, af_versions = [], set()
    by_acc = {}
    for a, p, c, cov in selected.itertuples(index=False):
        by_acc.setdefault(a, []).append((p, c, cov))
    for acc in accessions:
        for pdb, chain, cov in by_acc.get(acc, []):
            out = _chain_pdb(d, acc, pdb, chain)
            if not _valid_pdb(out):
                continue
            cif = d / acc / f"{pdb.lower()}.cif"
            hdr = headers.get((acc, pdb)) or (
                cif_header(cif) if _valid_cif(cif) else {"Resolution": "-", "Structure method": "-"})
            rows.append([acc, names[acc], pdb.upper(), chain, False, hdr["Resolution"],
                         hdr["Structure method"], cov, f"{STAGE_DIR}/{acc}/{out.name}"])
        af = _af_pdb(d, acc)
        if _valid_af(af):
            rows.append([acc, names[acc], "-", "A", True, "-", "-", "-",
                         f"{STAGE_DIR}/{acc}/{af.name}"])
            version = af_state.get(acc, ("", None))[1]
            if version is not None:
                af_versions.add(version)
    info = pd.DataFrame(rows, columns=INFO_COLUMNS)
    tmp = layout.protein_info.with_name(layout.protein_info.name + ".tmp")
    info.to_csv(tmp, index=False)
    os.replace(tmp, layout.protein_info)
    _atomic_write(af_log, "accession\tstatus\tversion\n" + "".join(
        f"{a}\t{s}\t{'' if v is None else v}\n" for a, (s, v) in sorted(af_state.items())))

    reasons = unmapped["reason"].value_counts()
    counts = {
        "accessions": len(accessions),
        "pdb_entries_listed": len(pairs),
        "chains_extracted": int((~info["AlphaFold"]).sum()) if len(info) else 0,
        "unmapped": {r: int(reasons[r]) for r in sorted(reasons.index)},
        "af_found": int(info["AlphaFold"].sum()) if len(info) else 0,
        "af_missing": sum(af_state.get(a, ("",))[0] == "missing" for a in accessions),
        "failed": failed,
    }
    versions = {"sifts_last_modified": last_modified,
                "alphafold_versions": sorted(af_versions)}
    release = _uniprot_release(layout)
    if release is not None:
        versions["uniprot_release"] = release
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE, args={"nprocs": args.nprocs}, inputs=[layout.uniprot_csv, sifts_path],
                   versions=versions, counts=counts, timing={})
    return counts
