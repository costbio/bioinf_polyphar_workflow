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
  Protein_info.csv               one row per extracted chain / AlphaFold model, rebuilt
                                 from the current uniprot.csv on every run; columns
                                 INFO_COLUMNS, incl. SIFTS segments / SIFTS other accessions
                                 (select_chains) for fusion-chain pocket filtering
  unmapped.tsv                   UniProt ID, PDB, reason
  alphafold.tsv                  accession, status (found|missing), version
  annotation_mismatch.tsv        extracted chains whose cif depositor annotation (_struct_ref
                                 UNP) names other accessions than SIFTS; a warning only
  rejected.tsv                   chains rejected for good (not in the cif, or no ATOM
                                 records); not retried on resume
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
                "Structure method", "SIFTS coverage", "File", "SIFTS segments",
                "SIFTS other accessions"]
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


SEL_COLUMNS = ["UniProt ID", "PDB", "Chain ID", "SP coverage", "SIFTS segments",
               "SIFTS other accessions"]


def select_chains(pairs, sifts):
    """Pick one author chain per (UniProt ID, PDB) pair.

    Returns (selected[UniProt ID, PDB, Chain ID, SP coverage, SIFTS segments,
    SIFTS other accessions], unmapped[UniProt ID, PDB, reason]). SP coverage is the number
    of UniProt residues the chain covers (sum of SP_END - SP_BEG + 1 over its segments);
    the largest wins, ties go to the smallest chain ID. SIFTS segments are this accession's
    segments on the chosen chain in author numbering, `PDB_BEG:PDB_END` joined by ';' in
    SIFTS order (a bound SIFTS leaves empty or None, e.g. an unobserved end, is kept as
    given). SIFTS other accessions are the other SP_PRIMARY accessions SIFTS maps to the
    same chain (';'-joined, sorted; empty for a single-protein chain): a fusion/chimeric
    construct whose pockets the pockets stage checks against the segments."""
    import pandas as pd

    accs = set(pairs["UniProt ID"])
    pdbs = {p.lower() for p in pairs["PDB"]}
    s = sifts[sifts["PDB"].isin(pdbs)].copy()
    present = set(s["PDB"])
    on_chain = {}   # (pdb, chain) -> every accession SIFTS maps to that chain
    for pdb, chain, acc in zip(s["PDB"], s["CHAIN"], s["SP_PRIMARY"]):
        on_chain.setdefault((pdb, chain), set()).add(acc)
    s = s[s["SP_PRIMARY"].isin(accs)]
    beg = pd.to_numeric(s["SP_BEG"], errors="coerce")
    end = pd.to_numeric(s["SP_END"], errors="coerce")
    s["cov"] = (end - beg + 1).fillna(0).clip(lower=0).astype(int)
    cov = s.groupby(["PDB", "SP_PRIMARY", "CHAIN"])["cov"].sum()
    chains = {}
    for (pdb, acc, chain), c in cov.items():
        chains.setdefault((pdb, acc), []).append((chain, int(c)))
    segments = {}
    for pdb, chain, acc, b, e in zip(s["PDB"], s["CHAIN"], s["SP_PRIMARY"], s["PDB_BEG"],
                                     s["PDB_END"]):
        segments.setdefault((pdb, chain, acc), []).append(f"{b.strip()}:{e.strip()}")

    selected, unmapped = [], []
    for acc, pdb in zip(pairs["UniProt ID"], pairs["PDB"]):
        key = (pdb.lower(), acc)
        if key[0] not in present:
            unmapped.append((acc, pdb, REASON_NO_PDB))
        elif key not in chains:
            unmapped.append((acc, pdb, REASON_NO_CHAIN))
        else:
            chain, c = min(chains[key], key=lambda t: (-t[1], t[0]))
            others = sorted(on_chain.get((key[0], chain), set()) - {acc})
            selected.append((acc, pdb, chain, c,
                             ";".join(segments.get((key[0], chain, acc), [])),
                             ";".join(others)))
    return (pd.DataFrame(selected, columns=SEL_COLUMNS),
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


class ChainRejected(ValueError):
    """The selected chain cannot be used, and retrying will not change that (it is
    not in the cif, or it has no ATOM records).
    Recorded in rejected.tsv and failed.tsv; not retried on resume."""


def _cif_dict(cif):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    return MMCIF2Dict(io.StringIO(_header_text(cif)))


def _header_from(d) -> dict:
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


def _base_acc(acc: str) -> str:
    return acc.split("-")[0]       # P12345-2 (isoform) -> P12345


def _chain_accessions(d) -> dict:
    """author chain -> set of UniProt accessions, from the cif's own reference
    annotation: _struct_ref (db_name UNP) joined to author chains through
    _struct_ref_seq.pdbx_strand_id, or, when that loop is absent, through
    _entity_poly.pdbx_strand_id. RCSB mmCIFs do not carry _pdbx_sifts_unp_segments."""
    ref = {rid: (db, acc) for rid, db, acc in zip(
        d.get("_struct_ref.id", []), d.get("_struct_ref.db_name", []),
        d.get("_struct_ref.pdbx_db_accession", []))}
    out = {}
    if "_struct_ref_seq.ref_id" in d:
        for rid, strand in zip(d["_struct_ref_seq.ref_id"],
                               d.get("_struct_ref_seq.pdbx_strand_id", [])):
            db, acc = ref.get(rid, (None, None))
            if db == "UNP":
                out.setdefault(strand, set()).add(_base_acc(acc))
    else:
        strands = {e: s.replace(" ", "").split(",") for e, s in zip(
            d.get("_entity_poly.entity_id", []), d.get("_entity_poly.pdbx_strand_id", []))}
        for rid, ent in zip(d.get("_struct_ref.id", []), d.get("_struct_ref.entity_id", [])):
            db, acc = ref.get(rid, (None, None))
            if db == "UNP":
                for strand in strands.get(ent, []):
                    out.setdefault(strand, set()).add(_base_acc(acc))
    return out


def _verify(chain_accs, chain, acc):
    accs = sorted(chain_accs.get(chain, ()))
    if not accs:
        return "unverified", accs   # no UniProt annotation for this chain (or the entry)
    return ("verified" if _base_acc(acc) in accs else "mismatch"), accs


def verify_chain(cif, chain, acc):
    """Cross-check a SIFTS-selected author chain against the cif's own depositor
    UniProt annotation. Returns (status, cif_accessions) with status "verified",
    "unverified" (no UniProt reference for that chain) or "mismatch" (it names
    other accessions). SIFTS stays the authority: a mismatch is only a warning."""
    return _verify(_chain_accessions(_cif_dict(cif)), chain, acc)


def cif_header(cif) -> dict:
    """Resolution (float, or "-") and Structure method (`_exptl.method`, ';'-joined)."""
    return _header_from(_cif_dict(cif))


def _write_chain(cif, chain, out_pdb) -> None:
    prody = _prody()
    out_pdb = Path(out_pdb)
    atoms = prody.parseMMCIF(str(cif), segment=chain)
    if atoms is None or atoms.numAtoms() == 0:
        raise ChainRejected(f"chain {chain} not found in {Path(cif).name}")
    if atoms.getFlags("hetatm").all():
        raise ChainRejected("chain has no ATOM records")
    atoms.setChids([chain if len(chain) == 1 else "A"] * atoms.numAtoms())
    tmp = out_pdb.with_name(out_pdb.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        prody.writePDBStream(fh, atoms)
    os.replace(tmp, out_pdb)


def extract_chain(cif, chain, out_pdb) -> dict:
    """Write author chain `chain` of `cif` (all models) to `out_pdb`.

    A multi-character chain ID is written as `A` (the PDB format has one column).
    Returns cif_header(cif). Raises ChainRejected (a ValueError), writing nothing,
    when the chain has no atoms or only HETATM records."""
    _write_chain(cif, chain, out_pdb)
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


def _read_rejected(d) -> dict:
    """(UniProt ID, PDB, Chain ID) -> reason for chains rejected in earlier runs."""
    path = d / "rejected.tsv"
    out = {}
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) == 4:
                out[tuple(parts[:3])] = parts[3]
    return out


def _af_done(d, acc, state) -> bool:
    return _valid_af(_af_pdb(d, acc)) or state.get(acc, ("",))[0] == "missing"


def units(layout, args) -> WorkUnits:
    accessions, _, pairs = _read_uniprot(layout)
    d = layout.structures
    sifts_path = _local_sifts(d)
    if sifts_path is not None:
        selected, _ = select_chains(pairs, load_sifts(sifts_path))
        rejected = _read_rejected(d)
        n_pairs = sum(not _valid_pdb(_chain_pdb(d, a, p, c)) and (a, p, c) not in rejected
                      for a, p, c in zip(selected["UniProt ID"], selected["PDB"],
                                         selected["Chain ID"]))
    else:   # not mapped yet: a pair is done if any extracted chain of it exists
        n_pairs = sum(not any(_valid_pdb(f) for f in (d / a).glob(f"{p.lower()}_*.pdb"))
                      for a, p in zip(pairs["UniProt ID"], pairs["PDB"]))
    state = _read_af_state(d)
    n_af = sum(not _af_done(d, a, state) for a in accessions)
    return WorkUnits(n_pairs + n_af, "download")


def _pdb_job(d, acc, pdb, chain, session) -> dict:
    """Download, cross-check against the cif's own annotation, extract."""
    cif = download_mmcif(pdb, d / acc, session)
    cd = _cif_dict(cif)
    status, accs = _verify(_chain_accessions(cd), chain, acc)
    _write_chain(cif, chain, _chain_pdb(d, acc, pdb, chain))
    return {**_header_from(cd), "cif_check": status, "cif_accessions": accs}


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
    rejected = _read_rejected(d)
    rejected_log = d / "rejected.tsv"
    if not rejected_log.is_file():
        _atomic_write(rejected_log, "UniProt ID\tPDB\tChain ID\treason\n")
    jobs = []
    for a, p, c, *_ in selected.itertuples(index=False):
        if _valid_pdb(_chain_pdb(d, a, p, c)):
            continue
        if (a, p, c) in rejected:      # deterministic: not retried, but still reported
            failed += 1
            append_failure(d, f"{a}/{p}", STAGE, f"ChainRejected: {rejected[(a, p, c)]}")
            continue
        jobs.append((("pdb", a, p, c), _pdb_job, (d, a, p, c, session)))
    jobs += [(("af", a, None, None), fetch_alphafold, (a, d / a, session))
             for a in accessions if not _af_done(d, a, af_state)]

    def on_result(key, result):
        kind, acc, pdb, _ = key
        if kind == "pdb":
            headers[(acc, pdb)] = result
        else:
            status, version = ("found", result[1]) if result[0] is not None else ("missing", None)
            af_state[acc] = (status, version)
            with open(af_log, "a", encoding="utf-8") as fh:
                fh.write(f"{acc}\t{status}\t{'' if version is None else version}\n")

    def on_error(key, exc):
        nonlocal failed
        kind, acc, pdb, chain = key
        failed += 1
        if isinstance(exc, ChainRejected):
            reason = str(exc).replace("\t", " ").replace("\n", " ")
            with open(rejected_log, "a", encoding="utf-8") as fh:
                fh.write(f"{acc}\t{pdb}\t{chain}\t{reason}\n")
        append_failure(d, f"{acc}/{pdb if kind == 'pdb' else 'AlphaFold'}", STAGE, _reason(exc))

    _run_pool(jobs, args.nprocs, progress, on_result, on_error)

    rows, af_versions, unverified, mismatches = [], set(), 0, []
    by_acc = {}
    for a, p, c, cov, segs, others in selected.itertuples(index=False):
        by_acc.setdefault(a, []).append((p, c, cov, segs, others))
    for acc in accessions:
        for pdb, chain, cov, segs, others in by_acc.get(acc, []):
            out = _chain_pdb(d, acc, pdb, chain)
            if not _valid_pdb(out):
                continue
            cif = d / acc / f"{pdb.lower()}.cif"
            hdr = headers.get((acc, pdb))
            if hdr is None:            # extracted in an earlier run: re-read from the cif
                hdr = {"Resolution": "-", "Structure method": "-", "cif_check": "unverified",
                       "cif_accessions": []}
                if _valid_cif(cif):
                    cd = _cif_dict(cif)
                    status, accs = _verify(_chain_accessions(cd), chain, acc)
                    hdr = {**_header_from(cd), "cif_check": status, "cif_accessions": accs}
            unverified += hdr["cif_check"] == "unverified"
            if hdr["cif_check"] == "mismatch":
                mismatches.append((acc, pdb.upper(), chain, ",".join(hdr["cif_accessions"])))
            rows.append([acc, names[acc], pdb.upper(), chain, False, hdr["Resolution"],
                         hdr["Structure method"], cov, f"{STAGE_DIR}/{acc}/{out.name}",
                         segs, others])
        af = _af_pdb(d, acc)
        if _valid_af(af):
            rows.append([acc, names[acc], "-", "A", True, "-", "-", "-",
                         f"{STAGE_DIR}/{acc}/{af.name}", "", ""])
            version = af_state.get(acc, ("", None))[1]
            if version is not None:
                af_versions.add(version)
    info = pd.DataFrame(rows, columns=INFO_COLUMNS)
    tmp = layout.protein_info.with_name(layout.protein_info.name + ".tmp")
    info.to_csv(tmp, index=False)
    os.replace(tmp, layout.protein_info)
    _atomic_write(d / "annotation_mismatch.tsv", "UniProt ID\tPDB\tChain ID\tcif_accessions\n"
                  + "".join("\t".join(m) + "\n" for m in mismatches))
    current = set(accessions)   # tables describe the current uniprot.csv only
    _atomic_write(af_log, "accession\tstatus\tversion\n" + "".join(
        f"{a}\t{s}\t{'' if v is None else v}\n" for a, (s, v) in sorted(af_state.items())
        if a in current))

    reasons = unmapped["reason"].value_counts()
    counts = {
        "accessions": len(accessions),
        "pdb_entries_listed": len(pairs),
        "chains_extracted": int((~info["AlphaFold"]).sum()) if len(info) else 0,
        "unmapped": {r: int(reasons[r]) for r in sorted(reasons.index)},
        "chimeric_chains": int((info["SIFTS other accessions"] != "").sum()) if len(info) else 0,
        "af_found": int(info["AlphaFold"].sum()) if len(info) else 0,
        "af_missing": sum(af_state.get(a, ("",))[0] == "missing" for a in accessions),
        "cif_unverified": unverified,
        "annotation_mismatch": len(mismatches),
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
