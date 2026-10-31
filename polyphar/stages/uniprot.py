"""Stage 01: fetch UniProt records for 00_ids/ids.txt via the UniProt REST API.

Outputs in 01_uniprot/:
  uniprot.csv     kept rows (reviewed, human), one per primary accession
  not_found.txt   requested accessions UniProt does not return at all
  excluded.tsv    accession, reason  (reason: unreviewed | non_human)
  renamed.tsv     requested, primary (secondary accession -> current primary)
  failed.tsv      header only: a failed batch aborts the stage instead
  manifest.json   inputs (ids.txt sha256), versions {"uniprot": release}, counts

Categories are distinct: an entry that is unreviewed is `unreviewed` even if it
is also non-human. `unreviewed` / `non_human` list primary accessions (the
`Entry` UniProt returned), so a renamed ID that resolves to an excluded entry
appears in both `renamed` and the excluded list.
"""

from __future__ import annotations

import io
import json
import math
import re
from pathlib import Path

import requests

from polyphar.common import WorkUnits, http_session, sha256_file, write_manifest

STAGE_DIR = "01_uniprot"
URL = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,id,gene_primary,reviewed,organism_id,xref_pdb"
BATCH = 100
TIMEOUT = 60          # seconds per request (http_session has no default timeout)
HUMAN_TAXID = 9606
COLUMNS = ["Entry", "Entry Name", "Gene Names (primary)", "Reviewed", "Organism (ID)", "PDB"]


# --- parsing ----------------------------------------------------------------

def parse_uniprot_tsv(text):
    """Parse a UniProt TSV response; strips the trailing ';' of the PDB column."""
    import pandas as pd

    if not text.strip():
        return pd.DataFrame(columns=COLUMNS)
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    if "PDB" in df.columns:
        df["PDB"] = df["PDB"].str.rstrip(";")
    if "Organism (ID)" in df.columns:
        df["Organism (ID)"] = df["Organism (ID)"].astype("int64")
    return df


# --- fetching ---------------------------------------------------------------

def _get(session, query, fields, what):
    """One GET; any failure becomes a RuntimeError naming `what`. Returns (text, release)."""
    try:
        resp = session.get(URL, params={"query": query, "fields": fields, "format": "tsv"},
                           timeout=TIMEOUT)
        resp.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        body = (getattr(exc.response, "text", "") or "").strip().replace("\n", " ")[:300]
        raise RuntimeError(f"UniProt request failed for {what}: {exc}. {body}") from exc
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(f"UniProt request failed for {what}: {exc!r}") from exc
    return resp.text, resp.headers.get("x-uniprot-release", "")


def _query(accessions):
    return " OR ".join(f"(accession:{a})" for a in accessions)


def fetch_uniprot(accessions, session, batch=BATCH, on_batch=None):
    """Fetch records for `accessions` in batches. Returns (DataFrame, release).

    `accession:` also matches secondary accessions and returns the primary entry,
    so the frame can contain entries under a different name than requested; rows
    are de-duplicated on Entry. Any failed batch raises RuntimeError."""
    import pandas as pd

    accessions = list(accessions)
    n_batches = math.ceil(len(accessions) / batch)
    frames, release = [], None
    for b in range(n_batches):
        chunk = accessions[b * batch:(b + 1) * batch]
        what = f"batch {b + 1}/{n_batches} ({chunk[0]}..{chunk[-1]})"
        text, rel = _get(session, _query(chunk), FIELDS, what)
        if release is None:
            release = rel
        elif rel != release:
            raise RuntimeError(
                f"UniProt release changed during the fetch ({release} -> {rel}) at {what}; "
                "rerun the stage")
        frames.append(parse_uniprot_tsv(text))
        if on_batch:
            on_batch(1)
    if not frames:
        return pd.DataFrame(columns=COLUMNS), ""
    df = pd.concat(frames, ignore_index=True).drop_duplicates("Entry").reset_index(drop=True)
    return df, release or ""


def requery_primary(session, accession):
    """Primary accession UniProt returns for `accession` (a secondary/merged ID),
    or None when UniProt returns nothing."""
    text, _ = _get(session, f"(accession:{accession})", "accession", f"accession {accession}")
    df = parse_uniprot_tsv(text)
    return str(df["Entry"].iloc[0]) if len(df) else None


# --- classification -----------------------------------------------------------

def classify(requested, df, requery=None):
    """Split `requested` accessions by outcome.

    Returns dict(kept: DataFrame, not_found: [acc], unreviewed: [acc],
    non_human: [acc], renamed: {requested: primary}). Requested IDs absent from
    df['Entry'] are re-queried one at a time with `requery(acc) -> primary|None`;
    without `requery` they are all not_found."""
    returned = set(df["Entry"])
    renamed, not_found = {}, []
    for acc in dict.fromkeys(requested):
        if acc in returned:
            continue
        primary = requery(acc) if requery is not None else None
        if primary is None:
            not_found.append(acc)
        else:
            renamed[acc] = primary
    reviewed = df["Reviewed"] == "reviewed"
    human = df["Organism (ID)"] == HUMAN_TAXID
    unreviewed = df.loc[~reviewed, "Entry"].tolist()
    non_human = df.loc[reviewed & ~human, "Entry"].tolist()
    kept = df[reviewed & human].reset_index(drop=True)
    return {"kept": kept, "not_found": not_found, "unreviewed": unreviewed,
            "non_human": non_human, "renamed": renamed}


# --- stage ------------------------------------------------------------------

def _read_ids(layout):
    if not layout.ids_txt.is_file():
        raise FileNotFoundError(
            f"{layout.ids_txt} not found; run `polyphar ids` first")
    ids = [l.strip() for l in layout.ids_txt.read_text(encoding="utf-8").splitlines()]
    return list(dict.fromkeys(i for i in ids if i))


def _previous_counts(layout):
    """Counts of a finished run on the same ids.txt, else None."""
    manifest = layout.stage_dir(STAGE_DIR) / "manifest.json"
    if not (manifest.is_file() and layout.uniprot_csv.is_file() and layout.ids_txt.is_file()):
        return None
    try:
        m = json.loads(manifest.read_text(encoding="utf-8"))
        if m["inputs"][0]["sha256"] == sha256_file(layout.ids_txt):
            return m["counts"]
    except (OSError, ValueError, KeyError, IndexError, TypeError):
        pass
    return None


def units(layout, args) -> WorkUnits:
    ids = _read_ids(layout)
    if _previous_counts(layout) is not None:
        return WorkUnits(0, "batch")
    return WorkUnits(math.ceil(len(ids) / BATCH), "batch")


def _atomic_write(path: Path, text: str):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def run(layout, args, progress) -> dict:
    ids = _read_ids(layout)
    previous = _previous_counts(layout)
    if previous is not None:
        return previous
    session = http_session()
    df, release = fetch_uniprot(ids, session, BATCH, on_batch=progress.advance)
    res = classify(ids, df, requery=lambda acc: requery_primary(session, acc))

    d = layout.stage_dir(STAGE_DIR)
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    _atomic_write(d / "not_found.txt", "".join(f"{a}\n" for a in res["not_found"]))
    _atomic_write(d / "renamed.tsv", "requested\tprimary\n" + "".join(
        f"{k}\t{v}\n" for k, v in res["renamed"].items()))
    _atomic_write(d / "excluded.tsv", "accession\treason\n"
                  + "".join(f"{a}\tunreviewed\n" for a in res["unreviewed"])
                  + "".join(f"{a}\tnon_human\n" for a in res["non_human"]))
    _atomic_write(d / "failed.tsv", "item\tstage\treason\n")
    tmp = layout.uniprot_csv.with_name("uniprot.csv.tmp")
    res["kept"].to_csv(tmp, index=False)
    tmp.replace(layout.uniprot_csv)

    counts = {"requested": len(ids), "kept": len(res["kept"]),
              "not_found": len(res["not_found"]), "unreviewed": len(res["unreviewed"]),
              "non_human": len(res["non_human"]), "renamed": len(res["renamed"]),
              "failed": 0}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, "uniprot", args={"batch": BATCH}, inputs=[layout.ids_txt],
                   versions={"uniprot": release}, counts=counts, timing={})
    return counts
