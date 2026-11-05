"""Stage 01: fetch UniProt records for 00_ids/ids.txt via the UniProt REST API.

Outputs in 01_uniprot/:
  uniprot.csv     kept rows (reviewed, human), one per primary accession
  not_found.txt   requested accessions UniProt does not return at all
  excluded.tsv    accession, reason  (reason: unreviewed | non_human | same workbook row as ...)
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
URL = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,gene_primary,reviewed,organism_id,xref_pdb"
BATCH = 100
PAGE_SIZE = 500       # rows per page; one BATCH is one page, so pagination rarely kicks in
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
        # `search` answers for an INACTIVE entry (deleted, or demerged into others) with a
        # stub row: accession and entry name only, every other field empty. int64 raised
        # ValueError on the empty string and took the whole stage down -- 14 of Cluster1's
        # 131 accessions are deleted TrEMBL entries. They parse as NA here and classify()
        # reports them as inactive.
        df["Organism (ID)"] = pd.to_numeric(df["Organism (ID)"], errors="coerce").astype("Int64")
    return df


# --- fetching ---------------------------------------------------------------

def _get(session, query, fields, what):
    """GET `query`, following UniProt's `next` page links; any failure becomes a
    RuntimeError naming `what`. Returns (concatenated TSV, release).

    The paginated `search` endpoint is used rather than `stream`: on 2026-11-02 `stream`
    answered these multi-accession OR-queries with intermittent HTTP 500s (a 100-accession
    batch failed from two different networks, while `search` returned all rows), and a 500
    is retried by http_session until the retries run out.
    """
    texts, release, seen = [], None, set()
    url, params = URL, {"query": query, "fields": fields, "format": "tsv", "size": PAGE_SIZE}
    while url and url not in seen:
        seen.add(url)
        try:
            resp = session.get(url, params=params, timeout=TIMEOUT)
            resp.raise_for_status()
        except requests.exceptions.HTTPError as exc:
            body = (getattr(exc.response, "text", "") or "").strip().replace("\n", " ")[:300]
            raise RuntimeError(f"UniProt request failed for {what}: {exc}. {body}") from exc
        except requests.exceptions.RequestException as exc:
            raise RuntimeError(f"UniProt request failed for {what}: {exc!r}") from exc
        texts.append(resp.text)
        if release is None:
            release = resp.headers.get("x-uniprot-release", "")
        url = (resp.links.get("next") or {}).get("url")
        params = None      # a next-link already carries query, fields, format and cursor
    return _join_tsv(texts), release or ""


def _join_tsv(texts):
    """One TSV from consecutive pages: the first page whole, later pages without
    their repeated header line (the header is whatever the first page's first line is,
    so this holds for any `fields`)."""
    if not texts:
        return ""
    header = texts[0].split("\n", 1)[0]
    out = []
    for i, text in enumerate(texts):
        if i and text.split("\n", 1)[0] == header:
            _, _, text = text.partition("\n")
        out.append(text if text.endswith("\n") or not text else text + "\n")
    return "".join(out)


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

    Returns dict(kept: DataFrame, not_found: [acc], inactive: [acc], unreviewed: [acc],
    non_human: [acc], renamed: {requested: primary}). Requested IDs absent from
    df['Entry'] are re-queried one at a time with `requery(acc) -> primary|None`;
    without `requery` they are all not_found.

    An inactive entry (deleted, or demerged into others) comes back as a stub row with no
    organism, so it is reported as inactive rather than re-queried: UniProt does return it
    for its own accession, so a re-query would "rename" it to itself."""
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
    live = df["Organism (ID)"].notna()
    reviewed = df["Reviewed"] == "reviewed"
    human = (df["Organism (ID)"] == HUMAN_TAXID).fillna(False)
    inactive = df.loc[~live, "Entry"].tolist()
    unreviewed = df.loc[live & ~reviewed, "Entry"].tolist()
    non_human = df.loc[live & reviewed & ~human, "Entry"].tolist()
    kept = df[live & reviewed & human].reset_index(drop=True)
    return {"kept": kept, "not_found": not_found, "inactive": inactive,
            "unreviewed": unreviewed, "non_human": non_human, "renamed": renamed}


# --- one protein per workbook row ----------------------------------------------

def read_rows(path):
    """Parse 00_ids/rows.tsv -> [(row, gene, [ids], {solo ids})]."""
    out = []
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:
        if not line.strip():
            continue
        f = line.split("\t")
        f += [""] * (4 - len(f))
        out.append((f[0], f[1], [i for i in f[2].split(";") if i], {i for i in f[3].split(";") if i}))
    return out


def collapse_rows(rows, kept, renamed):
    """Resolve multi-accession workbook rows to one protein each.

    `rows`: read_rows() output; `kept`: reviewed-human DataFrame; `renamed`: {requested: primary}.
    A row with several kept accessions keeps the one whose `Gene Names (primary)` equals the
    row's gene (case-insensitive); none or several matching raises ValueError. Returns
    {dropped accession: reason}. An accession also kept by another row (as that row's
    winner, or as a single-accession row) is not dropped."""
    genes = dict(zip(kept["Entry"], kept["Gene Names (primary)"]))
    winners, drops, protected = {}, {}, set()
    for row, gene, ids, solo in rows:
        gene = gene.strip()
        protected |= {renamed.get(i, i) for i in solo}
        cands = list(dict.fromkeys(renamed.get(i, i) for i in ids if renamed.get(i, i) in genes))
        if len(cands) < 2:
            if cands:
                winners[row] = cands[0]
            continue
        want = gene.strip().casefold()
        match = [a for a in cands if str(genes[a]).strip().casefold() == want] if want else []
        if len(match) != 1:
            listing = ", ".join(f"{a} ({genes[a]})" for a in cands)
            raise ValueError(
                f"workbook row {row} (gene {gene!r}): cannot pick one protein, "
                f"{len(match)} of its accessions match the gene; candidates: {listing}. "
                "Fix the workbook row or the gene column and rerun `polyphar ids`.")
        winners[row] = match[0]
        for a in cands:
            if a != match[0]:
                drops.setdefault(a, (match[0], gene))
    protected |= set(winners.values())
    return {a: f"same workbook row as {w} (gene {g})"
            for a, (w, g) in drops.items() if a not in protected}


# --- stage ------------------------------------------------------------------

def _read_ids(layout):
    if not layout.ids_txt.is_file():
        raise FileNotFoundError(
            f"{layout.ids_txt} not found; run `polyphar ids` first")
    ids = [l.strip() for l in layout.ids_txt.read_text(encoding="utf-8").splitlines()]
    return list(dict.fromkeys(i for i in ids if i))


def _input_files(layout):
    return [layout.ids_txt] + ([layout.rows_tsv] if layout.rows_tsv.is_file() else [])


def _previous_counts(layout):
    """Counts of a finished run on the same ids.txt (and rows.tsv), else None."""
    manifest = layout.stage_dir(STAGE_DIR) / "manifest.json"
    if not (manifest.is_file() and layout.uniprot_csv.is_file() and layout.ids_txt.is_file()):
        return None
    try:
        m = json.loads(manifest.read_text(encoding="utf-8"))
        now = [(p.name, sha256_file(p)) for p in _input_files(layout)]
        if [(Path(i["path"]).name, i["sha256"]) for i in m["inputs"]] == now:
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
    collapsed = {}
    if layout.rows_tsv.is_file():
        collapsed = collapse_rows(read_rows(layout.rows_tsv), res["kept"], res["renamed"])
        res["kept"] = res["kept"][~res["kept"]["Entry"].isin(collapsed)].reset_index(drop=True)

    d = layout.stage_dir(STAGE_DIR)
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    _atomic_write(d / "not_found.txt", "".join(f"{a}\n" for a in res["not_found"]))
    _atomic_write(d / "renamed.tsv", "requested\tprimary\n" + "".join(
        f"{k}\t{v}\n" for k, v in res["renamed"].items()))
    _atomic_write(d / "excluded.tsv", "accession\treason\n"
                  + "".join(f"{a}\tinactive\n" for a in res["inactive"])
                  + "".join(f"{a}\tunreviewed\n" for a in res["unreviewed"])
                  + "".join(f"{a}\tnon_human\n" for a in res["non_human"])
                  + "".join(f"{a}\t{why}\n" for a, why in collapsed.items()))
    _atomic_write(d / "failed.tsv", "item\tstage\treason\n")
    tmp = layout.uniprot_csv.with_name("uniprot.csv.tmp")
    res["kept"].to_csv(tmp, index=False)
    tmp.replace(layout.uniprot_csv)

    counts = {"requested": len(ids), "kept": len(res["kept"]),
              "not_found": len(res["not_found"]), "inactive": len(res["inactive"]),
              "unreviewed": len(res["unreviewed"]),
              "non_human": len(res["non_human"]), "renamed": len(res["renamed"]),
              "failed": 0}
    if layout.rows_tsv.is_file():
        counts["row_collapsed"] = len(collapsed)
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, "uniprot", args={"batch": BATCH}, inputs=_input_files(layout),
                   versions={"uniprot": release}, counts=counts, timing={})
    return counts
