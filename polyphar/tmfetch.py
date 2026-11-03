"""Transmembrane segment positions, fetched into a side file.

Why a side file and not new columns in `01_uniprot/uniprot.csv`: re-running the
`uniprot` stage refreshes its `PDB` cross-reference column from a current UniProt
release, and `structures` takes its whole work list from that column and rebuilds
`Protein_info.csv` on every run. New cross-references mean new chains and new models
that `pockets` would have to score with DoGSiteScorer — whose licence is host-bound and
expires 2026-11-16. Inside that window the cascade is unrecoverable, so the TM data is
fetched on its own and written beside uniprot.csv without touching it.

Output: 01_uniprot/tm.tsv — accession, tm_count, segments ("25..45;60..80").
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

from polyphar.common import http_session

TM_TSV = "tm.tsv"
BATCH = 80
FIELDS = "accession,ft_transmem"
URL = "https://rest.uniprot.org/uniprotkb/search"
COLUMNS = ["accession", "tm_count", "segments"]


def parse_transmem(text: str) -> list[tuple[int, int]]:
    """UniProt's ft_transmem text -> [(start, end)]. Only TRANSMEM features count:
    intramembrane and beta-barrel annotations are deliberately excluded, because this
    is the definition the panel's TM policy was measured against."""
    import re
    return [(int(a), int(b)) for a, b in
            re.findall(r"TRANSMEM (\d+)\.\.(\d+)", text or "")]


def fetch(accessions, session=None) -> dict[str, list[tuple[int, int]]]:
    """{accession: [(start, end)]} from UniProt, in batches."""
    session = session or http_session()
    out: dict[str, list[tuple[int, int]]] = {}
    missed: set[str] = set()
    accs = sorted(set(accessions))
    for i in range(0, len(accs), BATCH):
        batch = accs[i:i + BATCH]
        query = " OR ".join(f"accession:{a}" for a in batch)
        r = session.get(URL, params={"query": query, "fields": FIELDS,
                                     "format": "tsv", "size": "500"}, timeout=120)
        r.raise_for_status()
        got = set()
        for row in csv.DictReader(r.text.splitlines(), delimiter="\t"):
            out[row["Entry"]] = parse_transmem(row.get("Transmembrane", ""))
            got.add(row["Entry"])
        missed.update(a for a in batch if a not in got)
    if missed:
        # recording a miss as "no TM segments" would pass every TM policy, so a
        # multi-pass membrane protein could enter the panel silently
        raise RuntimeError(
            f"UniProt returned no record for {len(missed)} accession(s), so their "
            f"transmembrane status is unknown: {', '.join(sorted(missed)[:10])}"
            + (" ..." if len(missed) > 10 else "")
            + ". These may be demerged, secondary or renamed accessions.")
    return out


def write_tm_tsv(path, data: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tsv.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLUMNS)
        for acc in sorted(data):
            segs = data[acc]
            w.writerow([acc, len(segs), ";".join(f"{a}..{b}" for a, b in segs)])
    os.replace(tmp, path)
    return path


def read_tm_tsv(path) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            segs = []
            for part in (row.get("segments") or "").split(";"):
                if ".." in part:
                    a, b = part.split("..", 1)
                    segs.append((int(a), int(b)))
            out[row["accession"]] = segs
    return out


def ensure(layout, accessions, session=None) -> dict:
    """Read 01_uniprot/tm.tsv, fetching any accessions it is missing. Idempotent."""
    path = Path(layout.uniprot_csv).parent / TM_TSV
    have = read_tm_tsv(path) if path.is_file() else {}
    missing = sorted(set(accessions) - set(have))
    if missing:
        have.update(fetch(missing, session=session))
        write_tm_tsv(path, have)
    return have
