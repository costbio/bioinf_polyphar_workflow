"""Stage 00: build 00_ids/ids.txt from an xlsx workbook tab or a plain-text list."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from polyphar.common import WorkUnits, write_manifest

STAGE_DIR = "00_ids"
_SPLIT = re.compile(r"[;,\s]+")


def add_arguments(parser) -> None:
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--xlsx", help="workbook with the accession column (relative to --run-dir)")
    src.add_argument("--txt", help="text file, one or more accessions per line (relative to --run-dir)")
    parser.add_argument("--sheet", help="workbook tab (with --xlsx)")
    parser.add_argument("--id-column", help="column holding accessions (with --xlsx)")
    parser.add_argument("--reviewed-column",
                        help="optional column; only rows equal to 'reviewed' are kept (with --xlsx)")
    parser.add_argument("--gene-column",
                        help="optional gene-symbol column (with --xlsx). Writes 00_ids/rows.tsv for "
                             "rows whose cell holds several accessions, so the uniprot stage can keep "
                             "one protein per row (the one whose primary gene name equals the symbol)")


def _tokens(text) -> list[str]:
    if text is None or (not isinstance(text, str) and pd.isna(text)):
        return []
    return [t for t in _SPLIT.split(str(text).strip()) if t]


def _check_column(df, col, sheet, what):
    if col not in df.columns:
        raise ValueError(
            f"{what} {col!r} not found in sheet {sheet!r}; available columns: {list(df.columns)}")


def _collect(cells, reviewed_flags=None, genes=None):
    """Split, strip, de-duplicate (first-seen order). Returns (ids, stats).

    With `genes` (one per cell), stats["multi_rows"] lists (row, gene, [ids]) for every kept
    row with more than one ID (row = 1-based spreadsheet row, header is row 1), and
    stats["single_ids"] holds the IDs that are the sole ID of some kept row."""
    ids, seen = [], set()
    stats = {"rows_read": 0, "duplicates_dropped": 0, "unreviewed_dropped": 0,
             "multi_rows": [], "single_ids": set()}
    for i, cell in enumerate(cells):
        stats["rows_read"] += 1
        if reviewed_flags is not None:
            flag = reviewed_flags[i]
            if not isinstance(flag, str) or flag.strip().lower() != "reviewed":
                stats["unreviewed_dropped"] += 1
                continue
        row_ids = list(dict.fromkeys(_tokens(cell)))
        if genes is not None:
            gene = genes[i]
            gene = "" if gene is None or (not isinstance(gene, str) and pd.isna(gene)) else str(gene).strip()
            if len(row_ids) > 1:
                stats["multi_rows"].append((i + 2, gene, row_ids))
            elif row_ids:
                stats["single_ids"].add(row_ids[0])
        for tok in _tokens(cell):
            if tok in seen:
                stats["duplicates_dropped"] += 1
            else:
                seen.add(tok)
                ids.append(tok)
    stats["ids_kept"] = len(ids)
    return ids, stats


def _read_xlsx(xlsx, sheet, id_column, reviewed_column, gene_column=None):
    xlsx = Path(xlsx)
    with pd.ExcelFile(xlsx) as book:   # closed even when the sheet is missing
        if sheet not in book.sheet_names:
            raise ValueError(
                f"sheet {sheet!r} not found in {xlsx}; available sheets: {book.sheet_names}")
        df = pd.read_excel(book, sheet_name=sheet, dtype=str)
    _check_column(df, id_column, sheet, "id column")
    flags = None
    if reviewed_column:
        _check_column(df, reviewed_column, sheet, "reviewed column")
        flags = df[reviewed_column].tolist()
    genes = None
    if gene_column:
        _check_column(df, gene_column, sheet, "gene column")
        genes = df[gene_column].tolist()
    return _collect(df[id_column].tolist(), flags, genes)


def build_id_list(xlsx, sheet, id_column, reviewed_column=None) -> list[str]:
    return _read_xlsx(xlsx, sheet, id_column, reviewed_column)[0]


def _read_txt(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    return _collect(lines)


def _resolve(layout, p) -> Path:
    p = Path(p)
    return p if p.is_absolute() else layout.run_dir / p


def _validate(args):
    if args.xlsx and not (args.sheet and args.id_column):
        raise SystemExit("polyphar ids: --xlsx requires --sheet and --id-column")
    if args.txt and (args.sheet or args.id_column or args.reviewed_column or args.gene_column):
        raise SystemExit(
            "polyphar ids: --sheet/--id-column/--reviewed-column/--gene-column apply only to --xlsx")


def units(layout, args) -> WorkUnits:
    _validate(args)
    return WorkUnits(1, "id list")


def run(layout, args, progress) -> dict:
    _validate(args)
    if args.xlsx:
        src = _resolve(layout, args.xlsx)
        ids, stats = _read_xlsx(src, args.sheet, args.id_column, args.reviewed_column,
                                args.gene_column)
    else:
        src = _resolve(layout, args.txt)
        ids, stats = _read_txt(src)
    out = layout.ids_txt
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".txt.tmp")
    tmp.write_text("".join(f"{i}\n" for i in ids), encoding="utf-8")
    tmp.replace(out)
    rows_tsv = layout.rows_tsv
    if args.gene_column:
        text = "row\tgene\tids\tsolo\n" + "".join(
            f"{r}\t{g}\t{';'.join(rid)}\t{';'.join(i for i in rid if i in stats['single_ids'])}\n"
            for r, g, rid in stats["multi_rows"])
        tmp = rows_tsv.with_suffix(".tsv.tmp")
        tmp.write_text(text, encoding="utf-8")
        tmp.replace(rows_tsv)
    else:
        rows_tsv.unlink(missing_ok=True)   # a stale file would change the uniprot stage
    counts = {k: stats[k] for k in ("rows_read", "ids_kept", "duplicates_dropped",
                                    "unreviewed_dropped")}
    counts["failed"] = 0
    write_manifest(
        layout.stage_dir(STAGE_DIR), "ids",
        args={"xlsx": args.xlsx, "txt": args.txt, "sheet": args.sheet,
              "id_column": args.id_column, "reviewed_column": args.reviewed_column,
              "gene_column": args.gene_column},
        inputs=[src], versions={"pandas": pd.__version__}, counts=counts,
        timing={})
    progress.advance(1)
    return counts
