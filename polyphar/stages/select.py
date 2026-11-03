"""Stage 06: panel selection.

Inputs: 04_pockets/DrugScore.csv, 04_pockets/annotation/pocket_annotation.csv (the
annotate stage), 01_uniprot/tm.tsv (fetched on demand), 02_structures/Protein_info.csv
(for receptor paths).

The panel is defined by explicit pocket filters and nothing else — see `polyphar.panel`
for the filters and why each exists. Pocket similarity is **not** used here. It was,
under a hub-and-neighbours rule that took the highest-degree node of the similarity graph
plus its neighbours; that made the target count a property of how many structures
happened to be processed rather than of biology, and anchored the whole selection on one
pocket chosen by network topology. Similarity is now applied post hoc to the docking
results instead, which removes the circularity and needs no thresholds.

Outputs in 06_targets/: final_targets.csv (unchanged columns, so `dock` is unaffected),
panel_funnel.tsv, failed.tsv, manifest.json (written last; a rerun on the same inputs and
options changes nothing).

Pocket_Path and Receptor_Path are relative to the run folder. Receptor_Path comes from
Protein_info.csv: 03_pdbfixer/<UniProt>/<file name> for a PDB chain (looked up by UniProt
ID, PDB and chain, the chain case-insensitively because the Protein field is upper-cased,
e.g. 2MJW_B for author chain `b`; Chain then gives the real author chain), and the
AlphaFold row's File as listed. A target with no row there means the run folder is
inconsistent and is an error, as is one matching chains that differ only by case.

`Similarity` and `Is_Hub` are retained for schema compatibility only and are always empty
and False; `dock` reads neither.
"""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path, PurePosixPath

from polyphar import panel as panel_mod
from polyphar import tmfetch
from polyphar.common import WorkUnits, write_manifest

STAGE = "select"
# bumped whenever the filter SET changes, so a run folder records which rules made it
PANEL_VERSION = 4
DEFAULT_EXCLUDE = ("O15350",)
EXCLUDE_REASONS = {"O15350": "p73: problematic for docking (legacy decision)"}
OUTPUTS = ("final_targets.csv", "eligible_pockets.csv", "panel_funnel.tsv",
           "failed.tsv")
ELIGIBLE_COLUMNS = ["UniProt_ID", "Gene_Name", "Protein", "Pocket", "Pocket Filename",
                    "Pocket_Path", "AlphaFold", "drugScore", "radius_of_gyration"]
TARGET_COLUMNS = ["Node", "UniProt_ID", "Gene_Name", "Organism", "Structure_ID",
                  "Similarity", "Protein", "Chain", "Pocket", "Pocket_Path",
                  "Receptor_Path", "Is_Hub"]
AF_CHAIN = "A"   # Chain ID of AlphaFold rows in Protein_info.csv (structures stage)
NO_CHAIN = "-"   # a PDB-style Protein field without a chain


def add_arguments(parser) -> None:
    p = panel_mod
    parser.add_argument("--volume-min", type=float, default=p.VOLUME_MIN, metavar="V")
    parser.add_argument("--volume-max", type=float, default=p.VOLUME_MAX, metavar="V")
    parser.add_argument("--drugscore-min", type=float, default=p.DRUGSCORE_MIN, metavar="S")
    parser.add_argument("--plddt-min", type=float, default=p.PLDDT_MIN, metavar="P",
                        help=f"AlphaFold pocket pLDDT floor (default: {p.PLDDT_MIN}); a "
                             f"sanity check only, it does not detect oversized cavities")
    parser.add_argument("--tm-policy", choices=p.TM_POLICIES, default="C",
                        help="A: drop multi-pass TM proteins; B: drop any TM protein; "
                             "C (default): A plus drop pockets sitting in a TM segment")
    parser.add_argument("--tm-overlap", type=float, default=p.TM_OVERLAP_MAX, metavar="F")
    parser.add_argument("--no-fetch-tm", action="store_true",
                        help="do not contact UniProt for missing transmembrane data; "
                             "fail instead, so a run can be reproduced offline")
    parser.add_argument("--per-protein", type=int, default=p.PER_PROTEIN, metavar="N",
                        help=f"pockets docked per protein, best drugScore first "
                             f"(default: {p.PER_PROTEIN}); picking one is close to a "
                             f"coin flip against known ligand sites")
    parser.add_argument("--exclude", nargs="*", default=list(DEFAULT_EXCLUDE),
                        metavar="UNIPROT",
                        help=f"accessions to keep out of the panel "
                             f"(default: {' '.join(DEFAULT_EXCLUDE)}; give none for no "
                             f"exclusions)")


# --- inputs ---------------------------------------------------------------------------

def _options(args) -> panel_mod.Options:
    return panel_mod.Options(
        drugscore_min=getattr(args, "drugscore_min", panel_mod.DRUGSCORE_MIN),
        volume_min=getattr(args, "volume_min", panel_mod.VOLUME_MIN),
        volume_max=getattr(args, "volume_max", panel_mod.VOLUME_MAX),
        plddt_min=getattr(args, "plddt_min", panel_mod.PLDDT_MIN),
        tm_policy=getattr(args, "tm_policy", "C"),
        tm_overlap_max=getattr(args, "tm_overlap", panel_mod.TM_OVERLAP_MAX),
        per_protein=getattr(args, "per_protein", panel_mod.PER_PROTEIN),
    )


def _transmembrane(layout, accessions, no_fetch: bool) -> dict:
    """TM segments for every accession, from 01_uniprot/tm.tsv.

    The side file is fetched on demand because `uniprot` cannot be re-run safely: doing
    so refreshes its PDB cross-references, which cascades into new chains, new models
    and new DoGSiteScorer work. With --no-fetch-tm the stage stays offline and fails
    instead, so a run can be reproduced without network access.
    """
    path = Path(layout.uniprot_csv).parent / tmfetch.TM_TSV
    have = tmfetch.read_tm_tsv(path) if path.is_file() else {}
    missing = sorted(set(accessions) - set(have))
    if not missing:
        return have
    if no_fetch:
        raise FileNotFoundError(
            f"{path} is missing transmembrane data for {len(missing)} accession(s) "
            f"(e.g. {', '.join(missing[:5])}) and --no-fetch-tm was given. Drop the flag "
            f"to fetch them, or supply the file.")
    return tmfetch.ensure(layout, accessions)


def _exclusions(args) -> dict:
    ids = []
    for item in getattr(args, "exclude", list(DEFAULT_EXCLUDE)) or []:
        ids.extend(x.strip() for x in str(item).split(",") if x.strip())
    return {u: EXCLUDE_REASONS.get(u, "excluded with --exclude") for u in sorted(set(ids))}


def annotation_path(layout) -> Path:
    return layout.pockets / "annotation" / "pocket_annotation.csv"


def _read_csv(path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load(layout):
    """(selected DrugScore rows, annotations by pocket filename, parent volumes)."""
    if not layout.drugscore_csv.is_file():
        raise FileNotFoundError(f"{layout.drugscore_csv} not found: run `pockets` first")
    ann_path = annotation_path(layout)
    if not ann_path.is_file():
        raise FileNotFoundError(f"{ann_path} not found: run `annotate` first")

    all_rows = _read_csv(layout.drugscore_csv)
    rows = [r for r in all_rows
            if str(r.get("Selected", "")).strip().lower() in ("true", "1", "yes")]
    parent_volumes = {}
    for r in all_rows:
        try:
            parent_volumes[(r["UniProt ID"], r["Protein"], r["Pocket"])] = float(r["volume"])
        except (KeyError, TypeError, ValueError):
            continue

    from polyphar.stages.annotate import decode_residues
    annotations = {}
    for a in _read_csv(ann_path):
        a["_residues"] = decode_residues(a.get("residues", ""))
        # keyed on the PAIR, not the filename: a chimeric or shared chain can give the
        # same pocket file name under two accessions
        annotations[(a["UniProt ID"], a["Pocket Filename"])] = a
    return rows, annotations, parent_volumes


# --- receptor paths ---------------------------------------------------------------------

def split_entry_name(entry: str) -> tuple[str, str]:
    """`S100P_HUMAN` -> (`S100P`, `HUMAN`). Split at the last underscore; an entry name
    without one is the gene name."""
    gene, _, organism = entry.rpartition("_")
    return (gene, organism) if gene else (entry, "")


def split_protein(protein: str) -> tuple[str, str]:
    """`2MJW_B` -> (`2MJW`, `B`), `AF-P06493-F1` -> (`AF-P06493-F1`, `A`). PDB IDs have
    no underscore, so the first one ends the ID."""
    if protein.startswith("AF-"):
        return protein, AF_CHAIN
    pdb, _, chain = protein.partition("_")
    return pdb, chain or NO_CHAIN


def _receptor_index(protein_info):
    """PDB rows keyed by (UniProt ID, PDB, chain upper-cased) -> [(real chain, File)],
    and AlphaFold rows by UniProt ID. The Protein field upper-cases the chain, so the
    lookup ignores case; a key with several real chains (`a` and `A`) is ambiguous and
    is never resolved by guessing."""
    pdb, af = {}, {}
    for r in protein_info:
        acc, p, chain = r["UniProt ID"], str(r["PDB"]), str(r["Chain ID"])
        if str(r["AlphaFold"]).strip().lower() == "true":
            af.setdefault(acc, r["File"])
        else:
            entries = pdb.setdefault((acc, p.upper(), chain.upper()), [])
            if all(c != chain for c, _ in entries):
                entries.append((chain, r["File"]))
    return pdb, af


def _receptor_paths(panel, protein_info):
    """(Receptor_Path by accession, real author chain by accession)."""
    pdb, af = _receptor_index(protein_info)
    paths, chains, missing, ambiguous = {}, {}, [], []
    for acc, picks in panel.items():
        row, _ann = picks[0]
        structure_id, chain = split_protein(str(row.get("Protein", "")))
        if structure_id.startswith("AF-"):
            file = af.get(acc)
            if file is None:
                missing.append(f"{acc} AlphaFold model {structure_id}")
                continue
            paths[acc] = str(PurePosixPath(file))
            chains[acc] = AF_CHAIN
        else:
            entries = pdb.get((acc, structure_id.upper(), chain.upper()), [])
            if not entries:
                missing.append(f"{acc} PDB {structure_id} chain {chain}")
                continue
            if len(entries) > 1:
                ambiguous.append(f"{acc} PDB {structure_id} chains "
                                 + "/".join(c for c, _ in entries))
                continue
            real_chain, file = entries[0]
            paths[acc] = f"03_pdbfixer/{acc}/{PurePosixPath(file).name}"
            chains[acc] = real_chain
    if ambiguous:
        raise ValueError(
            f"{len(ambiguous)} target(s) are ambiguous: Protein_info.csv has chains that "
            f"differ only by case and the Protein field cannot tell them apart: "
            + "; ".join(ambiguous[:20]))
    if missing:
        raise ValueError(
            f"{len(missing)} target(s) have no row in Protein_info.csv (the run folder "
            f"is inconsistent): " + "; ".join(missing[:20])
            + ("; ..." if len(missing) > 20 else ""))
    return paths, chains


# --- outputs -----------------------------------------------------------------------------

def target_rows(panel, paths, chains) -> list[dict]:
    """final_targets.csv rows, sorted by accession so the output is reproducible."""
    out = []
    for acc in sorted(panel):
      for row, _ann in panel[acc]:
        entry = str(row.get("Gene Name", ""))
        gene, organism = split_entry_name(entry)
        structure_id, _ = split_protein(str(row.get("Protein", "")))
        pocket = str(row.get("Pocket", ""))
        out.append({
            "Node": f'{acc}_{entry}_{row.get("Protein", "")}_{pocket}',
            "UniProt_ID": acc,
            "Gene_Name": gene,
            "Organism": organism,
            "Structure_ID": structure_id,
            "Similarity": "",          # retained for schema compatibility only
            "Protein": row.get("Protein", ""),
            "Chain": chains.get(acc, ""),
            "Pocket": pocket,
            "Pocket_Path": f'04_pockets/{acc}/{row.get("Pocket Filename", "")}',
            "Receptor_Path": paths.get(acc, ""),
            "Is_Hub": False,           # retained for schema compatibility only
        })
    return out


def eligible_rows(eligible, excluded) -> list[dict]:
    """Every eligible pocket, not just the panel's one-per-protein pick. `similarity`
    runs on these: the post hoc test needs pairwise similarity among panel pockets,
    and this wider set also leaves a descriptive network. Sorted for reproducibility."""
    out = []
    for row, ann in eligible:
        acc = row["UniProt ID"]
        if acc in excluded:
            continue
        gene, _organism = split_entry_name(str(row.get("Gene Name", "")))
        fname = row.get("Pocket Filename", "")
        out.append({"UniProt_ID": acc, "Gene_Name": gene,
                    "Protein": row.get("Protein", ""), "Pocket": row.get("Pocket", ""),
                    "Pocket Filename": fname,
                    "Pocket_Path": f"04_pockets/{acc}/{fname}",
                    "AlphaFold": ann.get("AlphaFold", ""),
                    "drugScore": row.get("drugScore", ""),
                    "radius_of_gyration": ann.get("radius_of_gyration", "")})
    return sorted(out, key=lambda r: (r["UniProt_ID"], r["Pocket Filename"]))


def _atomic_rows(path: Path, columns, rows) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=columns)
            w.writeheader()
            w.writerows(rows)
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def units(layout, args) -> WorkUnits:
    return WorkUnits(1, "panel", [1.0])


def run(layout, args, progress) -> dict:
    opts = _options(args)
    excluded = _exclusions(args)
    d = layout.targets
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")

    rows, annotations, parent_volumes = _load(layout)
    accessions = sorted({r["UniProt ID"] for r in rows})
    tm = _transmembrane(layout, accessions, getattr(args, "no_fetch_tm", False))

    pan, eligible, funnel, rejections = panel_mod.select_panel(
        rows, annotations, tm, parent_volumes, opts)

    dropped = {a: excluded[a] for a in list(pan) if a in excluded}
    for acc in dropped:
        pan.pop(acc, None)

    protein_info = _read_csv(layout.protein_info)
    paths, chains = _receptor_paths(pan, protein_info)
    _atomic_rows(d / "final_targets.csv", TARGET_COLUMNS,
                 target_rows(pan, paths, chains))
    _atomic_rows(d / "eligible_pockets.csv", ELIGIBLE_COLUMNS,
                 eligible_rows(eligible, excluded))
    _atomic_rows(d / "panel_funnel.tsv", ["step", "pockets_remaining"],
                 [{"step": s, "pockets_remaining": n} for s, n in funnel])

    from polyphar.common import append_failure
    for item, step in sorted(rejections):
        append_failure(d, item, STAGE, panel_mod.STEP_REASON.get(step, step))
    for acc, why in sorted(dropped.items()):
        append_failure(d, acc, STAGE, why)
    progress.advance(1, size=1.0)

    af = sum(1 for picks in pan.values()
             if str(picks[0][1].get("AlphaFold", "")).strip().lower() == "true")
    targets = sum(len(v) for v in pan.values())
    # eligible_pockets must match eligible_pockets.csv, which drops excluded accessions
    eligible_kept = sum(1 for r, _a in eligible if r["UniProt ID"] not in excluded)
    counts = {"selected_pockets": len(rows), "eligible_pockets": eligible_kept,
              "panel": len(pan), "targets": targets, "alphafold": af,
              "experimental": len(pan) - af, "excluded": len(dropped),
              "rejected_pockets": len(rejections)}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE,
                   args={**opts.as_dict(), "exclude": sorted(excluded)},
                   inputs=[layout.drugscore_csv, annotation_path(layout),
                           layout.protein_info],
                   versions={"panel": "volume-drugscore", "panel_version": PANEL_VERSION}, counts=counts, timing={},
                   extra={"funnel": dict(funnel),
                          "outputs": {n: (d / n).stat().st_size for n in OUTPUTS}})
    return counts
