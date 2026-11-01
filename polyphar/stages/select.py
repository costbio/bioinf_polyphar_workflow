"""Stage 06: target selection from the pocket-similarity network.

Input: 05_similarity/Pocketsim.csv (and 02_structures/Protein_info.csv for receptor files).

Legacy semantics (Gene_Disease_Association/target_selection_analysis.py), made order
independent:
  1. one row per protein PAIR: the max-Pmax pocket pair, then Pmax >= --pmax. The pair is
     unordered ((A,B) and (B,A) are one pair; the legacy keyed on the ordered pair, so a
     reversed duplicate could add a second, weaker edge). Ties in Pmax go to the smallest
     (Node, Node_2), so the result does not depend on the row order;
  2. a graph over pocket nodes `<UniProt>_<EntryName>_<Protein>_<Pocket>`, edges weighted
     by Pmax; every node carries its fields as attributes (nothing is parsed back out of
     the node string);
  3. nodes are ranked by degree; the hub is the first non-excluded node, and the targets
     are the hub plus its neighbours (excluded neighbours dropped; hub Similarity 1.0,
     a neighbour's Similarity is its edge Pmax).
Ties at the top degree follow --hub-tie-break: `union` (default) takes every non-excluded
top-degree node as a hub (targets = hubs + the union of their neighbours; a neighbour of
several hubs gets its maximum Pmax to any hub); `weighted` and `betweenness` pick one hub
by Weighted_Degree or Betweenness (then Node).

Outputs in 06_targets/: final_targets.csv (legacy columns Node, UniProt_ID, Gene_Name,
Organism, Structure_ID, Similarity, plus Protein, Chain, Pocket, Pocket_Path,
Receptor_Path, Is_Hub), centrality.csv, unique_pairs.csv, failed.tsv, manifest.json
(written last; a rerun on the same inputs and options changes nothing).

Pocket_Path and Receptor_Path are relative to the run folder. Receptor_Path comes from
Protein_info.csv: 03_pdbfixer/<UniProt>/<file name> for a PDB chain (looked up by UniProt
ID, PDB, chain) and the AlphaFold row's File as listed. A target without a row there means
the run folder is inconsistent and is an error.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from pathlib import Path, PurePosixPath

from polyphar.common import WorkUnits, sha256_file, write_manifest

STAGE = "select"
DEFAULT_PMAX = 0.7
DEFAULT_EXCLUDE = ("O15350",)
EXCLUDE_REASONS = {"O15350": "p73: problematic for docking (legacy decision)"}
TIE_BREAKS = ("union", "weighted", "betweenness")
OUTPUTS = ("final_targets.csv", "centrality.csv", "unique_pairs.csv", "failed.tsv")
TEXT_COLUMNS = ["UniProt ID", "Gene Name", "Protein", "Pocket", "Pocket_Path",
                "UniProt ID_2", "Gene Name_2", "Protein_2", "Pocket_2", "Pocket_Path_2"]
TARGET_COLUMNS = ["Node", "UniProt_ID", "Gene_Name", "Organism", "Structure_ID", "Similarity",
                  "Protein", "Chain", "Pocket", "Pocket_Path", "Receptor_Path", "Is_Hub"]
RANK_COLUMNS = ["Node", "UniProt_ID", "Gene_Name", "Degree", "Weighted_Degree",
                "Degree_Centrality", "Betweenness_Centrality", "Closeness_Centrality"]
NO_CHAIN = "-"   # Chain for AlphaFold models, as in Protein_info.csv


def add_arguments(parser) -> None:
    parser.add_argument("--pmax", type=float, default=DEFAULT_PMAX,
                        help=f"minimum Pmax of a protein pair (default {DEFAULT_PMAX})")
    parser.add_argument("--exclude", nargs="*", default=list(DEFAULT_EXCLUDE), metavar="UNIPROT",
                        help="UniProt IDs never used as hub or target "
                             f"(default: {' '.join(DEFAULT_EXCLUDE)}; give none to exclude nothing)")
    parser.add_argument("--hub-tie-break", dest="hub_tie_break", choices=TIE_BREAKS,
                        default="union",
                        help="how a tie at the top degree is resolved: union = every tied "
                             "non-excluded node is a hub (default); weighted / betweenness = "
                             "pick one hub")


# --- reading ---------------------------------------------------------------------------

def read_pocketsim(path):
    """Pocketsim.csv without its index column. Text columns are categoricals (the file
    can have millions of rows)."""
    import pandas as pd

    dtype = {c: "category" for c in TEXT_COLUMNS}
    dtype.update({"Pmin": "float64", "Pmax": "float64"})
    return pd.read_csv(path, usecols=lambda c: c != "" and not c.startswith("Unnamed"),
                       dtype=dtype)


def _read_protein_info(path):
    import pandas as pd

    return pd.read_csv(path, dtype=str, keep_default_na=False)


def split_entry_name(entry: str) -> tuple[str, str]:
    """`S100P_HUMAN` -> (`S100P`, `HUMAN`); `2AAA_HUMAN` -> (`2AAA`, `HUMAN`). Split at the
    last underscore; an entry name without one is the gene name."""
    gene, _, organism = entry.rpartition("_")
    return (gene, organism) if gene else (entry, "")


def split_protein(protein: str) -> tuple[str, str]:
    """The Pocketsim `Protein` field: `2MJW_B` -> (`2MJW`, `B`), `AF-P06493-F1` ->
    (`AF-P06493-F1`, `-`). PDB IDs have no underscore, so the first one ends the ID."""
    if protein.startswith("AF-"):
        return protein, NO_CHAIN
    pdb, _, chain = protein.partition("_")
    return pdb, chain or NO_CHAIN


# --- pairs and graph -------------------------------------------------------------------

def _node(d, suffix: str):
    return (d["UniProt ID" + suffix] + "_" + d["Gene Name" + suffix] + "_"
            + d["Protein" + suffix] + "_" + d["Pocket" + suffix])


def unique_protein_pairs(df, pmax: float = DEFAULT_PMAX, ordered: bool = False):
    """One row per protein pair (the max-Pmax pocket pair), then Pmax >= `pmax`; Node and
    Node_2 columns added. The pair is the unordered {UniProt ID, UniProt ID_2}; Pmax ties
    go to the smallest (Node, Node_2) of the row read in sorted order, so the result does
    not depend on the order of `df`. Same-UniProt rows are dropped.

    `ordered=True` keys on the ordered pair like the 2023 notebook; it exists only to
    explain regression differences.

    Rows below `pmax` are dropped first: the retained row of a pair is its maximum, so the
    pair survives the final filter exactly when one of its rows passes it."""
    import pandas as pd

    d = df[df["Pmax"] >= pmax].drop(columns=["Node", "Node_2"], errors="ignore")
    d = d.reset_index(drop=True)
    for col in TEXT_COLUMNS:
        if col in d.columns:
            d[col] = d[col].astype(str)
    d = d[d["UniProt ID"] != d["UniProt ID_2"]].reset_index(drop=True)
    d["Node"], d["Node_2"] = _node(d, ""), _node(d, "_2")
    a, b = d["UniProt ID"], d["UniProt ID_2"]
    swap = (a > b) & (not ordered)
    d["_ka"], d["_kb"] = a.where(~swap, b), b.where(~swap, a)
    na, nb = d["Node"], d["Node_2"]
    swap = na > nb
    d["_lo"], d["_hi"] = na.where(~swap, nb), nb.where(~swap, na)
    d = d.sort_values(["Pmax", "_lo", "_hi", "Node", "Node_2"],
                      ascending=[False, True, True, True, True], kind="mergesort")
    d = d.drop_duplicates(["_ka", "_kb"]).drop(columns=["_ka", "_kb", "_lo", "_hi"])
    return d.reset_index(drop=True)


def build_graph(df):
    """Pocket-node graph; edge weight = Pmax. Node IDs `<UniProt>_<EntryName>_<Protein>_
    <Pocket>`; attributes uniprot, entry_name, protein, pocket, pocket_path (plus
    structure_id and chain, split from the Protein field)."""
    import networkx as nx

    g = nx.Graph()
    for suffix in ("", "_2"):
        node = _node(df, suffix) if "Node" not in df.columns else df["Node" + suffix]
        for n, u, e, p, pk, path in zip(node, df["UniProt ID" + suffix], df["Gene Name" + suffix],
                                        df["Protein" + suffix], df["Pocket" + suffix],
                                        df["Pocket_Path" + suffix]):
            if n not in g:
                structure_id, chain = split_protein(str(p))
                g.add_node(n, uniprot=str(u), entry_name=str(e), protein=str(p), pocket=str(pk),
                           pocket_path=str(path), structure_id=structure_id, chain=chain)
    node1 = df["Node"] if "Node" in df.columns else _node(df, "")
    node2 = df["Node_2"] if "Node_2" in df.columns else _node(df, "_2")
    for n1, n2, w in zip(node1, node2, df["Pmax"]):
        g.add_edge(n1, n2, weight=float(w))
    return g


# --- ranking and hub choice ------------------------------------------------------------

def rank_hubs(G, tie_break: str = "union"):
    """Centrality table for every node, best hub first: Degree desc, then Weighted_Degree
    desc (`weighted` and `union`) or Betweenness desc (`betweenness`), then Node asc. The
    float keys are rounded so summation order cannot change a tie."""
    import networkx as nx
    import pandas as pd

    if tie_break not in TIE_BREAKS:
        raise ValueError(f"unknown tie break {tie_break!r}; choose one of {', '.join(TIE_BREAKS)}")
    dc = nx.degree_centrality(G)
    bc = nx.betweenness_centrality(G)
    cc = nx.closeness_centrality(G)
    rows = []
    for n, attrs in G.nodes(data=True):
        rows.append({
            "Node": n, "UniProt_ID": attrs["uniprot"],
            "Gene_Name": split_entry_name(attrs["entry_name"])[0],
            "Degree": G.degree(n),
            "Weighted_Degree": math.fsum(w for _, _, w in G.edges(n, data="weight")),
            "Degree_Centrality": dc[n], "Betweenness_Centrality": bc[n],
            "Closeness_Centrality": cc[n]})
    df = pd.DataFrame(rows, columns=RANK_COLUMNS)
    second = "Betweenness_Centrality" if tie_break == "betweenness" else "Weighted_Degree"
    df["_key"] = df[second].round(10)
    df = df.sort_values(["Degree", "_key", "Node"], ascending=[False, False, True],
                        kind="mergesort").drop(columns="_key")
    return df.reset_index(drop=True)


def tied_top(ranking, exclude=()) -> list[str]:
    """The non-excluded nodes at the top degree, in ranking order (one node when there is
    no tie, none when every node is excluded)."""
    r = ranking[~ranking["UniProt_ID"].isin(set(exclude))]
    if r.empty:
        return []
    return r.loc[r["Degree"] == r["Degree"].max(), "Node"].tolist()


def _receptor_index(protein_info):
    pdb, af = {}, {}
    for acc, p, chain, is_af, file in zip(protein_info["UniProt ID"], protein_info["PDB"],
                                          protein_info["Chain ID"], protein_info["AlphaFold"],
                                          protein_info["File"]):
        if str(is_af).strip().lower() == "true":
            af.setdefault(acc, file)
        else:
            pdb.setdefault((acc, str(p).upper(), chain), file)
    return pdb, af


def _receptor_paths(G, nodes, protein_info) -> dict[str, str]:
    pdb, af = _receptor_index(protein_info)
    paths, missing = {}, []
    for n in nodes:
        a = G.nodes[n]
        if a["structure_id"].startswith("AF-"):
            file = af.get(a["uniprot"])
            if file is None:
                missing.append(f"{a['uniprot']} AlphaFold model {a['structure_id']}")
                continue
            paths[n] = str(PurePosixPath(file))
        else:
            file = pdb.get((a["uniprot"], a["structure_id"].upper(), a["chain"]))
            if file is None:
                missing.append(f"{a['uniprot']} PDB {a['structure_id']} chain {a['chain']}")
                continue
            paths[n] = f"03_pdbfixer/{a['uniprot']}/{PurePosixPath(file).name}"
    if missing:
        raise ValueError(
            f"{len(missing)} target(s) have no row in Protein_info.csv (the run folder is "
            f"inconsistent): " + "; ".join(missing[:20]) + ("; ..." if len(missing) > 20 else ""))
    return paths


def select_targets(G, ranking, exclude, tie_break: str = "union", hub: str | None = None,
                   protein_info=None):
    """The target table: the hub(s) with Similarity 1.0 and Is_Hub True, plus their
    non-excluded neighbours with Similarity = the largest Pmax to any hub, best first (ties
    by Node). `hub` overrides the hub choice (regression tests only). With `protein_info`
    (the Protein_info.csv table) Receptor_Path is filled in; without it, it is empty."""
    import pandas as pd

    exclude = set(exclude)
    if hub is not None:
        if hub not in G:
            raise ValueError(f"hub {hub} is not in the network")
        hubs = [hub]
    else:
        tied = tied_top(ranking, exclude)
        if not tied:
            raise ValueError("no node can be a hub: every node in the network is excluded "
                             f"({', '.join(sorted(exclude)) or 'none'}) or the network is empty")
        hubs = tied if tie_break == "union" else tied[:1]
    hub_set = set(hubs)
    sim = {h: 1.0 for h in hubs}
    for h in hubs:
        for nb, data in G[h].items():
            if nb in hub_set or G.nodes[nb]["uniprot"] in exclude:
                continue
            sim[nb] = max(sim.get(nb, 0.0), data["weight"])
    receptors = _receptor_paths(G, sim, protein_info) if protein_info is not None else {}
    rows = []
    for n, s in sim.items():
        a = G.nodes[n]
        gene, organism = split_entry_name(a["entry_name"])
        rows.append({
            "Node": n, "UniProt_ID": a["uniprot"], "Gene_Name": gene, "Organism": organism,
            "Structure_ID": a["structure_id"], "Similarity": s, "Protein": a["protein"],
            "Chain": a["chain"], "Pocket": a["pocket"], "Pocket_Path": a["pocket_path"],
            "Receptor_Path": receptors.get(n, ""), "Is_Hub": n in hub_set})
    df = pd.DataFrame(rows, columns=TARGET_COLUMNS)
    return df.sort_values(["Similarity", "Node"], ascending=[False, True],
                          kind="mergesort").reset_index(drop=True)


# --- stage -----------------------------------------------------------------------------

def _exclude(args) -> dict[str, str]:
    """{UniProt ID: reason} from --exclude (each value may also be comma separated)."""
    ids = []
    for item in getattr(args, "exclude", list(DEFAULT_EXCLUDE)) or []:
        ids += [x.strip() for x in str(item).split(",") if x.strip()]
    return {u: EXCLUDE_REASONS.get(u, "excluded with --exclude") for u in sorted(set(ids))}


def _options(args) -> dict:
    return {"pmax": float(getattr(args, "pmax", DEFAULT_PMAX)), "exclude": _exclude(args),
            "hub_tie_break": getattr(args, "hub_tie_break", "union")}


def _fingerprint(layout, options: dict) -> str:
    payload = {"options": options, "pocketsim_sha256": sha256_file(layout.pocketsim_csv),
               "protein_info_sha256": (sha256_file(layout.protein_info)
                                       if layout.protein_info.exists() else None)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _require_pocketsim(layout) -> None:
    if not layout.pocketsim_csv.exists():
        raise FileNotFoundError(
            f"{layout.pocketsim_csv} not found; run `polyphar similarity` first")


def _finished(layout, options) -> dict | None:
    d = layout.targets
    try:
        manifest = json.loads((d / "manifest.json").read_text())
        sizes = {name: (d / name).stat().st_size for name in OUTPUTS}
    except (OSError, ValueError):
        return None
    if manifest.get("outputs") != sizes or manifest.get("fingerprint") != _fingerprint(
            layout, options):
        return None
    return dict(manifest["counts"])


def units(layout, args) -> WorkUnits:
    """Trivial: one unit, or none when the outputs are finished on the current inputs."""
    _require_pocketsim(layout)
    if _finished(layout, _options(args)) is not None:
        return WorkUnits(0, "stage")
    return WorkUnits(1, "stage")


def _atomic_csv(df, path: Path) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            df.to_csv(fh, index=False)
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _announce(progress, ranking, tied: list[str], hubs: list[str], options: dict) -> None:
    """Name a tie at the top degree on stdout and in a notification."""
    degree = int(ranking.set_index("Node").loc[tied[0], "Degree"])
    how = ("all are hubs (union)" if options["hub_tie_break"] == "union"
           else f"{options['hub_tie_break']} picked {hubs[0]}")
    msg = f"{len(tied)} non-excluded nodes tie at degree {degree}: {', '.join(tied)}; {how}"
    print(f"select: tie at the top degree - {msg}", flush=True)
    notifier = getattr(progress, "notifier", None)
    if notifier is not None:
        try:
            notifier("select: tie at the top degree", msg)
        except Exception:
            pass


def run(layout, args, progress) -> dict:
    import networkx as nx
    import pandas as pd

    _require_pocketsim(layout)
    options = _options(args)
    done = _finished(layout, options)
    if done is not None:   # finished on exactly these inputs and options: touch nothing
        return done
    exclude = set(options["exclude"])
    if not layout.protein_info.exists():
        raise FileNotFoundError(f"{layout.protein_info} not found; run `polyphar structures` first")

    pairs = unique_protein_pairs(read_pocketsim(layout.pocketsim_csv), options["pmax"])
    if pairs.empty:
        raise ValueError(f"no protein pair has Pmax >= {options['pmax']} in {layout.pocketsim_csv}")
    g = build_graph(pairs)
    ranking = rank_hubs(g, options["hub_tie_break"])
    targets = select_targets(g, ranking, exclude, options["hub_tie_break"],
                             protein_info=_read_protein_info(layout.protein_info))
    tied = tied_top(ranking, exclude)
    hubs = targets.loc[targets["Is_Hub"], "Node"].tolist()

    d = layout.targets
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    _atomic_csv(targets, layout.final_targets_csv)
    _atomic_csv(ranking, d / "centrality.csv")
    _atomic_csv(pairs, d / "unique_pairs.csv")
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    if len(tied) > 1:
        _announce(progress, ranking, tied, hubs, options)
    progress.advance()

    counts = {"pairs": len(pairs), "nodes": g.number_of_nodes(), "edges": g.number_of_edges(),
              "hubs": len(hubs), "tied_top": len(tied), "targets": len(targets), "failed": 0}
    write_manifest(d, STAGE, args=options, inputs=[layout.pocketsim_csv, layout.protein_info],
                   versions={"networkx": nx.__version__, "pandas": pd.__version__},
                   counts=counts, timing={},
                   extra={"fingerprint": _fingerprint(layout, options),
                          "outputs": {name: (d / name).stat().st_size for name in OUTPUTS},
                          "top_tie": tied, "hubs": hubs})
    return counts
