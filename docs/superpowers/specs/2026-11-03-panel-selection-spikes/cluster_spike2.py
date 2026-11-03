"""THROWAWAY spike 2: why doesn't the mitotic group separate?
(a) is high Pmax confounded by pocket size? (b) does a quality/druggability
filter make the group separate? Induced subgraph == what filtering at
`similarity` would do, since filtering only removes pockets, never adds edges."""
import csv, sys, statistics, collections
import networkx as nx
from networkx.algorithms.community import louvain_communities

pairs_csv, drugscore_csv = sys.argv[1], sys.argv[2]

prop = {}
with open(drugscore_csv, newline="") as fh:
    for r in csv.DictReader(fh):
        node = f'{r["UniProt ID"]}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
        def f(k):
            try: return float(r[k])
            except (ValueError, KeyError): return None
        prop[node] = {"vol": f("volume"), "ds": f("drugScore"),
                      "atms": f("siteAtms"), "depth": f("depth"),
                      "af": r["Protein"].startswith("AF-"),
                      "sel": r.get("Selected", "")}
print(f"DrugScore.csv rows: {len(prop)}")

G = nx.Graph(); gene_of = {}
with open(pairs_csv, newline="") as fh:
    for r in csv.DictReader(fh):
        G.add_edge(r["Node"], r["Node_2"], weight=float(r["Pmax"]))
        gene_of[r["Node"]] = r["Gene Name"].replace("_HUMAN","")
        gene_of[r["Node_2"]] = r["Gene Name_2"].replace("_HUMAN","")
matched = sum(1 for n in G if n in prop)
print(f"graph nodes {G.number_of_nodes()}, matched to DrugScore: {matched}")
if matched < G.number_of_nodes():
    miss = [n for n in G if n not in prop][:3]
    print("  UNMATCHED examples:", miss)

# --- (a) size confound: pockets in the >=0.7 graph vs all selected pockets
sel_all = [p for p in prop.values() if p["sel"].lower() in ("true","1","yes")] or list(prop.values())
def med(vals, k):
    v = [x[k] for x in vals if x[k] is not None]
    return (round(statistics.median(v),1), round(statistics.mean(v),1), len(v)) if v else None
ingraph = [prop[n] for n in G if n in prop]
print("\n(a) SIZE CONFOUND  (median, mean, n)")
for k in ("vol","atms","ds","depth"):
    print(f"  {k:6s} all-selected {med(sel_all,k)}   in >=0.7 graph {med(ingraph,k)}")
# high-degree nodes specifically
deg = dict(G.degree())
hi = [prop[n] for n in sorted(deg, key=deg.get, reverse=True)[:200] if n in prop]
for k in ("vol","atms"):
    print(f"  {k:6s} top-200 degree nodes {med(hi,k)}")
print(f"  AF share: all-selected {sum(p['af'] for p in sel_all)/len(sel_all):.0%}, "
      f">=0.7 graph {sum(p['af'] for p in ingraph)/len(ingraph):.0%}, "
      f"top-200 degree {sum(p['af'] for p in hi)/len(hi):.0%}")

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}

def report(tag, H):
    if H.number_of_edges() == 0:
        print(f"\n{tag}: empty graph"); return
    comms = sorted(louvain_communities(H, weight="weight", resolution=1.0, seed=42),
                   key=len, reverse=True)
    best, score = None, -1
    for c in comms:
        g = {gene_of[n] for n in c} & MITOTIC
        if len(g) > score: best, score = c, len(g)
    sub = H.subgraph(best)
    genes = {gene_of[n] for n in best}
    mit = genes & MITOTIC
    iw = [d["weight"] for *_, d in sub.edges(data=True)]
    print(f"\n{tag}: {H.number_of_nodes()} nodes / {H.number_of_edges()} edges, "
          f"{len(comms)} communities")
    print(f"  best community: {len(best)} pockets / {len(genes)} proteins, "
          f"density {nx.density(sub):.3f}, mean internal Pmax {statistics.mean(iw):.3f}")
    print(f"  mitotic {len(mit)}/{len(MITOTIC & set(gene_of.values()))} present here; "
          f"purity (mitotic/all proteins) {len(mit)/len(genes):.2f}")
    print(f"  non-mitotic co-members: {sorted(genes - MITOTIC)[:25]}")

report("UNFILTERED", G)
keep = [n for n in G if n in prop and not prop[n]["af"]]
report("EXPERIMENTAL-ONLY (drop AF models)", G.subgraph(keep))
for vmin in (300, 600, 1000):
    keep = [n for n in G if n in prop and (prop[n]["vol"] or 0) >= vmin]
    report(f"volume >= {vmin} A^3", G.subgraph(keep))
for dsmin in (0.6, 0.8):
    keep = [n for n in G if n in prop and (prop[n]["ds"] or 0) >= dsmin]
    report(f"drugScore >= {dsmin}", G.subgraph(keep))
keep = [n for n in G if n in prop and not prop[n]["af"]
        and (prop[n]["vol"] or 0) >= 600 and (prop[n]["ds"] or 0) >= 0.6]
report("COMBINED (experimental + vol>=600 + ds>=0.6)", G.subgraph(keep))
