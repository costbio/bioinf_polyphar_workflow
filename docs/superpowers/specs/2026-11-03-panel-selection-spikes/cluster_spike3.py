"""THROWAWAY spike 3: size-conditional Pmax threshold.
Pmax = matched/max(n1,n2), so it is confounded by pocket size. Bin every pair by
(log2 size of n1, log2 size of n2) and keep the pairs that are extreme WITHIN their
size cell, at the same overall rate as the flat >=0.7 cutoff. Then ask the same
question: does the mitotic/replication group separate?
Legacy filters are matched to parse_scores: drop NULL, Pmin==1, Pmax==1, same-accession.
"""
import csv, sys, math, statistics, collections
import networkx as nx
from networkx.algorithms.community import louvain_communities

score_file, drugscore_csv, pairs_csv = sys.argv[1], sys.argv[2], sys.argv[3]
NULL = "NULL____"
NBINS, W = 501, 0.002          # Pmax histogram resolution
TARGET_EDGES = 11585           # match the flat-cutoff network for a fair comparison

# ---- name -> node / accession -------------------------------------------------------
node_of, acc_of, gene_of = {}, {}, {}
with open(drugscore_csv, newline="") as fh:
    for r in csv.DictReader(fh):
        fn = r["Pocket Filename"]
        node_of[fn] = f'{r["UniProt ID"]}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
        acc_of[fn] = r["UniProt ID"]
        gene_of[fn] = r["Gene Name"].replace("_HUMAN", "")
print(f"pocket names from DrugScore.csv: {len(node_of)}", flush=True)

def cell(n1, n2):
    b1, b2 = int(math.log2(max(n1, 1))), int(math.log2(max(n2, 1)))
    return (min(b1, b2), max(b1, b2))

# ---- pass 1: per-cell Pmax histogram ------------------------------------------------
hist = collections.defaultdict(lambda: [0] * NBINS)
kept = 0
with open(score_file, encoding="latin-1") as fh:
    for ln, line in enumerate(fh):
        if ln % 10_000_000 == 0 and ln:
            print(f"  pass1 {ln//1_000_000}M", flush=True)
        f = line.split()
        if f[2] == NULL or f[3] == NULL:
            continue
        pmin, pmax = float(f[2]), float(f[3])
        if pmin == 1.0 or pmax == 1.0:
            continue
        a, b = f[0].rstrip("_"), f[1].rstrip("_")
        if acc_of.get(a) == acc_of.get(b):
            continue
        n1, n2 = int(f[4]), int(f[5])
        hist[cell(n1, n2)][min(int(pmax / W), NBINS - 1)] += 1
        kept += 1
print(f"pass1 done: {kept} comparable pairs, {len(hist)} size cells", flush=True)

# ---- per-cell threshold at the global selection rate --------------------------------
rate = TARGET_EDGES / kept
thresh = {}
for c, h in hist.items():
    n = sum(h)
    want = n * rate
    acc = 0
    t = NBINS - 1
    for i in range(NBINS - 1, -1, -1):
        acc += h[i]
        if acc >= want:
            t = i
            break
    thresh[c] = t * W
print(f"\nglobal selection rate {rate:.3e}  (flat cutoff kept {TARGET_EDGES}/{kept})")
print("\nper-cell Pmax threshold (cell = log2 size of the two pockets):")
rows = sorted(hist, key=lambda c: -sum(hist[c]))[:18]
print(f"  {'cell':>10} {'pairs':>12} {'size range (n)':>22} {'thresh':>7}")
for c in rows:
    lo, hi = 2 ** c[0], 2 ** (c[1] + 1)
    print(f"  {str(c):>10} {sum(hist[c]):>12} {f'{lo}-{hi}':>22} {thresh[c]:>7.3f}")

# ---- pass 2: emit size-conditional edges --------------------------------------------
G = nx.Graph()
with open(score_file, encoding="latin-1") as fh:
    for ln, line in enumerate(fh):
        if ln % 20_000_000 == 0 and ln:
            print(f"  pass2 {ln//1_000_000}M", flush=True)
        f = line.split()
        if f[2] == NULL or f[3] == NULL:
            continue
        pmin, pmax = float(f[2]), float(f[3])
        if pmin == 1.0 or pmax == 1.0:
            continue
        a, b = f[0].rstrip("_"), f[1].rstrip("_")
        if acc_of.get(a) == acc_of.get(b):
            continue
        if pmax >= thresh[cell(int(f[4]), int(f[5]))]:
            G.add_edge(node_of[a], node_of[b], weight=pmax, n1=int(f[4]), n2=int(f[5]))
print(f"\nsize-conditional graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges,"
      f" {nx.number_connected_components(G)} components", flush=True)

gnode = {}
for fn, nd in node_of.items():
    gnode[nd] = gene_of[fn]

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}
pres = {gnode[n] for n in G} & MITOTIC
print(f"mitotic genes present: {len(pres)}/{len(MITOTIC)}  missing {sorted(MITOTIC-pres)}")

# size profile of the new graph, to confirm the confound is gone
sz = [max(d["n1"], d["n2"]) for *_, d in G.edges(data=True)]
print(f"edge max-pocket-size: median {statistics.median(sz):.0f}")
deg = dict(G.degree())
top = sorted(deg, key=deg.get, reverse=True)[:20]
print(f"top-20 degree genes: {[gnode[n] for n in top]}")
print(f"AF share of top-200 degree: "
      f"{sum('AF-' in n for n in sorted(deg,key=deg.get,reverse=True)[:200])/200:.0%}")

for res in (0.5, 1.0, 1.5):
    comms = sorted(louvain_communities(G, weight="weight", resolution=res, seed=42),
                   key=len, reverse=True)
    best, sc = None, -1
    for c in comms:
        k = len({gnode[n] for n in c} & MITOTIC)
        if k > sc: best, sc = c, k
    sub = G.subgraph(best)
    genes = {gnode[n] for n in best}
    mit = genes & MITOTIC
    iw = [d["weight"] for *_, d in sub.edges(data=True)]
    print(f"\n=== resolution {res}: {len(comms)} communities, "
          f"sizes top8 {[len(c) for c in comms[:8]]}")
    print(f"  best community: {len(best)} pockets / {len(genes)} proteins, "
          f"density {nx.density(sub):.3f}, mean internal Pmax {statistics.mean(iw):.3f}")
    print(f"  mitotic {len(mit)}/{len(pres)}, purity {len(mit)/len(genes):.2f}")
    print(f"  mitotic members: {sorted(mit)}")
    print(f"  non-mitotic co-members: {sorted(genes-MITOTIC)[:30]}")
