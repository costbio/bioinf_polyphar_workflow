"""THROWAWAY spike 4: corrected size-conditional threshold.
Pmax <= min(n1,n2)/max(n1,n2), so Pmax>=0.7 already implies size ratio >=0.7.
Spike 3 wrongly conditioned on disparity and admitted ~1%-overlap pairs.
Here: (1) keep only comparable-size pairs (ratio >= RATIO_MIN) so a high Pmax is
attainable in every bin, (2) bin by ABSOLUTE size, (3) take an EXACT per-bin top-k
by heap (no histogram bin-snapping), with sum(k) == TARGET_EDGES.
Readout is unchanged: does the mitotic/replication group separate?
"""
import csv, sys, math, heapq, statistics, collections
import networkx as nx
from networkx.algorithms.community import louvain_communities

score_file, drugscore_csv = sys.argv[1], sys.argv[2]
NULL = "NULL____"
RATIO_MIN = 0.7
TARGET_EDGES = 11585

node_of, acc_of, gene_of = {}, {}, {}
with open(drugscore_csv, newline="") as fh:
    for r in csv.DictReader(fh):
        fn = r["Pocket Filename"]
        node_of[fn] = f'{r["UniProt ID"]}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
        acc_of[fn] = r["UniProt ID"]
        gene_of[fn] = r["Gene Name"].replace("_HUMAN", "")

def sbin(n):            # absolute-size bin, half-log2 steps
    return int(2 * math.log2(max(n, 1)))

def rows():
    with open(score_file, encoding="latin-1") as fh:
        for ln, line in enumerate(fh):
            if ln % 20_000_000 == 0 and ln:
                print(f"    {ln//1_000_000}M", flush=True)
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
            lo, hi = min(n1, n2), max(n1, n2)
            if lo / hi < RATIO_MIN:
                continue
            yield a, b, pmax, hi

print("pass 1: counting comparable-size pairs per absolute-size bin", flush=True)
cnt = collections.Counter()
tot = 0
for _, _, _, hi in rows():
    cnt[sbin(hi)] += 1
    tot += 1
print(f"\ncomparable-size pairs (ratio>={RATIO_MIN}): {tot} of 64825170 "
      f"({tot/64825170:.2%}), {len(cnt)} size bins", flush=True)

rate = TARGET_EDGES / tot
k = {b: max(1, round(n * rate)) for b, n in cnt.items()}
print(f"selection rate {rate:.3e}; sum k = {sum(k.values())}")

print("\npass 2: exact per-bin top-k by Pmax", flush=True)
heaps = collections.defaultdict(list)
for a, b, pmax, hi in rows():
    bb = sbin(hi)
    h, kk = heaps[bb], k[bb]
    if len(h) < kk:
        heapq.heappush(h, (pmax, a, b))
    elif pmax > h[0][0]:
        heapq.heapreplace(h, (pmax, a, b))

print(f"\n{'bin':>5} {'size range':>16} {'pairs':>11} {'k':>6} {'Pmax cut':>9} {'Pmax med':>9}")
for bb in sorted(heaps):
    h = heaps[bb]
    lo, hi_ = int(2 ** (bb / 2)), int(2 ** ((bb + 1) / 2))
    ps = sorted(p for p, *_ in h)
    print(f"{bb:>5} {f'{lo}-{hi_}':>16} {cnt[bb]:>11} {len(h):>6} "
          f"{ps[0]:>9.3f} {statistics.median(ps):>9.3f}")

G = nx.Graph()
for h in heaps.values():
    for pmax, a, b in h:
        G.add_edge(node_of[a], node_of[b], weight=pmax)
print(f"\ngraph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
      f"{nx.number_connected_components(G)} components", flush=True)
gnode = {nd: gene_of[fn] for fn, nd in node_of.items()}

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}
pres = {gnode[n] for n in G} & MITOTIC
print(f"mitotic present: {len(pres)}/{len(MITOTIC)}  missing {sorted(MITOTIC-pres)}")
deg = dict(G.degree())
print(f"top-20 degree genes: {[gnode[n] for n in sorted(deg,key=deg.get,reverse=True)[:20]]}")
print(f"AF share top-200 degree: "
      f"{sum('AF-' in n for n in sorted(deg,key=deg.get,reverse=True)[:200])/200:.0%}")

for res in (0.5, 1.0, 1.5):
    comms = sorted(louvain_communities(G, weight="weight", resolution=res, seed=42),
                   key=len, reverse=True)
    best, sc = None, -1
    for c in comms:
        kk = len({gnode[n] for n in c} & MITOTIC)
        if kk > sc: best, sc = c, kk
    sub = G.subgraph(best); genes = {gnode[n] for n in best}; mit = genes & MITOTIC
    iw = [d["weight"] for *_, d in sub.edges(data=True)]
    print(f"\n=== resolution {res}: {len(comms)} communities, "
          f"sizes top8 {[len(c) for c in comms[:8]]}")
    print(f"  best community: {len(best)} pockets / {len(genes)} proteins, "
          f"density {nx.density(sub):.3f}, mean internal Pmax {statistics.mean(iw):.3f}")
    print(f"  mitotic {len(mit)}/{len(pres)}, purity {len(mit)/len(genes):.2f}")
    print(f"  mitotic members: {sorted(mit)}")
    print(f"  non-mitotic co-members: {sorted(genes-MITOTIC)[:30]}")
