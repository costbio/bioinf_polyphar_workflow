"""THROWAWAY spike 5: cluster the network built from FILTER-ELIGIBLE pockets only.
Every earlier clustering test used the 2594-node Pmax>=0.7 graph, which held only
2594 of 11599 pockets and was biased to large ones (median vol 1009 vs 269.7).
The compact eligible set (~791 pockets) has never been clustered. Pairs come from the
existing PocketMatch_score.txt, so PocketMatch is NOT re-run.
Also answers: how many eligible pockets are SUBpockets (P_n_m) of oversized parents?
"""
import csv, sys, json, os, statistics, collections
import networkx as nx
from networkx.algorithms.community import louvain_communities
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
NULL = "NULL____"

ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
tmr = json.load(open(f"{work}/tm_ranges.json"))
nseg = {a: len(v) for a, v in tmr.items()}

def num(r, k):
    try: return float(r[k])
    except (ValueError, KeyError): return 0.0

# volume of each pocket by (acc, pocket id), to look up a subpocket's parent
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r, "volume") for r in ds}

elig, subinfo = {}, []
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v = num(r, "volume")
    if num(r,"drugScore") < 0.5 or v < 300 or v > 2000: continue
    if g["segs"] > 3: continue
    isaf = r["Protein"].startswith("AF-")
    node = f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
    if isaf and pl.get(node, 0) < 70: continue
    t = nseg.get(acc, 0)
    if t >= 2: continue
    elig[fn] = {"acc": acc, "gene": r["Gene Name"].replace("_HUMAN",""),
                "node": node, "vol": v, "segs": g["segs"], "nres": g["nres"],
                "af": isaf, "pk": pk}
    parts = pk.split("_")               # P_n  or  P_n_m
    if len(parts) > 2:
        parent = "_".join(parts[:-1])
        subinfo.append((acc, pk, v, vol_of.get((acc, parent))))
print(f"eligible pockets: {len(elig)}, proteins: {len({v['acc'] for v in elig.values()})}")
print(f"  AF {sum(v['af'] for v in elig.values())}, "
      f"median volume {statistics.median([v['vol'] for v in elig.values()]):.0f}")

print(f"\nSUBPOCKETS among eligible: {len(subinfo)} of {len(elig)}")
withpar = [(a,p,v,pv) for a,p,v,pv in subinfo if pv]
if withpar:
    big = [x for x in withpar if x[3] > 2000]
    print(f"  with a known parent: {len(withpar)}; parent volume > 2000 A^3: {len(big)}")
    print(f"  parent volume: median {statistics.median([x[3] for x in withpar]):.0f}, "
          f"max {max(x[3] for x in withpar):.0f}")
    for a,p,v,pv in sorted(big, key=lambda x:-x[3])[:8]:
        print(f"    {a} {p}: sub {v:.0f} A^3, parent {pv:.0f} A^3")

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}

print("\nextracting pairs among eligible pockets from PocketMatch_score.txt", flush=True)
pairs = []
with open(f"{run}/05_similarity/PocketMatch_score.txt", encoding="latin-1") as fh:
    for ln, line in enumerate(fh):
        if ln % 20_000_000 == 0 and ln: print(f"  {ln//1_000_000}M", flush=True)
        f = line.split()
        if f[2] == NULL or f[3] == NULL: continue
        a, b = f[0].rstrip("_"), f[1].rstrip("_")
        ea, eb = elig.get(a), elig.get(b)
        if ea is None or eb is None: continue
        if ea["acc"] == eb["acc"]: continue
        pmin, pmax = float(f[2]), float(f[3])
        if pmin == 1.0 or pmax == 1.0: continue
        pairs.append((ea["node"], eb["node"], pmax))
print(f"pairs among eligible, different proteins: {len(pairs)}")
allp = sorted(p for *_, p in pairs)
def q(x): return allp[int(x*(len(allp)-1))]
print(f"Pmax over these pairs: median {q(.5):.3f} p90 {q(.9):.3f} p99 {q(.99):.3f} "
      f"p99.9 {q(.999):.3f} max {allp[-1]:.3f}")

gnode = {v["node"]: v["gene"] for v in elig.values()}
base = len({v["gene"] for v in elig.values()} & MITOTIC) / len({v["gene"] for v in elig.values()})
print(f"\nmitotic base rate among eligible proteins: {base:.3f}")

for thr in (0.5, 0.6, 0.7, 0.75, 0.8):
    G = nx.Graph()
    for a, b, p in pairs:
        if p >= thr: G.add_edge(a, b, weight=p)
    if G.number_of_edges() < 20:
        print(f"\n--- Pmax >= {thr}: only {G.number_of_edges()} edges, skipped"); continue
    print(f"\n--- Pmax >= {thr}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
          f"{nx.number_connected_components(G)} components, "
          f"mean degree {2*G.number_of_edges()/G.number_of_nodes():.1f}")
    for res in (1.0, 1.5):
        comms = sorted(louvain_communities(G, weight="weight", resolution=res, seed=42),
                       key=len, reverse=True)
        best, sc = None, -1
        for c in comms:
            k = len({gnode[n] for n in c} & MITOTIC)
            if k > sc: best, sc = c, k
        sub = G.subgraph(best); genes = {gnode[n] for n in best}
        mit = genes & MITOTIC
        iw = [d["weight"] for *_, d in sub.edges(data=True)]
        pur = len(mit)/len(genes)
        print(f"    res {res}: {len(comms)} communities "
              f"sizes top6 {[len(c) for c in comms[:6]]}")
        print(f"      best: {len(best)} pockets / {len(genes)} proteins, "
              f"density {nx.density(sub):.3f}, internal Pmax {statistics.mean(iw):.3f}, "
              f"mitotic {len(mit)}, purity {pur:.2f} ({pur/base:.1f}x base)")
        print(f"      members: {sorted(genes)[:28]}")
