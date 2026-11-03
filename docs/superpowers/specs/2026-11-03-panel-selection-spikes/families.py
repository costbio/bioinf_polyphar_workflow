"""THROWAWAY spike 6: list every pocket family on the filter-eligible network at
Pmax>=0.65/0.7, with coverage of the 241-protein panel. Also quantifies the subpocket
hole: eligible subpockets whose PARENT exceeds the 2000 A^3 cap."""
import csv, sys, json, os, statistics, collections
import networkx as nx
from networkx.algorithms.community import louvain_communities
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
NULL = "NULL____"
ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
nseg = {a: len(v) for a, v in json.load(open(f"{work}/tm_ranges.json")).items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError): return 0.0
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}

def eligible(exclude_sub_of_big):
    out = {}
    for r in ds:
        acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
        g = geo.get(f"{acc}/{fn}")
        if g is None: continue
        v = num(r,"volume")
        if num(r,"drugScore") < 0.5 or v < 300 or v > 2000: continue
        if g["segs"] > 3: continue
        isaf = r["Protein"].startswith("AF-")
        node = f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}'
        if isaf and pl.get(node,0) < 70: continue
        if nseg.get(acc,0) >= 2: continue
        if exclude_sub_of_big:
            parts = pk.split("_")
            if len(parts) > 2:
                pv = vol_of.get((acc, "_".join(parts[:-1])))
                if pv and pv > 2000: continue
        out[fn] = {"acc":acc, "gene":r["Gene Name"].replace("_HUMAN",""), "node":node}
    return out

for tag, excl in (("all eligible", False), ("subpockets of >2000 parents removed", True)):
    e = eligible(excl)
    print(f"{tag}: {len(e)} pockets, {len({v['acc'] for v in e.values()})} proteins")

elig = eligible(True)          # the corrected rule
gnode = {v["node"]: v["gene"] for v in elig.values()}
pairs = []
with open(f"{run}/05_similarity/PocketMatch_score.txt", encoding="latin-1") as fh:
    for ln, line in enumerate(fh):
        if ln % 25_000_000 == 0 and ln: print(f"  {ln//1_000_000}M", flush=True)
        f = line.split()
        if f[2] == NULL or f[3] == NULL: continue
        a, b = f[0].rstrip("_"), f[1].rstrip("_")
        ea, eb = elig.get(a), elig.get(b)
        if ea is None or eb is None or ea["acc"] == eb["acc"]: continue
        pmin, pmax = float(f[2]), float(f[3])
        if pmin == 1.0 or pmax == 1.0: continue
        pairs.append((ea["node"], eb["node"], pmax))
print(f"pairs among corrected-eligible: {len(pairs)}")

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}
npanel = len({v["acc"] for v in elig.values()})
for thr in (0.65, 0.70):
    G = nx.Graph()
    for a,b,p in pairs:
        if p >= thr: G.add_edge(a,b,weight=p)
    comms = sorted(louvain_communities(G, weight="weight", resolution=1.0, seed=42),
                   key=len, reverse=True)
    fams = []
    for c in comms:
        genes = sorted({gnode[n] for n in c})
        if len(genes) < 3: continue
        sub = G.subgraph(c)
        iw = [d["weight"] for *_,d in sub.edges(data=True)]
        fams.append((genes, len(c), nx.density(sub), statistics.mean(iw)))
    cov = len({g for f in fams for g in f[0]})
    print(f"\n=== Pmax >= {thr}: {G.number_of_nodes()} pockets in graph, "
          f"{G.number_of_edges()} edges; {len(fams)} families of >=3 proteins "
          f"covering {cov}/{npanel} panel proteins ({cov/npanel:.0%})")
    for genes, npk, dens, iw in sorted(fams, key=lambda x:-len(x[0])):
        mit = [g for g in genes if g in MITOTIC]
        print(f"  [{len(genes)} prot, {npk} pockets, density {dens:.2f}, "
              f"Pmax {iw:.3f}] {', '.join(genes)}"
              + (f"   <- mitotic: {len(mit)}" if mit else ""))
