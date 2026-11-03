"""REPAIR STEP 1: reconcile the panel to one number end-to-end, on the DOCKABLE set.

Applies every filter the design states, including the two that earlier spikes omitted:
  - TM policy C (not A): drop multi-pass proteins, AND drop AF pockets whose lining
    residues overlap a TM segment by >10%
  - box longest side <= 30 A  AND  Uni-Dock grid legality at 0.375 A spacing
Then re-derives clusters, coverage, control donor pool and the experimental-only subset
on that dockable set, so the downstream numbers stop being stale.
"""
import csv, json, math, os, statistics, sys, collections
import networkx as nx
from networkx.algorithms.community import louvain_communities
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
NULL, MAXAX, MAXTOT, SPACING, BOXCAP = "NULL____", 128, 531441, 0.375, 30.0

ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl  = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
tmr = json.load(open(f"{work}/tm_ranges.json"))
nseg = {a: len(v) for a, v in tmr.items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}

def pocket_res_and_box(acc, fn):
    res, lo, hi = set(), None, None
    try:
        for ln in open(f"{run}/04_pockets/{acc}/{fn}"):
            if ln.startswith(("ATOM","HETATM")):
                try: res.add(int(ln[22:26]))
                except ValueError: pass
                c = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                if lo is None: lo, hi = list(c), list(c)
                else:
                    for k in range(3): lo[k]=min(lo[k],c[k]); hi[k]=max(hi[k],c[k])
    except OSError: return None, None
    if lo is None: return None, None
    return res, [hi[k]-lo[k] for k in range(3)]

steps = collections.Counter()
elig = []
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    steps["0 selected pockets"] += 1
    v, d = num(r,"volume"), num(r,"drugScore")
    if v is None or d is None or d < 0.5: continue
    steps["1 drugScore>=0.5"] += 1
    if not (300 <= v <= 2000): continue
    steps["2 volume 300-2000"] += 1
    if g["segs"] > 3: continue
    steps["3 segments<=3"] += 1
    p = pk.split("_")
    if len(p) > 2:
        pv = vol_of.get((acc, "_".join(p[:-1])))
        if pv and pv > 2000: continue
    steps["4 subpocket parent ok"] += 1
    isaf = r["Protein"].startswith("AF-")
    if isaf and pl.get(f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}', 0) < 70: continue
    steps["5 AF pLDDT>=70"] += 1
    if nseg.get(acc, 0) >= 2: continue
    steps["6 TM: drop multi-pass"] += 1
    res, dims = pocket_res_and_box(acc, fn)
    if res is None: continue
    # TM policy C: pocket-level overlap, AF models only (PDB resnum == UniProt position)
    if isaf and nseg.get(acc, 0) >= 1:
        tmres = {x for a, b in tmr.get(acc, []) for x in range(a, b+1)}
        if res and len(res & tmres)/len(res) > 0.10: continue
    steps["7 TM policy C overlap"] += 1
    if max(dims) > BOXCAP: continue
    steps["8 box longest side<=30A"] += 1
    n = [math.floor(s/SPACING)+1 for s in dims]
    if max(n) > MAXAX or n[0]*n[1]*n[2] > MAXTOT: continue
    steps["9 Uni-Dock grid legal"] += 1
    elig.append({"acc":acc,"fn":fn,"pk":pk,"af":isaf,"ds":d,"vol":v,
                 "gene":r["Gene Name"].replace("_HUMAN",""),
                 "node":f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}'})

print("FUNNEL (pockets surviving each filter):")
for k in sorted(steps): print(f"  {k:<28} {steps[k]:>6}")
prot = {e["acc"] for e in elig}
print(f"\nDOCKABLE: {len(elig)} pockets, {len(prot)} proteins "
      f"(AF {len({e['acc'] for e in elig if e['af']})})")

# one pocket per protein, highest drugScore (ties: smallest pocket name)
best = {}
for e in sorted(elig, key=lambda x: (-x["ds"], x["pk"])):
    best.setdefault(e["acc"], e)
print(f"PANEL (one pocket per protein): {len(best)} proteins; "
      f"AF share {sum(b['af'] for b in best.values())/len(best):.0%}")
exp_only = [a for a,b in best.items() if not b["af"]]
print(f"experimental-only sensitivity subset: {len(exp_only)} proteins")

# clusters on the dockable eligible set
gnode = {e["node"]: e["gene"] for e in elig}
acc_of = {e["node"]: e["acc"] for e in elig}
ok = set(gnode)
pairs = []
with open(f"{run}/05_similarity/PocketMatch_score.txt", encoding="latin-1") as fh:
    name2node = {e["fn"]: e["node"] for e in elig}
    for ln in fh:
        f = ln.split()
        if f[2] == NULL or f[3] == NULL: continue
        a, b = f[0].rstrip("_"), f[1].rstrip("_")
        na, nb = name2node.get(a), name2node.get(b)
        if na is None or nb is None or acc_of[na] == acc_of[nb]: continue
        pmin, pmax = float(f[2]), float(f[3])
        if pmin == 1.0 or pmax == 1.0: continue
        if pmax >= 0.70: pairs.append((na, nb, pmax))
G = nx.Graph()
for a,b,p in pairs: G.add_edge(a,b,weight=p)
print(f"\nCLUSTER GRAPH on dockable set: {G.number_of_nodes()} pockets, "
      f"{G.number_of_edges()} edges")
comms = sorted(louvain_communities(G, weight="weight", resolution=1.0, seed=42),
               key=len, reverse=True) if G.number_of_edges() else []
fams = []
for c in comms:
    genes = sorted({gnode[n] for n in c}); accs = {acc_of[n] for n in c}
    if len(accs) < 3: continue
    sub = G.subgraph(c)
    iw = [d["weight"] for *_,d in sub.edges(data=True)]
    fams.append((genes, accs, len(c), nx.density(sub), statistics.mean(iw)))
cov = set().union(*[f[1] for f in fams]) if fams else set()
print(f"CLUSTERS (>=3 proteins): {len(fams)}, covering {len(cov)}/{len(best)} panel "
      f"proteins ({len(cov)/len(best):.0%})")
slots = sum(len(f[1]) for f in fams)
print(f"  protein-slots {slots} over {len(cov)} distinct proteins "
      f"(reuse factor {slots/max(len(cov),1):.2f})")
for genes, accs, npk, dens, iw in fams:
    print(f"  [{len(accs)} prot, {npk} pockets, density {dens:.2f}, Pmax {iw:.3f}] "
          f"{', '.join(genes)}")
print(f"\nCONTROL DONOR POOL (panel, not in any cluster): {len(set(best) - cov)}")
# how many panel proteins are in a cluster via a DIFFERENT pocket than the docked one
mism = sum(1 for a in cov if best[a]["node"] not in set().union(*[set(c) for c in comms]))
print(f"AUDIT FINDING 1 CHECK: panel proteins in a cluster whose DOCKED pocket is not "
      f"itself a clustered pocket: {mism} of {len(cov)}")
