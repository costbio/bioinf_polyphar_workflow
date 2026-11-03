"""THROWAWAY: mean pLDDT over pocket residues for every selected AlphaFold pocket.
pLDDT lives in the B-factor column of the AF model; DoGSiteScorer zeroes B-factors in
its _res_P_*.pdb files, so residues are matched back to 02_structures by (chain,resnum).
Also tests whether low pocket pLDDT marks the mega-cavity hubs of the flat network."""
import csv, sys, os, statistics, collections
run = sys.argv[1]
ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected", "").lower() in ("true", "1", "yes")]
af = [r for r in ds if r["Protein"].startswith("AF-")]
print(f"selected pockets {len(ds)}, of which AlphaFold {len(af)}")

def model_plddt(acc, protein):
    p = f"{run}/02_structures/{acc}/{protein}.pdb"
    out = {}
    try:
        with open(p) as fh:
            for ln in fh:
                if ln.startswith("ATOM"):
                    try: out[(ln[21], int(ln[22:26]))] = float(ln[60:66])
                    except ValueError: pass
    except OSError:
        return None
    return out

cache, rows, nofile, nomatch = {}, [], 0, 0
for i, r in enumerate(af):
    if i and i % 2000 == 0: print(f"  {i}/{len(af)}", flush=True)
    acc, prot = r["UniProt ID"], r["Protein"]
    if (acc, prot) not in cache:
        cache[(acc, prot)] = model_plddt(acc, prot)
    pl = cache[(acc, prot)]
    if pl is None:
        nofile += 1; continue
    fn = f'{run}/04_pockets/{acc}/{r["Pocket Filename"]}'
    res = set()
    try:
        with open(fn) as fh:
            for ln in fh:
                if ln.startswith("ATOM"):
                    try: res.add((ln[21], int(ln[22:26])))
                    except ValueError: pass
    except OSError:
        nofile += 1; continue
    vals = [pl[k] for k in res if k in pl]
    if not vals:
        nomatch += 1; continue
    def f(k):
        try: return float(r[k])
        except (ValueError, KeyError): return 0.0
    rows.append({"acc": acc, "gene": r["Gene Name"].replace("_HUMAN",""),
                 "node": f'{acc}_{r["Gene Name"]}_{prot}_{r["Pocket"]}',
                 "nres": len(res), "matched": len(vals),
                 "mean": statistics.mean(vals), "min": min(vals),
                 "frac70": sum(v >= 70 for v in vals)/len(vals),
                 "vol": f("volume"), "ds": f("drugScore")})
print(f"computed {len(rows)} AF pockets; missing files {nofile}, no residue match {nomatch}")

m = sorted(x["mean"] for x in rows)
def q(p): return m[int(p*(len(m)-1))]
print(f"\nmean pocket pLDDT: p5 {q(.05):.1f}  p25 {q(.25):.1f}  median {q(.5):.1f}  "
      f"p75 {q(.75):.1f}  p95 {q(.95):.1f}")
print("fraction of AF pockets with mean pLDDT >= X:")
for t in (50, 60, 70, 80, 90):
    k = sum(1 for x in rows if x["mean"] >= t)
    print(f"  >= {t}: {k:>6} pockets ({k/len(rows):.0%}), "
          f"{len({x['acc'] for x in rows if x['mean']>=t}):>4} proteins")

# does pLDDT separate the mega-cavity artifacts?
deg = collections.Counter()
try:
    for r in csv.DictReader(open(f"{run}/06_targets/unique_pairs.csv", newline="")):
        deg[r["Node"]] += 1; deg[r["Node_2"]] += 1
except OSError:
    pass
ing = [x for x in rows if deg.get(x["node"], 0) > 0]
print(f"\nAF pockets in the flat Pmax>=0.7 network: {len(ing)}")
if ing:
    hi = sorted(ing, key=lambda x: -deg[x["node"]])[:100]
    print(f"  all AF selected pockets : mean pLDDT {statistics.mean([x['mean'] for x in rows]):.1f}, "
          f"median volume {statistics.median([x['vol'] for x in rows]):.0f}")
    print(f"  AF pockets in network   : mean pLDDT {statistics.mean([x['mean'] for x in ing]):.1f}, "
          f"median volume {statistics.median([x['vol'] for x in ing]):.0f}")
    print(f"  top-100 degree AF nodes : mean pLDDT {statistics.mean([x['mean'] for x in hi]):.1f}, "
          f"median volume {statistics.median([x['vol'] for x in hi]):.0f}")
    big = [x for x in rows if x["vol"] >= 3000]
    sml = [x for x in rows if x["vol"] < 1000]
    if big and sml:
        print(f"  AF pockets volume >=3000 A^3 (n={len(big)}): mean pLDDT "
              f"{statistics.mean([x['mean'] for x in big]):.1f}")
        print(f"  AF pockets volume < 1000 A^3 (n={len(sml)}): mean pLDDT "
              f"{statistics.mean([x['mean'] for x in sml]):.1f}")

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}
print("\nbest AF pocket per protein, panel size by pLDDT cutoff:")
best = {}
for x in rows:
    if x["acc"] not in best or x["mean"] > best[x["acc"]]["mean"]:
        best[x["acc"]] = x
for t in (50, 60, 70, 80, 90):
    k = {a for a, x in best.items() if x["mean"] >= t}
    g = {best[a]["gene"] for a in k}
    print(f"  pLDDT >= {t}: {len(k):>4} AF-only-capable proteins, "
          f"mitotic {len(g & MITOTIC)}")
import json
json.dump(rows, open(os.path.expanduser("~/scratch_clust/plddt_rows.json"), "w"))
