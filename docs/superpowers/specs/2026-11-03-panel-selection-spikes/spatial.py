"""Can a SPATIAL compactness filter replace `segments <= 3` without the cross-class bias?

`segments` measures compactness in SEQUENCE, which depends on the residue-numbering
coordinate system and is confounded with structure provenance (AF 39.4% pass vs
experimental 17.5%). A spatial measure is identical across both by construction.

Candidates, all from the pocket PDB coordinates:
  rg          radius of gyration of lining atoms
  shape       rg / n_atoms^(1/3)   (dimensionless; compact objects have low values)
  asph        asphericity from the inertia tensor (0 = sphere, 1 = rod)
  maxdim      longest bounding-box side (already used by the 30 A cap)
Each is swept to the SAME stringency as `segments <= 3` so the comparison is fair.
"""
import csv, json, math, os, sys, statistics, collections
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl  = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
tmr = json.load(open(f"{work}/tm_ranges.json")); nseg={a:len(v) for a,v in tmr.items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}

def metrics(acc, fn):
    xs=[]
    try:
        for ln in open(f"{run}/04_pockets/{acc}/{fn}"):
            if ln.startswith(("ATOM","HETATM")):
                xs.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
    except OSError: return None
    n=len(xs)
    if n<4: return None
    cx=[sum(p[k] for p in xs)/n for k in range(3)]
    rg=math.sqrt(sum(sum((p[k]-cx[k])**2 for k in range(3)) for p in xs)/n)
    # inertia-tensor principal moments -> asphericity
    S=[[sum((p[i]-cx[i])*(p[j]-cx[j]) for p in xs)/n for j in range(3)] for i in range(3)]
    # symmetric 3x3 eigenvalues via characteristic polynomial (stable enough here)
    import numpy as np
    ev=sorted(np.linalg.eigvalsh(np.array(S)))
    l1,l2,l3=ev
    asph=(l3-0.5*(l1+l2))/max(l1+l2+l3,1e-9)
    lo=[min(p[k] for p in xs) for k in range(3)]; hi=[max(p[k] for p in xs) for k in range(3)]
    return {"n":n,"rg":rg,"shape":rg/(n**(1/3)),"asph":asph,
            "maxdim":max(hi[k]-lo[k] for k in range(3))}

pre=[]
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v,d = num(r,"volume"), num(r,"drugScore")
    if v is None or d is None or d<0.5 or not (300<=v<=2000): continue
    p=pk.split("_")
    if len(p)>2:
        pv=vol_of.get((acc,"_".join(p[:-1])))
        if pv and pv>2000: continue
    isaf=r["Protein"].startswith("AF-")
    if isaf and pl.get(f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}',0)<70: continue
    if nseg.get(acc,0)>=2: continue
    m=metrics(acc,fn)
    if m is None: continue
    pre.append({"acc":acc,"af":isaf,"segs":g["segs"],"ds":d,"pk":pk,**m})
af=[x for x in pre if x["af"]]; exp=[x for x in pre if not x["af"]]
print(f"population {len(pre)} pockets ({len(af)} AF, {len(exp)} experimental)")

base=[x for x in pre if x["segs"]<=3]
TARGET=len(base)
print(f"\nBASELINE  segments<=3 : keeps {TARGET} "
      f"(AF {sum(1 for x in base if x['af'])}, exp {sum(1 for x in base if not x['af'])})")
bg_af=sum(1 for x in af if x["segs"]<=3)/len(af); bg_ex=sum(1 for x in exp if x["segs"]<=3)/len(exp)
print(f"          pass rates   : AF {bg_af:.1%} vs exp {bg_ex:.1%}  gap {abs(bg_af-bg_ex)*100:.1f} pts")
print(f"          proteins kept: {len({x['acc'] for x in base})}")

print(f"\n{'metric':<10} {'cut':>9} {'kept':>6} {'AF pass':>8} {'exp pass':>9} {'gap(pts)':>9} {'proteins':>9}")
results={}
for key, lower_is_better in (("shape",True),("rg",True),("asph",True),("maxdim",True)):
    vals=sorted(x[key] for x in pre)
    cut=vals[min(TARGET-1, len(vals)-1)]          # keep the same number of pockets
    keep=[x for x in pre if x[key]<=cut]
    a=sum(1 for x in af if x[key]<=cut)/len(af); e=sum(1 for x in exp if x[key]<=cut)/len(exp)
    nprot=len({x["acc"] for x in keep})
    results[key]=(cut,keep,abs(a-e)*100)
    print(f"{key:<10} {cut:>9.3f} {len(keep):>6} {a:>7.1%} {e:>8.1%} "
          f"{abs(a-e)*100:>8.1f} {nprot:>9}")

best=min(results, key=lambda k: results[k][2])
cut,keep,gap=results[best]
print(f"\nSMALLEST CROSS-CLASS GAP: {best} <= {cut:.3f}  ({gap:.1f} pts, "
      f"against {abs(bg_af-bg_ex)*100:.1f} for segments<=3)")
print(f"  overlap with the segments<=3 set: "
      f"{len({(x['acc'],x['pk']) for x in keep} & {(x['acc'],x['pk']) for x in base})}"
      f" of {TARGET}")
# what the panel becomes with the spatial filter + the 30 A box cap
final=[x for x in keep if x["maxdim"]<=30.0]
bestp={}
for x in sorted(final,key=lambda y:(-y["ds"],y["pk"])): bestp.setdefault(x["acc"],x)
print(f"\nPANEL with {best}<= {cut:.3f} + box cap 30 A: {len(bestp)} proteins, "
      f"{len(final)} pockets, AF share {sum(b['af'] for b in bestp.values())/max(len(bestp),1):.0%}")
print(f"  (current segments<=3 + box cap gives 183 proteins / 492 pockets / 68% AF)")
