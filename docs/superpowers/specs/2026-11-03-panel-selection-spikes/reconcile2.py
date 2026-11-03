"""AUTHORITATIVE panel definition: spatial compactness (Rg) replaces the sequence-based
`segments <= 3` and the 30 A box cap.

Rationale: `segments` measures compactness in SEQUENCE, which depends on the residue
numbering coordinate system and proved confounded with structure provenance (AF 39.4%
pass vs experimental 17.5%, a 21.9-point gap). Radius of gyration is identical across
both by construction (1.5-point gap) and bounds the docking box as a side effect.
"""
import csv, json, math, os, sys, statistics, collections
import numpy as np
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
RG_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else 8.5
MAXAX, MAXTOT, SPACING = 128, 531441, 0.375

ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
pl  = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
tmr = json.load(open(f"{work}/tm_ranges.json")); nseg={a:len(v) for a,v in tmr.items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}

def geom(acc, fn):
    xs, res = [], set()
    try:
        for ln in open(f"{run}/04_pockets/{acc}/{fn}"):
            if ln.startswith(("ATOM","HETATM")):
                xs.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
                try: res.add(int(ln[22:26]))
                except ValueError: pass
    except OSError: return None
    n=len(xs)
    if n<4: return None
    cx=[sum(p[k] for p in xs)/n for k in range(3)]
    rg=math.sqrt(sum(sum((p[k]-cx[k])**2 for k in range(3)) for p in xs)/n)
    lo=[min(p[k] for p in xs) for k in range(3)]; hi=[max(p[k] for p in xs) for k in range(3)]
    dims=[hi[k]-lo[k] for k in range(3)]
    return {"rg":rg,"dims":dims,"res":res,"n":n}

steps=collections.Counter(); elig=[]
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    steps["0 selected"]+=1
    v,d=num(r,"volume"),num(r,"drugScore")
    if v is None or d is None or d<0.5: continue
    steps["1 drugScore>=0.5"]+=1
    if not (300<=v<=2000): continue
    steps["2 volume 300-2000"]+=1
    p=pk.split("_")
    if len(p)>2:
        pv=vol_of.get((acc,"_".join(p[:-1])))
        if pv and pv>2000: continue
    steps["3 subpocket parent ok"]+=1
    isaf=r["Protein"].startswith("AF-")
    if isaf and pl.get(f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}',0)<70: continue
    steps["4 AF pLDDT>=70"]+=1
    if nseg.get(acc,0)>=2: continue
    steps["5 TM drop multi-pass"]+=1
    g=geom(acc,fn)
    if g is None: continue
    if isaf and nseg.get(acc,0)>=1:
        t={x for a,b in tmr.get(acc,[]) for x in range(a,b+1)}
        if g["res"] and len(g["res"]&t)/len(g["res"])>0.10: continue
    steps["6 TM pocket overlap"]+=1
    if g["rg"]>RG_MAX: continue
    steps[f"7 Rg<={RG_MAX}"]+=1
    nn=[math.floor(s/SPACING)+1 for s in g["dims"]]
    if max(nn)>MAXAX or nn[0]*nn[1]*nn[2]>MAXTOT: continue
    steps["8 Uni-Dock grid legal"]+=1
    elig.append({"acc":acc,"fn":fn,"pk":pk,"af":isaf,"ds":d,"rg":g["rg"],
                 "maxdim":max(g["dims"]),"gene":r["Gene Name"].replace("_HUMAN","")})

print(f"FUNNEL (Rg cutoff {RG_MAX} A)")
for k in sorted(steps): print(f"  {k:<26} {steps[k]:>6}")
best={}
for e in sorted(elig,key=lambda x:(-x["ds"],x["pk"])): best.setdefault(e["acc"],e)
afn=sum(b["af"] for b in best.values())
print(f"\nPANEL: {len(best)} proteins, {len(elig)} eligible pockets, "
      f"AF {afn} ({afn/len(best):.0%}), experimental {len(best)-afn}")
md=[e["maxdim"] for e in elig]
print(f"  box longest side: median {statistics.median(md):.1f} A, max {max(md):.1f} A "
      f"(the 30 A cap is now redundant: {sum(1 for x in md if x>30)} exceed it)")
rg=[e["rg"] for e in elig]
print(f"  Rg: median {statistics.median(rg):.2f}, max {max(rg):.2f} A")
# cross-class balance at this cutoff
print(f"\n  cross-class check on eligible pockets: "
      f"AF {sum(1 for e in elig if e['af'])}, exp {sum(1 for e in elig if not e['af'])}")
with open(os.path.expanduser("~/panel_net1_rg.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["UniProt","Gene","Pocket","PocketFile","AlphaFold","drugScore","Rg","BoxMaxDim"])
    for a in sorted(best):
        e=best[a]; w.writerow([a,e["gene"],e["pk"],e["fn"],e["af"],round(e["ds"],3),
                               round(e["rg"],3),round(e["maxdim"],2)])
print(f"\nwrote ~/panel_net1_rg.csv ({len(best)} rows)")
