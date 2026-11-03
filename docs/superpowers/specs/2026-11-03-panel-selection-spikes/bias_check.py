"""Is the cross-class segment-count bias big enough to justify the SIFTS mapping?

Claim under test: pdbfixer models in missing loops, so experimental pocket residue
numbering is contiguous and segment counts come out LOWER than for AlphaFold pockets,
which are numbered in UniProt positions with real gaps. If true, `segments <= 3` admits
experimental pockets more readily than AF ones. Measure the distributions directly.
"""
import csv, json, math, os, sys, collections, statistics
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

# Pre-segment-filter population: everything that passes the OTHER pocket filters, so the
# segment cut is seen in isolation.
pre = []
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v, d = num(r,"volume"), num(r,"drugScore")
    if v is None or d is None or d < 0.5 or not (300 <= v <= 2000): continue
    p = pk.split("_")
    if len(p) > 2:
        pv = vol_of.get((acc, "_".join(p[:-1])))
        if pv and pv > 2000: continue
    isaf = r["Protein"].startswith("AF-")
    if isaf and pl.get(f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}', 0) < 70: continue
    if nseg.get(acc, 0) >= 2: continue
    pre.append({"acc":acc,"af":isaf,"segs":g["segs"],"nres":g["nres"],"span":g.get("span")})

af  = [x for x in pre if x["af"]]
exp = [x for x in pre if not x["af"]]
print(f"population before the segment cut: {len(pre)} pockets "
      f"({len(af)} AlphaFold, {len(exp)} experimental)\n")

def dist(label, xs):
    s = sorted(x["segs"] for x in xs)
    if not s: return
    c = collections.Counter(s)
    print(f"  {label:<14} n={len(s):>5}  median {statistics.median(s):>4.1f}  "
          f"mean {statistics.mean(s):>5.2f}  "
          f"segs<=3 {sum(1 for v in s if v<=3)/len(s):>5.1%}")
    print(f"    {'counts:':<12} " + "  ".join(f"{k}:{c[k]}" for k in sorted(c)[:10]))
print("SEGMENT COUNT DISTRIBUTION")
dist("AlphaFold", af); dist("experimental", exp)

a_pass = sum(1 for x in af if x["segs"] <= 3)/len(af) if af else 0
e_pass = sum(1 for x in exp if x["segs"] <= 3)/len(exp) if exp else 0
print(f"\n  pass rate at segments<=3: AF {a_pass:.1%} vs experimental {e_pass:.1%}  "
      f"-> gap {abs(e_pass-a_pass)*100:.1f} points "
      f"({'experimental favoured' if e_pass>a_pass else 'AF favoured'})")

print("\n  controlling for pocket size (lining residues), pass rate at segs<=3:")
print(f"    {'residue bin':<14} {'AF n':>6} {'AF pass':>8} {'exp n':>6} {'exp pass':>9} {'gap':>7}")
bins = [(0,10),(10,15),(15,20),(20,30),(30,50),(50,10**6)]
for lo,hi in bins:
    A=[x for x in af if lo<=x["nres"]<hi]; E=[x for x in exp if lo<=x["nres"]<hi]
    if len(A)<5 or len(E)<5: continue
    ap=sum(1 for x in A if x["segs"]<=3)/len(A); ep=sum(1 for x in E if x["segs"]<=3)/len(E)
    name = f"{lo}-{'inf' if hi >= 10**6 else hi}"
    print(f"    {name:<14} {len(A):>6} {ap:>7.1%} "
          f"{len(E):>6} {ep:>8.1%} {(ep-ap)*100:>+6.1f}")

# How many pockets would the segment cut decide differently if experimental counts were
# inflated to match AF? Upper bound on the fix's effect.
print("\n  UPPER BOUND on what the SIFTS fix could change:")
exp_at_boundary = sum(1 for x in exp if x["segs"] <= 3)
print(f"    experimental pockets currently passing segs<=3: {exp_at_boundary}")
print(f"    of which seqres-unmappable (~65%): ~{round(exp_at_boundary*0.65)}")
print(f"    i.e. the fix is unverifiable for ~{round(exp_at_boundary*0.65)} of the "
      f"{exp_at_boundary} pockets it would act on")
