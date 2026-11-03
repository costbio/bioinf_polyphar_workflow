"""THROWAWAY: docking-box volume vs pocket volume for the filter-eligible pockets.
smina search cost scales with box volume; an elongated groove has a huge box for a small
pocket, so it is both slow AND under-sampled at fixed exhaustiveness.
Box = pocket-atom bounding box + the padding dock.py uses."""
import csv, sys, json, os, statistics
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
PAD = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
nseg = {a: len(v) for a, v in json.load(open(f"{work}/tm_ranges.json")).items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}
rows = []
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v, dsc = num(r,"volume"), num(r,"drugScore")
    if v is None or dsc is None or dsc < 0.5 or v < 300 or v > 2000: continue
    if g["segs"] > 3: continue
    isaf = r["Protein"].startswith("AF-")
    if isaf and pl.get(f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}', 0) < 70: continue
    if nseg.get(acc, 0) >= 2: continue
    parts = pk.split("_")
    if len(parts) > 2:
        pv = vol_of.get((acc, "_".join(parts[:-1])))
        if pv and pv > 2000: continue
    xs, ys, zs = [], [], []
    try:
        with open(f"{run}/04_pockets/{acc}/{fn}") as fh:
            for ln in fh:
                if ln.startswith("ATOM"):
                    xs.append(float(ln[30:38])); ys.append(float(ln[38:46])); zs.append(float(ln[46:54]))
    except (OSError, ValueError):
        continue
    if not xs: continue
    d = [max(xs)-min(xs)+2*PAD, max(ys)-min(ys)+2*PAD, max(zs)-min(zs)+2*PAD]
    box = d[0]*d[1]*d[2]
    rows.append({"acc":acc,"gene":r["Gene Name"].replace("_HUMAN",""),"pk":pk,"vol":v,
                 "box":box,"eff":v/box,"maxdim":max(d),
                 "ellca":num(r,"ell c/a"),"ellba":num(r,"ell b/a")})
print(f"eligible pockets measured: {len(rows)}")
def q(a,p):
    a=sorted(a); return a[int(p*(len(a)-1))]
bx=[r["box"] for r in rows]; ef=[r["eff"] for r in rows]; md=[r["maxdim"] for r in rows]
print(f"box volume (A^3): median {q(bx,.5):.0f}  p75 {q(bx,.75):.0f}  p90 {q(bx,.9):.0f}  max {q(bx,1):.0f}")
print(f"box efficiency (pocket vol / box vol): median {q(ef,.5):.3f}  p10 {q(ef,.1):.3f}  min {q(ef,0):.3f}")
print(f"max box dimension (A): median {q(md,.5):.1f}  p90 {q(md,.9):.1f}  max {q(md,1):.1f}")
print("\nif search cost ~ box volume, relative cost vs the median pocket:")
for p in (0.5, 0.75, 0.9, 0.99, 1.0):
    print(f"  p{int(p*100):<3} box {q(bx,p):>8.0f} A^3  -> {q(bx,p)/q(bx,.5):>5.1f}x")
print("\ncandidate filters (panel = proteins keeping >=1 pocket):")
allp = {r["acc"] for r in rows}
for lab, f in (("none", lambda r: True),
               ("max dimension <= 30 A", lambda r: r["maxdim"] <= 30),
               ("max dimension <= 26 A", lambda r: r["maxdim"] <= 26),
               ("box <= 20000 A^3", lambda r: r["box"] <= 20000),
               ("box <= 12000 A^3", lambda r: r["box"] <= 12000),
               ("box eff >= 0.10", lambda r: r["eff"] >= 0.10),
               ("box eff >= 0.15", lambda r: r["eff"] >= 0.15)):
    k = [r for r in rows if f(r)]
    if not k: continue
    prot = {r["acc"] for r in k}
    print(f"  {lab:<24} {len(k):>4} pockets, {len(prot):>4} proteins "
          f"({len(prot)/len(allp):.0%} of panel), median box {statistics.median([r['box'] for r in k]):>7.0f}")
