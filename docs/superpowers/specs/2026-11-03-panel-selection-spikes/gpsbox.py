"""THROWAWAY: cavity-derived docking box vs the current residue-envelope box.
DoGSiteScorer writes three files per pocket: _res_ (complete lining residues, used today
by pocket_box and by PocketMatch), _atms_ (lining atoms), _gps_ (cavity GRID POINTS).
The cavity box is the _gps_ bounding box plus padding for ligand rotation."""
import csv, sys, json, os, statistics
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
nseg = {a: len(v) for a, v in json.load(open(f"{work}/tm_ranges.json")).items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}

def extent(path):
    xs=ys=zs=None
    try:
        with open(path) as fh:
            for ln in fh:
                if ln.startswith(("ATOM","HETATM")):
                    x,y,z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
                    if xs is None: xs=[x,x]; ys=[y,y]; zs=[z,z]
                    else:
                        xs[0]=min(xs[0],x); xs[1]=max(xs[1],x)
                        ys[0]=min(ys[0],y); ys[1]=max(ys[1],y)
                        zs[0]=min(zs[0],z); zs[1]=max(zs[1],z)
    except OSError: return None
    if xs is None: return None
    return [xs[1]-xs[0], ys[1]-ys[0], zs[1]-zs[0]], [sum(a)/2 for a in (xs,ys,zs)]

rows, missing = [], 0
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v, dsc = num(r,"volume"), num(r,"drugScore")
    if v is None or dsc is None or dsc < 0.5 or v < 300 or v > 2000: continue
    if g["segs"] > 3: continue
    if r["Protein"].startswith("AF-") and pl.get(f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}',0) < 70: continue
    if nseg.get(acc,0) >= 2: continue
    parts = pk.split("_")
    if len(parts) > 2:
        pv = vol_of.get((acc,"_".join(parts[:-1])))
        if pv and pv > 2000: continue
    base = f"{run}/04_pockets/{acc}/{fn}"
    res = extent(base)
    gps = extent(base.replace("_res_P", "_gps_P"))
    if res is None: continue
    if gps is None: missing += 1; continue
    npts = sum(1 for ln in open(base.replace("_res_P","_gps_P")) if ln.startswith(("ATOM","HETATM")))
    rows.append({"vol":v,"res_d":res[0],"gps_d":gps[0],"npts":npts,
                 "res_box":res[0][0]*res[0][1]*res[0][2]})
print(f"pockets measured: {len(rows)}   missing _gps_ file: {missing}")
def q(a,p):
    a=sorted(a); return a[int(p*(len(a)-1))]

print(f"\ncavity grid points per pocket: median {q([r['npts'] for r in rows],.5):.0f}")
gv=[r['npts']*0.5**3 for r in rows]
print(f"implied cavity volume at 0.5 A grid: median {q(gv,.5):.0f} A^3 "
      f"(DoGSiteScorer 'volume' median {q([r['vol'] for r in rows],.5):.0f} A^3)")

print(f"\n{'box source':<34} {'median':>9} {'p90':>9} {'max':>9}  {'median maxdim':>13}")
def show(lab, dims):
    bv=[d[0]*d[1]*d[2] for d in dims]; md=[max(d) for d in dims]
    print(f"{lab:<34} {q(bv,.5):>9.0f} {q(bv,.9):>9.0f} {q(bv,1):>9.0f}  {q(md,.5):>13.1f}")
show("CURRENT: _res_ envelope, no pad", [r["res_d"] for r in rows])
for pad in (0, 4, 6, 8):
    show(f"cavity _gps_ + {pad} A each side",
         [[d+2*pad for d in r["gps_d"]] for r in rows])

print("\nspeedup if smina cost ~ box volume:")
cur = q([r["res_box"] for r in rows], .5)
for pad in (4, 6, 8):
    nv = q([(r["gps_d"][0]+2*pad)*(r["gps_d"][1]+2*pad)*(r["gps_d"][2]+2*pad) for r in rows], .5)
    print(f"  cavity + {pad} A: median box {nv:>8.0f} A^3  -> {cur/nv:>4.1f}x faster")

print("\npockets whose CAVITY itself is elongated (max dim, no pad):")
md=[max(r["gps_d"]) for r in rows]
for t in (20, 25, 30):
    n=sum(1 for x in md if x > t)
    print(f"  cavity max dimension > {t} A: {n} of {len(rows)} ({n/len(rows):.0%})")
