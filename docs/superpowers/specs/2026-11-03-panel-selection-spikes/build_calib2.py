"""Rebuild the calibration run folder with a REPRESENTATIVE target sample.
build_calib.py used a maximin design over descriptors that correlate with box volume,
so its 48 pockets had a median box of 15979 A^3 against 9550 for the eligible set --
good for estimating slopes, useless for a throughput figure. Here: proportional
stratified random sample over box-volume quintiles, AF/experimental share proportional
within each quintile, fixed seed. Also writes calib_targets.csv with every descriptor
plus box volume, for the later score regression."""
import csv, sys, json, os, random, shutil, statistics, collections
from pathlib import Path
net1, calib, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(os.path.expanduser("~/scratch_clust"))
N = 48
rng = random.Random(42)

ds = [r for r in csv.DictReader(open(net1/"04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(work/"geo.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(work/"plddt_rows.json"))}
nseg = {a: len(v) for a, v in json.load(open(work/"tm_ranges.json")).items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}
DESC = ["volume","depth","lid/hull","hydrophobicity","siteAtms","surf/vol",
        "accept","donor","aromat","apolarAA","drugScore"]

def box(path):
    lo = hi = None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(("ATOM","HETATM")):
                c = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                if lo is None: lo, hi = list(c), list(c)
                else:
                    for k in range(3):
                        lo[k]=min(lo[k],c[k]); hi[k]=max(hi[k],c[k])
    if lo is None: return None
    return [hi[k]-lo[k] for k in range(3)]

elig = []
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v, dsc = num(r,"volume"), num(r,"drugScore")
    if v is None or dsc is None or dsc < 0.5 or v < 300 or v > 2000: continue
    if g["segs"] > 3: continue
    isaf = r["Protein"].startswith("AF-")
    node = f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}'
    if isaf and pl.get(node,0) < 70: continue
    if nseg.get(acc,0) >= 2: continue
    parts = pk.split("_")
    if len(parts) > 2:
        pv = vol_of.get((acc,"_".join(parts[:-1])))
        if pv and pv > 2000: continue
    base = fn.split("_res_P")[0]
    rec = (net1/"02_structures"/acc/f"{base}.pdb") if isaf else (net1/"03_pdbfixer"/acc/f"{base}.pdb")
    if not rec.is_file(): continue
    poc = net1/"04_pockets"/acc/fn
    d = box(poc)
    if d is None: continue
    elig.append({"r":r,"acc":acc,"fn":fn,"pk":pk,"af":isaf,"base":base,"rec":rec,"poc":poc,
                 "gene":r["Gene Name"],"protein":r["Protein"],
                 "boxv":d[0]*d[1]*d[2],"maxdim":max(d),"plddt":pl.get(node,""),
                 "desc":{k:num(r,k) for k in DESC}})
print(f"eligible: {len(elig)} pockets, {len({e['acc'] for e in elig})} proteins")
bv = sorted(e["boxv"] for e in elig)
print(f"eligible box volume: median {bv[len(bv)//2]:.0f}, p90 {bv[int(len(bv)*.9)]:.0f}")
af_share = sum(e["af"] for e in elig)/len(elig)
print(f"eligible AlphaFold share: {af_share:.0%}")

# proportional stratified sample over box-volume quintiles
elig.sort(key=lambda e: e["boxv"])
Q = 5
strata = [elig[i*len(elig)//Q:(i+1)*len(elig)//Q] for i in range(Q)]
sel = []
for s in strata:
    k = round(N/Q)
    n_af = round(k * sum(e["af"] for e in s)/len(s))
    afs = [e for e in s if e["af"]]; exps = [e for e in s if not e["af"]]
    pick = rng.sample(afs, min(n_af, len(afs))) + rng.sample(exps, min(k-n_af, len(exps)))
    while len(pick) < k:
        c = rng.choice(s)
        if c not in pick: pick.append(c)
    sel += pick[:k]
sel = sel[:N]
sv = sorted(e["boxv"] for e in sel)
print(f"\nselected {len(sel)}: AF {sum(e['af'] for e in sel)}, exp {sum(not e['af'] for e in sel)}")
print(f"selected box volume: median {sv[len(sv)//2]:.0f}, p90 {sv[int(len(sv)*.9)]:.0f}, "
      f"max {sv[-1]:.0f}   (eligible median {bv[len(bv)//2]:.0f})")
for k in ("volume","depth","hydrophobicity","siteAtms"):
    a=[e["desc"][k] for e in sel if e["desc"][k] is not None]
    b=[e["desc"][k] for e in elig if e["desc"][k] is not None]
    print(f"  {k:<15} selected {min(a):7.2f}-{max(a):7.2f} med {statistics.median(a):7.2f}"
          f" | eligible {min(b):7.2f}-{max(b):7.2f} med {statistics.median(b):7.2f}")

labels=[f'{e["protein"]}_{e["pk"]}' for e in sel]
dup=[l for l,c in collections.Counter(labels).items() if c>1]
if dup: sys.exit(f"duplicate labels: {dup}")

for sub in ("inputs","06_targets","logs"):
    (calib/sub).mkdir(parents=True, exist_ok=True)
rows, meta = [], []
for e in sel:
    acc, isaf = e["acc"], e["af"]
    rec_rel = (f"02_structures/{acc}/{e['base']}.pdb" if isaf
               else f"03_pdbfixer/{acc}/{e['base']}.pdb")
    poc_rel = f"04_pockets/{acc}/{e['fn']}"
    for rel, src in ((rec_rel, e["rec"]), (poc_rel, e["poc"])):
        dst = calib/rel; dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists(): shutil.copy2(src, dst)
    chain = e["protein"].split("_")[-1] if (not isaf and "_" in e["protein"]) else "A"
    rows.append({"Node": f'{acc}_{e["gene"]}_{e["protein"]}_{e["pk"]}', "UniProt_ID": acc,
                 "Gene_Name": e["gene"].replace("_HUMAN",""), "Organism":"HUMAN",
                 "Structure_ID": e["protein"] if isaf else e["protein"].split("_")[0],
                 "Similarity":"", "Protein": e["protein"], "Chain": chain,
                 "Pocket": e["pk"], "Pocket_Path": poc_rel,
                 "Receptor_Path": rec_rel, "Is_Hub": False})
    meta.append({"label": f'{e["protein"]}_{e["pk"]}', "UniProt": acc,
                 "Gene": e["gene"].replace("_HUMAN",""), "alphafold": isaf,
                 "box_volume": round(e["boxv"],1), "box_maxdim": round(e["maxdim"],2),
                 "pocket_plddt": e["plddt"], **{k:e["desc"][k] for k in DESC}})
cols=["Node","UniProt_ID","Gene_Name","Organism","Structure_ID","Similarity",
      "Protein","Chain","Pocket","Pocket_Path","Receptor_Path","Is_Hub"]
with open(calib/"06_targets/final_targets.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=cols); w.writeheader(); w.writerows(rows)
with open(calib/"calib_targets.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(meta[0].keys())); w.writeheader(); w.writerows(meta)
print(f"\nwrote {calib}/06_targets/final_targets.csv ({len(rows)} targets)")
print(f"wrote {calib}/calib_targets.csv (descriptors + box volume for the regression)")
