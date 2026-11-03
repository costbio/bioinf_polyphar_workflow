"""Build the pocket-geometry / ligand-size CALIBRATION run folder.
48 eligible net1 pockets by maximin (space-filling) design over standardised
[volume, depth, lid/hull, hydrophobicity, siteAtms], balanced 24 AF / 24 experimental,
plus a 250-ligand e-Drug3D sample stratified by heavy-atom count. Self-contained: only
the chosen receptors and pocket files are copied in, so no symlinks or extra mounts.
"""
import csv, sys, json, os, math, shutil, statistics, collections
from pathlib import Path
net1, calib, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(os.path.expanduser("~/scratch_clust"))
N_AF = N_EXP = 24
N_LIG = 250

ds = [r for r in csv.DictReader(open(net1/"04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(work/"geo.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(work/"plddt_rows.json"))}
nseg = {a: len(v) for a, v in json.load(open(work/"tm_ranges.json")).items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}

DESC = ["volume","depth","lid/hull","hydrophobicity","siteAtms"]
elig = []
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v = num(r,"volume")
    if v is None or num(r,"drugScore") is None: continue
    if num(r,"drugScore") < 0.5 or v < 300 or v > 2000: continue
    if g["segs"] > 3: continue
    isaf = r["Protein"].startswith("AF-")
    node = f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}'
    if isaf and pl.get(node,0) < 70: continue
    if nseg.get(acc,0) >= 2: continue
    parts = pk.split("_")
    if len(parts) > 2:
        pv = vol_of.get((acc,"_".join(parts[:-1])))
        if pv and pv > 2000: continue
    d = [num(r,k) for k in DESC]
    if any(x is None for x in d): continue
    base = fn.split("_res_P")[0]
    rec = (net1/"02_structures"/acc/f"{base}.pdb") if isaf else (net1/"03_pdbfixer"/acc/f"{base}.pdb")
    if not rec.is_file(): continue
    elig.append({"r":r,"acc":acc,"fn":fn,"pk":pk,"af":isaf,"desc":d,"base":base,
                 "rec":rec,"poc":net1/"04_pockets"/acc/fn,
                 "gene":r["Gene Name"],"protein":r["Protein"]})
print(f"eligible with all descriptors and a receptor on disk: {len(elig)}")

mu = [statistics.mean([e["desc"][i] for e in elig]) for i in range(len(DESC))]
sd = [statistics.pstdev([e["desc"][i] for e in elig]) or 1.0 for i in range(len(DESC))]
for e in elig:
    e["z"] = [(x-mu[i])/sd[i] for i,x in enumerate(e["desc"])]

def maximin(pool, k):
    if len(pool) <= k: return list(pool)
    start = max(pool, key=lambda e: sum(x*x for x in e["z"]))
    chosen = [start]
    while len(chosen) < k:
        best, bd = None, -1
        for e in pool:
            if e in chosen: continue
            d = min(math.dist(e["z"], c["z"]) for c in chosen)
            if d > bd: best, bd = e, d
        chosen.append(best)
    return chosen

sel = maximin([e for e in elig if e["af"]], N_AF) + \
      maximin([e for e in elig if not e["af"]], N_EXP)
print(f"selected {len(sel)}: AF {sum(e['af'] for e in sel)}, "
      f"experimental {sum(not e['af'] for e in sel)}")
for i,k in enumerate(DESC):
    vs = [e["desc"][i] for e in sel]
    al = [e["desc"][i] for e in elig]
    print(f"  {k:<15} selected {min(vs):8.2f}-{max(vs):8.2f} (median {statistics.median(vs):8.2f})"
          f"   eligible {min(al):8.2f}-{max(al):8.2f}")

labels = [f'{e["protein"]}_{e["pk"]}' for e in sel]
dup = [l for l,c in collections.Counter(labels).items() if c > 1]
if dup: sys.exit(f"duplicate dock labels: {dup}")

for sub in ("inputs","06_targets","logs"):
    (calib/sub).mkdir(parents=True, exist_ok=True)
rows = []
for e in sel:
    acc, isaf = e["acc"], e["af"]
    rec_rel = (f"02_structures/{acc}/{e['base']}.pdb" if isaf
               else f"03_pdbfixer/{acc}/{e['base']}.pdb")
    poc_rel = f"04_pockets/{acc}/{e['fn']}"
    for rel, src in ((rec_rel, e["rec"]), (poc_rel, e["poc"])):
        dst = calib/rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists(): shutil.copy2(src, dst)
    chain = e["protein"].split("_")[-1] if (not isaf and "_" in e["protein"]) else "A"
    rows.append({"Node": f'{acc}_{e["gene"]}_{e["protein"]}_{e["pk"]}',
                 "UniProt_ID": acc, "Gene_Name": e["gene"].replace("_HUMAN",""),
                 "Organism": "HUMAN",
                 "Structure_ID": e["protein"] if isaf else e["protein"].split("_")[0],
                 "Similarity": "", "Protein": e["protein"], "Chain": chain,
                 "Pocket": e["pk"], "Pocket_Path": poc_rel,
                 "Receptor_Path": rec_rel, "Is_Hub": False})
cols = ["Node","UniProt_ID","Gene_Name","Organism","Structure_ID","Similarity",
        "Protein","Chain","Pocket","Pocket_Path","Receptor_Path","Is_Hub"]
with open(calib/"06_targets/final_targets.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
print(f"wrote {calib}/06_targets/final_targets.csv with {len(rows)} targets")

# ---- ligand sample, stratified by heavy-atom count ----------------------------------
src_sdf = net1/"inputs/e-Drug3D_2083.sdf"
blocks, cur = [], []
for line in open(src_sdf, errors="replace"):
    if line.startswith("$$$$"):
        blocks.append(cur); cur = []
    else:
        cur.append(line)
if cur: blocks.append(cur)
mols = []
for b in blocks:
    if len(b) < 4: continue
    try: na = int(b[3][0:3])
    except ValueError: continue
    heavy = 0
    for ln in b[4:4+na]:
        el = ln[31:34].strip() if len(ln) > 34 else ""
        if el and el != "H": heavy += 1
    has_id = any(l.strip().startswith("> <ID>") for l in b)
    if heavy and has_id: mols.append((heavy, b))
mols.sort(key=lambda x: x[0])
print(f"library records usable: {len(mols)} of {len(blocks)}")
step = len(mols)/N_LIG
pick = [mols[min(int(i*step), len(mols)-1)] for i in range(N_LIG)]
seen, uniq = set(), []
for h, b in pick:
    key = id(b)
    if key in seen: continue
    seen.add(key); uniq.append((h,b))
i = 0
while len(uniq) < N_LIG and i < len(mols):
    if id(mols[i][1]) not in seen:
        seen.add(id(mols[i][1])); uniq.append(mols[i])
    i += 1
out = calib/"inputs/calib_250.sdf"
with open(out,"w") as fh:
    for h, b in uniq:
        fh.writelines(b); fh.write("$$$$\n")
hv = [h for h,_ in uniq]
print(f"wrote {out} with {len(uniq)} ligands; heavy atoms "
      f"{min(hv)}-{max(hv)} (median {statistics.median(hv)})")
est = len(rows)*len(uniq)*57/64/3600
print(f"\nestimated docking: {len(rows)} targets x {len(uniq)} ligands "
      f"= {len(rows)*len(uniq)} dockings, ~{est:.1f} h on 64 cores")
