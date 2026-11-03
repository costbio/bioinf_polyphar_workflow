"""Pre-stage-1 measurement: can experimental pocket residues actually be mapped to
UniProt? Two questions:
  A. Is the cached SIFTS flat file present, and does it carry SP_BEG/SP_END (the offset
     that Protein_info.csv discards)?
  B. For the panel's experimental chains, how many fall back to a seqres-derived range
     ('SIFTS segments source' = seqres), where the offset is unrecoverable by any route?
"""
import csv, gzip, json, math, os, sys, collections
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
MAXAX, MAXTOT, SPACING, BOXCAP = 128, 531441, 0.375, 30.0

# ---- A. the cached SIFTS file -------------------------------------------------------
sif = f"{run}/02_structures/sifts/pdb_chain_uniprot.tsv.gz"
print("A. CACHED SIFTS FLAT FILE")
print(f"   path exists: {os.path.exists(sif)}")
if os.path.exists(sif):
    print(f"   size: {os.path.getsize(sif)/1e6:.1f} MB")
    with gzip.open(sif, "rt") as fh:
        head = [next(fh) for _ in range(3)]
    for l in head: print("   " + l.rstrip()[:160])
    cols = head[1].rstrip().split("\t") if head[1].startswith("PDB") else head[-1].rstrip().split("\t")
    print(f"   has SP_BEG/SP_END: {('SP_BEG' in cols) and ('SP_END' in cols)}")

# ---- rebuild the panel (same filters as reconcile.py) -------------------------------
ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl  = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
tmr = json.load(open(f"{work}/tm_ranges.json")); nseg = {a: len(v) for a,v in tmr.items()}
def num(r,k):
    try: return float(r[k])
    except (ValueError,KeyError,TypeError): return None
vol_of = {(r["UniProt ID"], r["Pocket"]): num(r,"volume") for r in ds}
def res_box(acc, fn):
    res, lo, hi = set(), None, None
    try:
        for ln in open(f"{run}/04_pockets/{acc}/{fn}"):
            if ln.startswith(("ATOM","HETATM")):
                try: res.add(int(ln[22:26]))
                except ValueError: pass
                c=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
                if lo is None: lo,hi=list(c),list(c)
                else:
                    for k in range(3): lo[k]=min(lo[k],c[k]); hi[k]=max(hi[k],c[k])
    except OSError: return None,None
    return (res,[hi[k]-lo[k] for k in range(3)]) if lo else (None,None)

elig = []
for r in ds:
    acc, fn, pk = r["UniProt ID"], r["Pocket Filename"], r["Pocket"]
    g = geo.get(f"{acc}/{fn}")
    if g is None: continue
    v, d = num(r,"volume"), num(r,"drugScore")
    if v is None or d is None or d<0.5 or not (300<=v<=2000) or g["segs"]>3: continue
    p = pk.split("_")
    if len(p)>2:
        pv = vol_of.get((acc,"_".join(p[:-1])))
        if pv and pv>2000: continue
    isaf = r["Protein"].startswith("AF-")
    if isaf and pl.get(f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{pk}',0)<70: continue
    if nseg.get(acc,0)>=2: continue
    res, dims = res_box(acc, fn)
    if res is None: continue
    if isaf and nseg.get(acc,0)>=1:
        t={x for a,b in tmr.get(acc,[]) for x in range(a,b+1)}
        if res and len(res&t)/len(res)>0.10: continue
    if max(dims)>BOXCAP: continue
    n=[math.floor(s/SPACING)+1 for s in dims]
    if max(n)>MAXAX or n[0]*n[1]*n[2]>MAXTOT: continue
    elig.append({"acc":acc,"fn":fn,"pk":pk,"af":isaf,"ds":d,
                 "prot":r["Protein"],"gene":r["Gene Name"].replace("_HUMAN","")})
best = {}
for e in sorted(elig, key=lambda x:(-x["ds"], x["pk"])): best.setdefault(e["acc"], e)
print(f"\n   panel rebuilt: {len(best)} proteins, {len(elig)} eligible pockets")

# ---- B. segment source for the panel's experimental chains --------------------------
info = list(csv.DictReader(open(f"{run}/02_structures/Protein_info.csv", newline="")))
key = {}
for r in info:
    if r.get("AlphaFold","").lower() == "true": continue
    base = os.path.basename(r.get("File","")).replace(".pdb","")
    key[(r["UniProt ID"], base.lower())] = r
print("\nB. SEGMENT SOURCE FOR EXPERIMENTAL CHAINS")
print(f"   'SIFTS segments source' values over all {len(info)} Protein_info rows:")
for k,v in collections.Counter(r.get("SIFTS segments source","") for r in info).most_common():
    print(f"     {k or '(empty)':<20} {v}")

def lookup(e):
    base = e["fn"].split("_res_P")[0].lower()
    return key.get((e["acc"], base))

exp_pockets = [e for e in elig if not e["af"]]
exp_panel   = {a:b for a,b in best.items() if not b["af"]}
for label, items in (("eligible experimental pockets", exp_pockets),
                     ("panel rows that are experimental", list(exp_panel.values()))):
    srcs, miss = collections.Counter(), 0
    for e in items:
        row = lookup(e)
        if row is None: miss += 1; continue
        srcs[row.get("SIFTS segments source","") or "(empty)"] += 1
    tot = sum(srcs.values())
    unver = srcs.get("seqres", 0)
    print(f"\n   {label}: {len(items)} (unmatched to Protein_info: {miss})")
    for k,v in srcs.most_common(): print(f"     {k:<20} {v}")
    if tot:
        print(f"   -> UNMAPPABLE (seqres fallback): {unver}/{tot} = {unver/tot:.0%}")
