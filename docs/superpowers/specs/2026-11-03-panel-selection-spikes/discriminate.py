"""Does drugScore discriminate real drug-binding sites here, and does any geometric
descriptor add anything on top of it?

Evaluation set: experimental pockets labelled by whether a drug-like ligand is
co-crystallised in them. drugScore is DoGSiteScorer's trained druggability SVM, which is
the field-standard answer to this question; the geometric filters we invented sit on top
of it and have never been checked against anything.
"""
import csv, json, math, os, sys, statistics, collections
from pathlib import Path

run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
CACHE = Path(work) / "ligand_labels.json"
CONTACT, MIN_FRAC, MIN_HEAVY, MAX_HEAVY = 4.5, 0.70, 12, 60
SKIP = set("""HOH DOD UNX UNK NAG MAN BMA FUC GAL GLC NDG SIA XYS RAM BGC A2G
GOL EDO PEG PGE PG4 1PE 2PE P6G MPD DMS DMF TRS ACT ACY FMT OXL
SO4 PO4 NO3 CO3 EPE MES CIT TAR FLC IMD BME DTT TCE MLI
CL NA K MG CA ZN FE FE2 MN CU CU1 NI CD HG IOD BR CS SR BA CO PT AU AG
CLR POV PLM OLA MYR STE LDA C8E BOG DDQ LMT PC1 PEE PCW HEM HEC""".split())

ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if str(r.get("Selected","")).strip().lower() in ("true","1","yes")
      and not r["Protein"].startswith("AF-")]

def coords(path, het):
    out = collections.defaultdict(list)
    try:
        for ln in open(path, errors="replace"):
            if not ln.startswith("HETATM" if het else "ATOM"): continue
            if (ln[76:78].strip() or ln[12:14].strip()).upper().startswith("H"): continue
            resn = ln[17:20].strip()
            if het and resn in SKIP: continue
            try: c = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError: continue
            out[(resn, ln[21], ln[22:27].strip()) if het else None].append(c)
    except OSError: return {}
    return out

def rg(xyz):
    n=len(xyz); c=[sum(p[k] for p in xyz)/n for k in range(3)]
    return math.sqrt(sum(sum((p[k]-c[k])**2 for k in range(3)) for p in xyz)/n)

if CACHE.is_file():
    labels = json.load(open(CACHE))
    print(f"loaded {len(labels)} cached occupancy labels")
else:
    labels, cache = {}, {}
    for i, r in enumerate(ds):
        if i and i % 2000 == 0: print(f"  labelling {i}/{len(ds)}", flush=True)
        acc, fn = r["UniProt ID"], r["Pocket Filename"]
        stem = fn.split("_res_P")[0]
        if (acc, stem) not in cache:
            cache[(acc, stem)] = coords(Path(run)/"02_structures"/acc/f"{stem}.pdb", True)
        pocket = [c for v in coords(Path(run)/"04_pockets"/acc/fn, False).values() for c in v]
        if len(pocket) < 4: continue
        lig = None
        for (resn, _c, _n), lat in cache[(acc, stem)].items():
            if not (MIN_HEAVY <= len(lat) <= MAX_HEAVY): continue
            hit = sum(1 for a in lat if any(
                (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2 <= CONTACT*CONTACT for b in pocket))
            if hit/len(lat) >= MIN_FRAC: lig = resn; break
        labels[f"{acc}/{fn}"] = {"lig": lig, "rg": round(rg(pocket), 4)}
    json.dump(labels, open(CACHE, "w"))
    print(f"labelled {len(labels)} pockets")

def num(r, k):
    try: return float(r[k])
    except (TypeError, ValueError, KeyError): return None

DESC = ["drugScore", "simpleScore", "volume", "depth", "lid/hull", "surf/vol",
        "hydrophobicity", "apolarAA", "siteAtms", "accept", "donor", "aromat", "enclosure"]
rows = []
for r in ds:
    key = f'{r["UniProt ID"]}/{r["Pocket Filename"]}'
    lab = labels.get(key)
    if not lab: continue
    d = {k: num(r, k) for k in DESC if k in r}
    d["rg"] = lab["rg"]
    d["y"] = 1 if lab["lig"] else 0
    rows.append(d)
pos = sum(r["y"] for r in rows)
print(f"\nevaluation set: {len(rows)} pockets, {pos} ligand-bound ({pos/len(rows):.1%})")

def auc(vals_y):
    """Rank-based AUC (Mann-Whitney). >0.5 means higher values favour ligand-bound."""
    v = [(x, y) for x, y in vals_y if x is not None]
    if not v: return None
    v.sort(key=lambda t: t[0])
    ranks, i = {}, 0
    while i < len(v):
        j = i
        while j+1 < len(v) and v[j+1][0] == v[i][0]: j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j+1): ranks[k] = r
        i = j + 1
    npos = sum(y for _, y in v); nneg = len(v) - npos
    if not npos or not nneg: return None
    s = sum(ranks[k] for k, (_, y) in enumerate(v) if y)
    return (s - npos*(npos+1)/2) / (npos*nneg)

print(f"\n{'descriptor':<16} {'AUC (all)':>10} {'AUC | drugScore>=0.5':>22}")
hi = [r for r in rows if (r.get("drugScore") or 0) >= 0.5]
scored = []
for k in ["drugScore"] + [d for d in DESC[1:] if any(d in r for r in rows)] + ["rg"]:
    a = auc([(r.get(k), r["y"]) for r in rows])
    b = auc([(r.get(k), r["y"]) for r in hi])
    if a is None: continue
    scored.append((abs(a-0.5), k, a, b))
    print(f"{k:<16} {a:>10.3f} {('%.3f' % b) if b else '     n/a':>22}")
print(f"\n(AUC 0.5 = no discrimination; <0.5 means LOWER values favour ligand-bound)")
print(f"ranked by |AUC-0.5|: " + ", ".join(k for _, k, _, _ in sorted(scored, reverse=True)[:6]))

print(f"\ndrugScore thresholds on this set:")
print(f"{'cut':>6} {'sensitivity':>12} {'specificity':>12} {'precision':>10} {'pockets kept':>13}")
for cut in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
    tp = sum(1 for r in rows if (r.get("drugScore") or 0) >= cut and r["y"])
    fp = sum(1 for r in rows if (r.get("drugScore") or 0) >= cut and not r["y"])
    fn_ = pos - tp; tn = len(rows) - pos - fp
    print(f"{cut:>6.1f} {tp/max(pos,1):>11.1%} {tn/max(tn+fp,1):>11.1%} "
          f"{tp/max(tp+fp,1):>9.1%} {tp+fp:>13}")
