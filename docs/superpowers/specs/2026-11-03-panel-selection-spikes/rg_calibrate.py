"""Give the Rg cutoff an empirical basis: measure it on pockets that are KNOWN to bind
a drug-like ligand, because a co-crystallised ligand sits in them.

No published cutoff exists for pocket radius of gyration — Rg is used as a pocket size
descriptor in the literature, but with no standard definition (which atoms, which pocket
boundary, what normalisation) and no druggability threshold. So rather than borrow a
number, measure the Rg distribution of experimentally observed drug-binding sites in this
very dataset and set the cutoff to cover a stated fraction of them.

Important: this is computed over ALL selected experimental pockets, BEFORE any Rg
filtering, or we would be conditioning on the filter being calibrated.
"""
import csv, math, os, sys, statistics, collections
from pathlib import Path

run = sys.argv[1]
CONTACT = 4.5      # A, ligand atom to pocket lining atom
MIN_FRAC = 0.70    # fraction of ligand heavy atoms that must be in contact
MIN_HEAVY, MAX_HEAVY = 12, 60   # drug-like size range

# water, ions, buffers, cryoprotectants, glycans, lipids and detergents: present in
# crystals but not evidence of a drug-binding site
SKIP = set("""HOH DOD UNX UNK
NAG MAN BMA FUC GAL GLC NDG SIA XYS RAM BGC A2G
GOL EDO PEG PGE PG4 1PE 2PE P6G MPD DMS DMF TRS ACT ACY FMT OXL
SO4 PO4 NO3 CO3 EPE MES CIT TAR FLC IMD BME DTT TCE MLI
CL NA K MG CA ZN FE FE2 MN CU CU1 NI CD HG IOD BR CS SR BA CO PT AU AG
CLR POV PLM OLA MYR STE LDA C8E BOG DDQ LMT PC1 PEE PCW HEM HEC""".split())

ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if str(r.get("Selected","")).strip().lower() in ("true","1","yes")
      and not r["Protein"].startswith("AF-")]
print(f"selected experimental pockets: {len(ds)}")

def atoms(path, hetatm):
    out = collections.defaultdict(list)
    try:
        with open(path, errors="replace") as fh:
            for ln in fh:
                tag = "HETATM" if hetatm else "ATOM"
                if not ln.startswith(tag):
                    continue
                el = (ln[76:78].strip() or ln[12:14].strip())
                if el.upper().startswith("H"):
                    continue
                resn = ln[17:20].strip()
                if hetatm and resn in SKIP:
                    continue
                key = (resn, ln[21], ln[22:27].strip()) if hetatm else None
                try:
                    c = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                except ValueError:
                    continue
                out[key].append(c)
    except OSError:
        return {}
    return out

def rg(xyz):
    n = len(xyz)
    c = [sum(p[k] for p in xyz)/n for k in range(3)]
    return math.sqrt(sum(sum((p[k]-c[k])**2 for k in range(3)) for p in xyz)/n)

ligand_cache = {}
rows = []
for i, r in enumerate(ds):
    if i and i % 2000 == 0: print(f"  {i}/{len(ds)}", flush=True)
    acc, fn = r["UniProt ID"], r["Pocket Filename"]
    stem = fn.split("_res_P")[0]
    chain_pdb = Path(run)/"02_structures"/acc/f"{stem}.pdb"
    if (acc, stem) not in ligand_cache:
        ligand_cache[(acc, stem)] = atoms(chain_pdb, hetatm=True)
    ligs = ligand_cache[(acc, stem)]
    pocket = [c for v in atoms(Path(run)/"04_pockets"/acc/fn, hetatm=False).values() for c in v]
    if len(pocket) < 4:
        continue
    g = rg(pocket)
    occupied_by = None
    for (resn, _ch, _num), lat in ligs.items():
        if not (MIN_HEAVY <= len(lat) <= MAX_HEAVY):
            continue
        hit = 0
        for a in lat:
            for b in pocket:
                if ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) <= CONTACT*CONTACT:
                    hit += 1
                    break
        if hit / len(lat) >= MIN_FRAC:
            occupied_by = resn
            break
    rows.append({"acc": acc, "fn": fn, "rg": g, "ligand": occupied_by})

occ = [r for r in rows if r["ligand"]]
emp = [r for r in rows if not r["ligand"]]
print(f"\nmeasured {len(rows)} pockets: {len(occ)} hold a drug-like co-crystal ligand, "
      f"{len(emp)} do not")
print(f"distinct ligands: {len({r['ligand'] for r in occ})}; most common: "
      f"{collections.Counter(r['ligand'] for r in occ).most_common(6)}")

def q(vals, p):
    v = sorted(vals); return v[min(int(p*len(v)), len(v)-1)]
o = [r["rg"] for r in occ]
e = [r["rg"] for r in emp]
print(f"\nRg of LIGAND-OCCUPIED pockets (n={len(o)}):")
print(f"  median {statistics.median(o):.2f}  p75 {q(o,.75):.2f}  p90 {q(o,.90):.2f}  "
      f"p95 {q(o,.95):.2f}  p99 {q(o,.99):.2f}  max {max(o):.2f}")
print(f"Rg of pockets with NO such ligand (n={len(e)}):")
print(f"  median {statistics.median(e):.2f}  p75 {q(e,.75):.2f}  p95 {q(e,.95):.2f}")

print(f"\n{'cutoff':>8} {'keeps ligand-bound':>20} {'keeps unoccupied':>18} {'enrichment':>11}")
for cut in (6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.0):
    a = sum(1 for x in o if x <= cut)/len(o)
    b = sum(1 for x in e if x <= cut)/len(e)
    print(f"{cut:>8.1f} {a:>19.1%} {b:>17.1%} {a/b if b else float('nan'):>10.2f}x")
print(f"\ncutoff covering 95% of ligand-bound pockets: Rg <= {q(o,.95):.2f} A")
print(f"cutoff covering 90%: Rg <= {q(o,.90):.2f} A")
