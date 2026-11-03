"""Two audit findings that could invalidate the current panel.

B1: the AF-vs-experimental pass-rate check that condemned `segments<=3` (21.9-point gap)
    and vindicated Rg (1.5) was NEVER run on surf/vol or the volume range. The AF share
    fell 53% -> 32% under them, which is exactly what a class-differential filter looks
    like.
A7: nothing checks that the pocket the panel picks (highest drugScore) is the protein's
    real ligand site. The same failure killed the cluster layer. We have 687 pockets
    labelled by a co-crystal ligand, so this is directly measurable.
"""
import csv, json, os, sys, collections
sys.path.insert(0, "/media/ilayda/Elements/adep_csb/repos/bioinf_polyphar_workflow")
from polyphar.stages import annotate
from polyphar import panel

run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
allrows = list(csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline="")))
rows = [r for r in allrows
        if str(r.get("Selected","")).strip().lower() in ("true","1","yes")]
vols = {}
for r in allrows:
    try: vols[(r["UniProt ID"], r["Protein"], r["Pocket"])] = float(r["volume"])
    except (KeyError, TypeError, ValueError): pass
labels = json.load(open(f"{work}/ligand_labels.json"))
tm = {a: [tuple(s) for s in v] for a, v in json.load(open(f"{work}/tm_ranges.json")).items()}

anns = {}
for i, r in enumerate(rows):
    pv = annotate.parent_volume_of(r["UniProt ID"], r.get("Protein",""),
                                   r.get("Pocket",""), vols)
    try:
        a = annotate.annotate_one(r, f"{run}/04_pockets", f"{run}/02_structures", pv)
        a["_residues"] = annotate.decode_residues(a["residues"])
        anns[r["Pocket Filename"]] = a
    except Exception:
        pass

def is_af(r): return r["Protein"].startswith("AF-")
def num(r,k):
    try: return float(r[k])
    except (TypeError,ValueError,KeyError): return None

# --- B1: class-differential pass rates, each filter in isolation --------------------
print("B1. AF vs EXPERIMENTAL PASS RATE, each filter alone, over all selected pockets")
print(f"{'filter':<28} {'AF pass':>9} {'exp pass':>10} {'gap (pts)':>10}")
af_all = [r for r in rows if is_af(r)]; ex_all = [r for r in rows if not is_af(r)]
def rate(sub, pred): 
    ok = [r for r in sub if pred(r)]
    return len(ok)/len(sub) if sub else 0.0
tests = [
    ("volume 300-1500", lambda r: (num(r,"volume") or -1) >= 300 and (num(r,"volume") or 1e9) <= 1500),
    ("volume <= 1500 (ceiling)", lambda r: (num(r,"volume") or 1e9) <= 1500),
    ("volume >= 300 (floor)", lambda r: (num(r,"volume") or -1) >= 300),
    ("surf/vol <= 0.45", lambda r: (num(r,"surf/vol") or 9) <= 0.45),
    ("drugScore >= 0.6", lambda r: (num(r,"drugScore") or 0) >= 0.6),
]
for name, pred in tests:
    a, e = rate(af_all, pred), rate(ex_all, pred)
    print(f"{name:<28} {a:>8.1%} {e:>9.1%} {abs(a-e)*100:>9.1f}")
print(f"  (for reference: segments<=3 gave a 21.9-point gap and was rejected;")
print(f"   Rg<=8.5 gave 1.5 and was accepted on that basis)")

# --- A7: does the panel pick the pocket that actually binds a ligand? ---------------
pan, elig, funnel, rej = panel.select_panel(rows, anns, tm, vols, panel.Options())
by_acc = collections.defaultdict(list)
for r, a in elig: by_acc[r["UniProt ID"]].append(r)

checked = agree = 0
no_label = 0
for acc, (rep, _a) in pan.items():
    if is_af(rep): continue
    cands = by_acc[acc]
    lab = [r for r in cands
           if labels.get(f'{acc}/{r["Pocket Filename"]}', {}).get("lig")]
    if not lab:
        no_label += 1
        continue
    checked += 1
    if any(r["Pocket Filename"] == rep["Pocket Filename"] for r in lab):
        agree += 1
print(f"\nA7. IS THE PANEL'S CHOSEN POCKET THE KNOWN LIGAND SITE?")
print(f"  panel: {len(pan)} proteins ({sum(1 for r,_ in pan.values() if not is_af(r))} experimental)")
print(f"  experimental panel proteins with >=1 ligand-labelled eligible pocket: {checked}")
print(f"  of those, the drugScore pick IS a labelled pocket: {agree} "
      f"({agree/checked:.0%})" if checked else "  none to check")
print(f"  experimental panel proteins with no labelled eligible pocket: {no_label}")
# how often is a labelled pocket even eligible?
lab_elig = sum(1 for r, _ in elig
               if labels.get(f'{r["UniProt ID"]}/{r["Pocket Filename"]}', {}).get("lig"))
lab_total = sum(1 for k, v in labels.items() if v.get("lig"))
print(f"  ligand-labelled pockets surviving the filters: {lab_elig}/{lab_total} "
      f"({lab_elig/lab_total:.0%})")
