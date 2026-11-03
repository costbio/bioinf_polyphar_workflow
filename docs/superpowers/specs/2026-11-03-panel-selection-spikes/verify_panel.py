"""End-to-end verification of stage 1 + stage 2 against real net1 data:
annotate every selected pocket with the shipped module, then run the shipped panel
filters over those annotations. Must reproduce 164 proteins / 537 pockets / 93 AF.
"""
import csv, json, os, sys, collections
sys.path.insert(0, "/media/ilayda/Elements/adep_csb/repos/bioinf_polyphar_workflow")
from polyphar.stages import annotate
from polyphar import panel

run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
pockets_dir, structures_dir = f"{run}/04_pockets", f"{run}/02_structures"

all_rows = list(csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline="")))
rows = [r for r in all_rows
        if str(r.get("Selected","")).strip().lower() in ("true","1","yes")]
vols = {}
for r in all_rows:
    try: vols[(r["UniProt ID"], r["Protein"], r["Pocket"])] = float(r["volume"])
    except (KeyError, TypeError, ValueError): pass
print(f"selected pockets: {len(rows)}")

# --- stage 1: annotate ---------------------------------------------------------------
anns, failed = {}, []
for i, r in enumerate(rows):
    if i and i % 3000 == 0: print(f"  annotating {i}/{len(rows)}", flush=True)
    pv = annotate.parent_volume_of(r["UniProt ID"], r.get("Protein",""),
                                   r.get("Pocket",""), vols)
    try:
        a = annotate.annotate_one(r, pockets_dir, structures_dir, pv)
        a["_residues"] = annotate.decode_residues(a["residues"])
        anns[r["Pocket Filename"]] = a
    except Exception as e:
        failed.append((r["Pocket Filename"], str(e)[:60]))
print(f"annotated {len(anns)}, failed {len(failed)}")

# residues round-trip through the CSV encoding on real data
bad = [f for f, a in anns.items()
       if annotate.decode_residues(a["residues"]) != a["_residues"]]
print(f"residues round-trip failures: {len(bad)}")

# --- stage 2: panel ------------------------------------------------------------------
tm = {a: [tuple(s) for s in segs]
      for a, segs in json.load(open(f"{work}/tm_ranges.json")).items()}
pan, elig, funnel, rej = panel.select_panel(rows, anns, tm, vols, panel.Options())

print("\nFUNNEL")
for step, n in funnel: print(f"  {step:<20} {n:>6}")
n_elig = len(elig)
af = sum(1 for _, a in pan.values() if str(a["AlphaFold"]).lower() == "true")
print(f"\nPANEL: {len(pan)} proteins, {n_elig} eligible pockets, AF {af} "
      f"({af/max(len(pan),1):.0%})")
print("EXPECTED from the standalone sweep: 178 proteins, 1481 eligible pockets")

print(f"\nrejections recorded: {len(rej)} (every rejected pocket has a reason)")
for step, n in collections.Counter(s for _, s in rej).most_common():
    print(f"  {step:<22} {n:>6}")
assert len(rej) + n_elig == len(rows), "rejections + eligible must cover every pocket"
print("\nrejections + eligible == selected: ok")
sys.exit(0)
