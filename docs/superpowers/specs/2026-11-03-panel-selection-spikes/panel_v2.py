"""Panel under the common-sense definition: trust DoGSiteScorer's trained druggability
score, add an enclosure criterion (surf/vol) and a literature-anchored volume range, plus
the TM rule. No Rg, no sequence segments, no box cap — all three were invented thresholds
and one was measurably backwards.
"""
import csv, json, os, sys, statistics, collections
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
ds = list(csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline="")))
sel = [r for r in ds if str(r.get("Selected","")).strip().lower() in ("true","1","yes")]
tmr = json.load(open(f"{work}/tm_ranges.json")); nseg = {a: len(v) for a,v in tmr.items()}
labels = json.load(open(f"{work}/ligand_labels.json"))
def num(r,k):
    try: return float(r[k])
    except (TypeError,ValueError,KeyError): return None
vol_of = {(r["UniProt ID"], r["Protein"], r["Pocket"]): num(r,"volume") for r in ds}

deg = collections.Counter()
for r in csv.DictReader(open(f"{run}/06_targets/unique_pairs.csv", newline="")):
    deg[r["Node"]] += 1; deg[r["Node_2"]] += 1
top = {n for n,_ in deg.most_common(100)}
def node(r): return f'{r["UniProt ID"]}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'

def build(sv_max, vol_max, vol_min=300.0, dsmin=0.5):
    steps = collections.Counter(); keep = []
    for r in sel:
        steps["0 selected"] += 1
        d, v, sv = num(r,"drugScore"), num(r,"volume"), num(r,"surf/vol")
        if d is None or d < dsmin: continue
        steps["1 drugScore"] += 1
        if v is None or not (vol_min <= v <= vol_max): continue
        steps["2 volume"] += 1
        p = r["Pocket"].split("_")
        if len(p) > 2:
            pv = vol_of.get((r["UniProt ID"], r["Protein"], "_".join(p[:-1])))
            if pv and pv > vol_max: continue
        steps["3 subpocket parent"] += 1
        if sv is None or sv > sv_max: continue
        steps["4 surf/vol"] += 1
        if nseg.get(r["UniProt ID"], 0) >= 2: continue
        steps["5 TM multi-pass"] += 1
        keep.append(r)
    return steps, keep

print(f"{'surf/vol':>9} {'vol max':>8} {'pockets':>8} {'proteins':>9} {'AF%':>5} "
      f"{'artifacts':>10} {'ligand-bound kept':>18}")
occ_all = sum(1 for r in sel if labels.get(f'{r["UniProt ID"]}/{r["Pocket Filename"]}',{}).get("lig"))
for sv_max in (0.45, 0.50, 0.55):
    for vol_max in (1000.0, 1500.0, 2000.0):
        steps, keep = build(sv_max, vol_max)
        prot = {r["UniProt ID"] for r in keep}
        af = {r["UniProt ID"] for r in keep if r["Protein"].startswith("AF-")}
        art = sum(1 for r in keep if node(r) in top)
        occ = sum(1 for r in keep
                  if labels.get(f'{r["UniProt ID"]}/{r["Pocket Filename"]}',{}).get("lig"))
        print(f"{sv_max:>9.2f} {vol_max:>8.0f} {len(keep):>8} {len(prot):>9} "
              f"{len(af)/max(len(prot),1)*100:>4.0f}% {f'{art}/100':>10} "
              f"{f'{occ}/{occ_all} ({occ/occ_all:.0%})':>18}")

print("\nFUNNEL at surf/vol<=0.50, volume 300-2000:")
steps, keep = build(0.50, 2000.0)
for k in sorted(steps): print(f"  {k:<22} {steps[k]:>6}")
best = {}
for r in sorted(keep, key=lambda x: (-(num(x,"drugScore") or 0), x["Pocket"])):
    best.setdefault(r["UniProt ID"], r)
af = sum(1 for r in best.values() if r["Protein"].startswith("AF-"))
print(f"\nPANEL: {len(best)} proteins ({af} AF = {af/len(best):.0%}, "
      f"{len(best)-af} experimental), {len(keep)} eligible pockets")
print(f"  median volume {statistics.median([num(r,'volume') for r in keep]):.0f} A^3, "
      f"median surf/vol {statistics.median([num(r,'surf/vol') for r in keep]):.3f}, "
      f"median drugScore {statistics.median([num(r,'drugScore') for r in keep]):.2f}")
