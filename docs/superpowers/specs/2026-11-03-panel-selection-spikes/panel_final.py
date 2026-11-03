"""THROWAWAY: final panel definition + sweep of the volume ceiling and segment cap.
Base filters: drugScore>=0.5, volume>=300, drop multi-pass TM, AF pockets need
mean pocket pLDDT>=70. One best pocket per protein. Reports panel size/composition
and how many of the flat-network artifact hubs survive."""
import csv, sys, json, os, statistics, collections
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")

ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
print(f"selected pockets: {len(ds)}")

geo_path = f"{work}/geo.json"
if os.path.exists(geo_path):
    geo = json.load(open(geo_path))
else:
    geo = {}
    for i, r in enumerate(ds):
        if i and i % 2500 == 0: print(f"  geometry {i}/{len(ds)}", flush=True)
        acc = r["UniProt ID"]; key = f'{acc}/{r["Pocket Filename"]}'
        res = set()
        try:
            with open(f"{run}/04_pockets/{key}") as fh:
                for ln in fh:
                    if ln.startswith("ATOM"):
                        try: res.add(int(ln[22:26]))
                        except ValueError: pass
        except OSError:
            continue
        if not res: continue
        s = sorted(res); segs, prev = 1, s[0]
        for x in s[1:]:
            if x - prev > 5: segs += 1
            prev = x
        geo[key] = {"nres": len(s), "segs": segs}
    json.dump(geo, open(geo_path, "w"))
print(f"geometry for {len(geo)} pockets")

tm = json.load(open(f"{work}/tm_ann.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
deg = collections.Counter()
for r in csv.DictReader(open(f"{run}/06_targets/unique_pairs.csv", newline="")):
    deg[r["Node"]] += 1; deg[r["Node_2"]] += 1
artifact = {n for n, _ in deg.most_common(100)}
orig75 = {r["UniProt_ID"] for r in
          csv.DictReader(open(f"{run}/06_targets/final_targets.csv", newline=""))}
MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}

def num(r, k):
    try: return float(r[k])
    except (ValueError, KeyError): return 0.0

def build(volmax, segmax, plddt_min=70, dsmin=0.5, volmin=300, tmdrop=2):
    keep = []
    for r in ds:
        acc = r["UniProt ID"]; key = f'{acc}/{r["Pocket Filename"]}'
        g = geo.get(key)
        if g is None: continue
        if num(r,"drugScore") < dsmin or num(r,"volume") < volmin: continue
        if volmax and num(r,"volume") > volmax: continue
        if segmax and g["segs"] > segmax: continue
        if tm.get(acc, {}).get("tm", 0) >= tmdrop: continue
        node = f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
        if r["Protein"].startswith("AF-"):
            if pl.get(node, 0) < plddt_min: continue
        keep.append((r, node, g))
    best = {}
    for r, node, g in keep:
        a = r["UniProt ID"]
        if a not in best or num(r,"drugScore") > num(best[a][0],"drugScore"):
            best[a] = (r, node, g)
    return keep, best

print(f"\n{'volMax':>7} {'segMax':>7} {'pockets':>8} {'PANEL':>6} {'AF%':>5} {'TM%':>5} "
      f"{'mito':>5} {'of75':>5} {'artifacts':>10} {'dock(d)':>8}")
for volmax in (1500, 2000, 3000, None):
    for segmax in (3, 4, 5, None):
        keep, best = build(volmax, segmax)
        if not best: continue
        panel = set(best)
        af = sum(1 for a in panel if best[a][0]["Protein"].startswith("AF-"))
        tmc = sum(1 for a in panel if tm.get(a, {}).get("tm", 0) >= 1)
        genes = {best[a][0]["Gene Name"].replace("_HUMAN","") for a in panel}
        surv = sum(1 for _, n, _ in keep if n in artifact)
        days = len(panel) * 2083 * 57 / 64 / 86400
        print(f"{str(volmax):>7} {str(segmax):>7} {len(keep):>8} {len(panel):>6} "
              f"{af/len(panel)*100:>4.0f}% {tmc/len(panel)*100:>4.0f}% "
              f"{len(genes & MITOTIC):>5} {len(panel & orig75):>5} "
              f"{surv:>7}/100 {days:>8.1f}")

VOLMAX, SEGMAX = 2000, 3
keep, best = build(VOLMAX, SEGMAX)
panel = sorted(best)
print(f"\n=== RECOMMENDED: volume 300-{VOLMAX} A^3, segments <= {SEGMAX}, "
      f"drugScore >= 0.5, drop multi-pass TM, AF pLDDT >= 70")
print(f"panel: {len(panel)} proteins")
g = [(best[a][0]["Gene Name"].replace("_HUMAN",""), a,
      best[a][0]["Protein"], best[a][0]["Pocket"],
      round(num(best[a][0],"volume")), best[a][2]["nres"], best[a][2]["segs"],
      round(num(best[a][0],"drugScore"),2)) for a in panel]
print(f"median pocket: volume {statistics.median([x[4] for x in g]):.0f} A^3, "
      f"{statistics.median([x[5] for x in g]):.0f} residues, "
      f"{statistics.median([x[6] for x in g]):.0f} segments")
mit = sorted({x[0] for x in g} & MITOTIC)
print(f"mitotic/replication retained ({len(mit)}): {mit}")
with open(f"{work}/panel_net1.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Gene","UniProt","Structure","Pocket","Volume","Residues","Segments","drugScore"])
    w.writerows(sorted(g))
print(f"written: {work}/panel_net1.csv")
print("\nfirst 30 panel members:")
for x in sorted(g)[:30]:
    print("   ", "\t".join(str(y) for y in x))
