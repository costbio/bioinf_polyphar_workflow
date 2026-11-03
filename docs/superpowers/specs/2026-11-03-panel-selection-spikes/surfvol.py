"""Does surf/vol reject the oversized interdomain cavities that motivated the geometry
filters in the first place? If so, no volume ceiling is needed.

surf/vol is DoGSiteScorer's own surface-to-volume ratio: low means enclosed, high means
splayed open. It is a SHAPE ratio, so unlike Rg or a raw volume cut it is not simply a
size threshold in disguise.
"""
import csv, json, os, sys, statistics, collections
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")
ds = list(csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline="")))
sel = [r for r in ds if str(r.get("Selected","")).strip().lower() in ("true","1","yes")]
def num(r,k):
    try: return float(r[k])
    except (TypeError,ValueError,KeyError): return None

# the artifacts: top-degree nodes of the ORIGINAL unfiltered Pmax>=0.7 network
deg = collections.Counter()
for r in csv.DictReader(open(f"{run}/06_targets/unique_pairs.csv", newline="")):
    deg[r["Node"]] += 1; deg[r["Node_2"]] += 1
top = {n for n, _ in deg.most_common(100)}
def node(r): return f'{r["UniProt ID"]}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
art = [r for r in sel if node(r) in top]
print(f"known artifact pockets (top-100 degree in the old network): {len(art)}")

def dist(label, rows, key):
    v = sorted(x for x in (num(r,key) for r in rows) if x is not None)
    if not v: return
    def q(p): return v[min(int(p*len(v)), len(v)-1)]
    print(f"  {label:<34} n={len(v):>5}  median {statistics.median(v):>7.3f}  "
          f"p10 {q(.10):>7.3f}  p90 {q(.90):>7.3f}")

for key in ("surf/vol", "volume", "lid/hull"):
    print(f"\n{key}:")
    dist("all selected pockets", sel, key)
    dist("the 100 artifact pockets", art, key)

# ligand-bound reference, from the cached labels
labels = json.load(open(f"{work}/ligand_labels.json"))
occ = [r for r in sel if labels.get(f'{r["UniProt ID"]}/{r["Pocket Filename"]}', {}).get("lig")]
print("\nsurf/vol:")
dist("pockets with a co-crystal ligand", occ, "surf/vol")

print(f"\n{'surf/vol cut':>13} {'artifacts kept':>15} {'ligand-bound kept':>18} {'all kept':>10}")
for cut in (0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5):
    a = sum(1 for r in art if (num(r,"surf/vol") or 9) <= cut)
    l = sum(1 for r in occ if (num(r,"surf/vol") or 9) <= cut)
    s = sum(1 for r in sel if (num(r,"surf/vol") or 9) <= cut)
    print(f"{cut:>13.1f} {f'{a}/{len(art)}':>15} {f'{l}/{len(occ)}':>18} "
          f"{f'{s} ({s/len(sel):.0%})':>10}")
