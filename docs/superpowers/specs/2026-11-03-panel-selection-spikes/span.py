"""THROWAWAY: pLDDT did not flag the AF mega-cavity hubs. Hypothesis: they are
INTERDOMAIN voids between two confidently-predicted but wrongly-packed domains --
high pLDDT, bad PAE. Cheap proxy: do those pockets span far-apart sequence regions
in several disconnected segments? If so, a sequence-dispersion filter works today."""
import csv, sys, statistics, collections
run = sys.argv[1]
ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")
      and r["Protein"].startswith("AF-")]
deg = collections.Counter()
for r in csv.DictReader(open(f"{run}/06_targets/unique_pairs.csv", newline="")):
    deg[r["Node"]] += 1; deg[r["Node_2"]] += 1

rows = []
for i, r in enumerate(ds):
    if i and i % 1500 == 0: print(f"  {i}/{len(ds)}", flush=True)
    acc = r["UniProt ID"]
    res = set()
    try:
        with open(f'{run}/04_pockets/{acc}/{r["Pocket Filename"]}') as fh:
            for ln in fh:
                if ln.startswith("ATOM"):
                    try: res.add(int(ln[22:26]))
                    except ValueError: pass
    except OSError:
        continue
    if not res: continue
    s = sorted(res)
    segs, prev = 1, s[0]
    for x in s[1:]:
        if x - prev > 5: segs += 1      # gap >5 residues = new segment
        prev = x
    def f(k):
        try: return float(r[k])
        except (ValueError, KeyError): return 0.0
    node = f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
    rows.append({"nres": len(s), "span": s[-1]-s[0]+1, "segs": segs,
                 "disp": (s[-1]-s[0]+1)/len(s), "vol": f("volume"),
                 "deg": deg.get(node, 0),
                 "gene": r["Gene Name"].replace("_HUMAN","")})
print(f"\nAF selected pockets measured: {len(rows)}")

def show(tag, sub):
    if not sub: return
    print(f"  {tag:<34} n={len(sub):>5}  residues {statistics.median([x['nres'] for x in sub]):>5.0f}"
          f"  span {statistics.median([x['span'] for x in sub]):>6.0f}"
          f"  segments {statistics.median([x['segs'] for x in sub]):>4.0f}"
          f"  span/res {statistics.median([x['disp'] for x in sub]):>5.1f}"
          f"  vol {statistics.median([x['vol'] for x in sub]):>6.0f}")
print(f"\n{'group':<36} {'':>3}")
show("all AF selected pockets", rows)
show("in flat Pmax>=0.7 network", [x for x in rows if x["deg"] > 0])
ing = sorted([x for x in rows if x["deg"] > 0], key=lambda x: -x["deg"])
show("top-100 degree AF nodes", ing[:100])
show("top-25 degree AF nodes", ing[:25])
show("volume >= 3000 A^3", [x for x in rows if x["vol"] >= 3000])
show("volume < 1000 A^3", [x for x in rows if x["vol"] < 1000])

print("\npanel effect of a sequence-dispersion cap (best AF pocket per protein kept):")
by = collections.defaultdict(list)
for x in rows: by[x["gene"]].append(x)
for cap in (3, 5, 8, 12, 1e9):
    keep = {g for g, xs in by.items() if any(x["segs"] <= cap for x in xs)}
    hubs = sum(1 for x in ing[:100] if x["segs"] <= cap)
    print(f"  max segments <= {cap if cap<1e9 else 'inf':>4}: {len(keep):>4} proteins retained,"
          f" {hubs:>3}/100 top-degree artifact nodes still pass")
