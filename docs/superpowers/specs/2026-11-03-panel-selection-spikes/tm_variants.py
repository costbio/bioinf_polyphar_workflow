"""THROWAWAY: three TM policies compared on the 241-protein panel definition.
A) drop multi-pass TM only (current)  B) drop any TM protein  C) drop POCKETS that
overlap a TM segment (AF models: PDB resnum == UniProt position; experimental chains
of single-pass proteins are nearly always the soluble construct, so protein-level
multi-pass drop still applies there).
Needs TM segment POSITIONS, so it re-fetches ft_transmem."""
import csv, sys, json, os, re, urllib.request, urllib.parse, time, statistics, collections
run, work = sys.argv[1], os.path.expanduser("~/scratch_clust")

ds = [r for r in csv.DictReader(open(f"{run}/04_pockets/DrugScore.csv", newline=""))
      if r.get("Selected","").lower() in ("true","1","yes")]
geo = json.load(open(f"{work}/geo.json"))
pl = {x["node"]: x["mean"] for x in json.load(open(f"{work}/plddt_rows.json"))}
accs = sorted({r["UniProt ID"] for r in ds})

rng_path = f"{work}/tm_ranges.json"
if os.path.exists(rng_path):
    tmr = json.load(open(rng_path))
else:
    tmr = {}
    for i in range(0, len(accs), 80):
        q = " OR ".join(f"accession:{a}" for a in accs[i:i+80])
        url = ("https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
            {"query": q, "fields": "accession,ft_transmem", "format": "tsv", "size": "500"}))
        for att in range(4):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    txt = r.read().decode(); break
            except Exception as e:
                print("  retry", e, flush=True); time.sleep(3*(att+1))
        else: sys.exit("fetch failed")
        for row in csv.DictReader(txt.splitlines(), delimiter="\t"):
            s = row.get("Transmembrane","") or ""
            tmr[row["Entry"]] = [(int(a), int(b)) for a, b in
                                 re.findall(r"TRANSMEM (\d+)\.\.(\d+)", s)]
        print(f"  {len(tmr)}/{len(accs)}", flush=True)
    json.dump(tmr, open(rng_path, "w"))
nseg = {a: len(tmr.get(a, [])) for a in accs}
print(f"TM ranges for {len(tmr)} proteins; >=1 TM {sum(1 for a in accs if nseg[a]>=1)}, "
      f"multi-pass {sum(1 for a in accs if nseg[a]>=2)}")

resn = {}
for r in ds:
    key = f'{r["UniProt ID"]}/{r["Pocket Filename"]}'
    if key in geo: resn[key] = None   # residues reloaded lazily below only if needed

def pocket_res(acc, fn):
    out = set()
    try:
        with open(f"{run}/04_pockets/{acc}/{fn}") as fh:
            for ln in fh:
                if ln.startswith("ATOM"):
                    try: out.add(int(ln[22:26]))
                    except ValueError: pass
    except OSError: pass
    return out

def num(r, k):
    try: return float(r[k])
    except (ValueError, KeyError): return 0.0

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}
RTK = {"VGFR1","VGFR2","VGFR3","PGFRA","PGFRB","MET","EGFR","ERBB2","FGFR1","FGFR2",
       "FGFR3","FGFR4","KIT","RET","INSR","INSRR","IGF1R","NTRK1","NTRK2","TIE2","AXL"}

def build(policy, volmax=2000, segmax=3):
    keep = []
    for r in ds:
        acc, fn = r["UniProt ID"], r["Pocket Filename"]
        key = f"{acc}/{fn}"; g = geo.get(key)
        if g is None: continue
        if num(r,"drugScore") < 0.5: continue
        v = num(r,"volume")
        if v < 300 or v > volmax: continue
        if g["segs"] > segmax: continue
        isaf = r["Protein"].startswith("AF-")
        node = f'{acc}_{r["Gene Name"]}_{r["Protein"]}_{r["Pocket"]}'
        if isaf and pl.get(node, 0) < 70: continue
        t = nseg.get(acc, 0)
        if policy == "A" and t >= 2: continue
        if policy == "B" and t >= 1: continue
        if policy == "C":
            if t >= 2: continue
            if isaf and t >= 1:
                pr = pocket_res(acc, fn)
                tmres = {x for a, b in tmr.get(acc, []) for x in range(a, b+1)}
                if pr and len(pr & tmres)/len(pr) > 0.10: continue
        keep.append((r, acc))
    best = {}
    for r, acc in keep:
        if acc not in best or num(r,"drugScore") > num(best[acc],"drugScore"):
            best[acc] = r
    return best

print(f"\n{'policy':<42} {'panel':>6} {'TM%':>5} {'mito':>5} {'RTK kept':>9} {'dock(d)':>8}")
res = {}
for pol, lab in (("A","A: drop multi-pass TM only (current)"),
                 ("B","B: drop ANY TM protein"),
                 ("C","C: drop pockets overlapping TM (>10%)")):
    best = build(pol); res[pol] = best
    panel = set(best)
    genes = {best[a]["Gene Name"].replace("_HUMAN","") for a in panel}
    tmc = sum(1 for a in panel if nseg.get(a,0) >= 1)
    print(f"{lab:<42} {len(panel):>6} {tmc/len(panel)*100:>4.0f}% "
          f"{len(genes & MITOTIC):>5} {len(genes & RTK):>9} "
          f"{len(panel)*2083*57/64/86400:>8.1f}")

lost = set(res["A"]) - set(res["B"])
print(f"\nproteins lost by going A -> B ({len(lost)}):")
rows = sorted((res["A"][a]["Gene Name"].replace("_HUMAN",""), a, nseg.get(a,0))
              for a in lost)
for g, a, t in rows:
    print(f"   {g:<10} {a:<10} TM segments {t}")
