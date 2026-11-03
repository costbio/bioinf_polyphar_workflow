"""THROWAWAY: fetch TM annotation for net1 proteins, then sweep quality/druggability
filters and report the resulting PANEL (one best pocket per protein).
Uses /uniprotkb/search + Link pagination (stream 500s on multi-accession OR queries)."""
import csv, sys, json, time, urllib.request, urllib.parse, collections, os, re

drugscore_csv, protein_info, cache = sys.argv[1], sys.argv[2], sys.argv[3]

ds = list(csv.DictReader(open(drugscore_csv, newline="")))
info = list(csv.DictReader(open(protein_info, newline="")))
accs = sorted({r["UniProt ID"] for r in ds})
print(f"pockets {len(ds)}, proteins with pockets {len(accs)}")

# ---- UniProt TM / location annotation (cached) --------------------------------------
if os.path.exists(cache):
    ann = json.load(open(cache))
else:
    ann = {}
    F = "accession,ft_transmem,cc_subcellular_location,keyword,length"
    for i in range(0, len(accs), 80):
        batch = accs[i:i+80]
        q = " OR ".join(f"accession:{a}" for a in batch)
        url = ("https://rest.uniprot.org/uniprotkb/search?"
               + urllib.parse.urlencode({"query": q, "fields": F,
                                         "format": "tsv", "size": "500"}))
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    txt = r.read().decode()
                break
            except Exception as e:
                print(f"  retry {attempt}: {e}", flush=True)
                time.sleep(3 * (attempt + 1))
        else:
            sys.exit("UniProt fetch failed")
        rd = csv.DictReader(txt.splitlines(), delimiter="\t")
        for row in rd:
            a = row["Entry"]
            tm = row.get("Transmembrane", "") or ""
            ann[a] = {"tm": tm.count("TRANSMEM"),
                      "loc": (row.get("Subcellular location [CC]", "") or "")[:200],
                      "kw": (row.get("Keywords", "") or "")}
        print(f"  fetched {len(ann)}/{len(accs)}", flush=True)
    json.dump(ann, open(cache, "w"))
miss = [a for a in accs if a not in ann]
print(f"annotated {len(ann)}, missing {len(miss)}{(': '+str(miss[:5])) if miss else ''}")
tmn = sum(1 for a in accs if ann.get(a, {}).get("tm", 0) >= 1)
tmm = sum(1 for a in accs if ann.get(a, {}).get("tm", 0) >= 2)
print(f"proteins with >=1 TM segment: {tmn}; multi-pass (>=2): {tmm}")

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}
orig75 = set()
try:
    for r in csv.DictReader(open(sys.argv[4], newline="")):
        orig75.add(r["UniProt_ID"])
except Exception:
    pass

def f(r, k):
    try: return float(r[k])
    except (ValueError, KeyError, TypeError): return 0.0

print(f"\n{'dsMin':>6} {'volMin':>7} {'exp?':>5} {'TMrule':>8} "
      f"{'pockets':>8} {'PANEL':>6} {'AF%':>5} {'TM%':>5} {'mito':>5} {'of75':>5}")
for dsmin in (0.5, 0.6, 0.7, 0.8):
    for volmin in (0, 300, 600):
        for exp_only in (False, True):
            for tmrule in ("none", "drop>=2", "drop>=1"):
                keep = []
                for r in ds:
                    if f(r, "drugScore") < dsmin: continue
                    if f(r, "volume") < volmin: continue
                    isaf = r["Protein"].startswith("AF-")
                    if exp_only and isaf: continue
                    t = ann.get(r["UniProt ID"], {}).get("tm", 0)
                    if tmrule == "drop>=2" and t >= 2: continue
                    if tmrule == "drop>=1" and t >= 1: continue
                    keep.append(r)
                if not keep:
                    continue
                best = {}
                for r in keep:
                    a = r["UniProt ID"]
                    if a not in best or f(r, "drugScore") > f(best[a], "drugScore"):
                        best[a] = r
                panel = set(best)
                afs = sum(1 for a in panel if best[a]["Protein"].startswith("AF-"))
                tms = sum(1 for a in panel if ann.get(a, {}).get("tm", 0) >= 1)
                genes = {best[a]["Gene Name"].replace("_HUMAN","") for a in panel}
                print(f"{dsmin:>6} {volmin:>7} {str(exp_only):>5} {tmrule:>8} "
                      f"{len(keep):>8} {len(panel):>6} {afs/len(panel)*100:>4.0f}% "
                      f"{tms/len(panel)*100:>4.0f}% {len(genes & MITOTIC):>5} "
                      f"{len(panel & orig75):>5}")
