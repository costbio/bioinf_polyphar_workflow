"""THROWAWAY: heavy-atom count per eDrug3D ligand, and where the known HCC drugs
sit in that distribution. Tests whether the Fisher test on docking score is
confounded by ligand size."""
import sys, statistics
path = sys.argv[1]
mols = []
with open(path, errors="replace") as fh:
    block = []
    for line in fh:
        if line.startswith("$$$$"):
            mols.append(block); block = []
        else:
            block.append(line.rstrip("\n"))
    if block: mols.append(block)
out = {}
for b in mols:
    if len(b) < 4: continue
    try: na = int(b[3][0:3])
    except ValueError: continue
    heavy = 0
    for ln in b[4:4+na]:
        el = ln[31:34].strip() if len(ln) > 34 else ln.split()[3] if len(ln.split())>3 else ""
        if el and el != "H": heavy += 1
    name = ""
    for i, ln in enumerate(b):
        if ln.strip().startswith("> <name>"): name = b[i+1].strip(); break
    if name: out[name.upper()] = heavy
print(f"ligands parsed: {len(out)}")
hv = sorted(out.values())
print(f"heavy atoms: median {statistics.median(hv)}, mean {statistics.mean(hv):.1f}, "
      f"p10 {hv[len(hv)//10]}, p90 {hv[9*len(hv)//10]}")
HCC = ["SORAFENIB","LENVATINIB","REGORAFENIB","CABOZANTINIB","DONAFENIB","APATINIB",
       "CAPECITABINE","DOXORUBICIN","GEMCITABINE","OXALIPLATIN","CISPLATIN",
       "FLUOROURACIL","MITOXANTRONE","TIVANTINIB","BRIVANIB","SUNITINIB","ERLOTINIB"]
print(f"\n{'drug':>14} {'heavy':>6} {'percentile':>11}")
found = []
for d in HCC:
    hits = [k for k in out if d in k]
    if not hits: print(f"{d:>14} {'--':>6} {'not in lib':>11}"); continue
    h = out[hits[0]]
    pct = sum(1 for x in hv if x <= h) / len(hv) * 100
    found.append(h)
    print(f"{d:>14} {h:>6} {pct:>10.0f}%")
if found:
    print(f"\nfound {len(found)} HCC drugs: median heavy atoms {statistics.median(found)} "
          f"vs library median {statistics.median(hv)}")
