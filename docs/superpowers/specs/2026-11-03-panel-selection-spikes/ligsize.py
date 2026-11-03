"""THROWAWAY: size/flexibility profile of e-Drug3D, to pick a dockability cutoff.
Rotatable bonds = single, non-ring, non-terminal bonds between heavy atoms. Non-ring is
found by bridge detection on the heavy-atom bond graph (pure stdlib, no RDKit)."""
import sys, statistics, collections
MASS = {"C":12.011,"N":14.007,"O":15.999,"S":32.06,"P":30.974,"F":18.998,"CL":35.45,
        "BR":79.904,"I":126.90,"H":1.008,"NA":22.99,"K":39.098,"CA":40.078,"MG":24.305,
        "FE":55.845,"ZN":65.38,"SE":78.971,"B":10.81,"SI":28.085,"AS":74.922,"LI":6.94,
        "AL":26.982,"GD":157.25,"PT":195.08,"TC":98.0,"AU":196.97,"AG":107.87,"MN":54.938,
        "CU":63.546,"CR":51.996,"CO":58.933,"SR":87.62,"BA":137.33,"SM":150.36,"RA":226.0,
        "Y":88.906,"IN":114.82,"TL":204.38,"BI":208.98,"SB":121.76,"HG":200.59,"ST":0.0}
def bridges(n, edges):
    adj = collections.defaultdict(list)
    for i,(a,b) in enumerate(edges): adj[a].append((b,i)); adj[b].append((a,i))
    disc, low, out, t = {}, {}, set(), [0]
    for s in range(1, n+1):
        if s in disc: continue
        stack = [(s, -1, iter(adj[s]))]
        disc[s] = low[s] = t[0]; t[0] += 1
        while stack:
            u, pe, it = stack[-1]
            advanced = False
            for v, ei in it:
                if ei == pe: continue
                if v not in disc:
                    disc[v] = low[v] = t[0]; t[0] += 1
                    stack.append((v, ei, iter(adj[v]))); advanced = True; break
                low[u] = min(low[u], disc[v])
            if not advanced:
                stack.pop()
                if stack:
                    p, ppe, _ = stack[-1]
                    low[p] = min(low[p], low[u])
                    if low[u] > disc[p]: out.add(pe)
    return out

mols = []
cur = []
for line in open(sys.argv[1], errors="replace"):
    if line.startswith("$$$$"): mols.append(cur); cur = []
    else: cur.append(line)
if cur: mols.append(cur)
rows = []
for b in mols:
    if len(b) < 4: continue
    try:
        na, nb = int(b[3][0:3]), int(b[3][3:6])
    except ValueError: continue
    el, mw = [], 0.0
    for ln in b[4:4+na]:
        e = (ln[31:34].strip() or "C").upper()
        el.append(e); mw += MASS.get(e, 12.0)
    heavy_idx = {i+1 for i,e in enumerate(el) if e != "H"}
    heavy = len(heavy_idx)
    edges, orders = [], []
    for ln in b[4+na:4+na+nb]:
        try: a1, a2, o = int(ln[0:3]), int(ln[3:6]), int(ln[6:9])
        except ValueError: continue
        if a1 in heavy_idx and a2 in heavy_idx:
            edges.append((a1,a2)); orders.append(o)
    br = bridges(na, edges)
    deg = collections.Counter()
    for a,b2 in edges: deg[a]+=1; deg[b2]+=1
    rot = sum(1 for i,(a,b2) in enumerate(edges)
              if orders[i]==1 and i in br and deg[a]>1 and deg[b2]>1)
    name = ""
    for i,ln in enumerate(b):
        if ln.strip().startswith("> <name>"): name = b[i+1].strip(); break
    rows.append({"name":name,"heavy":heavy,"mw":mw,"rot":rot})
print(f"parsed {len(rows)} of {len(mols)} records\n")
hv=sorted(r["heavy"] for r in rows); mws=sorted(r["mw"] for r in rows); rt=sorted(r["rot"] for r in rows)
def q(a,p): return a[int(p*(len(a)-1))]
print(f"heavy atoms : median {q(hv,.5)}  p90 {q(hv,.9)}  p99 {q(hv,.99)}  max {hv[-1]}")
print(f"MW (Da)     : median {q(mws,.5):.0f}  p90 {q(mws,.9):.0f}  p99 {q(mws,.99):.0f}  max {mws[-1]:.0f}")
print(f"rot bonds   : median {q(rt,.5)}  p90 {q(rt,.9)}  p99 {q(rt,.99)}  max {rt[-1]}")
print("\n--- 25 largest by rotatable bonds (the ones that stall Vina-type search):")
for r in sorted(rows, key=lambda x:-x["rot"])[:25]:
    print(f"  {r['name'][:42]:<42} heavy {r['heavy']:>4}  MW {r['mw']:>7.0f}  rot {r['rot']:>3}")
print("\n--- how many excluded by candidate cutoffs:")
for lab, f in (("MW > 900", lambda r: r["mw"]>900),
               ("MW > 750", lambda r: r["mw"]>750),
               ("MW > 600", lambda r: r["mw"]>600),
               ("heavy > 60", lambda r: r["heavy"]>60),
               ("heavy > 50", lambda r: r["heavy"]>50),
               ("rot > 15", lambda r: r["rot"]>15),
               ("rot > 10", lambda r: r["rot"]>10),
               ("rot > 15 OR MW > 900", lambda r: r["rot"]>15 or r["mw"]>900),
               ("rot > 12 OR MW > 750", lambda r: r["rot"]>12 or r["mw"]>750)):
    n = sum(1 for r in rows if f(r))
    print(f"  {lab:<22} excludes {n:>4} of {len(rows)} ({n/len(rows):.1%}), keeps {len(rows)-n}")
