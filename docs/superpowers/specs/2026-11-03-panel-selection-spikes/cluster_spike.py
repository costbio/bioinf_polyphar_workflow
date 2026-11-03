"""THROWAWAY spike: does community detection on the net1 pocket-similarity network
recover the mitotic/replication group on its own? Not pipeline code."""
import csv, sys, statistics, collections
import networkx as nx
from networkx.algorithms.community import louvain_communities, modularity

pairs_csv, targets_csv = sys.argv[1], sys.argv[2]

G = nx.Graph()
gene_of, acc_of = {}, {}
with open(pairs_csv, newline="") as fh:
    for r in csv.DictReader(fh):
        a, b = r["Node"], r["Node_2"]
        w = float(r["Pmax"])
        G.add_edge(a, b, weight=w)
        for n, g, acc in ((a, r["Gene Name"], r["UniProt ID"]),
                          (b, r["Gene Name_2"], r["UniProt ID_2"])):
            gene_of[n] = g.replace("_HUMAN", "")
            acc_of[n] = acc

print(f"graph: {G.number_of_nodes()} pocket nodes, {G.number_of_edges()} edges, "
      f"{nx.number_connected_components(G)} components")
giant = max(nx.connected_components(G), key=len)
print(f"largest component: {len(giant)} nodes ({len(giant)/G.number_of_nodes():.0%})")

MITOTIC = {"PLK1","CDK1","TTK","NEK2","BUB1","MELK","CDC20","UBE2C","KI18A","NUF2",
           "CENPM","CENPI","CDC45","CDC6","ORC1","CCNE1","TOP2A","ECT2","CND3",
           "EXO1","RECQ4","TERT","PCH2"}
hub_genes = set()
with open(targets_csv, newline="") as fh:
    for r in csv.DictReader(fh):
        if r["Is_Hub"] == "True":
            hub_genes.add(r["Gene_Name"])

present = {g for g in gene_of.values() if g in MITOTIC}
print(f"mitotic/replication genes present in network: {len(present)}/{len(MITOTIC)}")
print(f"missing: {sorted(MITOTIC - present)}\n")

for res in (0.5, 1.0, 1.5, 2.0):
    comms = louvain_communities(G, weight="weight", resolution=res, seed=42)
    comms = sorted(comms, key=len, reverse=True)
    sizes = [len(c) for c in comms]
    q = modularity(G, comms, weight="weight", resolution=res)
    print(f"=== resolution {res}: {len(comms)} communities, Q={q:.3f}, "
          f"sizes top10 {sizes[:10]}, singletons {sizes.count(1)}")
    # where do the mitotic genes land?
    loc = collections.defaultdict(set)
    for ci, c in enumerate(comms):
        for n in c:
            if gene_of[n] in MITOTIC:
                loc[gene_of[n]].add(ci)
    spread = collections.Counter()
    for g, cs in loc.items():
        for c in cs:
            spread[c] += 1
    top = spread.most_common(5)
    print(f"    mitotic genes span {len(spread)} communities; "
          f"most populated: {top}")
    if top:
        ci = top[0][0]
        c = comms[ci]
        prot = collections.Counter(gene_of[n] for n in c)
        sub = G.subgraph(c)
        dens = nx.density(sub)
        iw = [d["weight"] for *_, d in sub.edges(data=True)]
        mit_in = sorted(g for g in prot if g in MITOTIC)
        print(f"    community {ci}: {len(c)} pockets / {len({acc_of[n] for n in c})} proteins,"
              f" density {dens:.3f}, mean internal Pmax {statistics.mean(iw):.3f}")
        print(f"    mitotic members ({len(mit_in)}): {mit_in}")
        other = sorted(g for g in prot if g not in MITOTIC)
        print(f"    other members ({len(other)}): {other[:40]}")
        print(f"    hub {sorted(hub_genes)} in this community: "
              f"{bool(hub_genes & set(prot))}")
    print()
