"""Are randomised controls favoured by short sensory->motor paths? Uses generation-0 stored genomes (W per island)."""
import numpy as np, glob
from collections import deque
B = np.load("data/brain.npz", allow_pickle=True)
K = B["W0"].shape[0]
orn = sorted(set(int(i) for i in B["in_orn"])); dn = set(int(i) for i in B["dn"])
loom = sorted(set(int(i) for i in B["in_loomL"]) | set(int(i) for i in B["in_loomR"]))
taste = sorted(set(int(i) for i in B["in_sugar"]) | set(int(i) for i in B["in_bitter"]))
def stats(W, src):
    # W[i, j] = weight from j to i (rows = targets)
    A = np.abs(W) > 0
    direct = np.abs(W[np.ix_(sorted(dn), src)]).sum() / (np.abs(W[:, src]).sum() + 1e-12)
    # BFS hop distance from any source to any DN
    dist = np.full(K, 99); q = deque()
    for s in src: dist[s] = 0; q.append(s)
    while q:
        u = q.popleft()
        for v in np.nonzero(A[:, u])[0]:
            if dist[v] == 99: dist[v] = dist[u] + 1; q.append(v)
    dd = [dist[d] for d in dn]
    # weight-based 2-hop influence: |W|^2 from sources to DNs, normalised
    M = np.abs(W); M = M / (M.sum(1, keepdims=True) + 1e-12)
    h1 = M[np.ix_(sorted(dn), src)].sum(); h2 = (M @ M)[np.ix_(sorted(dn), src)].sum(); h3 = (M @ M @ M)[np.ix_(sorted(dn), src)].sum()
    return direct, np.median(dd), h1, h2, h3
rows = {}
for f in sorted(glob.glob("runs/eco_P0.0_T0.0_s*/elite_g000000.npz")):
    Z = np.load(f)
    for c in ["A", "N1", "N2", "N3"]:
        W = Z[f"{c}_W"][0].astype(np.float32)
        for name, src in (("smell", orn), ("vision", loom), ("taste", taste)):
            rows.setdefault((c, name), []).append(stats(W, src))
print("per condition, median over 10 seeds. direct = share of sensory output weight landing on descending neurons;")
print("hops = median shortest path sensory->DN; h1/h2/h3 = normalised input influence reaching DNs in 1/2/3 steps")
for name in ("smell", "vision", "taste"):
    print(f"--- {name}")
    for c in ["A", "N1", "N2", "N3"]:
        a = np.array(rows[(c, name)])
        print(f"  {c}: direct {np.median(a[:,0]):.4f}  hops {np.median(a[:,1]):.1f}  h1 {np.median(a[:,2]):.4f}  h2 {np.median(a[:,3]):.4f}  h3 {np.median(a[:,4]):.4f}")
