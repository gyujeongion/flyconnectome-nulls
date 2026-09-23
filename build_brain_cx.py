"""FlyWire v783 -> compact cell-type-level brain (K groups) for evolution.
Output: data/brain_cx.npz  (brain.npz + neuron-level central complex layer)
"""
import numpy as np, pandas as pd, scipy.sparse as sp, json, sys, time
t0 = time.time()
K_TOTAL = int(sys.argv[1]) if len(sys.argv) > 1 else 512
rng = np.random.default_rng(0)
a = pd.read_csv("data/annotations.tsv", sep="\t", low_memory=False)
c = pd.read_parquet("data/Connectivity_783.parquet",
                    columns=["Presynaptic_ID", "Postsynaptic_ID", "Connectivity", "Excitatory x Connectivity"])
a = a[~a.super_class.isin(["optic", "visual_centrifugal"])].copy()
a = a[~((a.super_class == "sensory") & (a.cell_class.isin(["visual", "optic_lobes"])))]
side = a.side.fillna("center").map({"left": "L", "right": "R"}).fillna("C")

ORN_KEEP = ["ORN_V"] + [t for t in a[a.cell_class == "olfactory"].cell_type.value_counts().index if t != "ORN_V"][:15]
key = a.cell_type.fillna(a.hemibrain_type).fillna(a.cell_class).fillna("unk_" + a.super_class.astype(str))
key = key.astype(str)
is_orn = a.cell_class == "olfactory"
key[is_orn & ~a.cell_type.isin(ORN_KEEP)] = "ORN_other"
is_grn = a.cell_class == "gustatory"
key[is_grn & (a.cell_sub_class == "sugar/water")] = "GRN_sugar"
key[is_grn & (a.cell_sub_class == "bitter")] = "GRN_bitter"
key[is_grn & ~a.cell_sub_class.isin(["sugar/water", "bitter"])] = "GRN_other"
bilateral_merge = a.cell_class.isin(["MBON", "DAN", "MBIN"])
gkey = np.where(bilateral_merge, key, key + "|" + side)
a["g"] = gkey

# --- individual KCs: sample N_KC_IND real neurons, whole KC population = one bulk group for path scoring ---
N_KC_IND = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
kc = a[a.cell_class == "Kenyon_Cell"]
kc_sample = np.concatenate([rng.choice(kc[side.loc[kc.index] == s_].root_id.values, N_KC_IND // 2, replace=False) for s_ in "LR"])
a.loc[kc.index, "g"] = "KCBULK"
# --- central complex: navigation core + its sensory interface neurons at single-neuron resolution ---
CX_CORE = r"^(EPG|PEG|PEN|Delta7|PFL|PFN|hDeltaB|FC2)"
CX_IN = r"^(GLNO|LNO1|LNO2|LNOa|ER1|ER2|ER3|ER4|ER5|ER6)"
ctype = a.cell_type.fillna("")
is_core = (a.cell_class == "CX") & ctype.str.match(CX_CORE)
is_cxin = (a.cell_class == "CX") & ctype.str.match(CX_IN)
cxn = a[is_core | is_cxin].copy()
cxn["role"] = np.where(ctype[cxn.index].str.match(CX_IN), "in", "core")
cxn = cxn.sort_values(["role", "cell_type", "side", "root_id"])
cx_ids = cxn.root_id.values
a.loc[cxn.index, "g"] = "CXBULK"
print(f"CX layer neurons: {len(cx_ids)} (core {int((cxn.role=='core').sum())}, inputs {int((cxn.role=='in').sum())})")
# --- group-level matrix ---
gid_of = pd.Series(a.g.values, index=a.root_id.values)
names = np.array(sorted(a.g.unique())); name2i = {n: i for i, n in enumerate(names)}
nid = gid_of.map(name2i)
c = c[c.Presynaptic_ID.isin(nid.index) & c.Postsynaptic_ID.isin(nid.index)]
pre = nid.reindex(c.Presynaptic_ID).values; post = nid.reindex(c.Postsynaptic_ID).values
G = len(names)
Sabs = sp.coo_matrix((c.Connectivity.values.astype(np.float32), (post, pre)), shape=(G, G)).tocsr()
Ssgn = sp.coo_matrix((c["Excitatory x Connectivity"].values.astype(np.float32), (post, pre)), shape=(G, G)).tocsr()
nneur = np.bincount(nid.values, minlength=G)
print(f"groups={G} edges={len(c):,} t={time.time()-t0:.0f}s")

def members(pred):
    return [i for i, n in enumerate(names) if pred(n)]
base = lambda n: n.split("|")[0]
inp = {
    "orn": [name2i[f"{o}|{s}"] for o in ORN_KEEP for s in "LR" if f"{o}|{s}" in name2i],
    "sugar": members(lambda n: base(n) == "GRN_sugar"),
    "bitter": members(lambda n: base(n) == "GRN_bitter"),
    "loomL": members(lambda n: n in ("LC4|L", "LPLC2|L")),
    "loomR": members(lambda n: n in ("LC4|R", "LPLC2|R")),
}
kc_i = []
mbon_i = [name2i[n] for n in names if a[a.g == n].cell_class.eq("MBON").any()] if False else \
    sorted(set(nid[a[a.cell_class == "MBON"].root_id].values))
dan_i = sorted(set(nid[a[a.cell_class.isin(["DAN", "MBIN"])].root_id].values))
dn_all = sorted(set(nid[a[a.super_class == "descending"].root_id].values))
NAMED_DN = ["DNa01", "DNa02", "DNp09", "oDN1", "MDN", "DNp01", "DNg11", "DNg12_a"]
dn_named = [i for i in dn_all if base(names[i]) in NAMED_DN]
print("named DN found:", [names[i] for i in dn_named])

# path score: reachability from inputs x to outputs (3 hops, row-normalised transfer)
T = sp.diags(1 / (np.asarray(Sabs.sum(1)).ravel() + 1)) @ Sabs
src = np.zeros(G, np.float32); src[sum(inp.values(), [])] = 1
dst = np.zeros(G, np.float32); dst[dn_all] = 1
f = np.zeros(G, np.float32); v = src.copy()
for _ in range(4): v = T @ v; f += v
b = np.zeros(G, np.float32); v = dst.copy()
for _ in range(4): v = T.T @ v; b += v
score = f * b
# reinforcement pathways: GRN -> DAN and PN -> KC must survive compression
def reach(seed_idx, back=False, hops=4):
    v = np.zeros(G, np.float32); v[seed_idx] = 1; acc = np.zeros(G, np.float32)
    M = T.T if back else T
    for _ in range(hops): v = M @ v; acc += v
    return acc
grn = inp["sugar"] + inp["bitter"]
path_rein = reach(grn) * reach(dan_i, back=True)
path_rein[grn + dan_i] = 0
rein_keep = [int(i) for i in np.argsort(-path_rein)[:48] if path_rein[i] > 0 and names[i] not in ("KCBULK", "CXBULK")]
print("GRN->DAN path groups kept:", len(rein_keep), [names[i] for i in rein_keep[:12]])
# antennal lobe: every PN group getting >=5% of its input from kept ORNs, plus top LNs driven by kept ORNs
orn_idx = inp["orn"]
Sa_dense_cols = Sabs[:, orn_idx]
frac_from_orn = np.asarray(Sa_dense_cols.sum(1)).ravel() / (np.asarray(Sabs.sum(1)).ravel() + 1)
pn_groups = sorted(set(nid[a[a.cell_class == "ALPN"].root_id].values))
al_pn = [i for i in pn_groups if frac_from_orn[i] >= 0.05]
ln_groups = sorted(set(nid[a[a.cell_class.isin(["ALLN"])].root_id].values))
al_ln = [i for i in sorted(ln_groups, key=lambda i: -frac_from_orn[i])[:24] if frac_from_orn[i] > 0.01]
print("AL PN kept", len(al_pn), "LN kept", len(al_ln))
pfl3_ids = cxn[cxn.cell_type.str.startswith("PFL3")].root_id.values
cx_all_ids = set(cx_ids)
eo = c[c.Presynaptic_ID.isin(pfl3_ids) & ~c.Postsynaptic_ID.isin(cx_all_ids)]
tg = eo.assign(g=nid.reindex(eo.Postsynaptic_ID).values).groupby("g").Connectivity.sum().sort_values(ascending=False)
pfl3_targets = [int(i) for i in tg.index[:16]]
steer_dn = [i for i in dn_all if base(names[i]) in ("DNa01", "DNa02")]
path_steer = reach(pfl3_targets) * reach(steer_dn, back=True); path_steer[pfl3_targets + steer_dn] = 0
steer_keep = pfl3_targets + [int(i) for i in np.argsort(-path_steer)[:24] if path_steer[i] > 0 and names[i] not in ("KCBULK", "CXBULK")]
print("PFL3 targets:", [names[i] for i in pfl3_targets[:8]], "| steering path groups:", len(steer_keep))
mandatory = sorted(set(sum(inp.values(), []) + mbon_i + dan_i + dn_named + rein_keep + al_pn + al_ln + steer_keep
                       + members(lambda n: base(n) in ("APL", "DPM"))))
dn_extra = [i for i in np.argsort(-f[dn_all]) if dn_all[i] not in mandatory][:32]
mandatory = sorted(set(mandatory + [dn_all[i] for i in dn_extra]))
rest = [i for i in np.argsort(-score) if i not in set(mandatory) and not names[i].startswith(("ORN_other", "GRN_other", "KCBULK", "CXBULK"))]
keep = mandatory + rest[:max(0, K_TOTAL - len(mandatory))]
keep = np.array(sorted(keep)); K = len(keep)
remap = -np.ones(G, int); remap[keep] = np.arange(K)
print(f"mandatory={len(mandatory)} K={K}")

Ssub = Ssgn[keep][:, keep].toarray()
Abs_in_total = np.asarray(Sabs[keep].sum(1)).ravel()
W0 = Ssub / (Abs_in_total[:, None] + 1.0)          # fraction of each group's total input (incl. dropped)
sensor_rows = remap[sum(inp.values(), [])]; sensor_rows = sensor_rows[sensor_rows >= 0]
W0[sensor_rows, :] = 0.0                          # sensory neurons = transducers (no recurrent drive)
outputs = sorted(remap[[i for i in dn_all if remap[i] >= 0]])
Aabs = Sabs[keep][:, keep].toarray()
km = lambda L: np.array(sorted(remap[L]))
craw = pd.read_parquet("data/Connectivity_783.parquet", columns=["Presynaptic_ID", "Postsynaptic_ID", "Connectivity", "Excitatory x Connectivity"])
kc_pos = pd.Series(np.arange(len(kc_sample)), index=kc_sample)
grp_of_root = nid.map(lambda i: remap[i])                     # root -> kept group index or -1
tot_in_kc = craw[craw.Postsynaptic_ID.isin(kc_sample)].groupby("Postsynaptic_ID").Connectivity.sum().reindex(kc_sample).fillna(0).values
e_in = craw[craw.Postsynaptic_ID.isin(kc_sample) & craw.Presynaptic_ID.isin(grp_of_root.index)]
gi = grp_of_root.reindex(e_in.Presynaptic_ID).values; ok = gi >= 0
Wpk = np.zeros((len(kc_sample), K), np.float32)
np.add.at(Wpk, (kc_pos.reindex(e_in.Postsynaptic_ID).values[ok], gi[ok]), e_in["Excitatory x Connectivity"].values[ok])
Wpk /= (tot_in_kc[:, None] + 1.0)
e_out = craw[craw.Presynaptic_ID.isin(kc_sample) & craw.Postsynaptic_ID.isin(grp_of_root.index)]
go = grp_of_root.reindex(e_out.Postsynaptic_ID).values; ok = go >= 0
Wko = np.zeros((K, len(kc_sample)), np.float32)
np.add.at(Wko, (go[ok], kc_pos.reindex(e_out.Presynaptic_ID).values[ok]), e_out["Excitatory x Connectivity"].values[ok])
Wko *= (len(kc) / len(kc_sample)) / (Abs_in_total[:, None] + 1.0)
print(f"KC individual: n={len(kc_sample)} PN-group inputs per KC mean {(Wpk != 0).sum(1).mean():.1f} | KC input frac retained {np.abs(Wpk).sum(1).mean():.2f} | KC->kept targets per KC {(Wko != 0).sum(0).mean():.1f}")

NCX = len(cx_ids); cx_pos = pd.Series(np.arange(NCX), index=cx_ids)
tot_in_cx = craw[craw.Postsynaptic_ID.isin(cx_ids)].groupby("Postsynaptic_ID").Connectivity.sum().reindex(cx_ids).fillna(0).values
e = craw[craw.Postsynaptic_ID.isin(cx_ids) & craw.Presynaptic_ID.isin(cx_ids)]
Wcc = np.zeros((NCX, NCX), np.float32)
np.add.at(Wcc, (cx_pos.reindex(e.Postsynaptic_ID).values, cx_pos.reindex(e.Presynaptic_ID).values), e["Excitatory x Connectivity"].values)
Wcc /= (tot_in_cx[:, None] + 1.0)
e = craw[craw.Postsynaptic_ID.isin(cx_ids) & craw.Presynaptic_ID.isin(grp_of_root.index)]
gi = grp_of_root.reindex(e.Presynaptic_ID).values; ok = gi >= 0
Wgc = np.zeros((NCX, K), np.float32)
np.add.at(Wgc, (cx_pos.reindex(e.Postsynaptic_ID).values[ok], gi[ok]), e["Excitatory x Connectivity"].values[ok])
Wgc /= (tot_in_cx[:, None] + 1.0)
e = craw[craw.Presynaptic_ID.isin(cx_ids) & craw.Postsynaptic_ID.isin(grp_of_root.index)]
go = grp_of_root.reindex(e.Postsynaptic_ID).values; ok = go >= 0
Wcg = np.zeros((K, NCX), np.float32)
np.add.at(Wcg, (go[ok], cx_pos.reindex(e.Presynaptic_ID).values[ok]), e["Excitatory x Connectivity"].values[ok])
Wcg /= (Abs_in_total[:, None] + 1.0)
cx_role = cxn.role.values; cx_type = cxn.cell_type.values.astype(str); cx_side = cxn.side.fillna("c").values.astype(str)
cx_sensor = cx_role == "in"
Wcc[cx_sensor, :] = 0.0; Wgc[cx_sensor, :] = 0.0         # interface neurons are transducers of self-motion / visual cues
print(f"CX: Wcc nnz {(Wcc != 0).sum()} | input frac retained within CX (core) {np.abs(Wcc[~cx_sensor]).sum(1).mean():.2f} | from groups {np.abs(Wgc[~cx_sensor]).sum(1).mean():.2f} | CX->group edges {(Wcg != 0).sum()}")
mb = km(mbon_i); dn_ = km(dan_i)
Mdan = Aabs[np.ix_(mb, dn_)]
Mdan = Mdan / (Mdan.sum(1, keepdims=True) + 1e-6)
pn_k = np.array(sorted(remap[al_pn])); ln_k = np.array(sorted(remap[al_ln]))
out = dict(Wcc=Wcc, Wgc=Wgc, Wcg=Wcg, cx_type=cx_type, cx_side=cx_side, cx_sensor=cx_sensor, steer_groups=km(steer_keep),
           W0=W0.astype(np.float32), Wpk=Wpk, Wko=Wko, names=names[keep], pn=pn_k, ln=ln_k, nneur=nneur[keep],
           mbon=mb, dan=dn_, dan_to_mbon=Mdan.astype(np.float32), dn=np.array(outputs),
           dn_named=km(dn_named), orn_types=np.array(ORN_KEEP))
for k_, L in inp.items(): out["in_" + k_] = km(L)
np.savez_compressed("data/brain_cx.npz", **out)
print(json.dumps({k_: int(len(v_)) for k_, v_ in out.items() if hasattr(v_, "__len__") and k_ not in ("W0", "Wcc")}))
print("W0 nnz", int((W0 != 0).sum()), "neurons covered", int(nneur[keep].sum()), f"t={time.time()-t0:.0f}s")
print("DAN->MBON rows with input", int((Aabs[np.ix_(mb, dn_)].sum(1) > 0).sum()), "/", len(mb))
