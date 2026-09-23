"""evo_fix.py: evo_paper.py with defect 28 fixed (zero-mean turning noise); nothing else changed.
evo3 — connectome-seeded neuroevolution, clean engine (MLX / Apple Silicon).

Brain  : K cell-type groups (rate units, homeostatically calibrated) + N individual Kenyon cells
         with real PN->KC wiring, APL-like sparsening, plastic KC->MBON synapses gated by DANs.
World  : parallel arenas; 8 odor-coded food patches; toxic food tastes sweet and causes delayed malaise
         -> only odor->outcome memory avoids it. Adaptive predator (CO2 odor + loom).
Islands: A* = FlyWire connectome, B* = column-shuffled control (Dale sign + out-degree kept).
Modes  : --circuit (unit tests, seconds) | --diag | default evolution run.
"""
import argparse, json, math, os, time
import numpy as np
import mlx.core as mx

# ======================================================================== config
ap = argparse.ArgumentParser()
A_ = ap.add_argument
A_("--run", default="v3"); A_("--brain", default="data/brain.npz"); A_("--seed", type=int, default=0)
A_("--pop", type=int, default=2048); A_("--islands", type=int, default=4); A_("--T", type=int, default=400)
A_("--gens", type=int, default=10**6); A_("--hours", type=float, default=1e9); A_("--probe_every", type=int, default=25)
A_("--migrate_every", type=int, default=0)
A_("--npatch", type=int, default=24); A_("--nbad", type=int, default=8); A_("--bite", type=int, default=12)
A_("--E0", type=float, default=0.7); A_("--bite_bad", type=int, default=1); A_("--gain_good", type=float, default=0.04); A_("--gain_bad", type=float, default=0.4)
A_("--cost_scale", type=float, default=0.1); A_("--kc_sparsity", type=float, default=0.05)
A_("--homeo_target", type=float, default=0.15); A_("--eta0", type=float, default=0.0)
A_("--c_base", type=float, default=0.0025); A_("--sig", type=float, default=0.06); A_("--ant", type=float, default=0.05); A_("--sweep", action="store_true"); A_("--landscape", action="store_true"); A_("--contra", type=float, default=1.0); A_("--eat_scale", type=float, default=0.1); A_("--trace", type=float, default=0.0); A_("--flip_rate", type=float, default=0.0); A_("--nobrain", action="store_true"); A_("--regrow", type=int, default=0); A_("--regrow_time", type=int, default=100); A_("--ecosearch", type=int, default=0); A_("--mal_gain", type=float, default=0.25); A_("--wlearn", action="store_true"); A_("--rm0", type=float, default=0.5); A_("--behav", action="store_true"); A_("--ptox0", type=float, default=0.0); A_("--tox_curriculum", type=int, default=0); A_("--homeo", default="global"); A_("--homeo_init", type=float, default=1.0); A_("--conds", default="A,N1,N2,N3"); A_("--pred_fixed", type=float, default=-1.0); A_("--probe_ptox", type=float, default=0.5); A_("--probe_reps", type=int, default=8); A_("--gen0assay", action="store_true"); A_("--assay", default=""); A_("--assay_gen", type=int, default=250); A_("--circuit", action="store_true"); A_("--diag", action="store_true"); A_("--ladder", action="store_true")
args = ap.parse_args()
rng = np.random.default_rng(args.seed); mx.random.seed(args.seed)
OUT = f"runs/{args.run}"; os.makedirs(OUT, exist_ok=True)
DT = mx.float16

# ======================================================================== brain constants
B = np.load(args.brain, allow_pickle=True)
names = [str(n) for n in B["names"]]; K = len(names); n2i = {n: i for i, n in enumerate(names)}
W0, Wpk0, Wko0 = B["W0"].astype(np.float32), B["Wpk"].astype(np.float32), B["Wko"].astype(np.float32)
NKC = Wpk0.shape[0]
mb, dan, dn = B["mbon"], B["dan"], B["dn"]; nMB, nDN = len(mb), len(dn)
Mdan = B["dan_to_mbon"].astype(np.float32)
orn_types = [str(o) for o in B["orn_types"]]; NG = len(orn_types)      # [0] = ORN_V (CO2)
in_idx, in_kind = [], []
for s in "LR":
    for o in orn_types: in_idx.append(n2i[f"{o}|{s}"]); in_kind.append(0)
for key, kind in (("in_sugar", 1), ("in_bitter", 1), ("in_loomL", 2), ("in_loomR", 2)):
    for gi in B[key]: in_idx.append(int(gi)); in_kind.append(kind)
NS = len(in_idx); I_SUG, I_BIT, I_LOOM = 2 * NG, 2 * NG + 2, 2 * NG + 4
Bin = np.zeros((K, NS), np.float32); Bin[in_idx, np.arange(NS)] = 1
Gmap = np.zeros((NS, 3), np.float32); Gmap[np.arange(NS), in_kind] = 1
Pmb = np.zeros((K, nMB), np.float32); Pmb[mb, np.arange(nMB)] = 1
nneur = B["nneur"].astype(np.float32); N0 = float(nneur.sum()) + 5177.0
dsign = np.sign(W0.sum(0)); dsign[dsign == 0] = 1
sensor = np.zeros(K, bool); sensor[in_idx] = True
side_of = np.array([n.split("|")[1] if "|" in n else "C" for n in names])
if args.contra != 1.0:                                    # ipsilateral release bias (Gaudry et al. 2013)
    orn_cols = np.array(in_idx[:2 * NG])
    for j in orn_cols:
        sj = side_of[j]
        contra_rows = (side_of != sj) & (side_of != "C")
        W0[contra_rows, j] *= args.contra
    print(f"ORN contralateral output scaled x{args.contra}", flush=True)
Wko_fixed0 = Wko0.copy(); Wko_fixed0[mb] = 0.0            # KC->MBON handled by plastic matrix
Wkm0 = Wko0[mb]                                          # (nMB, NKC) initial plastic weights
dn_names = [names[i].split("|") for i in dn]
def dnw(base, side=None):
    return np.array([1.0 if (b == base and (side is None or s == side)) else 0 for b, s in dn_names], np.float32)
_types = sorted({b for b, s_ in dn_names if s_ == "L"} & {b for b, s_ in dn_names if s_ == "R"})
pL = np.array([[i for i, (b, s_) in enumerate(dn_names) if b == t and s_ == "L"][0] for t in _types])
pR = np.array([[i for i, (b, s_) in enumerate(dn_names) if b == t and s_ == "R"][0] for t in _types])
NP = len(_types)
def pw(t): return np.array([1.0 if x == t else 0.0 for x in _types], np.float32)
R_prior = np.zeros((3, NP), np.float32)                  # readout over bilateral DN pairs
R_prior[0] = 2 * pw("DNp09") - 2 * pw("MDN") + pw("DNp01")   # forward: (L+R)
R_prior[1] = 2 * pw("DNa01") + 2 * pw("DNa02")               # turn: (L-R), antisymmetric by construction
print(f"DN bilateral pairs: {NP}", flush=True)
def dn_readout(Rm, rdn, lrbar):
    """Rm (n,3,NP), rdn (n,nDN) -> (n,3) logits. forward/eat: L+R. turn: adapted (L-R) minus its running mean
    (removes fixed hemispheric bias; leaves bilateral + temporal contrast). returns logits, new lrbar."""
    DL, DR = mx.take(rdn, C(pL), axis=1), mx.take(rdn, C(pR), axis=1)
    lr = DL - DR; lrbar = 0.9 * lrbar + 0.1 * lr
    return mx.stack([mx.sum(Rm[:, 0] * (DL + DR), 1), mx.sum(Rm[:, 1] * (lr - lrbar), 1), mx.sum(Rm[:, 2] * (DL + DR), 1)], 1), lrbar
C = mx.array  # shorthand
CS = dict(Bin=C(Bin), Gmap=C(Gmap), Pmb=C(Pmb), Mdan=C(Mdan), mb=C(mb), dan=C(dan), dn=C(dn))

# ======================================================================== world constants
F, O = args.npatch, 8
ROLE = np.zeros(F, int); nneu = F // 4; ROLE[F - nneu - args.nbad:F - nneu] = 1; ROLE[F - nneu:] = 2
odor_vecs = np.zeros((O, NG), np.float32)
for o in range(O): odor_vecs[o, 1 + rng.choice(NG - 1, 4, replace=False)] = rng.uniform(0.5, 1.0, 4)
E = dict(vmax=0.012, wmax=0.35, ant=args.ant, phi=0.6, sig=args.sig, sigp=0.12, reat=0.055, rkill=0.035,
         c_base=args.c_base, c_move=0.0008, c_n=0.0006 * args.cost_scale, c_e=0.0004 * args.cost_scale)
CS["role"] = C(ROLE)
_starve_margin = args.E0 - E["c_base"] * args.T
print(f"energy loophole check: E0 - c_base*T = {_starve_margin:+.2f} ({'OK: must eat to survive' if _starve_margin < -0.1 else 'WARNING: can survive without eating'})", flush=True)
print(f"K={K} NKC={NKC} nMB={nMB} nDN={nDN} NS={NS} F={F} roles={ROLE.tolist()}", flush=True)

# ======================================================================== neural core (shared by all modes)
def kc_layer(rs, Wpk, kc_gain, kc_lambda):
    """rs (n,K) group rates -> sparse KC rates (n,NKC). APL: subtract mean + lambda*std, rectify."""
    u = (rs @ Wpk.T) * kc_gain
    thr = mx.mean(u, 1, keepdims=True) + kc_lambda * mx.sqrt(mx.var(u, 1, keepdims=True) + 1e-8)
    return mx.tanh(mx.maximum(u - thr, 0))

def brain_update(r, Wkm, dbar, S, g, W, eta_mask, etr=None):
    """one recurrent step. W: (n,K,K) or (K,K); Wkm: (n,nMB,NKC) plastic; returns r, Wkm, dbar, kc."""
    rs = r * g["csz"]
    x = (W @ rs[..., None].astype(W.dtype))[..., 0].astype(mx.float32) if W.ndim == 3 else rs @ W.T
    kc = kc_layer(rs, g["Wpk"], g["kc_gain"], g["kc_lambda"])
    x = x + kc @ g["Wko"].T + ((Wkm @ kc[..., None])[..., 0]) @ CS["Pmb"].T
    x = x + (S * (g["gin"] @ CS["Gmap"].T)) @ CS["Bin"].T + g["bias"]
    r = (1 - g["alph"]) * r + g["alph"] * mx.tanh(mx.maximum(x, 0))
    d = mx.take(r, CS["dan"], axis=1) @ CS["Mdan"].T
    dbar = dbar * 0.95 + 0.05 * d
    elig = kc if etr is None else etr
    Wkm = mx.clip(Wkm + (g["eta"] * eta_mask * (d - dbar))[:, :, None] * elig[:, None, :], 0, 2.0)
    return r, Wkm, dbar, kc

# ======================================================================== homeostatic calibration (A and B identical procedure)
def stim_ensemble():
    S = []
    for o in range(O):
        for c in (0.3, 1.0):
            v = np.zeros(NS, np.float32); v[:NG] = odor_vecs[o] * c; v[NG:2 * NG] = odor_vecs[o] * c * 0.8; S.append(v)
    for k in range(2 * NG, NS):
        v = np.zeros(NS, np.float32); v[k] = 1; S.append(v)
    v = np.zeros(NS, np.float32); v[:2 * NG] = 0.5; S.append(v)
    return np.stack(S)

CUR_BIAS = np.zeros(K, np.float32)
def static_genome(Wg, Wpk, n=1, bias=None):
    b_ = CUR_BIAS if bias is None else bias
    return dict(csz=mx.ones((n, K)), alph=mx.full((n, K), 0.27), bias=mx.broadcast_to(C(b_), (n, K)), gin=mx.full((n, 3), 2.0),
                eta=mx.zeros((n, nMB)), kc_gain=mx.ones((n, 1)), kc_lambda=mx.full((n, 1), 1.64),
                Wpk=C(Wpk), Wko=C(Wko_fixed0), W=C(Wg))

def run_static(Wg, Wpk, S, steps=60, kc_gain=1.0, Wkm=None, eta=0.0, eta_mask=1.0, bias=None):
    n = S.shape[0]; g = static_genome(Wg, Wpk, n, bias); g["kc_gain"] = mx.full((n, 1), kc_gain); g["eta"] = mx.full((n, nMB), eta)
    r = mx.zeros((n, K)); dbar = mx.zeros((n, nMB)); Wkm = mx.broadcast_to(C(Wkm0), (n, nMB, NKC)) if Wkm is None else Wkm
    accr = mx.zeros((n, K)); acck = mx.zeros((n, NKC)); Sm = C(S)
    for t in range(steps):
        r, Wkm, dbar, kc = brain_update(r, Wkm, dbar, Sm, g, g["W"], eta_mask)
        if t >= steps // 2: accr = accr + r; acck = acck + kc
        if t % 20 == 19: mx.eval(r, Wkm, dbar)
    h = steps - steps // 2
    return np.array(accr) / h, np.array(acck) / h, Wkm

def homeostatic_global(Wbase, Wpk, target, kc_target=0.05, normalise=True):
    """Identical, exactly-converging procedure for every wiring condition:
    (1) structural normalisation: each non-sensor group's total |input| weight = 1 (groups without inputs untouched);
    (2) ONE global gain g (bisection) so mean group activity over the stimulus ensemble = target;
    (3) ONE KC gain (bisection) so mean KC activity = kc_target. No per-group feedback control."""
    S = stim_ensemble()
    rs_ = np.abs(Wbase).sum(1, keepdims=True); rs_[rs_ == 0] = 1.0
    Wn = (Wbase / rs_ if normalise else Wbase.copy()); Wn[sensor] = 0.0
    zero = np.zeros(K, np.float32)
    def act(g, kg):
        R_, KCr, _ = run_static(Wn * g, Wpk, S, kc_gain=kg, bias=zero); return R_, KCr
    lo, hi = np.log(0.2), np.log(200.0); kg = 50.0
    for rnd in range(3):
        for _ in range(22):
            mid = 0.5 * (lo + hi); R_, KCr = act(np.exp(mid), kg)
            if R_[:, ~sensor].mean() < target: lo = mid
            else: hi = mid
        g = float(np.exp(0.5 * (lo + hi)))
        klo, khi = np.log(0.1), np.log(5000.0)
        for _ in range(22):
            kmid = 0.5 * (klo + khi); _, KCr = act(g, np.exp(kmid))
            if KCr.mean() < kc_target: klo = kmid
            else: khi = kmid
        kg = float(np.exp(0.5 * (klo + khi))); lo, hi = np.log(0.2), np.log(200.0)
    R_, KCr = act(g, kg); m = R_.mean(0)[~sensor]
    info = dict(g=round(g, 3), kc_gain=round(kg, 2), mean=round(float(m.mean()), 4), rel_err=round(abs(float(m.mean()) / target - 1), 4),
                sat_frac=round(float((R_[:, ~sensor] > 0.9).mean()), 4), silent_groups=int((m < 0.005).sum()), group_act_sd=round(float(m.std()), 4),
                kc_mean=round(float(KCr.mean()), 4), kc_active=round(float((KCr > 0.01).mean()), 4))
    return Wn * g, kg, zero.copy(), info

def homeostatic_distmatch(Wbase, Wpk, target, kc_target=0.05, iters=60, damp=0.5):
    """Activity-DISTRIBUTION matched control (reviewer request): after structural normalisation, per-group gains are
    tuned by quantile matching so that every condition has the SAME sorted activity profile (mean, SD, silent fraction),
    not just the same mean. Reference profile is fixed a priori (log-normal, matched mean=target, SD/mean=0.8)."""
    S = stim_ensemble()
    rs_ = np.abs(Wbase).sum(1, keepdims=True); rs_[rs_ == 0] = 1.0
    Wn = Wbase / rs_; Wn[sensor] = 0.0
    nlive = int((~sensor).sum())
    q = (np.arange(nlive) + 0.5) / nlive
    sig = np.sqrt(np.log(1 + 0.8 ** 2)); mu = np.log(target) - sig ** 2 / 2
    from scipy.stats import norm
    REF = np.sort(np.exp(mu + sig * norm.ppf(q))).astype(np.float32)          # reference sorted activity profile
    gr = np.full(K, 4.0, np.float32); kg = 50.0; zero = np.zeros(K, np.float32)
    for it in range(iters):
        R_, KCr, _ = run_static(Wn * gr[:, None], Wpk, S, kc_gain=kg, bias=zero)
        m = R_.mean(0); live = ~sensor
        order = np.argsort(m[live]); tgt = np.zeros(nlive, np.float32); tgt[order] = REF
        full_t = np.zeros(K, np.float32); full_t[live] = tgt
        upd = np.clip((full_t / (m + 1e-3)) ** damp, 0.7, 1.4); upd[sensor] = 1.0
        gr = np.clip(gr * upd, 0.05, 2000).astype(np.float32)
        kg = float(np.clip(kg * np.clip((kc_target / (KCr.mean() + 1e-4)) ** 0.3, 0.8, 1.3), 0.1, 200000))
    R_, KCr, _ = run_static(Wn * gr[:, None], Wpk, S, kc_gain=kg, bias=zero); m = R_.mean(0)[~sensor]
    info = dict(mean=round(float(m.mean()), 4), sd=round(float(m.std()), 4), silent_groups=int((m < 0.005).sum()),
                sat_frac=round(float((R_[:, ~sensor] > 0.9).mean()), 4), kc_mean=round(float(KCr.mean()), 4),
                kc_active=round(float((KCr > 0.01).mean()), 4), profile_err=round(float(np.abs(np.sort(m) - REF).mean()), 4), kc_gain=round(kg, 2))
    return Wn * gr[:, None], kg, zero.copy(), info

def homeostatic_pergroup(Wbase, Wpk, target, iters=150, damp=0.25, kc_target=0.05, tol=2e-3):
    """previous (exploratory-phase) scheme: per-group multiplicative gain feedback, no normalisation."""
    S = stim_ensemble(); gr = np.full(K, 4.0, np.float32); kg = 4.0; zero = np.zeros(K, np.float32); resid = None
    for it in range(iters):
        R_, KCr, _ = run_static(Wbase * gr[:, None], Wpk, S, kc_gain=kg, bias=zero)
        m = R_.mean(0)
        upd = np.clip((target / (m + 1e-3)) ** damp, 0.8, 1.25); upd[sensor] = 1.0
        gr = np.clip(gr * upd, 0.2, 200)
        kg = float(np.clip(kg * np.clip((kc_target / (KCr.mean() + 1e-4)) ** damp, 0.8, 1.25), 0.1, 500))
        live = (m > 1e-3) & ~sensor
        resid = float(np.abs(np.log((m[live] + 1e-3) / target)).mean())
        if it > 20 and np.max(np.abs(np.log(upd[live]))) < tol: break
    R_, KCr, _ = run_static(Wbase * gr[:, None], Wpk, S, kc_gain=kg, bias=zero); m = R_.mean(0)[~sensor]
    info = dict(resid=round(resid, 4), mean=round(float(m.mean()), 4), sat_frac=round(float((R_[:, ~sensor] > 0.9).mean()), 4),
                silent_groups=int((m < 0.005).sum()), group_act_sd=round(float(m.std()), 4), kc_mean=round(float(KCr.mean()), 4))
    return Wbase * gr[:, None], kg, zero.copy(), info

def homeostatic(Wbase, Wpk, target=args.homeo_target, iters=25):
    if args.homeo == "global":
        return homeostatic_global(Wbase, Wpk, target)
    if args.homeo == "rawglobal":
        return homeostatic_global(Wbase, Wpk, target, normalise=False)
    if args.homeo == "pergroup":
        return homeostatic_pergroup(Wbase, Wpk, target)
    if args.homeo == "distmatch":
        return homeostatic_distmatch(Wbase, Wpk, target)
    raise SystemExit("paper protocol requires --homeo global")
    S = stim_ensemble(); gr = np.full(K, 4.0 * args.homeo_init, np.float32); kg = 4.0
    for _ in range(iters):
        R, KCr, _ = run_static(Wbase * gr[:, None], Wpk, S, kc_gain=kg)
        m = R.mean(0); upd = np.clip((target / (m + 1e-3)) ** 0.5, 0.5, 2.0)
        upd[sensor] = 1.0; upd[m < 1e-4] = 1.3
        gr = np.clip(gr * upd, 0.2, 200)
        kact = (KCr > 0.01).mean()                        # KC gain: keep KCs responsive
        kg = float(np.clip(kg * (1.3 if KCr.mean() < 0.02 else 0.9 if KCr.mean() > 0.1 else 1.0), 0.1, 500))
    R, KCr, _ = run_static(Wbase * gr[:, None], Wpk, S, kc_gain=kg)
    info = dict(mean=float(R.mean()), sat=float((R > 0.9).mean()), silent=int((R.mean(0) < 0.01).sum()),
                kc_active=float((KCr > 0.01).mean()), kc_gain=kg)
    return Wbase * gr[:, None], kg, info

def shuffle_control(Wg, Wpk):
    Ws = np.zeros_like(Wg); rows = np.where(~sensor)[0]
    for j in range(K):
        nz = np.nonzero(Wg[:, j])[0]
        if len(nz): Ws[rng.choice(rows, len(nz), replace=False), j] = Wg[nz, j]
    Wps = np.zeros_like(Wpk)
    for j in range(K):                                   # PN->KC: shuffle which KCs receive each group's input
        nz = np.nonzero(Wpk[:, j])[0]
        if len(nz): Wps[rng.choice(NKC, len(nz), replace=False), j] = Wpk[nz, j]
    return Ws, Wps

# ---------------- null models ----------------
import pandas as _pd
_ann = _pd.read_csv("data/annotations.tsv", sep="\t", low_memory=False, usecols=["cell_type", "hemibrain_type", "cell_class", "super_class"])
_k = _ann.cell_type.fillna(_ann.hemibrain_type).fillna(_ann.cell_class)
_sc = _pd.Series(_ann.super_class.values, index=_k.values)
_sc = _sc[~_sc.index.duplicated()]
def _grp_block(nm):
    base, sd = (nm.split("|") + ["C"])[:2]
    b = base.split("_")[0] if base.startswith(("ORN", "GRN")) else base
    sup = "sensory" if base.startswith(("ORN", "GRN")) else str(_sc.get(base, "central"))
    return f"{sup}|{sd}"
GBLOCK = np.array([_grp_block(n_) for n_ in names]); KCSIDE = np.array(["L"] * (NKC // 2) + ["R"] * (NKC - NKC // 2))
print("null blocks:", len(set(GBLOCK)), "group blocks", flush=True)

def null_degswap(Wg, Wpk, nrng, rounds=10):
    """Maslov-Sneppen target swaps: keeps every group's in- and out-degree, each weight stays with its presynaptic group (Dale sign)."""
    def swap(M, rows_allowed):
        M = M.copy(); ii, jj = np.nonzero(M); w = M[ii, jj].copy(); occ = set(zip(ii.tolist(), jj.tolist())); ne = len(ii); ok_rows = set(rows_allowed.tolist())
        for _ in range(rounds * ne):
            a, b = nrng.integers(0, ne, 2)
            i1, j1, i2, j2 = ii[a], jj[a], ii[b], jj[b]
            if i1 == i2 or j1 == j2 or (i2, j1) in occ or (i1, j2) in occ: continue
            occ.discard((i1, j1)); occ.discard((i2, j2)); occ.add((i2, j1)); occ.add((i1, j2))
            ii[a], ii[b] = i2, i1
        out = np.zeros_like(M); out[ii, jj] = w; return out
    return swap(Wg, np.where(~sensor)[0]), swap(Wpk, np.arange(Wpk.shape[0]))

def null_block(Wg, Wpk, nrng):
    """permute weights within (post block, pre block) pairs; block = super_class x hemisphere. Keeps block-level wiring & laterality."""
    Wn = np.zeros_like(Wg); rows_ok = ~sensor
    for pb in set(GBLOCK):
        R_ = np.where((GBLOCK == pb) & rows_ok)[0]
        for qb in set(GBLOCK):
            Q = np.where(GBLOCK == qb)[0]
            if len(R_) == 0 or len(Q) == 0: continue
            sub = Wg[np.ix_(R_, Q)].ravel(); Wn[np.ix_(R_, Q)] = nrng.permutation(sub).reshape(len(R_), len(Q))
    Pn = np.zeros_like(Wpk)
    for ks in ("L", "R"):
        R_ = np.where(KCSIDE == ks)[0]
        for qb in set(GBLOCK):
            Q = np.where(GBLOCK == qb)[0]
            sub = Wpk[np.ix_(R_, Q)].ravel(); Pn[np.ix_(R_, Q)] = nrng.permutation(sub).reshape(len(R_), len(Q))
    return Wn, Pn

def _interior(Wg):
    """interior edges = source is not a sensory group AND target is not a descending-neuron group.
    Sensory->first-layer and last-layer->descending wiring (the agent's interface) stays exactly as in the connectome."""
    dnm = np.zeros(K, bool); dnm[np.asarray(dn, int)] = True
    return (Wg != 0) & (~sensor)[None, :] & (~dnm)[:, None], dnm

def null_interface_degswap(Wg, Wpk, nrng, rounds=10):
    """N4: Maslov-Sneppen swaps restricted to interior edges; interface edges untouched; PN->KC swapped as in N2."""
    msk, dnm = _interior(Wg)
    ii, jj = np.nonzero(msk); w = Wg[ii, jj].copy(); occ = set(zip(*np.nonzero(Wg != 0))); occ = {(int(a), int(b)) for a, b in occ}; ne = len(ii)
    for _ in range(rounds * ne):
        a, b = nrng.integers(0, ne, 2)
        i1, j1, i2, j2 = ii[a], jj[a], ii[b], jj[b]
        if i1 == i2 or j1 == j2 or (i2, j1) in occ or (i1, j2) in occ: continue
        occ.discard((i1, j1)); occ.discard((i2, j2)); occ.add((i2, j1)); occ.add((i1, j2))
        ii[a], ii[b] = i2, i1
    out = np.where(msk, 0.0, Wg).astype(Wg.dtype); out[ii, jj] = w
    return out, null_degswap(Wg, Wpk, nrng, rounds)[1]

def null_interface_colshuffle(Wg, Wpk, nrng):
    """N5: for every non-sensory source, move its interior outputs to random interior targets (non-sensory, non-descending rows);
    interface edges untouched; PN->KC swapped as in N2."""
    msk, dnm = _interior(Wg)
    out = np.where(msk, 0.0, Wg).astype(Wg.dtype); pool = np.where(~sensor & ~dnm)[0]
    for j in range(K):
        nzr = np.nonzero(msk[:, j])[0]
        if not len(nzr): continue
        free = pool[out[pool, j] == 0]
        tgt = nrng.choice(free, len(nzr), replace=False); out[tgt, j] = Wg[nzr, j]
    return out, null_degswap(Wg, Wpk, nrng)[1]

def make_condition(cond, nrng):
    if cond in ("N4", "N5"):
        Wn, Pn = (null_interface_degswap if cond == "N4" else null_interface_colshuffle)(W0, Wpk0, nrng)
        msk, dnm = _interior(W0)
        iface = np.ones_like(W0, bool); iface[msk] = False; iface &= ~_interior(Wn)[0]
        d_if = float(np.abs(Wn[iface] - W0[iface]).max())
        same_deg = bool(np.array_equal((W0 != 0).sum(0), (Wn != 0).sum(0))) if cond == "N4" else None
        print(f"{cond}: interface max|diff| vs connectome {d_if:.2e}, nnz {int((W0 != 0).sum())}->{int((Wn != 0).sum())}, "
              f"interior overlap {float(((W0 != 0) & (Wn != 0) & msk).sum() / max(msk.sum(), 1)):.3f}, out-degree kept {same_deg}", flush=True)
        assert d_if == 0.0, "interface changed"
        return Wn, Pn
    if cond == "A": return W0.copy(), Wpk0.copy()
    if cond == "N1": return shuffle_control(W0, Wpk0)
    if cond == "N2": return null_degswap(W0, Wpk0, nrng)
    if cond == "N3": return null_block(W0, Wpk0, nrng)
    raise ValueError(cond)

def odorS(o, c=1.0):
    v = np.zeros(NS, np.float32); v[:NG] = odor_vecs[o] * c; v[NG:2 * NG] = odor_vecs[o] * c; return v

# ======================================================================== circuit unit tests
def circuit_tests(label, Wg, Wpk, kg):
    cos = lambda u, v: float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-9))
    S = np.stack([odorS(o) for o in range(O)])
    R, KCr, _ = run_static(Wg, Wpk, S, kc_gain=kg)
    res = {}
    pn = B["pn"]; orn = np.array(in_idx[:2 * NG])
    for st, X in (("ORN", R[:, orn]), ("PN", R[:, pn]), ("KC", KCr), ("MBON", R[:, mb]), ("DN", R[:, dn])):
        res[f"cos_{st}"] = round(np.mean([cos(X[i], X[j]) for i in range(O) for j in range(i + 1, O)]), 3)
    res["kc_active_frac"] = round(float((KCr > 0.01).mean()), 3)
    Sd = np.zeros((3, NS), np.float32); Sd[1, I_SUG:I_SUG + 2] = 1; Sd[2, I_BIT:I_BIT + 2] = 1
    Rd, _, _ = run_static(Wg, Wpk, Sd, kc_gain=kg)
    res["DAN_sugar"] = round(float((Rd[1, dan] - Rd[0, dan]).max()), 3); res["DAN_malaise"] = round(float((Rd[2, dan] - Rd[0, dan]).max()), 3)
    # differential conditioning: odor0 + malaise, odor1 never trained; test both from fresh state
    for eta in (0.05, 0.2, -0.05, -0.2):
        pre, _, _ = run_static(Wg, Wpk, S[:2], kc_gain=kg)
        Sp = np.stack([np.concatenate([odorS(0)[:I_BIT], [1, 1], odorS(0)[I_BIT + 2:]]).astype(np.float32)])
        Sp[0, I_BIT:I_BIT + 2] = 1
        # phase: odor alone 40 -> odor+malaise 60 -> blank 40, plasticity on
        g1 = None; Wkm = None
        _, _, Wkm = run_static(Wg, Wpk, np.stack([odorS(0)]), steps=40, kc_gain=kg, eta=eta)
        _, _, Wkm = run_static(Wg, Wpk, Sp, steps=60, kc_gain=kg, Wkm=Wkm, eta=eta)
        _, _, Wkm = run_static(Wg, Wpk, np.zeros((1, NS), np.float32), steps=40, kc_gain=kg, Wkm=Wkm, eta=eta)
        post, _, _ = run_static(Wg, Wpk, S[:2], kc_gain=kg, Wkm=mx.concatenate([Wkm, Wkm]))
        c0 = np.abs(post[0, mb] - pre[0, mb]).sum(); c1 = np.abs(post[1, mb] - pre[1, mb]).sum()
        res[f"spec_eta{eta:+}"] = f"{c0 / (c1 + 1e-9):.2f} ({c0:.2f}/{c1:.2f})"
    # lateralization: does left-vs-right odor / loom change the steering readout (DNa01/02 L - R)?
    turnw = None
    SL = np.zeros((6, NS), np.float32)
    for i_, (l_, r_) in enumerate(((1.0, 0.0), (0.0, 1.0), (1.0, 0.8), (0.8, 1.0))):
        SL[i_, :NG] = odor_vecs[0] * l_; SL[i_, NG:2 * NG] = odor_vecs[0] * r_
    SL[4, I_LOOM:I_LOOM + 2] = 1; SL[5, I_LOOM + 2:I_LOOM + 4] = 1
    RL, _, _ = run_static(Wg, Wpk, SL, kc_gain=kg)
    t_ = (RL[:, dn][:, pL] - RL[:, dn][:, pR]) @ R_prior[1]
    dnL = [i for i, (b, s_) in enumerate(dn_names) if s_ == "L"]; dnR = [i for i, (b, s_) in enumerate(dn_names) if s_ == "R"]
    res["turn_odorL_minus_odorR"] = round(float(t_[0] - t_[1]), 3); res["turn_odor_20pct_asym"] = round(float(t_[2] - t_[3]), 3)
    res["turn_loomL_minus_loomR"] = round(float(t_[4] - t_[5]), 3)
    res["DN_L-R_act_odorL"] = round(float(RL[0, dn][dnL].mean() - RL[0, dn][dnR].mean()), 3)
    res["DN_L-R_act_odorR"] = round(float(RL[1, dn][dnL].mean() - RL[1, dn][dnR].mean()), 3)
    res["steerDN_act_odorL"] = [round(float(x), 3) for x in RL[0, dn][[i for i, (b, s_) in enumerate(dn_names) if b in ("DNa01", "DNa02")]]]
    res["steerDN_act_odorR"] = [round(float(x), 3) for x in RL[1, dn][[i for i, (b, s_) in enumerate(dn_names) if b in ("DNa01", "DNa02")]]]
    # best linear lateral readout available in DN layer (upper bound): correlation of (L-R) DN pattern
    dvec = RL[0, dn] - RL[1, dn]; res["DN_pattern_diff_norm_odor"] = round(float(np.linalg.norm(dvec)), 3)
    res["DN_pattern_diff_norm_20pct"] = round(float(np.linalg.norm(RL[2, dn] - RL[3, dn])), 3)
    print(f"CIRCUIT {label}: " + json.dumps(res), flush=True)
    return res

# ======================================================================== genome
def innate_eat_prior(Wg, Wpk, kg):
    """PER-like innate reflex: eat readout = DNs this brain's own sugar input excites (with odor background);
    bias set midway between odor-only and odor+sugar so the reflex fires on sweet contact only."""
    ens = stim_ensemble()[:2 * O]                              # odor backgrounds
    sug = ens.copy(); sug[:, I_SUG:I_SUG + 2] = 1
    bit = ens.copy(); bit[:, I_BIT:I_BIT + 2] = 1
    Ro, _, _ = run_static(Wg, Wpk, ens, kc_gain=kg); Rs, _, _ = run_static(Wg, Wpk, sug, kc_gain=kg); Rb_, _, _ = run_static(Wg, Wpk, bit, kc_gain=kg)
    sym = lambda X: X[:, dn][:, pL] + X[:, dn][:, pR]
    w = (sym(Rs) - sym(Ro)).mean(0) - 0.5 * (sym(Rb_) - sym(Ro)).mean(0)
    w = 8 * args.eat_scale * w / (np.abs(w).max() + 1e-6)
    zo, zs = sym(Ro) @ w, sym(Rs) @ w
    bias = -float((zo.mean() + zs.mean()) / 2)
    sep = float((zs.min() - zo.max()))
    print(f"  innate eat reflex: odor-only logit {zo.mean() + bias:+.2f}  sugar logit {zs.mean() + bias:+.2f}  worst-case separation {sep:+.2f}", flush=True)
    return w.astype(np.float32), bias

def init_island(Wg, Wpk, kg, bias, n):
    Rp = R_prior.copy() * 8; w, eb = innate_eat_prior(Wg, Wpk, kg); Rp[2] = w
    return dict(W=mx.broadcast_to(C(Wg), (n, K, K)).astype(DT),
                Wpk=C(Wpk), kc_gain0=kg,
                size=mx.zeros((n, K)), alpha=mx.zeros((n, K)), bias=mx.broadcast_to(C(bias), (n, K)) * 1.0,
                eta=mx.random.normal((n, nMB)) * 0.02 + args.eta0, gin=mx.full((n, 3), 2.0),
                kcg=mx.zeros((n, 1)), kcl=mx.full((n, 1), 1.64),
                R=C(np.broadcast_to(Rp, (n, 3, NP)).copy()) + mx.random.normal((n, 3, NP)) * 0.3,
                RM=mx.random.normal((n, 2, nMB)) * args.rm0,
                Rb=C(np.broadcast_to(np.array([0, 0, eb], np.float32), (n, 3)).copy()))

GENE_KEYS = ["W", "size", "alpha", "bias", "eta", "gin", "kcg", "kcl", "R", "Rb", "RM"]

def mutate(g, n):
    W = g["W"].astype(mx.float32).reshape(n, K * K)
    M1, M2, M3 = 3000, 3, 60
    i1 = mx.random.randint(0, K * K, (n, M1))
    W = mx.put_along_axis(W, i1, mx.take_along_axis(W, i1, 1) * mx.exp(mx.random.normal((n, M1)) * 0.25), 1)
    i2 = mx.random.randint(0, K * K, (n, M2)); cur = mx.take_along_axis(W, i2, 1)
    rowok = mx.take(C((~sensor).astype(np.float32)), i2 // K)
    W = mx.put_along_axis(W, i2, mx.where((cur == 0) & (rowok > 0), mx.take(C(dsign), i2 % K) * mx.abs(mx.random.normal((n, M2))) * 0.05, cur), 1)
    i3 = mx.random.randint(0, K * K, (n, M3)); cur = mx.take_along_axis(W, i3, 1)
    W = mx.put_along_axis(W, i3, mx.where(mx.random.uniform(shape=(n, M3)) < 0.5, 0.0, cur), 1)
    W = mx.where(mx.abs(W) < 1e-3, 0.0, W)
    def pert(x, sig, frac, lo, hi):
        return mx.clip(x + (mx.random.uniform(shape=x.shape) < frac) * mx.random.normal(x.shape) * sig, lo, hi)
    out = dict(g)
    out.update(W=W.reshape(n, K, K).astype(DT), size=pert(g["size"], 0.15, 0.05, -2, 2), alpha=pert(g["alpha"], 0.3, 0.05, -4, 4),
               bias=pert(g["bias"], 0.1, 0.05, -2, 2), eta=pert(g["eta"], 0.03, 0.2, -1, 1), gin=pert(g["gin"], 0.2, 0.3, 0, 8),
               kcg=pert(g["kcg"], 0.1, 0.3, -2, 2), kcl=pert(g["kcl"], 0.1, 0.3, 0.5, 3.0),
               R=pert(g["R"], 0.5, 0.2, -12, 12), RM=pert(g["RM"], 0.5, 0.2, -10, 10), Rb=pert(g["Rb"], 0.2, 0.3, -5, 5))
    return out

# ======================================================================== world step
def env_init(n, ptox=1.0):
    ro = np.stack([rng.choice(O, 3, replace=False) for _ in range(n)])
    return dict(pos=mx.random.uniform(0.2, 0.8, (n, 2)), th=mx.random.uniform(0, 2 * math.pi, (n,)),
                E=mx.full((n,), args.E0), alive=mx.ones((n,)), r=mx.zeros((n, K)), dbar=mx.zeros((n, nMB)),
                etr=mx.zeros((n, NKC)), memv=mx.zeros((n, NG - 1)), memf=mx.zeros((n, NG - 1)), memg=mx.zeros((n, NG - 1)), memb=mx.zeros((n, NG - 1)),
                ptim=mx.zeros((n, F)), lastg=mx.zeros((n, 2)), hasg=mx.zeros((n,)), places=mx.zeros((n, 4, 2)), pvalid=mx.zeros((n, 4)), ptr=mx.zeros((n,), dtype=mx.int32), tgt=mx.zeros((n,), dtype=mx.int32), lrbar=mx.zeros((n, NP)), ppos=mx.random.uniform(0.05, 0.95, (n, F, 2)), pres=mx.where(C(ROLE)[None] == 1, float(args.bite_bad), float(args.bite)) * mx.ones((n, F)),
                qpos=mx.random.uniform(0, 1, (n, 2)), taste=mx.zeros((n, 2)), mal=mx.zeros((n,)),
                good=mx.zeros((n, 2)), bad=mx.zeros((n, 2)), steps=mx.zeros((n,)), caught=mx.zeros((n,)),
                rvec=C(odor_vecs[ro]), role=mx.where((C(ROLE)[None] == 1) & (mx.random.uniform(shape=(n, F)) >= (C(ptox)[:, None] if isinstance(ptox, np.ndarray) else ptox)), 0, C(ROLE)[None]).astype(mx.int32), dist=mx.zeros((n,)), eatt=mx.zeros((n,)), walls=mx.zeros((n,)), onpatch=mx.zeros((n,)), dnact=mx.zeros((n,)))

def world_step(st, Wkm, g, nz, half, pred_v, mask, pol, cg):
    alive, pos, th = st["alive"], st["pos"], st["th"]
    ang = mx.stack([th + E["phi"], th - E["phi"]], 1)
    antp = pos[:, None, :] + E["ant"] * mx.stack([mx.cos(ang), mx.sin(ang)], -1)
    d2 = mx.sum((antp[:, :, None, :] - st["ppos"][:, None]) ** 2, -1)
    conc = mx.exp(-d2 / (2 * E["sig"] ** 2)) * (st["pres"] > 0)[:, None, :]
    orn = mx.sum(conc[..., None] * mx.take_along_axis(st["rvec"], mx.broadcast_to(st["role"][..., None], (st["role"].shape[0], F, NG)), axis=1)[:, None], 2)
    dqa = mx.sum((antp - st["qpos"][:, None, :]) ** 2, -1)
    co2 = mx.exp(-dqa / (2 * E["sigp"] ** 2)) * 1.5
    orn = mx.concatenate([orn[..., :1] + co2[..., None], orn[..., 1:]], -1)
    dq = st["qpos"] - pos; dist = mx.sqrt(mx.sum(dq ** 2, -1) + 1e-6); rel = mx.arctan2(dq[:, 1], dq[:, 0]) - th
    loom = mx.clip(0.04 / dist, 0, 1) * (mx.cos(rel) > -0.3)
    lL, lR = loom * (mx.sin(rel) > 0), loom * (mx.sin(rel) <= 0)
    S = mx.concatenate([orn[:, 0] * mask[:, :1], orn[:, 1] * mask[:, :1], st["taste"][:, :1], st["taste"][:, :1],
                        st["taste"][:, 1:], st["taste"][:, 1:], mx.stack([lL, lL, lR, lR], 1) * mask[:, 1:2]], 1)
    if args.nobrain:
        r, dbar, kcn = st["r"], st["dbar"], st["etr"]
    else:
        r, Wkm, dbar, kcn = brain_update(st["r"], Wkm, st["dbar"], S, g, g["W"], mask[:, 2:3], st["etr"])
    etr = args.trace * st["etr"] + kcn
    act, lrbar = dn_readout(g["R"], mx.take(r, CS["dn"], axis=1), st["lrbar"]); act = act + g["Rb"]
    mbo = (g["RM"] @ mx.take(r, CS["mb"], axis=1)[..., None])[..., 0]          # MBON valence -> forward, eat
    act = act + mx.stack([mbo[:, 0], mx.zeros_like(mbo[:, 0]), mbo[:, 1]], 1)
    # scripted reference policies (pol: 0 brain, 1 random walk, 2 reflex, 3 chemotaxis, 4 oracle memory)
    sweet = st["taste"][:, 0]; oL, oR = mx.sum(orn[:, 0, 1:], -1), mx.sum(orn[:, 1, 1:], -1)
    gv, bv = st["rvec"][:, 0, 1:], st["rvec"][:, 1, 1:]
    gL, gR = mx.sum(orn[:, 0, 1:] * gv, -1), mx.sum(orn[:, 1, 1:] * gv, -1)
    bL, bR = mx.sum(orn[:, 0, 1:] * bv, -1), mx.sum(orn[:, 1, 1:] * bv, -1)
    avoid = 4.0 * (lR - lL)
    eat_reflex = mx.where(sweet > 0, 3.0, -3.0)
    eat_oracle = mx.where((sweet > 0) & (gL + gR > bL + bR), 3.0, -3.0)
    a1 = mx.stack([mx.full(pol.shape, 1.0), nz[:, 0] * 2 - 1, mx.full(pol.shape, -3.0)], 1)
    a2 = mx.stack([mx.full(pol.shape, 1.5), nz[:, 0] * 2 - 1 + avoid, eat_reflex], 1)
    a3 = mx.stack([mx.full(pol.shape, 1.5), cg * (oL - oR) / (oL + oR + 1e-3) * (oL + oR > 0.02) + avoid, eat_reflex], 1)
    # pol 5: raw (un-normalised) bilateral difference, like a linear neuron: turn = cg * (oL - oR)
    a5 = mx.stack([mx.full(pol.shape, 1.5), cg * (oL - oR) + avoid, eat_reflex], 1)
    a4 = mx.stack([mx.full(pol.shape, 1.5), 30.0 * ((gL - bL) - (gR - bR)) / (gL + gR + bL + bR + 1e-3) * (gL + gR > 0.02) + avoid, eat_oracle], 1)
    P = pol[:, None]
    act = mx.where(P == 5, a5, act)
    # pol 6: orthokinesis (slow down in odor, symmetric - no left/right needed); pol 7: odor-specific kinesis with oracle memory
    otot = oL + oR
    a6 = mx.stack([1.5 - cg * mx.tanh(3 * otot), (nz[:, 0] * 2 - 1) * (1 + 2 * mx.tanh(3 * otot)) + avoid, eat_reflex], 1)
    gb = (gL + gR) - (bL + bR)
    a7 = mx.stack([1.5 - cg * mx.tanh(3 * mx.maximum(gb, 0)) + cg * mx.tanh(3 * mx.maximum(-gb, 0)), (nz[:, 0] * 2 - 1) * (1 + 2 * mx.tanh(3 * mx.maximum(gb, 0))) + avoid, eat_oracle], 1)
    act = mx.where(P == 6, a6, mx.where(P == 7, a7, act))
    # pol 8: kinesis (symmetric) + innate reflex eat; pol 9: kinesis + LEARNED eat decision (oracle odor value, no steering)
    gh, bh = gL + gR, bL + bR
    eat_mem = mx.where((sweet > 0) & (gh >= bh), 3.0, -3.0)
    kin_f = 1.5 - 1.0 * mx.tanh(3 * otot); kin_t = (nz[:, 0] * 2 - 1) * (1 + 2 * mx.tanh(3 * otot)) + avoid
    a8 = mx.stack([kin_f, kin_t, eat_reflex], 1); a9 = mx.stack([kin_f, kin_t, eat_mem], 1)
    act = mx.where(P == 8, a8, mx.where(P == 9, a9, act))
    eat_sick = mx.where((sweet > 0) & (st["mal"] < 0.05), 3.0, -3.0)
    # pol 11: realistic learner - remembers odor present when malaise started, refuses matching odor afterwards
    loc = orn[:, 0, 1:] + orn[:, 1, 1:]
    mv = st["memv"]; cosm = mx.sum(loc * mv, -1) / (mx.sqrt(mx.sum(loc * loc, -1) * mx.sum(mv * mv, -1)) + 1e-6)
    eat_real = mx.where((sweet > 0) & ~((mx.sum(mv * mv, -1) > 0) & (cosm > 0.85)), 3.0, -3.0)
    a10 = mx.stack([kin_f, kin_t, eat_sick], 1)
    act = mx.where(P == 10, a10, act)
    act = mx.where(P == 11, mx.stack([kin_f, kin_t, eat_real], 1), act)
    mf = st["memf"]; cosf = mx.sum(loc * mf, -1) / (mx.sqrt(mx.sum(loc * loc, -1) * mx.sum(mf * mf, -1)) + 1e-6)
    eat_rigid = mx.where((sweet > 0) & ~((mx.sum(mf * mf, -1) > 0) & (cosf > 0.85)), 3.0, -3.0)
    act = mx.where(P == 13, mx.stack([kin_f, kin_t, eat_rigid], 1), act)
    ln = loc / (mx.sqrt(mx.sum(loc * loc, -1, keepdims=True)) + 1e-6)
    sg = mx.sum(ln * st["memg"], -1) / (mx.sqrt(mx.sum(st["memg"] ** 2, -1)) + 1e-6)
    sb = mx.sum(ln * st["memb"], -1) / (mx.sqrt(mx.sum(st["memb"] ** 2, -1)) + 1e-6)
    hasb = mx.sum(st["memb"] ** 2, -1) > 0
    eat_two = mx.where((sweet > 0) & (~hasb | (sg >= sb)), 3.0, -3.0)
    act = mx.where(P == 14, mx.stack([kin_f, kin_t, eat_two], 1), act)
    # pol 15: fly local search (path integration back to last food spot, Kim & Dickinson 2017)
    def steer_to(xy):
        d_ = xy - pos; ang_ = mx.arctan2(d_[:, 1], d_[:, 0]); err = mx.arctan2(mx.sin(ang_ - th), mx.cos(ang_ - th))
        return 3.0 * err, mx.sqrt(mx.sum(d_ * d_, -1))
    tl, dl = steer_to(st["lastg"])
    go_back = (st["hasg"] > 0) & (sweet == 0)
    a15 = mx.stack([mx.where(go_back & (dl < 0.03), -1.0, kin_f), mx.where(go_back & (dl >= 0.03), tl + avoid, kin_t), eat_reflex], 1)
    act = mx.where(P == 15, a15, act)
    # pol 16: trapline - remembers up to 4 good places, visits them in turn
    tgt_xy = mx.take_along_axis(st["places"], mx.broadcast_to(st["tgt"][:, None, None], (st["tgt"].shape[0], 1, 2)), axis=1)[:, 0]
    tv = mx.take_along_axis(st["pvalid"], st["tgt"][:, None], axis=1)[:, 0]
    tt, dt_ = steer_to(tgt_xy)
    go_tl = (tv > 0) & (sweet == 0)
    a16 = mx.stack([kin_f, mx.where(go_tl, tt + avoid, kin_t), eat_reflex], 1)
    act = mx.where(P == 16, a16, act)
    act = mx.where(P == 1, a1, mx.where(P == 2, a2, mx.where(P == 3, a3, mx.where(P == 4, a4, act))))
    eat = act[:, 2] > 0
    v = E["vmax"] * mx.sigmoid(act[:, 0]) * mx.where(eat, 0.2, 1.0) * alive
    th = th + (E["wmax"] * mx.tanh(act[:, 1]) + 0.05 * (nz[:, 0] - 0.5)) * alive   # defect 28 fixed: zero-mean turning noise
    pos = mx.clip(pos + v[:, None] * mx.stack([mx.cos(th), mx.sin(th)], -1), 0, 1)
    dq = pos - st["qpos"]; dist = mx.sqrt(mx.sum(dq ** 2, -1) + 1e-6)
    wand = mx.stack([mx.cos(nz[:, 1] * 6.283), mx.sin(nz[:, 1] * 6.283)], -1)
    qpos = mx.clip(st["qpos"] + mx.where((dist < 0.45)[:, None], dq / dist[:, None], 0.3 * wand) * pred_v, 0, 1)
    caught = (mx.sqrt(mx.sum((pos - qpos) ** 2, -1)) < E["rkill"]) & (pred_v > 0)
    pd = mx.sqrt(mx.sum((st["ppos"] - pos[:, None]) ** 2, -1))
    ok = (pd < E["reat"]) & (st["pres"] > 0) & eat[:, None] & (alive > 0)[:, None]
    j = mx.argmin(mx.where(ok, pd, 9.0), 1)
    bite = (mx.arange(F)[None] == j[:, None]) & ok
    isg = mx.sum(bite * (st["role"] == 0), 1); isb = mx.sum(bite * (st["role"] == 1), 1)
    onp = (pd < E["reat"]) & (st["pres"] > 0)
    tg = mx.max(onp * (st["role"] <= 1), 1)                 # good and toxic both taste sweet
    mal = st["mal"] * 0.93 + isb * args.mal_gain
    memv = mx.where((isb > 0)[:, None], orn[:, 0, 1:] + orn[:, 1, 1:], st["memv"])
    # pol 14: two-memory learner (MB-like approach vs avoid): running average of odor at good bites and at toxic bites
    locn = orn[:, 0, 1:] + orn[:, 1, 1:]
    memg = mx.where((isg > 0)[:, None], 0.8 * st["memg"] + 0.2 * locn, st["memg"])
    memb = mx.where((isb > 0)[:, None], 0.5 * st["memb"] + 0.5 * locn, st["memb"])
    # pol 13: rigid learner - keeps FIRST toxic odor forever
    memf = mx.where(((isb > 0) & (mx.sum(st["memf"] * st["memf"], -1) == 0))[:, None], orn[:, 0, 1:] + orn[:, 1, 1:], st["memf"])
    # P1 unpredictable environment: good/bad odor roles swap at random moments
    flip = (nz[:, -1] < args.flip_rate)[:, None, None]
    rv = st["rvec"]
    rvec_new = mx.concatenate([mx.where(flip, rv[:, 1:2], rv[:, 0:1]), mx.where(flip, rv[:, 0:1], rv[:, 1:2]), rv[:, 2:]], 1)                       # delayed malaise
    pres = st["pres"] - bite
    full = mx.where(st["role"] == 1, float(args.bite_bad), float(args.bite))
    if args.regrow:
        newly = (pres <= 0) & (st["ptim"] <= 0)
        ptim = mx.where(newly, float(args.regrow_time), mx.maximum(st["ptim"] - 1, 0))
        refill = (st["ptim"] == 1)
        pres = mx.where(refill, full, mx.maximum(pres, 0))
        ppos = st["ppos"]
    else:
        resp = pres <= 0
        ppos = mx.where(resp[..., None], nz[:, 2:2 + 2 * F].reshape(-1, F, 2) * 0.9 + 0.05, st["ppos"])
        pres = mx.where(resp, full, pres); ptim = st["ptim"]
    # spatial memories (used only by scripted policies 15/16; brain agents must evolve their own)
    gb_ = (isg > 0)
    lastg = mx.where(gb_[:, None], pos, st["lastg"]); hasg = mx.maximum(st["hasg"], gb_.astype(mx.float32))
    dpl = mx.sqrt(mx.sum((st["places"] - pos[:, None, :]) ** 2, -1)) + (1 - st["pvalid"]) * 9
    newplace = gb_ & (mx.min(dpl, 1) > 0.1)
    slot = (mx.arange(4)[None] == st["ptr"][:, None]) & newplace[:, None]
    places = mx.where(slot[..., None], pos[:, None, :], st["places"]); pvalid = mx.maximum(st["pvalid"], slot.astype(mx.float32))
    ptr = mx.where(newplace, (st["ptr"] + 1) % 4, st["ptr"])
    arrived = (mx.sqrt(mx.sum((tgt_xy - pos) ** 2, -1)) < 0.03) | (tv <= 0)
    tgt = mx.where(arrived & (sweet == 0), (st["tgt"] + 1) % 4, st["tgt"])
    En = st["E"] + args.gain_good * isg - args.gain_bad * isb - (E["c_base"] + E["c_move"] * v / E["vmax"] + g["cbrain"]) * alive
    hv = mx.stack([1 - half, half])
    return dict(ptim=ptim, lastg=lastg, hasg=hasg, places=places, pvalid=pvalid, ptr=ptr, tgt=tgt, memg=memg, memb=memb, memf=memf, memv=memv, etr=etr, lrbar=lrbar, pos=pos, th=th, E=En, alive=alive * (En > 0) * (1 - caught), r=r, dbar=dbar, ppos=ppos, pres=pres, qpos=qpos,
                taste=mx.stack([tg, mx.clip(mal, 0, 1)], 1).astype(mx.float32), mal=mal,
                good=st["good"] + isg[:, None] * hv[None], bad=st["bad"] + isb[:, None] * hv[None],
                steps=st["steps"] + alive, caught=st["caught"] + alive * caught, rvec=rvec_new, role=st["role"],
                dist=st["dist"] + v, eatt=st["eatt"] + eat * alive, walls=st["walls"] + ((pos <= 0.001) | (pos >= 0.999)).any(1) * alive,
                onpatch=st["onpatch"] + mx.max(onp, 1) * alive, dnact=st["dnact"] + mx.mean(mx.take(r, CS["dn"], axis=1), 1) * alive), Wkm

cstep = mx.compile(world_step)

def express(g, n):
    """genome -> runtime parameters"""
    e = dict(RM=g["RM"], W=g["W"], Wpk=g["Wpk"], Wko=C(Wko_fixed0), csz=mx.exp(g["size"]), alph=mx.sigmoid(g["alpha"] - 1.0), bias=g["bias"],
             eta=g["eta"], gin=g["gin"], kc_gain=g["kc_gain0"] * mx.exp(g["kcg"]), kc_lambda=g["kcl"], R=g["R"], Rb=g["Rb"])
    nedge = mx.sum(mx.abs(g["W"]) > 1e-3, axis=(1, 2)).astype(mx.float32)
    neff = (e["csz"] * C(nneur)).sum(1) + 5177.0
    e["cbrain"] = E["c_n"] * neff / N0 + E["c_e"] * nedge / float((np.abs(W0) > 0).sum())
    return e, neff, nedge

def evaluate(g, n, pred_v, mask=None, swap=None, pol=None, ptox=None, cg=30.0):
    mask = mx.ones((n, 3)) if mask is None else mask
    pol = mx.zeros((n,)) if pol is None else pol
    e, neff, nedge = express(g, n)
    st = env_init(n, np.repeat(PTOX, isz).astype(np.float32) if ptox is None else ptox); Wkm = mx.broadcast_to(C(Wkm0), (n, nMB, NKC)); pv = C(pred_v); T = args.T
    for t in range(T):
        if swap is not None and t == T // 2:
            rv = st["rvec"]; sw = swap[:, None, None] > 0
            st["rvec"] = mx.concatenate([mx.where(sw, rv[:, 1:2], rv[:, 0:1]), mx.where(sw, rv[:, 0:1], rv[:, 1:2]), rv[:, 2:]], 1)
        st, Wkm = cstep(st, Wkm, e, mx.random.uniform(shape=(n, 3 + 2 * F)), C(float(t >= T // 2)), pv, mask, pol, C(cg))
        if t % 25 == 24: mx.eval(st, Wkm)
    mx.eval(st)
    net = args.gain_good * st["good"].sum(1) - args.gain_bad * st["bad"].sum(1)
    fit = net - e["cbrain"] * st["steps"] + 0.3 * st["alive"] + 0.0005 * st["steps"]
    return dict(st=st, fit=fit, alive=st["alive"], good=st["good"], bad=st["bad"], caught=st["caught"], neff=neff, nedge=nedge)

npf = lambda x: np.array(x.astype(mx.float32))

# ======================================================================== setup islands (paper protocol)
LAB = args.conds.split(","); NI = len(LAB); isz = args.pop // NI
PTOX = np.full(NI, args.ptox0, np.float32)
t0 = time.time(); calib = {}; CALIB_INFO = {}
nrng = np.random.default_rng(10_000 + args.seed)            # null realisation differs per seed
for cnd in LAB:
    Wb, Wp = make_condition(cnd, nrng)
    Wg, kg, bias, info = homeostatic(Wb, Wp)
    calib[cnd] = (Wg, Wp, kg, bias); CALIB_INFO[cnd] = info
    print(f"homeostasis {cnd}: {info} | W nnz {(Wb != 0).sum()} Wpk nnz {(Wp != 0).sum()}  ({time.time() - t0:.0f}s)", flush=True)
json.dump({c_: dict(zip(["resid"], [None])) for c_ in LAB}, open(f"{OUT}/calib_placeholder.json", "w"))

if args.circuit:
    cf = open(f"{OUT}/circuit.jsonl", "a")
    for cnd in LAB:
        CUR_BIAS = calib[cnd][3]; res_ = circuit_tests(cnd, *calib[cnd][:3])
        cf.write(json.dumps(dict(seed=args.seed, homeo=args.homeo, cond=cnd, calib=CALIB_INFO[cnd], res=res_)) + "\n")
    raise SystemExit

pops = []
for cnd in LAB:
    CUR_BIAS = calib[cnd][3]
    pops.append(init_island(*calib[cnd], isz))
pop = {k: (mx.concatenate([p_[k] for p_ in pops]) if k in GENE_KEYS else pops[0][k]) for k in pops[0]}
pop["Wpk_idx"] = C(np.repeat(np.arange(NI), isz))
# sparse PN->KC (padded neighbour lists) per condition: avoids NI dense (n,K,NKC) products
_deg = max(int((calib[c_][1] != 0).sum(1).max()) for c_ in LAB)
KCIDX = np.zeros((NI, NKC, _deg), np.int32); KCVAL = np.zeros((NI, NKC, _deg), np.float32)
for ci, c_ in enumerate(LAB):
    Wp = calib[c_][1]
    for k_ in range(NKC):
        nz = np.nonzero(Wp[k_])[0]; KCIDX[ci, k_, :len(nz)] = nz; KCVAL[ci, k_, :len(nz)] = Wp[k_, nz]
KCIDX, KCVAL = C(KCIDX), C(KCVAL)
KG = C(np.array([calib[c_][2] for c_ in LAB], np.float32))
print(f"sparse PN->KC max degree {_deg}", flush=True)
_rs = mx.random.uniform(shape=(8, K)); _idx = C(np.repeat(np.arange(NI), 2))
for _ci, _c in enumerate(LAB):
    _dense = np.array(_rs @ C(calib[_c][1]).T)
    _ix = mx.take(KCIDX, mx.full((8,), _ci, dtype=mx.int32), axis=0).reshape(8, -1); _vv = mx.take(KCVAL, mx.full((8,), _ci, dtype=mx.int32), axis=0).reshape(8, -1)
    _sp = np.array(mx.sum((mx.take_along_axis(_rs, _ix, axis=1) * _vv).reshape(8, NKC, -1), -1))
    print(f"sparse==dense check {_c}: max abs diff {np.abs(_dense - _sp).max():.2e}", flush=True)
del pops
pop.pop("Wpk", None); pop.pop("kc_gain0", None)
GENE_KEYS += ["Wpk_idx"]

def with_shared(g):
    g = dict(g); n_ = g["Wpk_idx"].shape[0]
    g["Wpk"] = ("sparse", mx.take(KCIDX, g["Wpk_idx"], axis=0).reshape(n_, -1), mx.take(KCVAL, g["Wpk_idx"], axis=0).reshape(n_, -1))
    g["kc_gain0"] = mx.take(KG, g["Wpk_idx"])[:, None]
    return g

def kc_layer(rs, Wpk, kc_gain, kc_lambda):
    if isinstance(Wpk, tuple):
        _, ix, vv = Wpk
        u = mx.sum((mx.take_along_axis(rs, ix, axis=1) * vv).reshape(rs.shape[0], NKC, -1), -1) * kc_gain
    else:
        u = (rs @ Wpk.T) * kc_gain
    thr = mx.mean(u, 1, keepdims=True) + kc_lambda * mx.sqrt(mx.var(u, 1, keepdims=True) + 1e-8)
    return mx.tanh(mx.maximum(u - thr, 0))
cstep = mx.compile(world_step)

def take(g, idx):
    i = C(np.asarray(idx)); return {k: mx.take(v, i, axis=0) for k, v in g.items()}

if args.gen0assay:
    # Is the generation-0 advantage due to wiring, or to the per-brain innate feeding readout / to vision?
    out = open(f"{OUT}/gen0assay.jsonl", "a"); reps = 6
    def block(tag, zero_innate=False, mask_cols=(), rw=False):
        sel = np.tile(np.arange(args.pop), reps); n = len(sel)
        gp = with_shared(take(pop, sel))
        if zero_innate:
            Rr = np.array(gp["R"]); Rr[:, 2, :] = 0.0; gp["R"] = C(Rr)
            Rb = np.array(gp["Rb"]); Rb[:, 2] = 0.0; gp["Rb"] = C(Rb)     # eat readout removed -> eats whenever logit noise > 0
        if rw:
            Rr = np.array(gp["R"]); Rr[:] = 0; gp["R"] = C(Rr)
            RMr = np.array(gp["RM"]); RMr[:] = 0; gp["RM"] = C(RMr)
        m = np.ones((n, 3), np.float32)
        for c_ in mask_cols: m[:, c_] = 0
        mx.random.seed(909)
        r = evaluate(gp, n, max(args.pred_fixed, 0.0) * E["vmax"], C(m), None, ptox=np.full(n, args.ptox0, np.float32))
        mx.random.seed(int(time.time() * 1000) % 1000003)
        res = {}
        for i, isl in enumerate(LAB):
            ix = np.concatenate([np.arange(j * args.pop + i * isz, j * args.pop + (i + 1) * isz) for j in range(reps)])
            res[isl] = dict(fit=float(npf(r["fit"])[ix].mean()), surv=float(npf(r["alive"])[ix].mean()),
                            good=float(npf(r["good"])[ix].sum(1).mean()), caught=float(npf(r["caught"])[ix].mean()))
        out.write(json.dumps(dict(homeo=args.homeo, seed=args.seed, assay=tag, res=res)) + "\n"); out.flush()
        print(f"GEN0 {tag}: " + " ".join(f"{k}:{v['fit']:.3f}" for k, v in res.items()), flush=True)
    block("normal")
    block("no_innate_readout", zero_innate=True)
    block("blind", mask_cols=(1,))
    block("anosmic", mask_cols=(0,))
    block("blind_and_no_innate", zero_innate=True, mask_cols=(1,))
    block("random_walk", rw=True)
    raise SystemExit


if args.assay:
    # post-hoc assays on evolved elites: (i) odor-reversal x plasticity (is the nulls' olfaction use associative learning?),
    # (ii) transfer to shifted worlds never seen during evolution.
    import glob
    f_ = f"{args.assay}/elite_g{args.assay_gen:06d}.npz"
    Z = np.load(f_)
    nel = Z[f"{LAB[0]}_W"].shape[0]
    genome = {}
    for k in GENE_KEYS:
        if k == "Wpk_idx": continue
        genome[k] = mx.concatenate([C(Z[f"{isl}_{k}"]) for isl in LAB])
    genome["W"] = genome["W"].astype(DT)
    genome["Wpk_idx"] = C(np.repeat(np.arange(NI), nel))
    n0 = nel * NI
    out = open(f"{args.assay}/assay.jsonl", "a")
    def run_block(tag, reps, ptox, swap_flag, eta_off, mask_cols=(), world=None):
        if world:
            for k_, v_ in world.items():
                if k_ in ("sig",): E["sig"] = v_
                elif k_ == "ptox": pass
                else: setattr(args, k_, v_)
        sel = np.tile(np.arange(n0), reps); n = len(sel)
        gp = with_shared(take(genome, sel))
        m = np.ones((n, 3), np.float32)
        for c_ in mask_cols: m[:, c_] = 0
        if eta_off: m[:, 2] = 0
        mx.random.seed(4242)
        r = evaluate(gp, n, max(args.pred_fixed, 0.0) * E["vmax"], C(m), C(np.full(n, float(swap_flag), np.float32)), ptox=np.full(n, ptox, np.float32))
        mx.random.seed(int(time.time() * 1000) % 1000003)
        res = {}
        for i, isl in enumerate(LAB):
            ix = np.concatenate([np.arange(j * n0 + i * nel, j * n0 + (i + 1) * nel) for j in range(reps)])
            gd, bd = npf(r["good"])[ix], npf(r["bad"])[ix]
            res[isl] = dict(fit=float(npf(r["fit"])[ix].mean()), surv=float(npf(r["alive"])[ix].mean()),
                            bad_h2=float(bd[:, 1].mean()), good_h2=float(gd[:, 1].mean()), bad=float(bd.sum(1).mean()), good=float(gd.sum(1).mean()))
        out.write(json.dumps(dict(run=os.path.basename(args.assay), gen=args.assay_gen, homeo=args.homeo, seed=args.seed, assay=tag, res=res)) + "\n"); out.flush()
        print(f"ASSAY {tag}: " + " ".join(f"{k}:{v['fit']:.3f}/bad_h2 {v['bad_h2']:.2f}" for k, v in res.items()), flush=True)
    R_ = 12
    run_block("base", R_, 0.5, 0, False)
    run_block("base_noplast", R_, 0.5, 0, True)
    run_block("reversal", R_, 0.5, 1, False)
    run_block("reversal_noplast", R_, 0.5, 1, True)
    run_block("transfer_toxic0.9", R_, 0.9, 0, False)
    # ceiling-effect stress grid: harder predator / higher toxicity, no further evolution
    _pf = args.pred_fixed
    for pv_ in (0.12, 0.15, 0.18):
        args.pred_fixed = pv_; run_block(f"stress_pred{pv_}", R_, 0.5, 0, False)
    args.pred_fixed = _pf
    for pt_ in (0.7, 1.0):
        run_block(f"stress_toxic{pt_}", R_, pt_, 0, False)
    # cross-ecology common garden: evaluate every evolved population in all four ecologies
    for pv_ in (0.0, 0.2):
        for pt_ in (0.0, 1.0):
            args.pred_fixed = pv_; run_block(f"eco_P{pv_}_T{pt_}", R_, pt_, 0, False)
    args.pred_fixed = _pf
    _sig0 = E["sig"]; run_block("transfer_plume0.12", R_, 0.5, 0, False, world=dict(sig=0.12)); E["sig"] = _sig0
    _np0 = args.npatch
    raise SystemExit


pred_v = (args.pred_fixed if args.pred_fixed >= 0 else 0.25) * E["vmax"]
json.dump(dict(args=vars(args), E=E, K=K, NKC=NKC, islands=LAB, roles=ROLE.tolist(), protocol="paper-v1"), open(f"{OUT}/config.json", "w"), indent=1)
flog = open(f"{OUT}/gen.jsonl", "a"); plog = open(f"{OUT}/probe.jsonl", "a")
CONDS = ["normal", "no_plasticity", "anosmic", "blind", "random_walk"]

def probe(gen, orders):
    top, reps = 16, args.probe_reps
    ids = np.concatenate([o[:top] for o in orders]); per = top * reps; blk = per * NI
    sel = np.tile(np.repeat(ids, reps), len(CONDS)); n = len(sel)
    gp = with_shared(take(pop, sel))
    m = np.ones((n, 3), np.float32); m[blk:2 * blk, 2] = 0; m[2 * blk:3 * blk, 0] = 0; m[3 * blk:4 * blk, 1] = 0
    swp = np.zeros(n, np.float32)
    pt = np.full(n, args.probe_ptox, np.float32)
    Rr = np.array(gp["R"]); Rr[4 * blk:5 * blk] = 0; gp["R"] = C(Rr)
    Rb = np.array(gp["Rb"]); Rb[4 * blk:5 * blk, :2] = 0; gp["Rb"] = C(Rb)
    RMr = np.array(gp["RM"]); RMr[4 * blk:5 * blk] = 0; gp["RM"] = C(RMr)      # random_walk: no MBON-driven behaviour either
    mx.random.seed(777); r = evaluate(gp, n, pred_v, C(m), C(swp), ptox=pt); mx.random.seed(int(time.time() * 1000) % 1000003)
    out = dict(gen=gen, pred_v=pred_v / E["vmax"])
    for ci, cn in enumerate(CONDS):
        for ii, isl in enumerate(LAB):
            s0 = ci * blk + ii * per; s1 = s0 + per; gd, bd = npf(r["good"][s0:s1]), npf(r["bad"][s0:s1])
            out[f"{isl}.{cn}"] = dict(fit=float(npf(r["fit"][s0:s1]).mean()), surv=float(npf(r["alive"][s0:s1]).mean()),
                                      good=float(gd.sum(1).mean()), bad=float(bd.sum(1).mean()),
                                      good_h2=float(gd[:, 1].mean()), bad_h2=float(bd[:, 1].mean()))
    plog.write(json.dumps(out) + "\n"); plog.flush()

def select(fit, lo, hi, ne):
    order = np.argsort(-fit[lo:hi]) + lo
    cand = rng.integers(lo, hi, (hi - lo - ne, 4))
    return order[:ne], cand[np.arange(len(cand)), np.argmax(fit[cand], 1)], order


if args.ecosearch:
    # env-only random search over world parameters using scripted reference policies (brain off)
    import itertools
    rs = np.random.default_rng(1234); n = 256
    gp = with_shared(take(pop, np.arange(n)))
    out = open(f"{OUT}/ecosearch.jsonl", "a")
    pols = {"noEat": 1, "reflex": 8, "sickStop": 10, "twoMem": 14, "localSearch": 15, "trapline": 16, "oracle": 9}
    for it in range(args.ecosearch):
        cfg = dict(npatch=int(rs.choice([16, 24, 32, 48, 64])), badfrac=float(rs.choice([0.25, 0.4, 0.5])),
                   bite=int(rs.choice([1, 3, 6, 12])), bite_bad=int(rs.choice([1, 2])),
                   gain_good=float(rs.choice([0.03, 0.05, 0.08, 0.12])), gain_bad=float(rs.choice([0.1, 0.2, 0.3, 0.5])),
                   c_base=float(rs.choice([0.0008, 0.0012, 0.0016])), sig=float(rs.choice([0.04, 0.06, 0.1])),
                   T=int(rs.choice([400, 800])), E0=float(rs.choice([0.3, 0.5, 0.7])),
                   regrow=int(rs.choice([0, 1])), regrow_time=int(rs.choice([30, 100, 300])))
        if cfg["E0"] - cfg["c_base"] * cfg["T"] > -0.1:
            continue                                   # energy loophole: could survive without eating
        # apply world params (globals read inside world_step at trace time -> recompile)
        F = cfg["npatch"]; ROLE = np.zeros(F, int); nneu = F // 4; nb = int(round(cfg["badfrac"] * F))
        ROLE[F - nneu - nb:F - nneu] = 1; ROLE[F - nneu:] = 2; CS["role"] = C(ROLE)
        args.npatch, args.nbad, args.bite, args.bite_bad = F, nb, cfg["bite"], cfg["bite_bad"]
        args.gain_good, args.gain_bad, args.T = cfg["gain_good"], cfg["gain_bad"], cfg["T"]
        E["c_base"] = cfg["c_base"]; E["sig"] = cfg["sig"]; args.E0 = cfg["E0"]; args.regrow = cfg["regrow"]; args.regrow_time = cfg["regrow_time"]
        cstep = mx.compile(world_step)
        res = {}
        for fr in (0.0,):
            args.flip_rate = fr
            cstep = mx.compile(world_step)
            for nm, k in pols.items():
                mx.random.seed(7); r = evaluate(gp, n, 0.1 * E["vmax"], pol=mx.full((n,), float(k)), ptox=np.ones(n, np.float32))
                res[f"{nm}@{fr}"] = (round(float(npf(r["fit"]).mean()), 3), round(float(npf(r["alive"]).mean()), 3), round(float(npf(r["bad"]).sum(1).mean()), 2), round(float(npf(r["good"]).sum(1).mean()), 1))
        nonlearn = max(res["noEat@0.0"][0], res["reflex@0.0"][0], res["sickStop@0.0"][0])
        learn = max(res["twoMem@0.0"][0], res["localSearch@0.0"][0], res["trapline@0.0"][0])
        score = dict(learn_margin=round(learn - nonlearn, 3), learn_ratio=round(learn / max(nonlearn, 0.05), 3),
                     learner_surv=max(res["twoMem@0.0"][1], res["localSearch@0.0"][1], res["trapline@0.0"][1]), twomem_ratio=round(res["twoMem@0.0"][0] / max(nonlearn, 0.05), 3),
                     local_ratio=round(res["localSearch@0.0"][0] / max(nonlearn, 0.05), 3), trap_ratio=round(res["trapline@0.0"][0] / max(nonlearn, 0.05), 3), nonlearn_best=round(nonlearn, 3),
                     oracle_ratio=round(res["oracle@0.0"][0] / max(nonlearn, 0.05), 3))
        out.write(json.dumps(dict(cfg=cfg, score=score, res=res)) + "\n"); out.flush()
        print(f"ECO {it} {json.dumps(cfg)} -> {json.dumps(score)}", flush=True)
    raise SystemExit

if args.wlearn:
    # can an agent with a hand-ALIGNED MBON->eat readout learn odor->toxin in this world? (substrate + world learnability)
    n = 512; i0 = 0
    Wg, Wp, kg = calib["A"]
    pre, _, _ = run_static(Wg, Wp, odorS(0)[None], kc_gain=kg)
    Sp = odorS(0)[None].copy(); Sp[0, I_SUG:I_SUG + 2] = 1; Sp[0, I_BIT:I_BIT + 2] = 1
    _, _, Wkm = run_static(Wg, Wp, Sp, steps=80, kc_gain=kg, eta=-0.05)
    post, _, _ = run_static(Wg, Wp, odorS(0)[None], kc_gain=kg, Wkm=Wkm)
    dmb = post[0, mb] - pre[0, mb]; wal = dmb / (np.abs(dmb).max() + 1e-9)
    base = take(pop, np.arange(n))
    for eta in (0.0, -0.02, -0.05, -0.2):
        for scale in (5.0, 20.0):
            g2 = dict(base); RM = np.zeros((n, 2, nMB), np.float32); RM[:, 1] = wal * scale; g2["RM"] = C(RM)
            g2["eta"] = mx.full((n, nMB), eta)
            mx.random.seed(99); r = evaluate(with_shared(g2), n, 0.1 * E["vmax"], ptox=np.ones(n, np.float32))
            bd = npf(r["bad"]); gd = npf(r["good"])
            print(f"WLEARN trace={args.trace} mal={args.mal_gain} eta={eta:+.2f} scale={scale}: fit {float(npf(r['fit']).mean()):.3f} bad h1 {bd[:, 0].mean():.2f} h2 {bd[:, 1].mean():.2f} | good h1 {gd[:, 0].mean():.1f} h2 {gd[:, 1].mean():.1f}", flush=True)
    raise SystemExit

if args.behav:
    # behavioral conditioning: does pairing odor0 with malaise reduce the EAT logit for odor0+sweet but not odor1+sweet?
    for typ in "AB":
        Wg, Wp, kg = calib[typ]
        w, eb = innate_eat_prior(Wg, Wp, kg)
        wm = np.random.default_rng(7).normal(size=nMB).astype(np.float32) * args.rm0
        def eat_logit(o, Wkm=None):
            S = odorS(o)[None].copy(); S[0, I_SUG:I_SUG + 2] = 1
            R_, _, _ = run_static(Wg, Wp, S, kc_gain=kg, Wkm=Wkm)
            return float((R_[0, dn][pL] + R_[0, dn][pR]) @ w + eb + R_[0, mb] @ wm)
        rows = []
        for eta in (-0.02, -0.05, -0.2, -0.5, 0.05, 0.2):
            pre0, pre1 = eat_logit(0), eat_logit(1)
            Sp = odorS(0)[None].copy(); Sp[0, I_SUG:I_SUG + 2] = 1; Sp[0, I_BIT:I_BIT + 2] = 1
            _, _, Wkm = run_static(Wg, Wp, Sp, steps=80, kc_gain=kg, eta=eta)
            _, _, Wkm = run_static(Wg, Wp, np.zeros((1, NS), np.float32), steps=30, kc_gain=kg, Wkm=Wkm, eta=eta)
            p0, p1 = eat_logit(0, Wkm), eat_logit(1, Wkm)
            rows.append(f"eta{eta:+}: odor0 {pre0:+.1f}->{p0:+.1f} (d{p0-pre0:+.2f}) odor1 {pre1:+.1f}->{p1:+.1f} (d{p1-pre1:+.2f})")
        print(f"BEHAV {typ} eat_scale={args.eat_scale}: " + " | ".join(rows), flush=True)
    raise SystemExit

if args.landscape:
    n = 300
    for typ in "AB":
        Wg, Wp, kg = calib[typ]
        SL = np.zeros((2, NS), np.float32); SL[0, :NG] = odor_vecs[0]; SL[0, NG:2*NG] = odor_vecs[0] * 0.8
        SL[1, :NG] = odor_vecs[0] * 0.8; SL[1, NG:2*NG] = odor_vecs[0]
        # average direction over all odors
        dv = np.zeros(NP, np.float32)
        for o in range(O):
            SL[0, :NG] = odor_vecs[o]; SL[0, NG:2*NG] = odor_vecs[o] * 0.8; SL[1, :NG] = odor_vecs[o] * 0.8; SL[1, NG:2*NG] = odor_vecs[o]
            RL, _, _ = run_static(Wg, Wp, SL, kc_gain=kg); lr = RL[:, dn][:, pL] - RL[:, dn][:, pR]; dv += lr[0] - lr[1]
        dv /= (np.linalg.norm(dv) + 1e-9)
        i0 = LAB.index(typ + "1") * isz
        base = take(pop, np.arange(i0, i0 + n)); out = []
        for c in (-100, -30, -10, 0, 10, 30, 100, 300):
            g2 = dict(base); Rn = np.array(base["R"]); Rn[:, 1, :] += c * dv; g2["R"] = C(Rn)
            mx.random.seed(123); r = evaluate(with_shared(g2), n, 0.1 * E["vmax"], ptox=np.zeros(n, np.float32))
            st = r["st"]; out.append(f"c{c}: fit {float(npf(r['fit']).mean()):.2f} good {float(npf(st['good']).sum(1).mean()):.1f}")
        print(f"LANDSCAPE {typ} (pair L-R direction): " + " | ".join(out), flush=True)
    raise SystemExit

if args.sweep:
    n = 400; gp = with_shared(take(pop, np.arange(n)))
    out = []
    for k, lab in ((2, "reflex"),):
        mx.random.seed(123); r = evaluate(gp, n, 0.1 * E["vmax"], pol=mx.full((n,), 2.0), ptox=0.0); out.append(f"reflex {float(npf(r['fit']).mean()):.2f}")
    for cg in ():
        row = []
        for pt in (0.0, 1.0):
            mx.random.seed(123); r6 = evaluate(gp, n, 0.1 * E["vmax"], pol=mx.full((n,), 6.0), ptox=pt, cg=cg)
            mx.random.seed(123); r7 = evaluate(gp, n, 0.1 * E["vmax"], pol=mx.full((n,), 7.0), ptox=pt, cg=cg)
            row.append(f"ptox{pt:.0f} kinesis {float(npf(r6['fit']).mean()):.2f} mem-kinesis {float(npf(r7['fit']).mean()):.2f}")
        out.append(f"k{cg}: " + ", ".join(row))
    for pt in (0.5, 1.0):
        vals = []
        for k, nm in ((8, "kin+reflex"), (10, "kin+sickStop"), (13, "kin+RIGIDlearner"), (11, "kin+FLEXlearner"), (9, "kin+oracleEat")):
            mx.random.seed(123); rr = evaluate(gp, n, 0.1 * E["vmax"], pol=mx.full((n,), float(k)), ptox=pt)
            vals.append(f"{nm} {float(npf(rr['fit']).mean()):.2f}/s{float(npf(rr['alive']).mean()):.2f}/b{float(npf(rr['bad']).sum(1).mean()):.1f}")
        out.append(f"ptox{pt}: " + ", ".join(vals))
    print(f"SWEEP sig={args.sig} ant={args.ant}: " + " | ".join(out), flush=True)
    raise SystemExit

if args.ladder:
    n = 400; gp = with_shared(take(pop, np.arange(n)))
    names_p = ["brain_init", "random_walk", "reflex", "chemotaxis", "oracle_memory"]
    for pv, pt in ((0.15, 0.0), (0.15, 0.5), (0.15, 1.0)):
        rows = []
        for k in range(5):
            mx.random.seed(123)
            r = evaluate(gp, n, pv * E["vmax"], pol=mx.full((n,), float(k)), ptox=pt); st = r["st"]
            rows.append((names_p[k], float(npf(r["fit"]).mean()), float(npf(st["alive"]).mean()), float(npf(st["caught"]).mean()),
                         float(npf(st["good"]).sum(1).mean()), float(npf(st["bad"]).sum(1).mean()), float(npf(st["steps"]).mean())))
        print(f"LADDER pred={pv} ptox={pt}: " + " | ".join(f"{a} fit {b:.3f} surv {c:.2f} caught {d:.2f} good {e:.1f} bad {f_:.1f} life {h:.0f}" for a, b, c, d, e, f_, h in rows), flush=True)
    raise SystemExit

if args.diag:
    for lab in ("A1", "B1"):
        i = LAB.index(lab); n = 512
        gp = with_shared(take(pop, np.arange(i * isz, i * isz + n)))
        for pv in (0.0, 0.25):
            r = evaluate(gp, n, pv * E["vmax"]); st = r["st"]
            f = lambda k: float(npf(st[k]).mean())
            steps = npf(st["steps"])
            caught = npf(st["caught"]) > 0; alive = npf(st["alive"]) > 0; starved = ~alive & ~caught
            print(f"DIAG {lab} pred={pv}: alive {alive.mean():.2f} caught {caught.mean():.2f} starved {starved.mean():.2f} | mean life {steps.mean():.0f} | "
                  f"speed/vmax {f('dist')/max(steps.mean(),1)/E['vmax']:.2f} eat-frac {f('eatt')/max(steps.mean(),1):.2f} wall-frac {f('walls')/max(steps.mean(),1):.2f} "
                  f"onpatch-frac {f('onpatch')/max(steps.mean(),1):.3f} DN act {f('dnact')/max(steps.mean(),1):.3f} | good {float(npf(st['good']).sum(1).mean()):.2f} bad {float(npf(st['bad']).sum(1).mean()):.2f} "
                  f"E_end {float(npf(st['E']).mean()):.2f}", flush=True)
    raise SystemExit

tstart = time.time()
for gen in range(args.gens):
    tg = time.time()
    r = evaluate(with_shared(pop), args.pop, pred_v)
    fit = npf(r["fit"]); ne = max(2, isz // 25)
    sels = [select(fit, i * isz, (i + 1) * isz, ne) for i in range(NI)]
    alive, good, bad, neff, nedge = (npf(r[k]) for k in ("alive", "good", "bad", "neff", "nedge"))
    rec = dict(gen=gen, t=round(time.time() - tstart), sec=round(time.time() - tg, 2), pred_v=round(pred_v / E["vmax"], 3), ptox=[round(float(x), 3) for x in PTOX],
               caught=round(float(npf(r["caught"]).mean()), 3))
    for i, isl in enumerate(LAB):
        s = slice(i * isz, (i + 1) * isz)
        eta_mean = float(npf(pop["eta"][i * isz:(i + 1) * isz]).mean())
        rec[isl] = dict(fit=float(fit[s].mean()), fit_max=float(fit[s].max()), surv=float(alive[s].mean()), good=float(good[s].sum(1).mean()),
                        bad=float(bad[s].sum(1).mean()), neff=float(neff[s].mean()), nedge=float(nedge[s].mean()), eta=eta_mean)
    if args.probe_every > 0 and gen % args.probe_every == 0: probe(gen, [o for _, _, o in sels])
    flog.write(json.dumps(rec) + "\n"); flog.flush()
    if gen % 10 == 0: print(json.dumps(rec)[:260], flush=True)
    if args.tox_curriculum:
        for i, l in enumerate(LAB):                      # per-island curriculum
            if rec[l]["good"] > 10 and rec[l]["surv"] > 0.3: PTOX[i] = min(1.0, PTOX[i] + 0.01)
            elif rec[l]["surv"] < 0.1: PTOX[i] = max(0.0, PTOX[i] - 0.01)
    if args.pred_fixed < 0:
        if rec["caught"] < 0.15: pred_v = min(pred_v + 0.01 * E["vmax"], 1.2 * E["vmax"])
        elif rec["caught"] > 0.35: pred_v = max(pred_v - 0.01 * E["vmax"], 0.05 * E["vmax"])
    if args.migrate_every and gen and gen % args.migrate_every == 0:
        for i in range(NI):
            sis = i + 2 if i + 2 < NI else i - 2
            if 0 <= sis < NI and sis != i: sels[i][1][-2:] = sels[sis][2][:2]
    idx = np.concatenate([np.concatenate([e_, p_]) for e_, p_, _ in sels])
    flag = C(np.concatenate([np.concatenate([np.zeros(len(e_)), np.ones(len(p_))]) for e_, p_, _ in sels]).astype(np.float32))
    ga = take(pop, idx); gm = mutate(ga, len(idx))
    pop = {k: (mx.where(flag.reshape((-1,) + (1,) * (ga[k].ndim - 1)) > 0, gm[k], ga[k]) if k != "Wpk_idx" else ga[k]) for k in ga}
    mx.eval(pop)
    if gen % 250 == 0:
        np.savez_compressed(f"{OUT}/elite_g{gen:06d}.npz", pred_v=pred_v,
                            **{f"{isl}_{k}": npf(mx.take(v, C(sels[i][2][:8]), axis=0)) for i, isl in enumerate(LAB) for k, v in pop.items() if k != "Wpk_idx"})
    if (time.time() - tstart) / 3600 > args.hours: break
print("done", flush=True)
