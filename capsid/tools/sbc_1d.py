"""Localise M2's calibration failure, one parameter at a time, with no MCMC.

For each replicate: draw theta from the prior, simulate a dataset, then free ONE
coordinate and hold the other three at their generating values. The one-dimensional
posterior is evaluated exactly on a grid and normalised by trapezoid, so the rank of
the truth and the coverage of the 94% interval carry no sampling error at all. If
these ranks are uniform the likelihood and the noise model are right and the problem
lies in the joint geometry; if they are not, the offending coordinate is named.
"""
import json, io, contextlib, time, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
from multiprocessing import Pool

ns = {"__name__": "__main__", "display": lambda *a, **k: None}
_nb = json.load(open("BruSafe_03_M2_soluble_partition.ipynb"))
with contextlib.redirect_stdout(io.StringIO()):
    for _i, _c in enumerate(_nb["cells"]):
        if _c["cell_type"] != "code" or _i > 25: continue
        _src = "".join(_c["source"])
        if _i == 25: _src = _src.split("import time; _t0")[0]
        exec(_src, ns)

PARAMS, PRIOR, N_BOUNDS = ns["PARAMS"], ns["PRIOR"], ns["N_BOUNDS"]
predict_grid, beta_logpdf = ns["predict_grid"], ns["beta_logpdf"]
beta_loglik, censor = ns["beta_loglik"], ns["censor"]
TEMPS_C, DOSES, T_OBS_MIN, TGRID = ns["TEMPS_C"], ns["DOSES"], ns["T_OBS_MIN"], ns["_TGRID"]
TI = [int(np.argmin(np.abs(TGRID - t))) for t in T_OBS_MIN]
GRID_N = 401


def draw_prior(rng):
    v = np.array([rng.normal(*PRIOR[k]) for k in PARAMS])
    while not (N_BOUNDS[0] <= v[2] <= N_BOUNDS[1]):
        v[2] = rng.normal(*PRIOR["n"])
    return v


def simulate_dataset(v, rng):
    P = predict_grid(v); nu = np.exp(v[3])
    iT, iD, iTi, y = [], [], [], []
    for k in range(len(TEMPS_C)):
        for i in range(len(DOSES)):
            for tt in TI:
                p = float(np.clip(P[k, i, tt], 1e-6, 1-1e-6))
                for _ in range(2):
                    iT.append(k); iD.append(i); iTi.append(tt)
                    y.append(float(censor(rng.beta(p*nu, (1-p)*nu))))
    return dict(T_idx=np.array(iT), dose_idx=np.array(iD), t_idx=np.array(iTi), y=np.array(y))


def loglik(v, data):
    P = predict_grid(v)
    p = np.clip(P[data["T_idx"], data["dose_idx"], data["t_idx"]], 1e-6, 1-1e-6)
    return float(np.sum(beta_loglik(data["y"], p, np.exp(v[3]))))


def grid_for(j, truth_j):
    m, sd = PRIOR[PARAMS[j]]
    lo, hi = m - 4*sd, m + 4*sd
    if j == 2: lo, hi = max(lo, N_BOUNDS[0]), min(hi, N_BOUNDS[1])
    lo, hi = min(lo, truth_j - 1e-9), max(hi, truth_j + 1e-9)   # never exclude the truth
    return np.linspace(lo, hi, GRID_N)


def _logpost_on(g, j, th, data, m, sd):
    v = th.copy(); out = np.empty(len(g))
    for a, x in enumerate(g):
        v[j] = x
        out[a] = loglik(v, data) - 0.5*((x-m)/sd)**2
    return out


def one_replicate(args):
    r, j = args
    rng = np.random.default_rng(500000 + 1000*j + r)
    th = draw_prior(rng)
    data = simulate_dataset(th, rng)
    m, sd = PRIOR[PARAMS[j]]
    # The posterior is 3-15% of the prior width, so a grid spanning the PRIOR
    # resolves it with ~2 points and its trapezoid CDF is meaningless. Locate the
    # mode on a coarse scan, get a curvature-based width, then grid the posterior.
    g0 = grid_for(j, th[j])
    lp0 = _logpost_on(g0, j, th, data, m, sd)
    a0 = int(np.argmax(lp0)); mode = g0[a0]
    hw = (g0[1]-g0[0])
    gq = np.array([mode-hw, mode, mode+hw])
    lq = _logpost_on(gq, j, th, data, m, sd)
    curv = (lq[0] - 2*lq[1] + lq[2])/hw**2
    sd_hat = float(np.sqrt(-1.0/curv)) if curv < 0 else sd
    sd_hat = float(np.clip(sd_hat, 1e-6*sd, sd))
    lo, hi = mode - 12*sd_hat, mode + 12*sd_hat
    if j == 2: lo, hi = max(lo, N_BOUNDS[0]), min(hi, N_BOUNDS[1])
    g = np.linspace(lo, hi, GRID_N)
    lp = _logpost_on(g, j, th, data, m, sd)
    lp -= lp.max(); w = np.exp(lp)
    cdf = np.concatenate([[0.0], np.cumsum(0.5*(w[1:]+w[:-1])*np.diff(g))]); cdf /= cdf[-1]
    rank = (0.0 if th[j] < g[0] else 1.0 if th[j] > g[-1] else float(np.interp(th[j], g, cdf)))
    q = lambda a: float(np.interp(a, cdf, g))
    return {"rep": r, "param": PARAMS[j], "rank": rank,
            "inside94": bool(q(0.03) <= th[j] <= q(0.97)),
            "inside50": bool(q(0.25) <= th[j] <= q(0.75)),
            "post_sd_over_prior_sd": float((q(0.84)-q(0.16))/2.0/sd),
            "grid_covers_truth": bool(g[0] <= th[j] <= g[-1]),
            "tail_mass_at_edges": float(min(cdf[1], 1-cdf[-2]))}


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    jobs = [(r, j) for j in range(4) for r in range(R)]
    t0 = time.time()
    with Pool(8) as pool:
        out = pool.map(one_replicate, jobs)
    summ = {"R": R, "wall_s": round(time.time()-t0, 1), "grid_n": GRID_N, "params": {}}
    for j, k in enumerate(PARAMS):
        rr = [o for o in out if o["param"] == k]
        ranks = np.array([o["rank"] for o in rr])
        x = np.sort(ranks); n = len(x)
        D = float(np.max(np.abs(np.arange(1, n+1)/n - x)))
        c94 = float(np.mean([o["inside94"] for o in rr]))
        summ["params"][k] = dict(coverage94=c94, se94=float(np.sqrt(c94*(1-c94)/n)),
                                 coverage50=float(np.mean([o["inside50"] for o in rr])),
                                 mean_rank=float(ranks.mean()), ks_D=D,
                                 ks_p_bound=float(min(1.0, 2*np.exp(-2*n*D*D))),
                                 shrink=float(np.median([o["post_sd_over_prior_sd"] for o in rr])))
    json.dump({"summary": summ, "reps": out}, open("handoff/sbc_1d.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))
