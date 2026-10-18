"""Simulation-based calibration for M2.

One synthetic dataset gives one Bernoulli draw per parameter, so "2/4 recovered"
cannot distinguish a biased posterior from an unlucky dataset. SBC repeats the
experiment: draw theta from the prior, simulate data, refit, and record where the
truth falls in the posterior. If the model+likelihood+sampler triple is correct the
ranks are uniform and the 94% interval covers 94% of the time, whatever the dataset.

Sampler: the adaptive random-walk MH validated against the notebook's SMC in
diag_m2.py (posterior means and SDs agree to 3 decimals).
"""
import json, io, contextlib, time, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
from multiprocessing import Pool

NB = "BruSafe_03_M2_soluble_partition.ipynb"
ns = {"__name__": "__main__", "display": lambda *a, **k: None}
_nb = json.load(open(NB))
with contextlib.redirect_stdout(io.StringIO()):
    for _i, _c in enumerate(_nb["cells"]):
        if _c["cell_type"] != "code" or _i > 25: continue
        _src = "".join(_c["source"])
        if _i == 25: _src = _src.split("import time; _t0")[0]
        exec(_src, ns)

PARAMS, PRIOR, N_BOUNDS = ns["PARAMS"], ns["PRIOR"], ns["N_BOUNDS"]
predict_grid, beta_logpdf = ns["predict_grid"], ns["beta_logpdf"]
log_prior, DATA0 = ns["log_prior"], ns["DATA"]
beta_loglik, censor = ns["beta_loglik"], ns["censor"]
TEMPS_C, DOSES, T_OBS_MIN, TGRID = ns["TEMPS_C"], ns["DOSES"], ns["T_OBS_MIN"], ns["_TGRID"]
TI = [int(np.argmin(np.abs(TGRID - t))) for t in T_OBS_MIN]


def draw_prior(rng):
    v = np.array([rng.normal(*PRIOR[k]) for k in PARAMS])
    while not (N_BOUNDS[0] <= v[2] <= N_BOUNDS[1]):
        v[2] = rng.normal(*PRIOR["n"])
    return v


def simulate_dataset(v, rng):
    P = predict_grid(v)                                   # (T, dose, t)
    nu = np.exp(v[3])
    idx_T, idx_d, idx_t, y = [], [], [], []
    for k in range(len(TEMPS_C)):
        for i in range(len(DOSES)):
            for j, tt in enumerate(TI):
                p = float(np.clip(P[k, i, tt], 1e-6, 1-1e-6))
                for _ in range(2):                        # two replicates, as in the notebook
                    idx_T.append(k); idx_d.append(i); idx_t.append(tt)
                    y.append(float(censor(rng.beta(p*nu, (1-p)*nu))))
    return dict(T_idx=np.array(idx_T), dose_idx=np.array(idx_d),
                t_idx=np.array(idx_t), y=np.array(y))


def log_post_data(v, data):
    lp = log_prior(v, True)
    if not np.isfinite(lp): return -np.inf
    P = predict_grid(v)
    p = np.clip(P[data["T_idx"], data["dose_idx"], data["t_idx"]], 1e-6, 1-1e-6)
    return lp + float(np.sum(beta_loglik(data["y"], p, np.exp(v[3]))))


BASE_SCALE = np.array([0.3, 2000., 0.3, 0.3])


def mh(data, start, seed, n_keep=8000, n_tune=5000, thin=4):
    """Adaptive random-walk MH with two guards that the first version lacked:
    the adapted covariance keeps a ridge proportional to the prior scale, and the
    step size cannot collapse. Without them a chain that stalls early adapts to a
    near-zero covariance and stays stalled for the rest of the run, which shows up
    downstream as a posterior of zero width rather than as an obvious failure."""
    rng = np.random.default_rng(seed); d = len(start)
    x = np.array(start, float); lp = log_post_data(x, data)
    ridge = np.diag((0.05*BASE_SCALE)**2)
    L = np.linalg.cholesky(np.diag(BASE_SCALE**2)); s = 2.38/np.sqrt(d)
    keep, acc, acc_total, hist = [], 0, 0, []
    for it in range(n_tune + n_keep):
        prop = x + s*(L @ rng.standard_normal(d))
        lpp = log_post_data(prop, data)
        if np.log(rng.random()) < lpp - lp:
            x, lp, acc, acc_total = prop, lpp, acc + 1, acc_total + 1
        if it < n_tune:
            hist.append(x.copy())
            if (it+1) % 200 == 0:
                s = float(np.clip(s*np.exp(acc/200.0 - 0.234), 1e-3, 20.0)); acc = 0
                if (it+1) >= 1000 and (it+1) % 500 == 0:
                    c = np.cov(np.array(hist[-1000:]).T) + ridge
                    try: L = np.linalg.cholesky(c)
                    except np.linalg.LinAlgError: pass
            if it+1 == n_tune: acc_total = 0
        elif it % thin == 0:
            keep.append(x.copy())
    return np.array(keep), acc_total/float(n_keep)


def ess_1d(v):
    """Effective sample size from the initial-positive-sequence autocorrelations."""
    v = v - v.mean(); n = len(v)
    var = float(np.dot(v, v)/n)
    if var <= 0: return 0.0
    tot, k = 0.0, 1
    while k < min(n//2, 500):
        r = float(np.dot(v[:-k], v[k:])/(n*var))
        if r < 0.05: break
        tot += r; k += 1
    return float(n/(1.0 + 2.0*tot))


def one_replicate(r):
    rng = np.random.default_rng(90000 + r)
    th = draw_prior(rng)
    data = simulate_dataset(th, rng)
    P, acc = mh(data, th + 0.05*BASE_SCALE*rng.standard_normal(4), 70000 + r)
    ess = [ess_1d(P[:, j]) for j in range(P.shape[1])]
    out = {"rep": r, "theta": th.tolist(), "n_draws": int(len(P)),
           "accept": float(acc), "ess_min": float(min(ess)),
           "ok": bool(acc > 0.02 and min(ess) > 100)}
    for j, k in enumerate(PARAMS):
        v = P[:, j]
        out[k] = dict(rank=float((v < th[j]).mean()),
                      inside94=bool(np.percentile(v, 3) <= th[j] <= np.percentile(v, 97)),
                      inside50=bool(np.percentile(v, 25) <= th[j] <= np.percentile(v, 75)),
                      z=float((th[j]-v.mean())/v.std(ddof=1)))
    return out


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    t0 = time.time()
    with Pool(8) as pool:
        reps = pool.map(one_replicate, range(R))
    wall = round(time.time()-t0, 1)
    good = [r for r in reps if r["ok"]]
    summ = {"R": R, "R_used": len(good), "wall_s": wall,
            "accept_min": float(min(r["accept"] for r in reps)),
            "ess_min": float(min(r["ess_min"] for r in reps)), "params": {}}
    reps_all, reps = reps, good
    R = len(reps)
    for k in PARAMS:
        ranks = np.array([r[k]["rank"] for r in reps])
        c94 = float(np.mean([r[k]["inside94"] for r in reps]))
        c50 = float(np.mean([r[k]["inside50"] for r in reps]))
        # KS test of rank uniformity, computed without scipy's dependency on ties
        x = np.sort(ranks); n = len(x)
        D = float(np.max(np.abs(np.arange(1, n+1)/n - x)).item())
        ks_p = float(2*np.exp(-2*n*D*D))                  # asymptotic two-sided bound
        summ["params"][k] = dict(coverage94=c94, coverage50=c50,
                                 se94=float(np.sqrt(c94*(1-c94)/R)),
                                 mean_rank=float(ranks.mean()), ks_D=D, ks_p_bound=min(1.0, ks_p),
                                 mean_z=float(np.mean([r[k]["z"] for r in reps])))
    json.dump({"summary": summ, "reps": reps_all}, open("handoff/sbc_m2.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))
