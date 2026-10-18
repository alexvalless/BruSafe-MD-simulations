"""Is M2's recovery miss a sampler failure or a property of this dataset?

Builds nb03's namespace (definitions + the same seeded synthetic dataset),
then compares three posteriors over the same likelihood:
  REF   long adaptive random-walk MH, 4 independent chains  (the reference)
  SMC   the notebook's sampler as written, 400 particles
  SMC2k the notebook's sampler at 2000 particles (what the full run used)
Writes handoff/diag_m2.json.
"""
import json, time, sys, io, contextlib
import numpy as np
import matplotlib
matplotlib.use("Agg")

NB = "BruSafe_03_M2_soluble_partition.ipynb"
ns = {"__name__": "__main__", "display": lambda *a, **k: None}
nb = json.load(open(NB))
code = [i for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code"]
with contextlib.redirect_stdout(io.StringIO()):
    for i in code:
        if i > 25: break                      # stop after smc() is defined and run
        if i == 25:                           # define smc() but skip its own run
            src = "".join(nb["cells"][i]["source"])
            src = src.split("import time; _t0")[0]
            exec(src, ns); continue
        exec("".join(nb["cells"][i]["source"]), ns)

log_post, PARAMS, PRIOR = ns["log_post"], ns["PARAMS"], ns["PRIOR"]
log_lik, log_prior = ns["log_lik"], ns["log_prior"]
TRUTH = np.array([np.log(ns["S50_TRUE"]), ns["DE_TRUE"], ns["N_TRUE"], np.log(ns["NU_TRUE"])])

t0 = time.time(); [log_post(TRUTH) for _ in range(20)]
ms = (time.time() - t0) / 20 * 1000
print(f"log_post: {ms:.1f} ms", flush=True)


def mh_chain(seed, n_iter, n_tune, start):
    """Adaptive random-walk MH. Scale tuned to ~0.234 acceptance during tuning."""
    rng = np.random.default_rng(seed)
    d = len(start)
    x = np.array(start, float); lp = log_post(x)
    cov = np.diag((np.array([0.3, 2000., 0.3, 0.3]))**2)
    L = np.linalg.cholesky(cov); s = 2.38/np.sqrt(d)
    keep, acc_w, hist = [], 0, []
    for it in range(n_tune + n_iter):
        prop = x + s*(L @ rng.standard_normal(d))
        lpp = log_post(prop)
        if np.log(rng.random()) < lpp - lp:
            x, lp = prop, lpp; acc_w += 1
        if it < n_tune:
            hist.append(x.copy())
            if (it+1) % 200 == 0:
                a = acc_w/200.0
                s *= np.exp(a - 0.234)
                if (it+1) >= 1000 and (it+1) % 1000 == 0:
                    c = np.cov(np.array(hist[-1000:]).T) + 1e-12*np.eye(d)
                    try: L = np.linalg.cholesky(c)
                    except np.linalg.LinAlgError: pass
                acc_w = 0
        else:
            keep.append(x.copy())
            if (it+1) % 1000 == 0: acc_w = 0
    return np.array(keep)


def summarise(P, tag):
    out = {"tag": tag, "n": int(len(P))}
    for j, k in enumerate(PARAMS):
        v = P[:, j]
        out[k] = dict(mean=float(v.mean()), sd=float(v.std(ddof=1)),
                      lo94=float(np.percentile(v, 3)), hi94=float(np.percentile(v, 97)),
                      truth=float(TRUTH[j]),
                      inside=bool(np.percentile(v, 3) <= TRUTH[j] <= np.percentile(v, 97)))
    return out


if __name__ == "__main__":
    from multiprocessing import Pool
    res = {"log_post_ms": ms}

    N_ITER, N_TUNE = int(sys.argv[1]) if len(sys.argv) > 1 else 15000, 4000
    starts = [TRUTH + np.array([0.2, 1500., 0.2, 0.2])*np.random.default_rng(s).standard_normal(4)
              for s in (11, 12, 13, 14)]
    t0 = time.time()
    with Pool(4) as pool:
        chains = pool.starmap(mh_chain, [(s, N_ITER, N_TUNE, st)
                                         for s, st in zip((101, 102, 103, 104), starts)])
    res["ref_wall_s"] = round(time.time()-t0, 1)
    # between/within chain check on the reference itself
    M = np.stack(chains)                                   # (chain, iter, param)
    W = M.var(axis=1, ddof=1).mean(axis=0); B = M.mean(axis=1).var(axis=0, ddof=1)*M.shape[1]
    n = M.shape[1]
    res["ref_rhat"] = dict(zip(PARAMS, np.sqrt(((n-1)/n*W + B/n)/W).round(4).tolist()))
    REF = np.concatenate(chains)
    res["REF"] = summarise(REF, "reference MH")
    np.save("handoff/m2_reference_posterior.npy", REF)

    for npart in (400, 2000):
        t0 = time.time()
        with contextlib.redirect_stdout(io.StringIO()):
            part, logZ = ns["smc"](n_particles=npart, rng=np.random.default_rng(2026))
        r = summarise(part, f"SMC {npart}")
        r["wall_s"] = round(time.time()-t0, 1); r["logZ"] = float(logZ)
        res[f"SMC{npart}"] = r
        print(f"SMC {npart}: {r['wall_s']} s", flush=True)

    json.dump(res, open("handoff/diag_m2.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k.startswith(("REF", "SMC", "ref_", "log_"))},
                     indent=1)[:4000])
