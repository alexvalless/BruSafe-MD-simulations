#!/usr/bin/env python3
"""
Convergence and sanity analysis for the BruSafe MD campaign. numpy only.

Subcommands
-----------
  gate       pass/fail check on an equilibration (density, temperature)
  block      block-averaged standard error vs block size for one .xvg series
  cosine     cosine content of a series (Hess 2000/2002)
  replicas   agreement between per-replica RMSF profiles

Why these three, specifically
-----------------------------
* Block averaging: the naive standard error over MD frames is meaningless
  because frames are correlated. The SEM only becomes honest once the block
  size exceeds the correlation time, which is where the curve plateaus.

* Cosine content: a PC1 projection drawn from pure diffusion has cosine
  content approaching 1. If your headline "conformational transition" scores
  0.9, it is random drift, not a transition. Better to find that yourself
  than to have a judge find it.

* Replica agreement: the error bar you report should be the spread across
  independent replicas, not the spread across correlated frames within one.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np


# --------------------------------------------------------------------------

def read_xvg(path: str) -> np.ndarray:
    """Read a GROMACS .xvg into an (N, ncol) float array, skipping headers."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#@&":
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    if not rows:
        sys.exit(f"no numeric data in {path}")
    width = min(len(r) for r in rows)
    return np.array([r[:width] for r in rows], dtype=float)


# --------------------------------------------------------------------------

def block_average(y: np.ndarray, min_blocks: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """SEM of the mean as a function of block size.

    Returns (block_sizes, sem). Read off the plateau; that is the real
    uncertainty. If it never plateaus, the run is too short to quote an
    uncertainty at all.
    """
    n = len(y)
    sizes, sems = [], []
    size = 1
    while n // size >= min_blocks:
        nb = n // size
        means = y[: nb * size].reshape(nb, size).mean(axis=1)
        sems.append(means.std(ddof=1) / np.sqrt(nb))
        sizes.append(size)
        size = max(size + 1, int(size * 1.4))
    return np.array(sizes), np.array(sems)


def cosine_content(y: np.ndarray, index: int = 1) -> float:
    """Cosine content of a time series (Hess, PRE 2000; PRE 2002).

    c_i = (2/T) * (integral cos(i*pi*t/T) p(t) dt)^2 / integral p(t)^2 dt

    ~0  -> the motion carries real signal
    ~1  -> indistinguishable from free diffusion; do not interpret it
    """
    n = len(y)
    if n < 4:
        sys.exit("series too short for cosine content")
    y = y - y.mean()
    t = np.arange(n)
    cos = np.cos(index * np.pi * t / (n - 1))
    num = (np.trapezoid(cos * y, t)) ** 2 if hasattr(np, "trapezoid") \
        else (np.trapz(cos * y, t)) ** 2
    den = np.trapezoid(y * y, t) if hasattr(np, "trapezoid") else np.trapz(y * y, t)
    if den == 0:
        return 0.0
    return float(2.0 / (n - 1) * num / den)


# --------------------------------------------------------------------------

def cmd_gate(a: argparse.Namespace) -> None:
    """Hard pass/fail before production. Called by 02_equilibrate.sh."""
    ok = True
    tail = lambda arr: arr[len(arr) // 2:]        # noqa: E731 — second half only

    d = tail(read_xvg(a.density)[:, 1])
    dm, ds = d.mean(), d.std()
    print(f"density      {dm:8.1f} +/- {ds:5.1f} kg/m^3", end="  ")
    if 960.0 <= dm <= 1060.0:
        print("PASS")
    else:
        print("FAIL  (expected ~1000 for TIP3P at 310 K / 1 bar)")
        ok = False

    t = tail(read_xvg(a.temperature)[:, 1])
    tm, ts = t.mean(), t.std()
    print(f"temperature  {tm:8.2f} +/- {ts:5.2f} K", end="  ")
    if abs(tm - a.ref_t) <= 2.0:
        print("PASS")
    else:
        print(f"FAIL  (expected {a.ref_t} K)")
        ok = False

    # A density that is still drifting means the box has not settled, whatever
    # the mean says.
    half = len(d) // 2
    drift = abs(d[half:].mean() - d[:half].mean())
    print(f"density drift{drift:8.2f} kg/m^3 across second half", end="  ")
    if drift < 2.0:
        print("PASS")
    else:
        print("FAIL  (still equilibrating -- extend npt_free)")
        ok = False

    sys.exit(0 if ok else 1)


def cmd_block(a: argparse.Namespace) -> None:
    y = read_xvg(a.xvg)[:, a.col]
    if a.skip:
        y = y[a.skip:]
    sizes, sems = block_average(y)
    print(f"n = {len(y)}   mean = {y.mean():.4f}")
    print(f"{'block':>8}{'SEM':>12}")
    for s, e in zip(sizes, sems):
        print(f"{s:>8}{e:>12.5f}")
    plateau = sems[len(sems) * 2 // 3:].mean()
    naive = y.std(ddof=1) / np.sqrt(len(y))
    print(f"\nnaive SEM (WRONG, assumes independent frames): {naive:.5f}")
    print(f"block-averaged SEM (use this):                 {plateau:.5f}")
    print(f"underestimate factor:                          {plateau / naive:.1f}x")


def cmd_cosine(a: argparse.Namespace) -> None:
    """Cosine content of a PC projection.

    Calibration (verified against analytic limits):
      exact half-cosine  -> 1.00   the pathological diffusive shape
      white noise        -> 0.00

    Caveat that matters: a single 1D random walk gives a highly variable
    cosine content (ensemble mean ~0.43, 10th-90th percentile 0.02-0.85), so
    ONE number on ONE replica is not proof of anything. PCA additionally
    biases PC1 upward, because selecting the largest-variance mode favours
    cosine-like shapes. Treat a high value as a prompt to check the other
    replicas, not as a verdict on its own.
    """
    y = read_xvg(a.xvg)[:, a.col]
    if a.skip:
        y = y[a.skip:]
    for i in (1, 2, 3):
        c = cosine_content(y, i)
        verdict = ("likely diffusive -- check other replicas before interpreting"
                   if c > 0.5 else
                   "inconclusive on its own" if c > 0.2 else
                   "not diffusive; safe to interpret")
        print(f"cosine content PC{i}: {c:.3f}   {verdict}")
    print("\nRun this on every replica. Consistently high PC1 across replicas "
          "means your 'transition' is drift.")


def cmd_replicas(a: argparse.Namespace) -> None:
    """Pairwise agreement between per-replica RMSF profiles."""
    profiles = [read_xvg(p) for p in a.xvg]
    n = min(len(p) for p in profiles)
    mat = np.array([p[:n, 1] for p in profiles])
    print(f"{len(mat)} replicas, {n} residues\n")
    print("pairwise Pearson r:")
    for i in range(len(mat)):
        for j in range(i + 1, len(mat)):
            r = np.corrcoef(mat[i], mat[j])[0, 1]
            flag = "" if r > 0.8 else "   <-- replicas disagree"
            print(f"  rep{i+1} vs rep{j+1}:  r = {r:.3f}{flag}")
    print("\nper-residue mean +/- spread across replicas is your error bar")
    print(f"mean RMSF        : {mat.mean():.4f} nm")
    print(f"mean inter-replica spread: {mat.std(axis=0).mean():.4f} nm")
    worst = int(np.argmax(mat.std(axis=0)))
    print(f"least reproducible residue index: {worst} "
          f"(spread {mat.std(axis=0)[worst]:.4f} nm)")


# --------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("gate")
    p.add_argument("--density", required=True)
    p.add_argument("--temperature", required=True)
    p.add_argument("--ref-t", type=float, default=310.0)
    p.set_defaults(fn=cmd_gate)

    p = sub.add_parser("block")
    p.add_argument("xvg"); p.add_argument("--col", type=int, default=1)
    p.add_argument("--skip", type=int, default=0)
    p.set_defaults(fn=cmd_block)

    p = sub.add_parser("cosine")
    p.add_argument("xvg"); p.add_argument("--col", type=int, default=1)
    p.add_argument("--skip", type=int, default=0)
    p.set_defaults(fn=cmd_cosine)

    p = sub.add_parser("replicas")
    p.add_argument("xvg", nargs="+", help="one rmsf .xvg per replica")
    p.set_defaults(fn=cmd_replicas)

    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
