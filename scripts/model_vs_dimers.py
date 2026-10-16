#!/usr/bin/env python3
"""CA RMSD of a predicted CP dimer against the dimers of this project.

    python scripts/model_vs_dimers.py <model.pdb> [chainX chainY]

The model's two CP chains (default A and C; the 129-residue ones) are compared,
in both orders, with chains A+B of S1 (C/C dimer), S6 (A/B dimer) and S4 (A/B,
RNA-bound), after a least-squares superposition. The lowest value says which
dimer the model resembles. Above about 2 A the S9 graft will warn.
"""
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent


def ca_by_chain(path):
    out = {}
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.setdefault(line[21], {})[int(line[22:26])] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return out


def rmsd(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(((P @ R - Q) ** 2).sum(1).mean()))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    model = ca_by_chain(sys.argv[1])
    x, y = (sys.argv[2], sys.argv[3]) if len(sys.argv) >= 4 else ("A", "C")
    for ch in (x, y):
        print(f"model chain {ch}: {len(model[ch])} CA")
    res = range(1, 130)
    for name in ("S1_wt_cc_apo", "S6_wt_ab_apo", "S4_cp_pacdesign"):
        p = REPO / "input" / name / f"{name}.pdb"
        if not p.exists():
            print(f"{name}: no {p.relative_to(REPO)} here (run build first)"); continue
        ref = ca_by_chain(p)
        R = np.array([ref["A"][i] for i in res] + [ref["B"][i] for i in res])
        for a, b in ((x, y), (y, x)):
            M = np.array([model[a][i] for i in res] + [model[b][i] for i in res])
            print(f"{name:18s} model {a}+{b} vs ref A+B: {rmsd(M, R):5.2f} A")


main()
