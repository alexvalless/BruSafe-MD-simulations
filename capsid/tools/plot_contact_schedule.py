"""Figure for the structure-derived contact schedule of the MS2 T = 3 capsid.

Reads capsid/data/structure/contact_schedule.csv (written by enumerate_contacts.py from PDB 2MS2)
and compares it with the placeholder used by notebook 04 (`contact_schedule`: contacts rising
linearly from 1 to 3 and rescaled to the total, `degeneracy`: constant 2).

    python3 capsid/tools/plot_contact_schedule.py

Output: capsid/figures/fig29_M3_contact_schedule.png (150 dpi, team palette).
"""
import csv
import json
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", message="findfont")
HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "structure"
OUT = HERE.parent / "figures" / "fig29_M3_contact_schedule.png"

# team palette, as in notebook 04
CREAM, PALE, STEEL = "#F4F3EF", "#CCD8E2", "#89A1C2"
SKY, BLUE, MID = "#4CA3DF", "#3E77E1", "#2C5BA5"
NAVY, INK = "#173779", "#122047"
plt.rcParams.update({
    "font.family": ["SF Pro Text", "SF Pro Display", "Helvetica Neue", "DejaVu Sans"],
    "figure.facecolor": CREAM, "axes.facecolor": CREAM, "savefig.facecolor": CREAM,
    "text.color": INK, "axes.labelcolor": INK, "axes.titlecolor": INK, "axes.edgecolor": INK,
    "xtick.color": INK, "ytick.color": INK, "grid.color": STEEL, "grid.alpha": 0.35,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
    "legend.facecolor": CREAM, "legend.edgecolor": STEEL, "legend.framealpha": 0.9,
})

rows = list(csv.DictReader(open(DATA / "contact_schedule.csv")))
col = lambda name: np.array([float(r[name]) for r in rows])  # noqa: E731
k = col("k_dimers_after_addition")
meta = json.load(open(DATA / "contact_schedule.json"))

# placeholder of notebook 04
N, total = 90, 180.0
c_ph = 1.0 + (3.0 - 1.0) * (k - 2) / (N - 2)
c_ph *= total / c_ph.sum()
g_ph = np.full_like(k, 2.0)

fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.7))
m, s = col("c_random_mean"), col("c_random_sd")
a.fill_between(k, np.maximum(m - s, 0), m + s, color=PALE, alpha=0.8, lw=0, label="± SD (random growth)")
a.plot(k, m, color=NAVY, lw=2, label="random growth, mean of 4000 paths")
a.plot(k, col("c_max"), color=BLUE, lw=0.8, ls="--", label="greedy, most contacts first")
a.plot(k, col("c_min"), color=SKY, lw=0.8, ls="--", label="greedy, fewest contacts first")
a.plot(k, c_ph, color=INK, lw=1.2, ls=":", label="placeholder of notebook 04")
a.set_xlabel("dimers in the cluster after the addition, k")
a.set_ylabel("new dimer–dimer interfaces formed")
a.set_title("(a) Contacts formed by the k-th dimer", loc="left", fontsize=10)
a.set_ylim(0, 5)
a.grid(True)
a.legend(fontsize=7, loc="upper left")

g = col("g_boundary_random_mean")
b.plot(k, g, color=NAVY, lw=2, label="random growth")
b.plot(k, col("g_boundary_max"), color=BLUE, lw=0.8, ls="--", label="greedy, most contacts first")
b.plot(k, col("g_boundary_min"), color=SKY, lw=0.8, ls="--", label="greedy, fewest contacts first")
b.plot(k, g_ph, color=INK, lw=1.2, ls=":", label="placeholder (constant 2)")
b.set_xlabel("dimers in the cluster after the addition, k")
b.set_ylabel("dimers that can be added at that step")
b.set_title("(b) Boundary sites", loc="left", fontsize=10)
b.grid(True)
b.legend(fontsize=7, loc="lower center")

fig.suptitle(f"MS2 T = 3 capsid, PDB {meta['structure']}: {meta['total_interfaces']} interfaces between 90 dimers "
             f"({meta['cutoff_A']} Å, ≥{meta['edge_min_residue_pairs']} residue pairs)", fontsize=10, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(OUT, dpi=150)
print("wrote", OUT, "| sum of c (random mean):", round(float(m.sum()), 2))
