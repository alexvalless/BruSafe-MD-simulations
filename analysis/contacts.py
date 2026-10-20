#!/usr/bin/env python3
"""
Protein-RNA contact enumeration for one system, all replicas pooled.

Reads what 04_postprocess.sh already wrote (analysis/clean.gro + clean.xtc,
solute only, PBC-clean) and reports, over the production window:

  * contacts      residue pairs (protein residue x nucleotide) with any heavy-atom
                  distance < --cutoff, as an occupancy = fraction of frames
  * hydrogen bonds  protein<->RNA, occupancy per residue pair and atom pair
  * interface     residue pairs between protein chains (dimer interface)
  * time series   number of contacts and minimum protein-RNA distance per frame
  * key residues  optional: do the known operator-contact residues show up
                  among the top contacts? (the hot-spot gate in TUTORIAL 10)

CPU only -- no GPU needed. Needs MDAnalysis and matplotlib.

  python3 analysis/contacts.py S3_cp_operator
  python3 analysis/contacts.py S4_cp_pacdesign --key-residues 45,47,59,61

Output: analysis_out/contacts/<system>/
  contacts_pairs.csv  contacts_by_residue.csv  contacts_by_nucleotide.csv
  hbonds_pairs.csv  interface_pairs.csv  summary_contacts.csv/.md
  figures/*.png,svg   index.html

Read before quoting numbers
---------------------------
* Occupancy is a property of a REPLICA. The +/- is the SD across replicas, and
  "n_reps>=" counts in how many replicas the contact is persistent. A contact
  seen in one replica only is not a result.
* clean.xtc holds the solute only. Water-mediated bridges and Mg2+ contacts
  need prod.xtc, so they are NOT computed here.
* Contacts are geometry (a cutoff), not energy. Do not call them binding energy.
* Chains are detected where residue numbering restarts. The detected chains are
  printed: check them, and use --chain-sizes if they are wrong.
"""

from __future__ import annotations

import argparse
import csv
import html
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
from md_report import (GRID, INK, INK_2, PAGE_CSS, SERIES, find_replicas,  # noqa: E402
                       need_modules, registry_row)

RNA_RESNAMES = ("A G C U RA RG RC RU ADE GUA CYT URA RA5 RA3 RG5 RG3 RC5 RC3 RU5 RU3 "
                "A5 A3 G5 G3 C5 C3 U5 U3").split()
NT_LETTER = {"ADE": "A", "GUA": "G", "CYT": "C", "URA": "U"}


def log(msg: str) -> None:
    print(f"[contacts] {msg}", file=sys.stderr)


# ==========================================================================
# residues and chains

def chain_labels(resids: np.ndarray, sizes: list[int] | None, alphabet: str) -> list[str]:
    """One label per residue. A new chain starts where numbering does not increase."""
    n = len(resids)
    lab = np.zeros(n, dtype=int)
    if sizes:
        if sum(sizes) != n:
            sys.exit(f"--chain-sizes adds up to {sum(sizes)} but there are {n} protein residues")
        lab = np.repeat(np.arange(len(sizes)), sizes)
    else:
        c = 0
        for i in range(1, n):
            if resids[i] <= resids[i - 1]:
                c += 1
            lab[i] = c
    return [alphabet[i % len(alphabet)] for i in lab]


def residue_table(ag, kind: str, sizes: list[int] | None = None, rna_first: int | None = None):
    """Per-residue labels for an atom group + for every atom its local residue index."""
    uniq, inv = np.unique(ag.resindices, return_inverse=True)
    res = ag.universe.residues[uniq]
    resids, names = res.resids, res.resnames
    if kind == "prot":
        chains = chain_labels(resids, sizes, "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        labels = [f"{c}:{n}{r}" for c, n, r in zip(chains, names, resids)]
    else:
        chains = chain_labels(resids, None, "RSTUVWXYZ")
        nums = list(resids)
        if rna_first is not None:        # label with the paper's numbering (no position 0)
            nums = [int(r) - int(resids[0]) + rna_first for r in resids]
            nums = [p + 1 if (rna_first < 0 and p >= 0) else p for p in nums]
        labels = [f"{c}:{NT_LETTER.get(n, n)}{r}" for c, n, r in zip(chains, names, nums)]
    return dict(labels=labels, chains=np.array(chains), resids=np.array(resids),
                names=list(names), atom_to_res=inv, resindex=uniq)


# ==========================================================================
# one replica

def analyse_replica(rep, a) -> dict | None:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance, distance_array

    top, trj = rep.path / "clean.gro", rep.path / "clean.xtc"
    if not (top.is_file() and trj.is_file()):
        log(f"{rep.name}: no clean.gro/clean.xtc -- run 04_postprocess.sh first, skipped")
        return None
    u = mda.Universe(str(top), str(trj))
    rna_sel = a.rna_sel or "resname " + " ".join(RNA_RESNAMES)
    prot = u.select_atoms("protein and not name H*")
    rna = u.select_atoms(f"({rna_sel}) and not name H*")
    has_rna = len(rna) > 0
    if len(prot) == 0:
        log(f"{rep.name}: no protein atoms in clean.gro, skipped")
        return None
    P = residue_table(prot, "prot", a.chain_sizes)
    R = residue_table(rna, "rna", rna_first=a.rna_first) if has_rna else None
    nP, nR = len(P["labels"]), (len(R["labels"]) if has_rna else 0)
    sizes = Counter(P["chains"].tolist())
    log(f"{rep.name}: {nP} protein residues in chains "
        + ", ".join(f"{c}={n}" for c, n in sorted(sizes.items()))
        + (f"; {nR} nucleotides" if has_rna else "; no RNA found"))

    t = np.array([ts.time for ts in u.trajectory])
    frames = np.flatnonzero(t >= a.equil_ns * 1000.0)
    fallback = len(frames) < 10
    if fallback:
        log(f"{rep.name}: trajectory ends at {t[-1] / 1000:.1f} ns, fewer than 10 frames after "
            f"{a.equil_ns:g} ns -- analysing ALL frames, relaxation included")
        frames = np.arange(len(t))
    if len(frames) > a.max_frames:
        frames = frames[np.linspace(0, len(frames) - 1, a.max_frames).astype(int)]

    multi = len(sizes) > 1
    occ_pr = np.zeros((nP, nR))
    occ_pp = np.zeros((nP, nP)) if multi else None
    n_pairs, n_prot, n_nt, min_d, n_iface = [], [], [], [], []
    for fi in frames:
        ts = u.trajectory[fi]
        box = ts.dimensions
        if has_rna:
            D = distance_array(prot.positions, rna.positions, box=box)
            min_d.append(D.min())
            ai, aj = np.nonzero(D < a.cutoff)
            M = np.zeros((nP, nR), dtype=bool)
            M[P["atom_to_res"][ai], R["atom_to_res"][aj]] = True
            occ_pr += M
            n_pairs.append(int(M.sum()))
            n_prot.append(int(M.any(1).sum()))
            n_nt.append(int(M.any(0).sum()))
        if multi:
            pairs = capped_distance(prot.positions, prot.positions, a.cutoff, box=box,
                                    return_distances=False)
            ci, cj = P["atom_to_res"][pairs[:, 0]], P["atom_to_res"][pairs[:, 1]]
            keep = P["chains"][ci] != P["chains"][cj]
            Mi = np.zeros((nP, nP), dtype=bool)
            Mi[np.minimum(ci[keep], cj[keep]), np.maximum(ci[keep], cj[keep])] = True
            occ_pp += Mi
            n_iface.append(int(Mi.sum()))
    nf = len(frames)
    out = dict(name=rep.name, color=rep.color, t=t[frames] / 1000.0, nf=nf, P=P, R=R,
               has_rna=has_rna, occ_pr=occ_pr / nf, occ_pp=(occ_pp / nf if multi else None),
               n_pairs=np.array(n_pairs), n_prot=np.array(n_prot), n_nt=np.array(n_nt),
               min_d=np.array(min_d), n_iface=np.array(n_iface), hb=None, hb_count=None,
               window=(float(t[frames[0]]) / 1000.0, float(t[frames[-1]]) / 1000.0), fallback=fallback)
    if has_rna and not a.no_hbonds:
        out["hb"], out["hb_count"] = hbonds(u, frames, P, R, rna_sel, a)
    return out


def hbonds(u, frames, P, R, rna_sel, a):
    """Protein<->RNA hydrogen bonds. Hydrogens must be in clean.gro (they are)."""
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

    both = f"(protein or ({rna_sel}))"
    hba = HydrogenBondAnalysis(
        u, donors_sel=f"{both} and name N* O*", hydrogens_sel=f"{both} and name H*",
        acceptors_sel=f"{both} and name N* O*", between=["protein", f"({rna_sel})"],
        d_h_cutoff=1.2, d_a_cutoff=a.hb_dist, d_h_a_angle_cutoff=a.hb_angle,
        update_selections=False)
    hba.run(frames=frames, verbose=False)
    hb = hba.results.hbonds
    n = len(frames)
    if len(hb) == 0:
        return {}, np.zeros(n, dtype=int)
    rix = u.atoms.resindices
    p_of = np.full(len(u.residues), -1)
    r_of = np.full(len(u.residues), -1)
    p_of[P["resindex"]] = np.arange(len(P["resindex"]))
    r_of[R["resindex"]] = np.arange(len(R["resindex"]))
    names = u.atoms.names
    resn = u.atoms.resnames
    frames_of = defaultdict(set)
    atoms_of = defaultdict(Counter)
    pos = {int(f): i for i, f in enumerate(frames)}
    per_frame = np.zeros(n, dtype=int)
    for f, d, _, acc, _, _ in hb:
        f, d, acc = int(f), int(d), int(acc)
        dr, ar = rix[d], rix[acc]
        if p_of[dr] >= 0 and r_of[ar] >= 0:
            key, desc = (p_of[dr], r_of[ar]), f"{P['names'][p_of[dr]]}:{names[d]} -> {resn[acc]}:{names[acc]}"
        elif r_of[dr] >= 0 and p_of[ar] >= 0:
            key, desc = (p_of[ar], r_of[dr]), f"{resn[d]}:{names[d]} -> {P['names'][p_of[ar]]}:{names[acc]}"
        else:
            continue
        frames_of[key].add(f)
        atoms_of[key][desc] += 1
        per_frame[pos[f]] += 1
    out = {k: (len(v) / n, atoms_of[k].most_common(1)[0][0]) for k, v in frames_of.items()}
    return out, per_frame


# ==========================================================================
# pooling and tables

def stack(reps, key):
    arr = np.array([r[key] for r in reps])
    sd = arr.std(axis=0, ddof=1) if len(arr) > 1 else np.zeros(arr.shape[1:])
    return arr, arr.mean(0), sd


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow([f"{x:.3f}" if isinstance(x, float) else x for x in r])


def pool(reps, a, outdir: Path) -> dict:
    ref = reps[0]
    for r in reps[1:]:
        if r["P"]["labels"] != ref["P"]["labels"] or (
                r["has_rna"] and r["R"]["labels"] != ref["R"]["labels"]):
            sys.exit(f"{r['name']} has different residues than {ref['name']}: the replicas "
                     "were not built from the same input, they cannot be pooled")
    names = [r["name"] for r in reps]
    nrep = len(reps)
    res: dict = {"names": names}
    pl = ref["P"]["labels"]

    if ref["has_rna"]:
        rl = ref["R"]["labels"]
        arr, mean, sd = stack(reps, "occ_pr")
        res.update(arr_pr=arr, mean_pr=mean, sd_pr=sd)
        rows = []
        for i, j in zip(*np.nonzero(mean >= a.list_min)):
            hbo = [r["hb"].get((i, j), (0.0, ""))[0] if r["hb"] is not None else 0.0 for r in reps]
            top = max((r["hb"].get((i, j), (0.0, "")) for r in reps if r["hb"] is not None),
                      default=(0.0, ""))[1]
            rows.append([pl[i], rl[j]] + [float(arr[k, i, j]) for k in range(nrep)]
                        + [float(mean[i, j]), float(sd[i, j]), int((arr[:, i, j] >= a.min_occ).sum()),
                           float(np.mean(hbo)), top])
        rows.sort(key=lambda x: -x[2 + nrep])
        write_csv(outdir / "contacts_pairs.csv",
                  ["protein_residue", "nucleotide"] + [f"occ_{n}" for n in names]
                  + ["occ_mean", "occ_sd", f"n_reps>={a.min_occ:g}", "hbond_occ_mean", "top_hbond_atoms"],
                  rows)
        res["pairs_rows"] = rows

        # per residue / nucleotide: fraction of frames with any contact
        # (computed per replica from the pair matrix would overcount: use the max over partners)
        po = np.array([r["occ_pr"].max(1) for r in reps])     # lower bound on "any contact"
        no = np.array([r["occ_pr"].max(0) for r in reps])
        res.update(p_any=po, n_any=no)
        for fname, labels, arr_ in (("contacts_by_residue.csv", pl, po),
                                    ("contacts_by_nucleotide.csv", rl, no)):
            m = arr_.mean(0)
            s = arr_.std(0, ddof=1) if nrep > 1 else np.zeros_like(m)
            rws = [[labels[i], float(m[i]), float(s[i]), int((arr_[:, i] >= a.min_occ).sum())]
                   for i in np.argsort(-m) if m[i] > 0]
            write_csv(outdir / fname, ["residue", "max_pair_occ_mean", "sd", f"n_reps>={a.min_occ:g}"],
                      rws)

        hb_rows = []
        keys = set()
        for r in reps:
            keys |= set((r["hb"] or {}).keys())
        for k in keys:
            occ = [(r["hb"] or {}).get(k, (0.0, ""))[0] for r in reps]
            top = max(((r["hb"] or {}).get(k, (0.0, "")) for r in reps))[1]
            hb_rows.append([pl[k[0]], rl[k[1]]] + occ + [float(np.mean(occ)),
                            float(np.std(occ, ddof=1)) if nrep > 1 else 0.0,
                            int(sum(o >= a.min_occ for o in occ)), top])
        hb_rows.sort(key=lambda x: -x[2 + nrep])
        write_csv(outdir / "hbonds_pairs.csv",
                  ["protein_residue", "nucleotide"] + [f"occ_{n}" for n in names]
                  + ["occ_mean", "occ_sd", f"n_reps>={a.min_occ:g}", "top_atom_pair"], hb_rows)
        res["hb_rows"] = hb_rows

    if ref["occ_pp"] is not None:
        arr, mean, sd = stack(reps, "occ_pp")
        rows = []
        for i, j in zip(*np.nonzero(mean >= a.list_min)):
            rows.append([pl[i], pl[j]] + [float(arr[k, i, j]) for k in range(nrep)]
                        + [float(mean[i, j]), float(sd[i, j]), int((arr[:, i, j] >= a.min_occ).sum())])
        rows.sort(key=lambda x: -x[2 + nrep])
        write_csv(outdir / "interface_pairs.csv",
                  ["residue_1", "residue_2"] + [f"occ_{n}" for n in names]
                  + ["occ_mean", "occ_sd", f"n_reps>={a.min_occ:g}"], rows)
        res["iface_rows"] = rows
    return res


def window_note(reps, a) -> str:
    """What was really analysed, per replica: never claim a window that was not used."""
    parts = [f"{r['name']} {r['window'][0]:.0f}-{r['window'][1]:.0f} ns"
             + (" (ALL frames: shorter than the relaxation window)" if r["fallback"] else "")
             for r in reps]
    head = (f"WARNING: some trajectories are shorter than {a.equil_ns:g} ns, so the relaxation was "
            "not excluded. " if any(r["fallback"] for r in reps) else
            f"Frames t >= {a.equil_ns:g} ns. ")
    return head + "Window analysed: " + "; ".join(parts) + "."


def summary(reps, a, outdir: Path) -> tuple[list[str], list[list]]:
    def m(x): return float(np.mean(x)) if len(x) else float("nan")
    header = ["replica", "frames", "contact pairs", "protein residues in contact",
              "nucleotides in contact", "bound fraction", "min distance (A)", "H-bonds per frame",
              "interface pairs"]
    rows = []
    for r in reps:
        rows.append([r["name"], r["nf"], m(r["n_pairs"]), m(r["n_prot"]), m(r["n_nt"]),
                     float((r["min_d"] < a.cutoff).mean()) if len(r["min_d"]) else float("nan"),
                     m(r["min_d"]), m(r["hb_count"]) if r["hb_count"] is not None else float("nan"),
                     m(r["n_iface"])])
    if len(rows) > 1:
        arr = np.array([[x for x in row[1:]] for row in rows], dtype=float)
        rows.append(["mean"] + [float(v) for v in np.nanmean(arr, 0)])
        rows.append(["SD across replicas"] + [float(v) for v in np.nanstd(arr, 0, ddof=1)])
    write_csv(outdir / "summary_contacts.csv", header, rows)
    fmt = lambda v: "" if (isinstance(v, float) and np.isnan(v)) else (f"{v:.4g}" if isinstance(v, float) else str(v))  # noqa: E731
    md = ["| " + " | ".join(header) + " |", "|---" * len(header) + "|"]
    md += ["| " + " | ".join(fmt(v) for v in r) + " |" for r in rows]
    note = (f"Cutoff {a.cutoff:g} A heavy atoms. {window_note(reps, a)} Report mean +/- SD "
            "across replicas. Geometry only, not energy; no water bridges or Mg2+ (solute-only trajectory).")
    (outdir / "summary_contacts.md").write_text("\n".join(md) + f"\n\n{note}\n", encoding="utf-8")
    return header, rows


def key_report(res, reps, a) -> list[str]:
    if not a.key_residues:
        return []
    ref = reps[0]
    pl = ref["P"]["labels"]
    resids = ref["P"]["resids"]
    score = res["p_any"].mean(0)
    top = np.argsort(-score)[:a.top_n]
    top_ids = {int(resids[i]) for i in top}
    found = [k for k in a.key_residues if k in top_ids]
    lines = [f"hot-spot recovery: {len(found)}/{len(a.key_residues)} key residues are among the "
             f"top {a.top_n} contacting residues ({', '.join(map(str, found)) or 'none'})"]
    for k in a.key_residues:
        idx = [i for i in range(len(pl)) if int(resids[i]) == k]
        for i in idx:
            lines.append(f"  {pl[i]:<12} occupancy {score[i]:.2f} +/- "
                         f"{res['p_any'][:, i].std(ddof=1) if len(reps) > 1 else 0:.2f}")
        if not idx:
            lines.append(f"  residue {k}: not in the system")
    return lines


# ==========================================================================
# figures

def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK_2, "axes.labelcolor": INK,
                         "xtick.color": INK_2, "ytick.color": INK_2, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                         "grid.linewidth": 0.6, "legend.frameon": False, "svg.fonttype": "none"})
    return plt


def save(plt, fig, figdir: Path, name: str) -> None:
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(figdir / f"{name}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def timeseries(plt, reps, key, ylabel, title, figdir, name):
    series = [(r, r[key]) for r in reps if len(r[key]) == len(r["t"])]
    if not series:
        return False
    fig, ax = plt.subplots(figsize=(7.2, 4))
    for r, y in series:
        ax.plot(r["t"], y, color=r["color"], lw=1.2, label=r["name"])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=11, pad=12)
    if len(series) > 1:
        ax.legend(ncol=4, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    save(plt, fig, figdir, name)
    return True


def occupancy_bars(plt, labels, arr, key_ids, resids, title, figdir, name, ylabel):
    """arr: (n_reps, n). Bars = mean over replicas, error bars = SD across replicas."""
    mean = arr.mean(0)
    sd = arr.std(0, ddof=1) if len(arr) > 1 else np.zeros_like(mean)
    idx = np.flatnonzero(mean > 0.02)
    if len(idx) == 0:
        return False
    fig, ax = plt.subplots(figsize=(max(7.2, 0.22 * len(idx)), 4))
    colors = [SERIES[1] if (resids is not None and int(resids[i]) in key_ids) else SERIES[0] for i in idx]
    ax.bar(range(len(idx)), mean[idx], yerr=sd[idx], color=colors, error_kw=dict(lw=0.8, ecolor=INK_2))
    ax.set_xticks(range(len(idx)), [labels[i] for i in idx], rotation=90, fontsize=7)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1.05)
    ax.set_title(title, loc="left", fontsize=11, pad=12)
    if key_ids:
        ax.text(0.99, 0.97, "orange = key residues", transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color=INK_2)
    save(plt, fig, figdir, name)
    return True


def contact_map(plt, pl, rl, mean, a, figdir):
    rows = np.flatnonzero(mean.max(1) >= a.map_min)
    if len(rows) == 0:
        rows = np.argsort(-mean.max(1))[:25]
        rows = np.sort(rows)
    cols = np.flatnonzero(mean.max(0) > 0)
    if len(cols) == 0:
        return False
    fig, ax = plt.subplots(figsize=(max(5, 0.4 * len(cols) + 2), max(4, 0.22 * len(rows) + 1.5)))
    im = ax.imshow(mean[np.ix_(rows, cols)], aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)), [rl[j] for j in cols], rotation=90, fontsize=7)
    ax.set_yticks(range(len(rows)), [pl[i] for i in rows], fontsize=7)
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="occupancy (mean of replicas)", pad=0.01)
    ax.set_title("Protein residue x nucleotide contact map", loc="left", fontsize=11, pad=12)
    save(plt, fig, figdir, "contact_map")
    return True


def hbond_bars(plt, rows, nrep, figdir):
    top = rows[:20]
    if not top:
        return False
    fig, ax = plt.subplots(figsize=(7.2, 0.28 * len(top) + 1.5))
    y = np.arange(len(top))[::-1]
    ax.barh(y, [r[2 + nrep] for r in top], xerr=[r[3 + nrep] for r in top], color=SERIES[0],
            error_kw=dict(lw=0.8, ecolor=INK_2))
    ax.set_yticks(y, [f"{r[0]}  -  {r[1]}" for r in top], fontsize=7)
    ax.set_xlabel("H-bond occupancy (mean +/- SD across replicas)")
    ax.set_xlim(0, 1.05)
    ax.set_title("Most persistent protein-RNA hydrogen bonds", loc="left", fontsize=11, pad=12)
    save(plt, fig, figdir, "hbond_occupancy")
    return True


# ==========================================================================

def table_html(header, rows, limit=None) -> str:
    fmt = lambda v: "" if (isinstance(v, float) and np.isnan(v)) else (f"{v:.4g}" if isinstance(v, float) else html.escape(str(v)))  # noqa: E731
    body = "".join("<tr>" + "".join(f"<td>{fmt(v)}</td>" for v in r) + "</tr>" for r in rows[:limit])
    head = "".join(f"<th>{html.escape(h)}</th>" for h in header)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def cmd_run(a: argparse.Namespace) -> None:
    need_modules("MDAnalysis", "matplotlib")
    warnings.filterwarnings("ignore")
    runs = Path(a.runs)
    reps = find_replicas(runs / a.system, a.reps)
    results = [r for r in (analyse_replica(rep, a) for rep in reps) if r is not None]
    if not results:
        sys.exit(f"nothing to analyse under {runs / a.system}: every replica needs "
                 f"analysis/clean.gro and analysis/clean.xtc (./scripts/04_postprocess.sh {a.system} <n>)")
    outdir = Path(a.out or REPO / "analysis_out" / "contacts" / a.system)
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    if any(r["fallback"] for r in results):
        log(window_note(results, a))
    res = pool(results, a, outdir)
    header, rows = summary(results, a, outdir)
    key_lines = key_report(res, results, a) if results[0]["has_rna"] else []
    if key_lines:
        (outdir / "key_residues.txt").write_text("\n".join(key_lines) + "\n", encoding="utf-8")
        print("\n".join(key_lines))

    plt = style()
    nrep = len(results)
    ref = results[0]
    figs = []
    key_ids = set(a.key_residues or [])
    if ref["has_rna"]:
        if timeseries(plt, results, "n_pairs", "Residue-pair contacts", "Protein-RNA contacts over time",
                      figdir, "contacts_time"):
            figs.append(("contacts_time", "Number of protein-RNA residue pairs closer than "
                         f"{a.cutoff:g} A in each frame. A drop to zero means the RNA left."))
        if timeseries(plt, results, "min_d", "Minimum distance (A)", "Closest protein-RNA approach",
                      figdir, "min_distance"):
            figs.append(("min_distance", "Smallest heavy-atom distance between protein and RNA per frame."))
        if occupancy_bars(plt, ref["P"]["labels"], res["p_any"], key_ids, ref["P"]["resids"],
                          "Protein residues contacting RNA", figdir, "protein_occupancy",
                          "occupancy (max over partners)"):
            figs.append(("protein_occupancy", "Per protein residue: strongest contact with any "
                         "nucleotide, mean over replicas, error bars = SD across replicas."))
        occupancy_bars(plt, ref["R"]["labels"], res["n_any"], set(), None,
                       "Nucleotides contacting protein", figdir, "nucleotide_occupancy",
                       "occupancy (max over partners)") and figs.append(
            ("nucleotide_occupancy", "Per nucleotide, same definition."))
        if contact_map(plt, ref["P"]["labels"], ref["R"]["labels"], res["mean_pr"], a, figdir):
            figs.append(("contact_map", "Occupancy of each protein residue - nucleotide pair, mean of replicas."))
        if res.get("hb_rows") and hbond_bars(plt, res["hb_rows"], nrep, figdir):
            figs.append(("hbond_occupancy", f"Distance < {a.hb_dist:g} A and angle > {a.hb_angle:g} deg."))
    if ref["occ_pp"] is not None and timeseries(
            plt, results, "n_iface", "Residue-pair contacts", "Dimer interface contacts over time",
            figdir, "interface_time"):
        figs.append(("interface_time", "Residue pairs between different protein chains."))

    sections = "".join(
        f"<section><h2>{html.escape(n.replace('_', ' '))}</h2><img src='figures/{n}.svg' "
        f"alt='{n}' style='max-width:100%'><p class='cap'>{html.escape(c)}</p></section>"
        for n, c in figs)
    top_rows = res.get("pairs_rows", [])
    pair_hdr = ["protein residue", "nucleotide"] + res["names"] + ["mean", "SD", f"n_reps>={a.min_occ:g}",
                                                                   "H-bond occ.", "top H-bond atoms"]
    tail = (f"<h2>Top contacts</h2><p class='cap'>The 25 most persistent pairs; the full list is "
            f"contacts_pairs.csv.</p>{table_html(pair_hdr, top_rows, 25)}") if top_rows else ""
    keyhtml = f"<pre>{html.escape(chr(10).join(key_lines))}</pre>" if key_lines else ""
    meta = registry_row(a.system)
    page = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{html.escape(a.system)} contacts</title><style>{PAGE_CSS}</style></head><body><main>"
            f"<h1>{html.escape(a.system)} - contacts</h1>"
            f"<p class='lede'>{html.escape(meta.get('consumer', ''))} - {nrep} replica(s)</p>"
            f"<h2>Summary</h2>{table_html(header, rows)}"
            f"<p class='cap'>Cutoff {a.cutoff:g} A. {html.escape(window_note(results, a))} Mean &plusmn; SD across "
            f"replicas. Geometry only; solute-only trajectory (no water bridges, no Mg&sup2;&#8314;).</p>"
            f"{keyhtml}{sections}{tail}<footer>Generated by analysis/contacts.py</footer></main></body></html>")
    (outdir / "index.html").write_text(page, encoding="utf-8")
    log(f"done: {outdir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system")
    ap.add_argument("--runs", default=None, help="runs directory (default: $BRUSAFE_RUNS, else $REPO/runs)")
    ap.add_argument("--out", help="output directory (default: analysis_out/contacts/<system>)")
    ap.add_argument("--reps", type=int, nargs="+", help="replica numbers (default: all)")
    ap.add_argument("--equil-ns", type=float, default=50.0, help="relaxation window to skip (default 50)")
    ap.add_argument("--cutoff", type=float, default=4.5, help="contact distance, heavy atoms, A (default 4.5)")
    ap.add_argument("--hb-dist", type=float, default=3.0, help="H-bond donor-acceptor distance, A")
    ap.add_argument("--hb-angle", type=float, default=150.0, help="H-bond D-H...A angle, degrees")
    ap.add_argument("--min-occ", type=float, default=0.5,
                    help="occupancy counted as 'persistent' in a replica (default 0.5)")
    ap.add_argument("--list-min", type=float, default=0.05, help="smallest mean occupancy written to CSV")
    ap.add_argument("--map-min", type=float, default=0.2, help="smallest occupancy shown in the contact map")
    ap.add_argument("--max-frames", type=int, default=1000, help="frames analysed per replica (default 1000)")
    ap.add_argument("--no-hbonds", action="store_true", help="skip the hydrogen-bond analysis")
    ap.add_argument("--chain-sizes", type=lambda s: [int(x) for x in s.split(",")],
                    help="protein residues per chain, e.g. 129,129 (default: split at numbering restarts)")
    ap.add_argument("--rna-first", type=int,
                    help="number of the first nucleotide in the paper's numbering, e.g. -15 for the "
                         "19-nt operator hairpin (-15..-1, +1..+4); default: topology numbering")
    ap.add_argument("--rna-sel", help="MDAnalysis selection for the RNA (default: common nucleotide names)")
    ap.add_argument("--key-residues", type=lambda s: [int(x) for x in s.split(",")],
                    help="residue numbers of the known operator-contact residues, e.g. 45,47,59")
    ap.add_argument("--top-n", type=int, default=15, help="top contacting residues used for the hot-spot check")
    a = ap.parse_args()
    import os
    a.runs = a.runs or os.environ.get("BRUSAFE_RUNS", str(REPO / "runs"))
    cmd_run(a)


if __name__ == "__main__":
    main()
