#!/usr/bin/env python3
"""
Figures, tables and trajectory animations for the BruSafe MD campaign.

Reads what 04_postprocess.sh already wrote (the .xvg/.dat files under
runs/<system>/rep*/analysis/) and turns it into wiki-ready output:

  * static figures   PNG (300 dpi) + SVG, one per analysis
  * interactive      Plotly HTML, one per analysis, plus index.html (dashboard)
  * animations       3Dmol.js trajectory viewer per replica (needs MDAnalysis)
  * tables           summary.csv / summary.md with replica means and spreads

Subcommands
-----------
  report  <system>              everything for one system, all replicas overlaid
  compare <sysA> <sysB> [...]   cross-system overlays (RMSF, RMSD, distributions)
  animate --top T --traj X -o F one standalone trajectory viewer

Quick start
-----------
  pip install numpy matplotlib plotly MDAnalysis
  python3 analysis/md_report.py report  S1_wt_cc_apo
  python3 analysis/md_report.py compare S1_wt_cc_apo S2_sccp_apo

Output goes to analysis_out/<system>/ (or analysis_out/compare_A__B/):

  figures/*.png, figures/*.svg     static figures for slides / the wiki
  html/index.html                  the dashboard -- open this first
  html/<plot>.html                 one interactive plot per file, for <iframe>
  html/traj_rep<N>.html            trajectory animation per replica
  summary.csv, summary.md          numbers that go into the wiki tables

iGEM wiki note
--------------
The iGEM wiki does not allow scripts from external CDNs. The default
--js local therefore copies plotly.min.js and 3Dmol-min.js next to the HTML
(upload the whole html/ folder). --js inline makes every page a single
self-contained file (bigger). --js cdn is only for quick local previews.

What the numbers mean -- read before quoting them
-------------------------------------------------
* Everything before --equil-ns (default 50 ns, same as 04_postprocess.sh) is
  treated as relaxation and excluded from means, distributions and RMSF.
* The +/- in the summary table is the spread ACROSS REPLICAS. The per-replica
  "block SEM" column is shown only as convergence evidence; it is not an
  error bar (README, reporting rules).
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
from convergence import block_average, cosine_content  # noqa: E402

VENDOR_3DMOL = HERE / "vendor" / "3Dmol-min.js"
CDN_3DMOL = "https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"

# Colour-blind-checked categorical order (validated light/dark). Replicas and
# systems take slots in fixed order; never cycled.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
MEAN_COLOR = "#2b2b2a"

SOLVENT = {"SOL", "TIP3", "TIP3P", "HOH", "WAT", "POT", "CLA", "SOD",
           "K", "CL", "NA"}
IONS = {"MG": "Mg", "MGA": "Mg", "ZN2": "Zn", "ZN": "Zn", "CAL": "Ca",
        "POT": "K", "CLA": "Cl", "SOD": "Na", "K": "K", "CL": "Cl", "NA": "Na"}
NUC_MAP = {"ADE": "A", "RA": "A", "RA5": "A", "RA3": "A", "A": "A",
           "GUA": "G", "RG": "G", "RG5": "G", "RG3": "G", "G": "G",
           "CYT": "C", "RC": "C", "RC5": "C", "RC3": "C", "C": "C",
           "URA": "U", "RU": "U", "RU5": "U", "RU3": "U", "U": "U"}

# gmx dssp one-letter codes -> (label, colour). Loops/breaks recede in grey.
DSSP = [("H", "α-helix", "#2a78d6"), ("G", "3₁₀-helix", "#6250d6"),
        ("I", "π-helix", "#e87ba4"), ("P", "PPII", "#008300"),
        ("E", "β-strand", "#eb6834"), ("B", "β-bridge", "#eda100"),
        ("T", "Turn", "#1baf7a"), ("S", "Bend", "#b9b8b1"),
        ("~", "Loop", "#f2f1ee"), ("=", "Break", "#ffffff")]


def log(msg: str) -> None:
    print(f"[md_report] {msg}", file=sys.stderr)


# ==========================================================================
# .xvg input

@dataclass
class Xvg:
    data: np.ndarray
    title: str = ""
    xlabel: str = ""
    ylabel: str = ""
    legends: list[str] = field(default_factory=list)

    @property
    def x(self) -> np.ndarray:
        return self.data[:, 0]

    def col(self, i: int = 1) -> np.ndarray:
        return self.data[:, i]

    @property
    def t_ns(self) -> np.ndarray:
        """x axis in ns, whatever unit gmx wrote it in."""
        lab = self.xlabel.lower()
        if "(ps)" in lab or lab.endswith("ps"):
            return self.x / 1000.0
        if "(fs)" in lab:
            return self.x / 1e6
        if "(us)" in lab or "(µs)" in lab:
            return self.x * 1000.0
        return self.x


def _clean_label(s: str) -> str:
    """Strip xmgrace markup: nm\\S2\\N -> nm²."""
    s = s.replace("\\S2\\N", "²").replace("\\S3\\N", "³")
    return re.sub(r"\\[A-Za-z]", "", s)


def read_xvg(path: Path) -> Xvg | None:
    if not path.is_file():
        return None
    rows, meta = [], {"title": "", "xlabel": "", "ylabel": ""}
    legends: dict[int, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] == "#":
                continue
            if line[0] == "@":
                m = re.match(r'@\s+(title|xaxis\s+label|yaxis\s+label)\s+"(.*)"', line)
                if m:
                    key = {"title": "title"}.get(m.group(1), m.group(1)[0] + "label")
                    meta[key] = _clean_label(m.group(2))
                m = re.match(r'@\s+s(\d+)\s+legend\s+"(.*)"', line)
                if m:
                    legends[int(m.group(1))] = _clean_label(m.group(2))
                continue
            if line[0] == "&":
                break                       # second data set: not used here
            try:
                rows.append([float(v) for v in line.split()])
            except ValueError:
                continue
    if not rows:
        log(f"no numeric data in {path}")
        return None
    width = min(len(r) for r in rows)
    data = np.array([r[:width] for r in rows], dtype=float)
    leg = [legends.get(i, f"s{i}") for i in range(width - 1)]
    return Xvg(data, meta["title"], meta["xlabel"], meta["ylabel"], leg)


def read_dssp(path: Path) -> list[str] | None:
    """gmx dssp -o dssp.dat: one line per frame, one character per residue."""
    if not path.is_file():
        return None
    lines = [ln.rstrip("\n") for ln in open(path) if ln.strip() and ln[0] not in "#@"]
    return lines or None


def chain_segments(resnr: np.ndarray) -> list[tuple[int, int]]:
    """[start, stop) index ranges of each chain: numbering restarts = new chain."""
    breaks = [0] + [i for i in range(1, len(resnr)) if resnr[i] <= resnr[i - 1]] + [len(resnr)]
    return [(breaks[i], breaks[i + 1]) for i in range(len(breaks) - 1)]


def running_mean(y: np.ndarray, n: int) -> np.ndarray:
    if n <= 1 or len(y) < n:
        return y
    k = np.ones(n) / n
    pad = np.pad(y, (n // 2, n - 1 - n // 2), mode="edge")
    return np.convolve(pad, k, mode="valid")


# ==========================================================================
# A tiny plot spec, rendered by both matplotlib (static) and plotly (html),
# so every figure is defined exactly once.

@dataclass
class Trace:
    x: np.ndarray
    y: np.ndarray
    name: str = ""
    color: str = SERIES[0]
    kind: str = "line"          # line | band | scatter | bar
    lower: np.ndarray | None = None
    upper: np.ndarray | None = None
    width: float = 2.0
    opacity: float = 1.0
    dash: str | None = None     # None | "dash" | "dot"
    legend: bool = True
    hover: list[str] | None = None
    cvalues: np.ndarray | None = None   # scatter coloured by value (e.g. time)


@dataclass
class Heatmap:
    z: np.ndarray               # (n_rows, n_cols) integer category codes
    x: np.ndarray
    y_labels: list[str]
    categories: list[tuple[str, str, str]]


@dataclass
class Plot:
    key: str
    title: str
    xlabel: str
    ylabel: str
    traces: list[Trace] = field(default_factory=list)
    caption: str = ""
    vlines: list[tuple[float, str]] = field(default_factory=list)
    vbands: list[tuple[float, float, str]] = field(default_factory=list)
    heatmap: Heatmap | None = None
    xticks: tuple[list[float], list[str]] | None = None
    logx: bool = False
    logy: bool = False
    height: int = 420


def render_mpl(p: Plot, outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    plt.rcParams.update({
        "font.size": 10, "axes.edgecolor": INK_2, "axes.labelcolor": INK,
        "xtick.color": INK_2, "ytick.color": INK_2, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "legend.frameon": False, "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(7.2, p.height / 70), dpi=100)
    dash = {None: "-", "dash": "--", "dot": ":"}

    if p.heatmap is not None:
        h = p.heatmap
        cmap = ListedColormap([c for _, _, c in h.categories])
        extent = (h.x[0], h.x[-1], -0.5, h.z.shape[0] - 0.5)
        ax.imshow(h.z, aspect="auto", origin="lower", cmap=cmap, vmin=-0.5,
                  vmax=len(h.categories) - 0.5, interpolation="nearest", extent=extent)
        ax.grid(False)
        present = sorted(set(np.unique(h.z).tolist()))
        handles = [plt.Rectangle((0, 0), 1, 1, fc=h.categories[i][2], ec=GRID)
                   for i in present]
        ax.legend(handles, [h.categories[i][1] for i in present], ncol=min(5, len(present)),
                  loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=8)

    for t in p.traces:
        lab = t.name if (t.legend and t.name) else None
        if t.kind == "band":
            ax.fill_between(t.x, t.lower, t.upper, color=t.color, alpha=0.18,
                            linewidth=0, label=lab)
        elif t.kind == "scatter":
            if t.cvalues is not None:
                sc = ax.scatter(t.x, t.y, c=t.cvalues, cmap="Blues", s=10,
                                edgecolors="none", label=lab)
                fig.colorbar(sc, ax=ax, label="time (ns)", pad=0.01)
            else:
                ax.scatter(t.x, t.y, color=t.color, s=10, alpha=t.opacity, label=lab)
        elif t.kind == "bar":
            ax.bar(t.x, t.y, color=t.color, width=0.8, label=lab)
        else:
            ax.plot(t.x, t.y, color=t.color, lw=t.width * 0.6, alpha=t.opacity,
                    ls=dash[t.dash], label=lab)

    for x0, x1, lab in p.vbands:
        ax.axvspan(x0, x1, color="#8a8984", alpha=0.10, lw=0)
        ax.text(x0 + (x1 - x0) * 0.02, 0.98, lab, transform=ax.get_xaxis_transform(),
                va="top", fontsize=8, color=INK_2)
    for x0, lab in p.vlines:
        ax.axvline(x0, color=INK_2, lw=0.7, ls=":")
        if lab:
            ax.text(x0, 1.0, f" {lab}", transform=ax.get_xaxis_transform(),
                    va="bottom", fontsize=8, color=INK_2)
    if p.xticks:
        ax.set_xticks(p.xticks[0], p.xticks[1])
    if p.heatmap is not None and len(p.heatmap.y_labels) <= 40:
        ax.set_yticks(range(len(p.heatmap.y_labels)), p.heatmap.y_labels, fontsize=7)
    if p.logx:
        ax.set_xscale("log")
    if p.logy:
        ax.set_yscale("log")
    ax.set_title(p.title, loc="left", fontsize=11, color=INK, pad=14)
    ax.set_xlabel(p.xlabel)
    ax.set_ylabel(p.ylabel)
    if p.heatmap is None and sum(1 for t in p.traces if t.legend and t.name) >= 2:
        ax.legend(fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(outdir / f"{p.key}.{ext}", dpi=300 if ext == "png" else None,
                    bbox_inches="tight")
    plt.close(fig)


# Interactive pages carry at most this many points per line. The static
# figures keep every frame; the html is for looking, and the wiki has a size limit.
HTML_MAX_POINTS = 1500


def _thin(t: Trace) -> tuple[np.ndarray, np.ndarray, list[str] | None]:
    if t.kind != "line" or len(t.x) <= HTML_MAX_POINTS:
        return t.x, t.y, t.hover
    k = int(np.ceil(len(t.x) / HTML_MAX_POINTS))
    return t.x[::k], t.y[::k], (t.hover[::k] if t.hover else None)


def to_plotly(p: Plot):
    import plotly.graph_objects as go

    fig = go.Figure()
    if p.heatmap is not None:
        h = p.heatmap
        n = len(h.categories)
        scale = []
        for i, (_, _, c) in enumerate(h.categories):
            scale += [[i / n, c], [(i + 1) / n, c]]
        z, x = h.z, h.x
        if z.shape[1] > 300:                       # keep the page light
            k = int(np.ceil(z.shape[1] / 300))
            z, x = z[:, ::k], x[::k]
        labels = np.array([c[0] for c in h.categories], dtype=object)[z]
        fig.add_trace(go.Heatmap(
            z=z, x=x, y=h.y_labels, colorscale=scale, zmin=-0.5, zmax=n - 0.5,
            customdata=labels, showscale=False,
            hovertemplate="t = %{x:.1f} ns<br>%{y}<br>DSSP %{customdata}<extra></extra>"))
        for i in sorted(set(np.unique(h.z).tolist())):   # legend entries
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=h.categories[i][1],
                                     marker=dict(size=10, symbol="square",
                                                 color=h.categories[i][2],
                                                 line=dict(color=GRID, width=1))))

    for t in p.traces:
        common = dict(name=t.name, showlegend=t.legend and bool(t.name), legendgroup=t.name)
        if t.kind == "band":
            fig.add_trace(go.Scatter(
                x=np.concatenate([t.x, t.x[::-1]]),
                y=np.concatenate([t.upper, t.lower[::-1]]),
                fill="toself", fillcolor=_rgba(t.color, 0.18), line=dict(width=0),
                hoverinfo="skip", **common))
        elif t.kind == "scatter":
            marker = dict(size=6, color=t.color, opacity=t.opacity)
            if t.cvalues is not None:
                marker.update(color=t.cvalues, colorscale="Blues",
                              colorbar=dict(title="ns", thickness=12))
            fig.add_trace(go.Scatter(x=t.x, y=t.y, mode="markers", marker=marker,
                                     text=t.hover, **common))
        elif t.kind == "bar":
            fig.add_trace(go.Bar(x=t.x, y=t.y, marker_color=t.color, text=t.hover,
                                 **common))
        else:
            x, y, hover = _thin(t)
            fig.add_trace(go.Scatter(
                x=x, y=y, mode="lines", opacity=t.opacity, text=hover,
                line=dict(color=t.color, width=t.width, dash=t.dash or "solid"),
                hovertemplate=("%{text}<br>" if hover else "") +
                f"{t.name}: %{{y:.3g}}<extra></extra>", **common))

    for x0, x1, lab in p.vbands:
        fig.add_vrect(x0=x0, x1=x1, fillcolor="#8a8984", opacity=0.10, line_width=0,
                      annotation_text=lab, annotation_position="top left",
                      annotation_font=dict(size=11, color=INK_2))
    for x0, lab in p.vlines:
        fig.add_vline(x=x0, line=dict(color=INK_2, width=1, dash="dot"),
                      annotation_text=lab, annotation_font=dict(size=11, color=INK_2))
    axis = dict(gridcolor=GRID, zeroline=False, linecolor=INK_2, ticks="outside",
                tickcolor=GRID)
    fig.update_layout(
        title=dict(text=p.title, x=0.01, font=dict(size=16, color=INK)),
        xaxis=dict(title=p.xlabel, type="log" if p.logx else None, **axis),
        yaxis=dict(title=p.ylabel, type="log" if p.logy else None, **axis),
        template="plotly_white", height=p.height, hovermode="x unified"
        if p.heatmap is None and not any(t.kind == "scatter" for t in p.traces) else "closest",
        margin=dict(l=60, r=20, t=60, b=50), font=dict(family="system-ui, sans-serif",
                                                       color=INK, size=12),
        legend=dict(orientation="h", y=-0.2, x=0),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="#fcfcfb")
    if p.xticks:
        fig.update_xaxes(tickvals=p.xticks[0], ticktext=p.xticks[1])
    return fig


def _rgba(hexcol: str, a: float) -> str:
    h = hexcol.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{a})"


# ==========================================================================
# HTML output

PAGE_CSS = """
:root{--bg:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--rule:#e4e3df;--card:#ffffff;
      --accent:#2a78d6}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
     font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:1.6rem;margin:0 0 4px}
h2{font-size:1.15rem;margin:40px 0 8px;padding-top:12px;border-top:1px solid var(--rule)}
p.lede,p.cap{color:var(--ink2);margin:4px 0 12px}
p.cap{font-size:.9rem;max-width:75ch}
table{border-collapse:collapse;font-size:.9rem;font-variant-numeric:tabular-nums;
      width:100%;display:block;overflow-x:auto}
th,td{padding:6px 10px;border-bottom:1px solid var(--rule);text-align:right;
      white-space:nowrap}
th:first-child,td:first-child{text-align:left}
th{color:var(--ink2);font-weight:600}
.card{background:var(--card);border:1px solid var(--rule);border-radius:8px;
      padding:8px;margin:8px 0}
iframe{width:100%;height:620px;border:1px solid var(--rule);border-radius:8px;
       background:var(--card)}
nav a{color:var(--accent);margin-right:12px;font-size:.9rem}
footer{color:var(--ink2);font-size:.8rem;margin-top:48px}
"""


class JsMode:
    """Where plotly.js / 3Dmol.js come from: local file, inline, or CDN."""

    def __init__(self, mode: str, outdir: Path):
        self.mode, self.outdir = mode, outdir

    def plotly_tag(self) -> str:
        import plotly.offline
        if self.mode == "cdn":
            return '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
        src = plotly.offline.get_plotlyjs()
        if self.mode == "inline":
            return f"<script>{src}</script>"
        dest = self.outdir / "plotly.min.js"
        if not dest.exists():
            dest.write_text(src, encoding="utf-8")
        return '<script src="plotly.min.js"></script>'

    def mol_tag(self) -> str:
        if self.mode == "cdn" or not VENDOR_3DMOL.is_file():
            if self.mode != "cdn":
                log(f"{VENDOR_3DMOL} missing -- falling back to the CDN copy of 3Dmol.js")
            return f'<script src="{CDN_3DMOL}"></script>'
        if self.mode == "inline":
            return "<script>" + VENDOR_3DMOL.read_text(encoding="utf-8") + "</script>"
        dest = self.outdir / "3Dmol-min.js"
        if not dest.exists():
            shutil.copy(VENDOR_3DMOL, dest)
            lic = VENDOR_3DMOL.parent / "3Dmol.LICENSE.txt"
            if lic.is_file():
                shutil.copy(lic, self.outdir / "3Dmol.LICENSE.txt")
        return '<script src="3Dmol-min.js"></script>'


def write_plot_page(p: Plot, html_dir: Path, js: JsMode) -> str:
    """Standalone page for one plot; returns the dashboard fragment."""
    fig = to_plotly(p)
    cfg = {"displaylogo": False, "responsive": True,
           "toImageButtonOptions": {"format": "svg", "filename": p.key}}
    div = fig.to_html(full_html=False, include_plotlyjs=False, config=cfg, div_id=p.key)
    page = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{html.escape(p.title)}</title>{js.plotly_tag()}"
            f"<style>{PAGE_CSS}main{{padding:8px}}</style></head><body><main>"
            f"{div}<p class='cap'>{html.escape(p.caption)}</p></main></body></html>")
    (html_dir / f"{p.key}.html").write_text(page, encoding="utf-8")
    return (f"<section id='sec-{p.key}'><div class='card'>{div}</div>"
            f"<p class='cap'>{html.escape(p.caption)}</p></section>")


def write_dashboard(html_dir: Path, title: str, lede: str, table_html: str,
                    fragments: list[tuple[str, str]], animations: list[str],
                    js: JsMode) -> None:
    nav = "".join(f"<a href='#sec-{k}'>{html.escape(h)}</a>" for k, h in
                  [("summary", "Summary")] + [(k, k.replace('_', ' ')) for k, _ in fragments]
                  + ([("traj", "Trajectories")] if animations else []))
    body = [f"<h1>{html.escape(title)}</h1><p class='lede'>{html.escape(lede)}</p>",
            f"<nav>{nav}</nav>",
            f"<h2 id='sec-summary'>Summary</h2>{table_html}"]
    for key, frag in fragments:
        body.append(frag)
    if animations:
        body.append("<h2 id='sec-traj'>Trajectories</h2>"
                    "<p class='cap'>Protein/RNA heavy atoms, aligned on the backbone. "
                    "Colour by RMSF to see which regions move most.</p>")
        for a in animations:
            body.append(f"<iframe src='{a}' loading='lazy' title='{a}'></iframe>")
    page = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{html.escape(title)}</title>{js.plotly_tag()}"
            f"<style>{PAGE_CSS}</style></head><body><main>{''.join(body)}"
            f"<footer>Generated by analysis/md_report.py</footer></main></body></html>")
    (html_dir / "index.html").write_text(page, encoding="utf-8")


def emit(plots: list[Plot], outdir: Path, js_mode: str) -> tuple[Path, JsMode, list]:
    fig_dir, html_dir = outdir / "figures", outdir / "html"
    fig_dir.mkdir(parents=True, exist_ok=True)
    html_dir.mkdir(parents=True, exist_ok=True)
    js = JsMode(js_mode, html_dir)
    frags = []
    for p in plots:
        try:
            render_mpl(p, fig_dir)
        except Exception as e:      # never lose the html because a png failed
            log(f"static figure {p.key} failed: {e}")
        frags.append((p.key, write_plot_page(p, html_dir, js)))
        log(f"  {p.key}")
    return html_dir, js, frags


# ==========================================================================
# Per-system analysis

@dataclass
class Replica:
    name: str
    path: Path
    color: str
    files: dict = field(default_factory=dict)

    def x(self, name: str) -> Xvg | None:
        if name not in self.files:
            self.files[name] = read_xvg(self.path / name)
        return self.files[name]


def find_replicas(sys_dir: Path, wanted: list[int] | None) -> list[Replica]:
    reps = []
    for d in sorted(sys_dir.glob("rep*"), key=lambda p: int(re.sub(r"\D", "", p.name) or 0)):
        n = int(re.sub(r"\D", "", d.name) or 0)
        if wanted and n not in wanted:
            continue
        if (d / "analysis").is_dir():
            reps.append(Replica(d.name, d / "analysis", SERIES[len(reps) % len(SERIES)]))
    if len(reps) > len(SERIES):
        log("more than 8 replicas: colours would repeat -- plot fewer with --reps")
    return reps


def registry_row(system: str) -> dict:
    tsv = REPO / "config" / "systems.tsv"
    if not tsv.is_file():
        return {}
    lines = [ln for ln in open(tsv) if ln.strip() and not ln.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        if row.get("name") == system:
            return row
    return {}


def timeseries_plot(reps: list[Replica], fname: str, key: str, title: str, ylabel: str,
                    caption: str, equil: float, col: int = 1, scale: float = 1.0,
                    smooth_ns: float = 1.0) -> Plot | None:
    p = Plot(key, title, "Time (ns)", ylabel, caption=caption)
    for r in reps:
        x = r.x(fname)
        if x is None or x.data.shape[1] <= col:
            continue
        t, y = x.t_ns, x.col(col) * scale
        dt = np.median(np.diff(t)) if len(t) > 1 else 1.0
        n = max(1, int(round(smooth_ns / dt))) if dt > 0 else 1
        if n > 1:
            p.traces.append(Trace(t, y, r.name + " raw", r.color, width=1, opacity=0.25,
                                  legend=False))
            p.traces.append(Trace(t, running_mean(y, n), r.name, r.color))
        else:
            p.traces.append(Trace(t, y, r.name, r.color))
    if not p.traces:
        return None
    if equil > 0:
        p.vbands.append((0, equil, "relaxation (excluded)"))
    return p


def distribution_plot(series: list[tuple[str, str, np.ndarray]], key: str, title: str,
                      xlabel: str, caption: str) -> Plot | None:
    series = [s for s in series if len(s[2]) > 5]
    if not series:
        return None
    lo = min(s[2].min() for s in series)
    hi = max(s[2].max() for s in series)
    edges = np.linspace(lo, hi, 61)
    mid = 0.5 * (edges[1:] + edges[:-1])
    p = Plot(key, title, xlabel, "Probability density", caption=caption, height=360)
    for name, color, v in series:
        dens, _ = np.histogram(v, bins=edges, density=True)
        p.traces.append(Trace(mid, dens, name, color))
        p.vlines.append((float(v.mean()), ""))
    return p


def rmsf_plot(reps: list[Replica], key: str = "rmsf",
              title: str = "Cα RMSF per residue") -> tuple[Plot | None, np.ndarray | None]:
    prof = [(r, r.x("rmsf_ca.xvg")) for r in reps]
    prof = [(r, x) for r, x in prof if x is not None]
    if not prof:
        return None, None
    n = min(len(x.x) for _, x in prof)
    if any(len(x.x) != n for _, x in prof):
        log("RMSF profiles differ in length between replicas -- truncating to shortest")
    resnr = prof[0][1].x[:n].astype(int)
    segs = chain_segments(resnr)
    idx = np.arange(n)
    letters = [chr(ord("A") + i) for i in range(len(segs))]
    chain_of = np.empty(n, dtype=object)
    for (a, b), c in zip(segs, letters):
        chain_of[a:b] = c
    hover = [f"chain {c} · res {r}" for c, r in zip(chain_of, resnr)]
    p = Plot(key, title, "Residue", "RMSF (Å)",
             caption="Cα root-mean-square fluctuation after the relaxation window "
                     "(04_postprocess.sh uses -b 50 ns). Thin lines: replicas; "
                     "dark line and band: mean ± SD across replicas — that spread is "
                     "the error bar. Peaks are flexible loops (e.g. the FG loop).")
    mat = np.array([x.col(1)[:n] * 10 for _, x in prof])
    for r, x in prof:
        p.traces.append(Trace(idx, x.col(1)[:n] * 10, r.name, r.color, width=1.4,
                              opacity=0.8, hover=hover))
    if len(prof) > 1:
        m, s = mat.mean(0), mat.std(0, ddof=1)
        p.traces.append(Trace(idx, m, "± SD", MEAN_COLOR, kind="band", lower=m - s, upper=m + s))
        p.traces.append(Trace(idx, m, "mean", MEAN_COLOR, width=2.4, hover=hover))
    # ticks at round residue numbers, labelled with the real numbering
    ticks, labels = [], []
    step = next((s for s in (5, 10, 20, 25, 50, 100, 200, 500) if n / s <= 14), 1000)
    for (a, b), c in zip(segs, letters):
        for i in range(a, b):
            if resnr[i] % step == 0:
                ticks.append(i)
                labels.append(f"{resnr[i]}" if len(segs) == 1 else f"{c}{resnr[i]}")
        if a > 0:
            p.vlines.append((a - 0.5, f"chain {c}"))
    p.xticks = (ticks, labels)
    return p, mat


def dssp_plot(reps: list[Replica], t_ref: dict[str, np.ndarray]) -> list[Plot]:
    plots = []
    codes = {c: i for i, (c, _, _) in enumerate(DSSP)}
    for r in reps:
        frames = read_dssp(r.path / "dssp.dat")
        if not frames:
            continue
        nres = min(len(f) for f in frames)
        stride = max(1, len(frames) // 800)
        frames = frames[::stride]
        z = np.array([[codes.get(ch, codes["~"]) for ch in f[:nres]] for f in frames]).T
        t = t_ref.get(r.name)
        t = (np.linspace(t[0], t[-1], len(frames)) if t is not None and len(t)
             else np.arange(len(frames), dtype=float))
        ylab = [f"res {i + 1}" for i in range(nres)]
        p = Plot(f"dssp_{r.name}", f"Secondary structure over time — {r.name}",
                 "Time (ns)", "Residue (sequential index)",
                 caption="gmx dssp assignment per residue and frame. Persistent "
                         "colour bands are stable elements; flicker in loop regions is "
                         "normal; a band that disappears mid-run is a real unfolding event.",
                 heatmap=Heatmap(z, t, ylab, DSSP), height=520)
        plots.append(p)
    return plots


def summarise(reps: list[Replica], equil: float) -> tuple[list[dict], list[str]]:
    """Per-replica numbers over the production window, + across-replica mean/SD."""
    metrics = [  # file, column, label, scale, unit
        ("rmsd_backbone.xvg", 1, "Backbone RMSD", 10.0, "Å"),
        ("rmsd_rna.xvg", 1, "RNA RMSD", 10.0, "Å"),
        ("gyrate.xvg", 1, "Rg", 10.0, "Å"),
        ("sasa.xvg", 1, "SASA", 1.0, "nm²"),
        ("hbond_num.xvg", 1, "H-bonds", 1.0, ""),
        ("mindist_mg.xvg", 1, "Mg²⁺ min. distance", 10.0, "Å"),
    ]
    rows, cols = [], []
    for r in reps:
        row = {"replica": r.name}
        for fname, c, lab, sc, unit in metrics:
            x = r.x(fname)
            if x is None:
                continue
            t, y = x.t_ns, x.col(c) * sc
            y = y[t >= equil] if (t >= equil).sum() > 10 else y
            sizes, sems = block_average(y) if len(y) >= 16 else (None, np.array([np.nan]))
            head = f"{lab} ({unit})" if unit else lab
            row[head] = y.mean()
            row[f"{head} block SEM"] = sems[len(sems) * 2 // 3:].mean()
            if head not in cols:
                cols += [head, f"{head} block SEM"]
            if fname == "rmsd_backbone.xvg":
                row["length (ns)"] = t[-1]
        rmsf = r.x("rmsf_ca.xvg")
        if rmsf is not None:
            row["mean Cα RMSF (Å)"] = rmsf.col(1).mean() * 10
        pc1 = r.x("pc1.xvg")
        if pc1 is not None and len(pc1.x) > 4:
            row["PC1 cosine content"] = cosine_content(pc1.col(1), 1)
        rows.append(row)
    for extra in ("mean Cα RMSF (Å)", "PC1 cosine content", "length (ns)"):
        if any(extra in r for r in rows):
            cols.append(extra)
    if len(rows) > 1:
        agg_m, agg_s = {"replica": "mean"}, {"replica": "SD across replicas"}
        for c in cols:
            if "block SEM" in c:
                continue
            v = np.array([r[c] for r in rows if c in r and np.isfinite(r[c])])
            if len(v):
                agg_m[c] = v.mean()
            if len(v) > 1:
                agg_s[c] = v.std(ddof=1)
        rows += [agg_m, agg_s]
    return rows, cols


def table_outputs(rows: list[dict], cols: list[str], outdir: Path, note: str) -> str:
    fmt = lambda v: "" if v is None else (f"{v:.3g}" if isinstance(v, float) else str(v))  # noqa: E731
    with open(outdir / "summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["replica"] + cols)
        for r in rows:
            w.writerow([r["replica"]] + [fmt(r.get(c)) for c in cols])
    md = ["| replica | " + " | ".join(cols) + " |", "|---" * (len(cols) + 1) + "|"]
    md += ["| " + r["replica"] + " | " + " | ".join(fmt(r.get(c)) for c in cols) + " |"
           for r in rows]
    (outdir / "summary.md").write_text("\n".join(md) + f"\n\n{note}\n", encoding="utf-8")
    th = "".join(f"<th>{html.escape(c)}</th>" for c in ["replica"] + cols)
    tr = "".join("<tr>" + "".join(f"<td>{html.escape(fmt(r.get(c)) if c != 'replica' else r[c])}</td>"
                                  for c in ["replica"] + cols) + "</tr>" for r in rows)
    return f"<table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table><p class='cap'>{html.escape(note)}</p>"


def need_modules(*names: str) -> None:
    """Fail up front, with the install line, instead of a traceback halfway through."""
    import importlib.util
    missing = [n for n in names if importlib.util.find_spec(n) is None]
    if missing:
        sys.exit(f"missing Python package(s): {', '.join(missing)}\n"
                 f"install with:  python3 -m pip install --user {' '.join(missing)}")


def cmd_report(a: argparse.Namespace) -> None:
    need_modules("matplotlib", "plotly")
    runs = Path(a.runs)
    sys_dir = runs / a.system
    reps = find_replicas(sys_dir, a.reps)
    # 04_postprocess.sh creates analysis/ before it checks for prod.xtc, so an
    # empty analysis/ is common: only count replicas that have actual output
    empty = [r.name for r in reps if not any(r.path.glob("*.xvg"))]
    reps = [r for r in reps if r.name not in empty]
    for name in empty:
        log(f"{name}: analysis/ has no .xvg files -- skipped")
    if not reps:
        sys.exit(f"nothing to plot under {sys_dir}: run "
                 f"./scripts/04_postprocess.sh {a.system} <replica> first "
                 f"(it needs prod.xtc, prod.tpr and index.ndx in rep<N>/)")
    outdir = Path(a.out or REPO / "analysis_out" / a.system)
    log(f"{a.system}: {', '.join(r.name for r in reps)} -> {outdir}")
    eq = a.equil_ns
    plots: list[Plot] = []

    rmsd = timeseries_plot(reps, "rmsd_backbone.xvg", "rmsd_backbone", "Backbone RMSD",
                           "RMSD (Å)", "Backbone RMSD to the production starting structure, "
                           "least-squares fitted on the backbone. Faint: raw; solid: 1 ns "
                           "running mean. A plateau after the shaded window means the "
                           "structure has relaxed; a steady climb means it has not.",
                           eq, scale=10.0)
    plots += [rmsd] if rmsd else []
    p = timeseries_plot(reps, "rmsd_rna.xvg", "rmsd_rna", "RNA RMSD", "RMSD (Å)",
                        "RMSD of the RNA, fitted on the RNA itself. Large values with a "
                        "stable protein usually mean the hairpin loop or the termini are "
                        "fraying.", eq, scale=10.0)
    plots += [p] if p else []

    dist = []
    for fname, key, lab, sc, unit in [("rmsd_backbone.xvg", "dist_rmsd", "Backbone RMSD", 10, "Å"),
                                      ("gyrate.xvg", "dist_rg", "Radius of gyration", 10, "Å"),
                                      ("sasa.xvg", "dist_sasa", "SASA", 1, "nm²")]:
        series = []
        for r in reps:
            x = r.x(fname)
            if x is not None:
                series.append((r.name, r.color, x.col(1)[x.t_ns >= eq] * sc))
        d = distribution_plot(series, key, f"{lab} distribution (t ≥ {eq:g} ns)",
                              f"{lab} ({unit})", "Distribution over the production window, "
                              "one curve per replica; dotted lines mark replica means. "
                              "Curves that do not overlap mean the replicas sample different "
                              "states — report that, do not average it away.")
        if d:
            dist.append(d)

    rmsf, _ = rmsf_plot(reps)
    plots += [rmsf] if rmsf else []

    p = timeseries_plot(reps, "gyrate.xvg", "rg", "Radius of gyration (backbone)", "Rg (Å)",
                        "Compactness of the solute. Drift indicates opening/closing of the "
                        "dimer or partial unfolding.", eq, scale=10.0)
    plots += [p] if p else []
    p = timeseries_plot(reps, "sasa.xvg", "sasa", "Solvent-accessible surface area",
                        "SASA (nm²)", "Total SASA of the solute (gmx sasa). Rising SASA "
                        "together with rising Rg is the signature of unfolding or "
                        "dissociation.", eq)
    plots += [p] if p else []
    p = timeseries_plot(reps, "hbond_num.xvg", "hbonds", "Hydrogen bonds within the solute",
                        "Number of H-bonds", "Count of solute–solute hydrogen bonds "
                        "(gmx hbond, default geometric criteria).", eq)
    plots += [p] if p else []
    p = timeseries_plot(reps, "mindist_mg.xvg", "mindist_mg", "Mg²⁺ – solute minimum distance",
                        "Distance (Å)", "Closest approach of any Mg²⁺ to the solute. An ion "
                        "sitting below ~3 Å inside the protein–RNA interface contaminates "
                        "every MM/PBSA energy from that replica (TUTORIAL §7).", eq,
                        scale=10.0, smooth_ns=0)
    plots += [p] if p else []

    # DSSP: counts over time + per-residue map
    t_ref = {r.name: r.x("rmsd_backbone.xvg").t_ns for r in reps if r.x("rmsd_backbone.xvg")}
    for r in reps:
        x = r.x("dssp_num.xvg")
        if x is None:
            continue
        t = t_ref.get(r.name, x.x)
        t = np.linspace(t[0], t[-1], len(x.x)) if len(t) != len(x.x) else t
        p = Plot(f"dssp_counts_{r.name}", f"Secondary-structure content — {r.name}",
                 "Time (ns)", "Residues", caption="Number of residues in each DSSP class.")
        lut = {lab.lower(): col for _, lab, col in DSSP}
        for i, leg in enumerate(x.legends):
            col = next((c for k, c in lut.items() if k.split("-")[0] in leg.lower()), SERIES[i % 8])
            p.traces.append(Trace(t, running_mean(x.col(i + 1), 5), leg, col))
        plots.append(p)
    plots += dssp_plot(reps, t_ref)

    # PCA
    ev = [(r, r.x("eigenval.xvg")) for r in reps]
    ev = [(r, x) for r, x in ev if x is not None]
    if ev:
        p = Plot("pca_eigenvalues", "PCA: cumulative variance explained",
                 "Eigenvector index", "Cumulative variance (%)",
                 caption="Cα covariance analysis (gmx covar). If the first ~10 "
                         "eigenvectors capture most of the variance, the dynamics are "
                         "dominated by a few collective motions worth showing.")
        for r, x in ev:
            k = min(30, len(x.x))
            p.traces.append(Trace(x.x[:k], 100 * np.cumsum(x.col(1))[:k] / x.col(1).sum(),
                                  r.name, r.color))
        plots.append(p)
    p = Plot("pca_pc1", "Projection on PC1", "Time (ns)", "PC1 (nm)",
             caption="Projection of each replica onto its own first eigenvector. "
                     "Cosine content near 1 means the 'motion' is indistinguishable "
                     "from random diffusion — do not interpret it (Hess 2002).")
    for r in reps:
        x = r.x("pc1.xvg")
        if x is not None:
            cc = cosine_content(x.col(1), 1)
            p.traces.append(Trace(x.t_ns, x.col(1), f"{r.name} (cos = {cc:.2f})", r.color,
                                  width=1.2))
    plots += [p] if p.traces else []
    for r in reps:
        x = r.x("pc12.xvg")
        if x is not None and x.data.shape[1] >= 2:     # gmx anaeig -2d: PC1, PC2
            # -2d has no time column; pc1.xvg covers the same frames (-b 50 ns)
            pc1 = r.x("pc1.xvg")
            tt = pc1.t_ns if pc1 is not None and len(pc1.x) == len(x.x) \
                else np.arange(len(x.x), dtype=float)
            plots.append(Plot(f"pca_2d_{r.name}", f"Essential subspace PC1 × PC2 — {r.name}",
                              "PC1 (nm)", "PC2 (nm)",
                              [Trace(x.x, x.col(1), r.name, r.color, kind="scatter",
                                     cvalues=tt, legend=False,
                                     hover=[f"t = {v:.1f} ns" for v in tt])],
                              caption="Each point is a frame, coloured by time. Separate "
                                      "clusters are distinct conformational states.",
                              height=480))

    # Convergence: block-averaged SEM of the backbone RMSD
    p = Plot("block_sem", "Convergence: block-averaged SEM of backbone RMSD",
             "Block length (ns)", "SEM of the mean (Å)", logx=True,
             caption="The naive SEM (block length = one frame, left end) assumes "
                     "independent frames and is wrong. The plateau on the right is the "
                     "honest uncertainty of one replica; no plateau = run too short. "
                     "Quote the ratio in the methods (convergence.py block).")
    for r in reps:
        x = r.x("rmsd_backbone.xvg")
        if x is None:
            continue
        m = x.t_ns >= eq
        y = x.col(1)[m] * 10
        if len(y) < 32:
            continue
        sizes, sems = block_average(y)
        dt = np.median(np.diff(x.t_ns))
        p.traces.append(Trace(sizes * dt, sems, r.name, r.color))
    plots += [p] if p.traces else []
    plots += dist

    html_dir, js, frags = emit(plots, outdir, a.js)

    rows, cols = summarise(reps, eq)
    note = (f"Means over t ≥ {eq:g} ns. Report mean ± SD across replicas; the per-replica "
            "block SEM is convergence evidence, not the error bar.")
    table = table_outputs(rows, cols, outdir, note)

    anims = []
    if not a.no_anim:
        for r in reps:
            top, trj = r.path / "clean.gro", r.path / "clean.xtc"
            if not (top.is_file() and trj.is_file()):
                log(f"{r.name}: no clean.gro/clean.xtc -- skipping animation")
                continue
            out = html_dir / f"traj_{r.name}.html"
            try:
                build_animation(top, trj, out, f"{a.system} · {r.name}", a.frames, a.sel,
                                eq, a.js, a.smooth)
                anims.append(out.name)
            except ImportError:
                log("MDAnalysis not installed -- skipping animations (pip install MDAnalysis)")
                break

    meta = registry_row(a.system)
    lede = " · ".join(v for v in [meta.get("consumer", ""),
                                  f"{len(reps)} replica(s)",
                                  f"source {meta['source_pdb']}" if meta.get("source_pdb") else ""] if v)
    write_dashboard(html_dir, f"{a.system} — MD analysis", lede, table, frags, anims, js)
    log(f"done: open {html_dir / 'index.html'}")


# ==========================================================================
# Cross-system comparison

def cmd_compare(a: argparse.Namespace) -> None:
    need_modules("matplotlib", "plotly")
    runs = Path(a.runs)
    systems = [(s, find_replicas(runs / s, None)) for s in a.systems]
    for s, reps in systems:
        if not reps:
            sys.exit(f"no replicas with analysis/ for {s}")
    outdir = Path(a.out or REPO / "analysis_out" / ("compare_" + "__".join(a.systems)))
    eq = a.equil_ns
    log(f"comparing {', '.join(a.systems)} -> {outdir}")
    plots = []

    # RMSF: mean ± SD across replicas per system
    p = Plot("cmp_rmsf", "Cα RMSF — mean ± SD across replicas", "Residue (sequential index)",
             "RMSF (Å)", caption="A difference between systems only means something where "
             "the bands separate, i.e. where it exceeds the inter-replica spread. Systems "
             "with extra residues (e.g. the scCP linker) are not residue-aligned here: "
             "compare by sequence position, and use the eigenvector overlap from "
             "compare_systems.sh for the dynamics claim.")
    for i, (s, reps) in enumerate(systems):
        _, mat = rmsf_plot(reps)
        if mat is None:
            continue
        idx = np.arange(mat.shape[1])
        m = mat.mean(0)
        sd = mat.std(0, ddof=1) if len(mat) > 1 else np.zeros_like(m)
        p.traces.append(Trace(idx, m, s, SERIES[i], kind="band", lower=m - sd, upper=m + sd,
                              legend=False))
        p.traces.append(Trace(idx, m, s, SERIES[i]))
    plots += [p] if p.traces else []

    for fname, key, lab, sc, unit in [("rmsd_backbone.xvg", "cmp_rmsd", "Backbone RMSD", 10, "Å"),
                                      ("gyrate.xvg", "cmp_rg", "Radius of gyration", 10, "Å"),
                                      ("sasa.xvg", "cmp_sasa", "SASA", 1, "nm²"),
                                      ("hbond_num.xvg", "cmp_hbond", "H-bonds", 1, "")]:
        p = Plot(key, f"{lab} — mean ± SD across replicas", "Time (ns)",
                 f"{lab} ({unit})" if unit else lab,
                 caption="Replicas are averaged frame by frame (truncated to the shortest); "
                         "band = SD across replicas. Curves smoothed with a 1 ns window.")
        series = []
        for i, (s, reps) in enumerate(systems):
            xs = [r.x(fname) for r in reps]
            xs = [x for x in xs if x is not None]
            if not xs:
                continue
            n = min(len(x.x) for x in xs)
            t = xs[0].t_ns[:n]
            mat = np.array([x.col(1)[:n] * sc for x in xs])
            dt = np.median(np.diff(t)) if n > 1 else 1
            w = max(1, int(round(1.0 / dt))) if dt > 0 else 1
            sm = np.array([running_mean(row, w) for row in mat])
            m = sm.mean(0)
            sd = sm.std(0, ddof=1) if len(sm) > 1 else np.zeros_like(m)
            p.traces.append(Trace(t, m, s, SERIES[i], kind="band", lower=m - sd, upper=m + sd,
                                  legend=False))
            p.traces.append(Trace(t, m, s, SERIES[i]))
            series.append((s, SERIES[i], mat[:, t >= eq].ravel()))
        if p.traces:
            p.vbands.append((0, eq, "relaxation"))
            plots.append(p)
            d = distribution_plot(series, f"{key}_dist", f"{lab} distribution (t ≥ {eq:g} ns)",
                                  f"{lab} ({unit})" if unit else lab,
                                  "All replicas pooled per system. Dotted lines: system means.")
            if d:
                plots.append(d)

    html_dir, js, frags = emit(plots, outdir, a.js)
    rows = []
    for s, reps in systems:
        r, cols = summarise(reps, eq)
        for row in r:
            row["replica"] = f"{s} · {row['replica']}"
        rows += r
    allcols = []
    for row in rows:
        allcols += [c for c in row if c != "replica" and c not in allcols and "block SEM" not in c]
    table = table_outputs(rows, allcols, outdir,
                          f"Means over t ≥ {eq:g} ns; ± is the SD across replicas.")
    write_dashboard(html_dir, " vs ".join(a.systems), "Cross-system comparison", table,
                    frags, [], js)
    log(f"done: open {html_dir / 'index.html'}")


# ==========================================================================
# Trajectory animation (3Dmol.js)

def _kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rotation that best maps centred P onto centred Q (both (n,3))."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def _element(name: str, resname: str) -> str:
    if resname in IONS:
        return IONS[resname]
    letters = re.sub(r"[^A-Za-z]", "", name)
    return letters[:1].upper() if letters else "X"


def build_animation(top: Path, traj: Path, out: Path, title: str, nframes: int,
                    sel: str, equil_ns: float, js_mode: str, smooth: int = 1) -> None:
    import warnings
    import MDAnalysis as mda

    warnings.filterwarnings("ignore", module="MDAnalysis")
    u = mda.Universe(str(top), str(traj))
    ag = u.select_atoms(sel)
    if len(ag) == 0:
        raise SystemExit(f"selection '{sel}' matched no atoms in {top}")
    fit = ag.select_atoms("name CA") or ag.select_atoms("name P") or ag
    fit_idx = np.searchsorted(ag.indices, fit.indices)
    nt = len(u.trajectory)
    log(f"animating {top.parent}: {len(ag)} atoms, {nt} frames in trajectory")
    if len(ag) * nframes > 4_000_000:
        log("large animation: consider --frames 50 or --sel 'name CA P' for the wiki")

    # Read up to ~1500 evenly spaced frames: enough for RMSF, cheap to hold.
    read = np.unique(np.linspace(0, nt - 1, min(nt, max(nframes * smooth, 1500))).astype(int))
    X = np.empty((len(read), len(ag), 3), dtype=np.float32)
    times = np.empty(len(read))
    for k, ts in enumerate(u.trajectory[read]):
        X[k] = ag.positions
        times[k] = ts.time / 1000.0
    # a molecule split across the periodic box shows up as Cα–Cα "bonds" of
    # tens of Å; the movie would be garbage, so say why instead of guessing
    ca = np.flatnonzero(np.isin(ag.names, ["CA"]))
    if len(ca) > 1:
        same = np.diff(ag.resids[ca]) == 1
        gap = np.linalg.norm(np.diff(X[0, ca], axis=0), axis=1)[same]
        if len(gap) and gap.max() > 8.0:
            log(f"WARNING: Cα–Cα distance of {gap.max():.0f} Å in frame 0 -- the trajectory "
                "is not PBC-whole. Animate clean.xtc from 04_postprocess.sh.")
    # fit every frame onto frame 0 on the fit atoms; centre on the origin
    ref = X[0, fit_idx] - X[0, fit_idx].mean(0)
    rmsd = np.empty(len(read))
    for k in range(len(read)):
        c = X[k, fit_idx].mean(0)
        R = _kabsch(X[k, fit_idx] - c, ref)
        X[k] = (X[k] - c) @ R.T
        rmsd[k] = np.sqrt(((X[k, fit_idx] - ref) ** 2).sum(1).mean())

    # per-residue RMSF over the production window, stored in the B-factor column
    prod = times >= equil_ns if (times >= equil_ns).sum() >= 10 else np.ones(len(times), bool)
    mean = X[prod].mean(0)
    rmsf_atom = np.sqrt(((X[prod] - mean) ** 2).sum(2).mean(0))
    resix = ag.resindices
    bfac = np.zeros(len(ag))
    for ri in np.unique(resix):
        m = resix == ri
        bfac[m] = rmsf_atom[m].mean()

    # frames for the movie, optionally window-averaged to remove thermal jitter
    pick = np.linspace(0, len(read) - 1, min(nframes, len(read))).astype(int)
    half = smooth // 2
    movie = np.array([X[max(0, i - half):i + half + 1].mean(0) for i in pick]) if smooth > 1 \
        else X[pick]
    lim = float(np.abs(movie).max())
    scale = 0.01 if lim < 320 else lim / 32000.0
    q = np.round(movie / scale).astype("<i2")

    # chain IDs: from the topology when it has them (PDB/TPR), otherwise a new
    # chain wherever residue numbering restarts or protein turns into nucleic.
    chains = []
    has_chain = hasattr(ag.atoms, "chainIDs") and any(c.strip() for c in ag.atoms.chainIDs)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    cur, prev_resid, prev_kind, prev_res = 0, None, None, None
    for at in ag:
        if has_chain:
            chains.append(at.chainID.strip()[:1] or "A")
            continue
        kind = "n" if at.resname in NUC_MAP else ("i" if at.resname in IONS else "p")
        if at.resindex != prev_res:
            if prev_resid is not None and (at.resid <= prev_resid or kind != prev_kind):
                cur += 1
            prev_resid, prev_kind, prev_res = at.resid, kind, at.resindex
        chains.append(letters[cur % len(letters)])

    resn = [NUC_MAP.get(r, r)[:3] for r in ag.resnames]
    data = {
        "title": title, "n": len(ag), "nf": len(pick), "scale": scale,
        "name": [str(n)[:4] for n in ag.names], "resn": resn,
        "resi": [int(r) for r in ag.resids], "chain": chains,
        "elem": [_element(n, r) for n, r in zip(ag.names, ag.resnames)],
        "het": [r in IONS for r in ag.resnames],
        "b": [round(float(v), 2) for v in bfac],
        "bmin": round(float(np.percentile(bfac, 2)), 2),
        "bmax": round(float(np.percentile(bfac, 98)), 2),
        "t": [round(float(v), 2) for v in times[pick]],
        "rmsd": [round(float(v), 2) for v in rmsd[pick]],
        "xyz": base64.b64encode(q.tobytes()).decode("ascii"),
        "nucleic": any(r in NUC_MAP for r in ag.resnames),
    }
    js = JsMode(js_mode, out.parent)
    page = ANIM_HTML.replace("__TITLE__", html.escape(title)) \
                    .replace("__MOLJS__", js.mol_tag()) \
                    .replace("__DATA__", json.dumps(data, separators=(",", ":")))
    out.write_text(page, encoding="utf-8")
    log(f"  {out.name}: {len(pick)} frames, {out.stat().st_size / 1e6:.1f} MB")


ANIM_HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
__MOLJS__
<style>
:root{--bg:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--rule:#e4e3df;--btn:#ffffff;--accent:#2a78d6}
@media (prefers-color-scheme:dark){:root{--bg:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--rule:#3a3a38;--btn:#262624;--accent:#3987e5}}
*{box-sizing:border-box}
html,body{margin:0;height:100%;background:var(--bg);color:var(--ink);
  font:14px/1.4 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
.wrap{display:flex;flex-direction:column;height:100%;min-height:420px;padding:10px 12px}
h1{font-size:15px;margin:0 0 6px;font-weight:600}
#viewer{position:relative;flex:1;min-height:300px;border:1px solid var(--rule);border-radius:8px;overflow:hidden}
.bar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-top:8px}
button,select{font:inherit;color:var(--ink);background:var(--btn);border:1px solid var(--rule);
  border-radius:6px;padding:4px 10px;cursor:pointer}
button.primary{background:var(--accent);color:#fff;border-color:var(--accent);min-width:72px}
input[type=range]{flex:1;min-width:140px;accent-color:var(--accent)}
.info{font-variant-numeric:tabular-nums;color:var(--ink2);min-width:210px}
.legend{display:none;align-items:center;gap:6px;color:var(--ink2);font-size:12px}
.legend .grad{width:120px;height:10px;border-radius:3px;background:linear-gradient(90deg,#2166ac,#f7f7f7,#b2182b)}
label{color:var(--ink2);font-size:12px}
</style></head>
<body><div class="wrap">
<h1>__TITLE__</h1>
<div id="viewer"></div>
<div class="bar">
  <button id="play" class="primary">▶ Play</button>
  <input id="slider" type="range" min="0" value="0" aria-label="frame">
  <span class="info" id="info"></span>
</div>
<div class="bar">
  <label>Style <select id="style">
    <option value="cartoon">Cartoon</option>
    <option value="cartoon+sticks">Cartoon + side chains</option>
    <option value="sticks">Sticks</option>
    <option value="spheres">Spheres</option>
    <option value="trace">Backbone trace</option>
  </select></label>
  <label>Colour <select id="color">
    <option value="chain">Chain</option>
    <option value="spectrum">N → C (spectrum)</option>
    <option value="rmsf">RMSF</option>
    <option value="ss">Secondary structure</option>
  </select></label>
  <label>Speed <select id="speed">
    <option value="200">slow</option><option value="80" selected>normal</option><option value="30">fast</option>
  </select></label>
  <label><input type="checkbox" id="loop" checked> loop</label>
  <label><input type="checkbox" id="rock"> rock</label>
  <button id="png">Save PNG</button>
  <span class="legend" id="legend">RMSF <span id="bmin"></span><span class="grad"></span><span id="bmax"></span> Å</span>
</div>
</div>
<script>
const D = __DATA__;
(function () {
  // ---- decode int16 coordinates and build a multi-model PDB -----------------
  const bin = atob(D.xyz), buf = new ArrayBuffer(bin.length), u8 = new Uint8Array(buf);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  const q = new Int16Array(buf);
  const f3 = v => (v * D.scale).toFixed(3).padStart(8);
  const head = [];
  for (let i = 0; i < D.n; i++) {
    const nm = D.name[i], el = D.elem[i];
    const name = (nm.length >= 4 || el.length === 2) ? nm.padEnd(4) : (" " + nm).padEnd(4);
    head.push((D.het[i] ? "HETATM" : "ATOM  ") + String((i + 1) % 100000).padStart(5) + " " +
      name + " " + D.resn[i].padStart(3) + " " + D.chain[i] + String(D.resi[i] % 10000).padStart(4) + "    ");
  }
  const tail = D.b.map((b, i) => "  1.00" + b.toFixed(2).padStart(6) + "          " + D.elem[i].toUpperCase().padStart(2));
  const parts = [];
  for (let f = 0; f < D.nf; f++) {
    parts.push("MODEL     " + String(f + 1).padStart(4));
    const o = f * D.n * 3;
    for (let i = 0; i < D.n; i++) {
      const k = o + 3 * i;
      parts.push(head[i] + f3(q[k]) + f3(q[k + 1]) + f3(q[k + 2]) + tail[i]);
    }
    parts.push("ENDMDL");
  }
  const pdb = parts.join("\n");

  // ---- viewer ----------------------------------------------------------------
  const dark = matchMedia("(prefers-color-scheme: dark)").matches;
  const viewer = $3Dmol.createViewer("viewer", {backgroundColor: dark ? "#1a1a19" : "#fcfcfb"});
  viewer.addModelsAsFrames(pdb, "pdb");
  const $ = id => document.getElementById(id);
  const slider = $("slider"); slider.max = D.nf - 1;
  $("bmin").textContent = D.bmin.toFixed(1); $("bmax").textContent = D.bmax.toFixed(1);

  function colorSpec() {
    switch ($("color").value) {
      case "spectrum": return {color: "spectrum"};
      case "rmsf": return {colorscheme: {prop: "b", gradient: "rwb", min: D.bmax, max: D.bmin}};
      case "ss": return {colorscheme: "ssJmol"};
      default: return {colorscheme: "chainHetatm"};
    }
  }
  function applyStyle() {
    const c = colorSpec(), s = $("style").value;
    viewer.setStyle({}, {});
    if (s === "cartoon" || s === "cartoon+sticks") {
      viewer.setStyle({}, {cartoon: Object.assign({arrows: true}, c)});
      if (D.nucleic) viewer.addStyle({resn: ["A", "G", "C", "U"]}, {stick: Object.assign({radius: 0.15}, c)});
      if (s === "cartoon+sticks")
        viewer.addStyle({not: {atom: ["N", "C", "O"]}, hetflag: false},
                        {stick: Object.assign({radius: 0.12}, c)});
    } else if (s === "sticks") viewer.setStyle({}, {stick: Object.assign({radius: 0.18}, c)});
    else if (s === "spheres") viewer.setStyle({}, {sphere: c});
    else viewer.setStyle({atom: ["CA", "P"]}, {sphere: Object.assign({radius: 0.5}, c)});
    viewer.setStyle({hetflag: true}, {sphere: {radius: 1.0, color: "#1baf7a"}});  // ions
    $("legend").style.display = $("color").value === "rmsf" ? "inline-flex" : "none";
    viewer.render();
  }

  let frame = 0, timer = null;
  async function show(f) {
    frame = f; slider.value = f;
    await viewer.setFrame(f);
    viewer.render();
    $("info").textContent = "t = " + D.t[f].toFixed(1) + " ns · RMSD " + D.rmsd[f].toFixed(1) +
      " Å · frame " + (f + 1) + "/" + D.nf;
  }
  function step() {
    let f = frame + 1;
    if (f >= D.nf) { if (!$("loop").checked) { stop(); return; } f = 0; }
    show(f).then(() => { if (timer !== null) timer = setTimeout(step, +$("speed").value); });
  }
  function play() { timer = 0; $("play").textContent = "❚❚ Pause"; step(); }
  function stop() { clearTimeout(timer); timer = null; $("play").textContent = "▶ Play"; }

  $("play").onclick = () => (timer === null ? play() : stop());
  slider.oninput = () => { stop(); show(+slider.value); };
  $("style").onchange = applyStyle;
  $("color").onchange = applyStyle;
  $("rock").onchange = e => viewer.spin(e.target.checked ? "y" : false);
  $("png").onclick = () => {
    const a = document.createElement("a");
    a.href = viewer.pngURI(); a.download = "frame_" + (frame + 1) + ".png"; a.click();
  };
  applyStyle(); viewer.zoomTo(); show(0);
})();
</script>
</body></html>
"""


def cmd_animate(a: argparse.Namespace) -> None:
    need_modules("MDAnalysis")
    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    build_animation(Path(a.top), Path(a.traj), out, a.title or out.stem, a.frames, a.sel,
                    a.equil_ns, a.js, a.smooth)


# ==========================================================================

DEFAULT_SEL = "(protein or nucleic or resname MG MGA) and not name H*"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, anim=True):
        p.add_argument("--runs", default=os.environ.get("BRUSAFE_RUNS", str(REPO / "runs")),
                       help="runs directory (default: $BRUSAFE_RUNS, else $REPO/runs)")
        p.add_argument("--out", help="output directory (default: analysis_out/...)")
        p.add_argument("--equil-ns", type=float, default=50.0,
                       help="relaxation window excluded from statistics (default 50)")
        p.add_argument("--js", choices=["local", "inline", "cdn"], default="local",
                       help="local: copy .js next to the html (iGEM wiki); inline: "
                            "self-contained files; cdn: quick preview only")
        if anim:
            p.add_argument("--frames", type=int, default=100,
                           help="frames in each animation (default 100)")
            p.add_argument("--smooth", type=int, default=1,
                           help="average each movie frame over N neighbours to damp "
                                "thermal jitter (default 1 = off)")
            p.add_argument("--sel", default=DEFAULT_SEL,
                           help=f"MDAnalysis selection to animate (default: {DEFAULT_SEL})")

    p = sub.add_parser("report", help="all figures + animations for one system")
    p.add_argument("system")
    p.add_argument("--reps", type=int, nargs="+", help="replica numbers (default: all)")
    p.add_argument("--no-anim", action="store_true", help="skip trajectory animations")
    common(p)
    p.set_defaults(fn=cmd_report)

    p = sub.add_parser("compare", help="overlay two or more systems")
    p.add_argument("systems", nargs="+")
    common(p, anim=False)
    p.set_defaults(fn=cmd_compare)

    p = sub.add_parser("animate", help="one standalone trajectory viewer")
    p.add_argument("--top", required=True, help="topology matching the trajectory "
                                                "(clean.gro, or a .pdb with chain IDs)")
    p.add_argument("--traj", required=True, help="trajectory (clean.xtc)")
    p.add_argument("-o", "--output", required=True, help="output .html")
    p.add_argument("--title")
    p.add_argument("--equil-ns", type=float, default=50.0)
    p.add_argument("--js", choices=["local", "inline", "cdn"], default="inline")
    p.add_argument("--frames", type=int, default=100)
    p.add_argument("--smooth", type=int, default=1)
    p.add_argument("--sel", default=DEFAULT_SEL)
    p.set_defaults(fn=cmd_animate)

    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
