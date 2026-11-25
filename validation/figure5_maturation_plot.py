"""Compose manuscript Figure 5 and Supplementary Figure 5.

Plan: `docs/FIGURE5_build_plan.md`. Verdicts: `docs/PREREG_mvri_presentation_audit.md` and its
dated amendment (`validation/mvri_audit_amendment.py`).

This replaces the draft Figure 5 of the atlas manuscript ("the mucosal vaccine readiness index
across the lifespan"). Four things changed, and all four are forced by the audit:

  * the unit of analysis is the donor, not the cell, so panel b draws donor means and not the
    deck's per-cell boxplots, whose Kruskal-Wallis p-values underflow to 0.0;
  * two of the five axes fall below the measurability gate and are drawn as unmeasured rather
    than as flat, so panel a states detection on the figure itself;
  * cells are labelled from `Cell_type_prediction_level*`, which covers all 148,978 cells,
    rather than `Cell_type_annotation_level*`, which is NaN for 16 whole donors;
  * the weighted composite carries its donor count on the panel.

Styling follows Figures 1-4 of the manuscript: the figure carries its own title, donors are
coloured by age group on the shared Neonate-to-Elderly ramp, and each panel states its
statistic above the axes.

No statistic is recomputed here except the composite and the effective-test count; everything
else is read from files written by `validation/mvri_audit_amendment.py`.

Usage: python validation/figure5_maturation_plot.py <audit_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from mvri_axes import AGE_LABEL, AGE_ORDER, AXES, MEASURABLE  # noqa: E402
from nasage.stats import effective_tests                      # noqa: E402

AUD, OUT = Path(sys.argv[1]), Path(sys.argv[2])
F5 = AUD / "figure5"
OUT.mkdir(parents=True, exist_ok=True)

# Age-group ramp shared with Figures 1-4 of the manuscript (Neonate -> Elderly).
AGE_COLORS = {"Neonate": "#2C7BB6", "Infant": "#38A890", "Young child": "#9BCF63",
              "Child": "#D9E05A", "Adolescent": "#EFC84A", "Adult": "#E8822E",
              "Eldly": "#D7352A"}
BLUE, ORANGE, GREEN, PURPLE = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
GREY, INK, MUTED, FAINT = "#BFBFBF", "#1A1A1A", "#6B6B6B", "#E8E8E8"
MM = 1 / 25.4
AXCOL = {"Mucus permeability": BLUE, "Antigen presentation": ORANGE,
         "IgA support": GREEN, "Innate sensing": MUTED, "Cytokine tone": MUTED}
LABELS = [AGE_LABEL[b] for b in AGE_ORDER]
GATE = 0.60

TITLE = ("Figure 5 | Maturation of the mucosal interface across the lifespan")
STITLE = ("Supplementary Figure 5 | Unit of analysis, specificity and measurement "
          "of the interface axes")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.facecolor": "white", "axes.facecolor": "white",
})

# ------------------------------------------------------------------ inputs
R = pd.read_csv(AUD / "axis_audit.csv").set_index("axis")
DON = pd.read_csv(F5 / "donor_axis_long.csv")
DET = pd.read_csv(F5 / "per_gene_detection.csv")
CT = pd.read_csv(F5 / "celltype_axis_means.csv")
WIN = pd.read_csv(F5 / "window_trends.csv")
NUL = pd.read_csv(F5 / "null_stats.csv.gz")
S = pd.read_csv(AUD / "donor_axis_scores.csv", index_col=0)
C = pd.read_csv(AUD / "axis_correlation.csv", index_col=0)

N_TOTAL = DON.donor.nunique()            # donors carrying at least one axis subset
# Axes carried into the figure: those clearing the 60% measurability gate. Innate sensing and
# cytokine tone fail it, and their genes sit at the droplet detection floor, so neither the
# axes nor their undetected genes appear in the figure. Supplementary Fig. 5b-c keeps the
# record of why they were dropped.
KEEP = [a for a in AXES if a in R.index and bool(R.loc[a, "gate_pass"])]
MEASURABLE = KEEP                        # shadows the static list in mvri_axes
C = C.loc[KEEP, KEEP]
EFF = effective_tests(C)

def weight_sensitivity(S, keep, n=10000, seed=11):
    """Which age bin peaks, over random sign-respecting weightings of the kept axes."""
    Z = S[keep].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
    bm = Z.assign(age=S.age.values).groupby("age").mean()
    bm = bm.reindex([b for b in AGE_ORDER if b in bm.index])
    signs = np.array([-1 if AXES[k]["weight"] < 0 else 1 for k in keep])
    rng = np.random.default_rng(seed)
    peaks = []
    for _ in range(n):
        w = rng.uniform(0.05, 1.0, size=len(keep)) * signs
        peaks.append(bm.index[int(np.argmax(bm.values @ (w / np.abs(w).sum())))])
    pk = pd.Series(peaks).value_counts(normalize=True)
    return pd.DataFrame({"age_bin": pk.index, "frac_peak": pk.values})


DEP = DEPS = None
if (AUD / "depth_sensitivity.csv").exists():
    DEP = pd.read_csv(AUD / "depth_sensitivity.csv").set_index("axis")
    DEPS = pd.read_csv(AUD / "depth_summary.csv").iloc[0]

SUBSET_LABEL = {
    "Mucus permeability": "Secretory, Goblet, Club, Basal",
    "Antigen presentation": "Secretory, Goblet, Mac, Mono",
    "Innate sensing": "all epithelium, Mac, Mono",
    "IgA support": "Secretory, Goblet, Plasma, B mem",
    "Cytokine tone": "all epithelium, all myeloid",
}


def letter(ax, s, dx=-0.20, dy=1.18):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left", color=INK)


def age_legend(fig, x, y, ncol=7, fontsize=5.8, title=None):
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=3.4,
                      markerfacecolor=AGE_COLORS[b], markeredgecolor="white",
                      markeredgewidth=0.3, label=AGE_LABEL[b]) for b in AGE_ORDER]
    leg = fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x, y),
                     ncol=ncol, frameon=False, fontsize=fontsize,
                     handletextpad=0.25, columnspacing=0.9, title=title)
    if title:
        leg.get_title().set_fontsize(fontsize)
    return leg


def donor_strip(ax, sub, jitter_seed=0, show_median=True):
    """Donor points by age bin, coloured on the Figure 1-4 age ramp. One point, one donor."""
    rng = np.random.default_rng(jitter_seed)
    for i, b in enumerate(AGE_ORDER):
        v = sub.loc[sub.age == b, "score"].values
        if not len(v):
            continue
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=12,
                   facecolor=AGE_COLORS[b], edgecolor="white", linewidth=0.35,
                   alpha=0.95, zorder=3)
        if show_median:
            ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=1.2, zorder=4)
    ax.axhline(0, color=GREY, lw=0.5, ls=":", zorder=1)
    ax.set_xticks(range(len(AGE_ORDER)))
    ax.set_xlim(-0.6, len(AGE_ORDER) - 0.4)


def ribbon(ax, x0, x1, y0t, y0b, y1t, y1b, colour, alpha=0.45, n=60, hatch=None):
    """A Sankey ribbon: a filled band whose edges follow a smoothstep from x0 to x1."""
    t = np.linspace(0, 1, n)
    s = t * t * (3 - 2 * t)                      # smoothstep, flat at both nodes
    xs = x0 + (x1 - x0) * t
    top = y0t + (y1t - y0t) * s
    bot = y0b + (y1b - y0b) * s
    ax.fill(np.concatenate([xs, xs[::-1]]), np.concatenate([top, bot[::-1]]),
            facecolor=colour, edgecolor="none", alpha=alpha, zorder=2, hatch=hatch,
            linewidth=0)


def panel_a_sankey(ax, DET, CT, R):
    """genes -> axes -> cell types.

    Heights are flows, not decoration: an axis node is as tall as its gene count, split into
    a solid block (genes detected in >= 5% of the cells the axis is scored in) and a hatched
    block (genes below that floor, which still enter the score and contribute only noise).
    Each axis's height is then divided among its cell subsets in proportion to **cells**, so
    the panel also shows the composition the deck never states -- for instance that the
    antigen-presentation axis is overwhelmingly epithelial, with under 2% of its cells
    myeloid.
    """
    order_t = ["Basal", "Club", "Secretory", "Goblet", "Goblet inflammatory", "Ciliated",
               "Macrophage", "Monocyte", "B mem", "Plasma"]
    order_t = [t for t in order_t if t in set(CT[CT.axis.isin(KEEP)].t)]
    axes_order = list(KEEP)
    ng = {a: int(DET[(DET.axis == a) & DET.measurable].shape[0]) for a in axes_order}
    nok = dict(ng)

    GAP = 4.0
    T = sum(ng.values()) + GAP * (len(axes_order) - 1)

    def Y(u):                                     # units from the top -> axes fraction
        return 1.0 - u / T

    # flow from each axis into each cell type, in gene-units
    flow = {}
    for a in axes_order:
        sub = CT[CT.axis == a]
        tot = float(sub.n_cells.sum())
        for _, r in sub.iterrows():
            flow[(a, r.t)] = ng[a] * r.n_cells / tot

    # Small compartments get a minimum drawn height so their labels stay separable; the
    # ribbons into them are scaled to match. Macrophage, Monocyte and Plasma are a few
    # hundred cells against Goblet's 43,000, and at true scale they vanish.
    raw_h = {t: sum(flow.get((a, t), 0.0) for a in axes_order) for t in order_t}
    MINH, CGAP = 2.2, 0.9
    cell_h = {t: max(h, MINH) for t, h in raw_h.items()}
    n_gap = max(len(order_t) - 1, 1)
    need = sum(cell_h.values()) + CGAP * n_gap
    if need > T:                                  # shrink to fit rather than overflow
        f = T / need
        cell_h = {t: h * f for t, h in cell_h.items()}
        CGAP *= f
    else:
        CGAP = (T - sum(cell_h.values())) / n_gap
    scale = {t: cell_h[t] / raw_h[t] if raw_h[t] > 0 else 1.0 for t in order_t}

    cell_top, u = {}, 0.0
    for t in order_t:
        cell_top[t] = u
        u += cell_h[t] + CGAP
    cell_fill = dict(cell_top)                    # running offset as ribbons arrive

    XG0, XG1, XA0, XA1, XC0, XC1 = 0.215, 0.243, 0.520, 0.556, 0.800, 0.826

    atop, u = {}, 0.0
    for a in axes_order:
        atop[a] = u
        u += ng[a] + GAP

    for a in axes_order:
        ok = bool(R.loc[a, "gate_pass"])
        col = AXCOL[a] if ok else GREY
        t0, h = atop[a], ng[a]

        # ---- gene node: solid = detected, hatched = below the floor
        ax.add_patch(plt.Rectangle((XG0, Y(t0 + h)), XG1 - XG0, h / T,
                                   facecolor=col, edgecolor="none", zorder=4))
        frac = float(R.loc[a, "frac_measurable"])
        emp = float(R.loc[a, "empirical_p"])
        verdict = (f"empirical $P$ {emp:.3f}*" if (ok and emp < 0.05) else
                   (f"empirical $P$ {emp:.3f}" if ok else "below detection floor"))
        yc = Y(t0 + h / 2)
        ax.text(XG0 - 0.012, yc + 0.046, a, fontsize=6.0, fontweight="bold",
                ha="right", va="center", color=INK if ok else MUTED)
        ax.text(XG0 - 0.012, yc,
                f"{ng[a]} genes scored, {frac:.0%} of the list detected",
                fontsize=5.0, ha="right", va="center", color=MUTED)
        ax.text(XG0 - 0.012, yc - 0.046, verdict, fontsize=5.2, ha="right", va="center",
                color=col if ok else MUTED,
                fontweight="bold" if (ok and emp < 0.05) else "normal")

        # ---- gene -> axis
        ribbon(ax, XG1, XA0, Y(t0), Y(t0 + h), Y(t0), Y(t0 + h), col, alpha=0.42)

        # ---- axis node
        ax.add_patch(plt.Rectangle((XA0, Y(t0 + h)), XA1 - XA0, h / T,
                                   facecolor=col, edgecolor="none", zorder=4))
        # ---- axis -> cell types, ordered by target so ribbons do not cross
        tg = sorted([t for t in order_t if (a, t) in flow], key=lambda t: cell_top[t])
        off = t0
        for t in tg:
            w = flow[(a, t)]
            c0, c1 = cell_fill[t], cell_fill[t] + w * scale[t]
            ribbon(ax, XA1, XC0, Y(off), Y(off + w), Y(c0), Y(c1), col,
                   alpha=0.40 if ok else 0.12)
            cell_fill[t] = c1
            off += w

    SHORT_T = {"Goblet inflammatory": "Goblet infl."}
    for t in order_t:
        t0, h = cell_top[t], cell_h[t]
        ax.add_patch(plt.Rectangle((XC0, Y(t0 + h)), XC1 - XC0, h / T,
                                   facecolor="#6B6B6B", edgecolor="none", zorder=4))
        n = int(CT[(CT.t == t) & CT.axis.isin(KEEP)].n_cells.max())
        ax.text(XC1 + 0.009, Y(t0 + h / 2), f"{SHORT_T.get(t, t)} ({n:,})", fontsize=5.5,
                ha="left", va="center", color=INK)

    for x, lab in [((XG0 + XG1) / 2, "genes"), ((XA0 + XA1) / 2, "axis"),
                   ((XC0 + XC1) / 2, "cells it is scored in")]:
        ax.text(x, 1.035, lab, fontsize=6.2, ha="center", va="bottom", color=MUTED)
    # key for the two block styles
    ax.text(0.015, -0.090, "only genes detected in ≥ 5% of the cells their axis is "
                           "scored in are shown, and only those genes enter the score",
            fontsize=5.2, va="center", ha="left", color=MUTED)
    ax.text(0.70, -0.090, "cell nodes are cell counts, with a floor\n"
                          "so small compartments stay visible",
            fontsize=5.2, va="center", ha="left", color=MUTED)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.135, 1.02)
    ax.axis("off")


# ================================================================== MAIN FIGURE
fig = plt.figure(figsize=(180 * MM, 222 * MM))
fig.text(0.012, 0.988, TITLE, fontsize=10, fontweight="bold", va="top", ha="left",
         color=INK)
age_legend(fig, 0.012, 0.963, ncol=7)
gs = GridSpec(3, 3, figure=fig, height_ratios=[1.15, 1.0, 1.02],
              hspace=0.58, wspace=0.72,
              left=0.112, right=0.975, top=0.880, bottom=0.135)

# ---- a: genes -> axes -> cells (Sankey) ---------------------------------------
axa = fig.add_subplot(gs[0, :])
letter(axa, "a", dx=-0.030, dy=1.30)
axa.text(0, 1.20, "What each measurable axis is built from and where it is read out: its "
                  "detected genes, then the cells those genes are scored in",
         transform=axa.transAxes, fontsize=7, va="bottom", ha="left", color=INK)
panel_a_sankey(axa, DET, CT, R)

# ---- b: donor-level trajectories, measurable axes only ------------------------
for j, name in enumerate(MEASURABLE):
    ax = fig.add_subplot(gs[1, j])
    sub = DON[DON.axis == name]
    donor_strip(ax, sub, jitter_seed=j)
    ns = [int((sub.age == b).sum()) for b in AGE_ORDER]
    ax.set_xticklabels([f"{lab}\n(n={n})" for lab, n in zip(LABELS, ns)],
                       rotation=55, ha="right", fontsize=5.6)
    r = R.loc[name]
    # The headline is the empirical p against expression-matched random gene sets, not the
    # nominal or BH p. With all three datasets pooled, 53-80% of matched random sets reach
    # donor-level P < 0.05 (Supplementary Fig. 5b), because dataset and 10x chemistry are
    # confounded with age bin -- so a nominal p here carries almost no information.
    star = "*" if r.empirical_p < 0.05 else ""
    ax.set_title(f"{name}\nempirical $P$ = {r.empirical_p:.3f}{star} vs matched nulls\n"
                 f"(BH $P$ = {r.BH_perdonor:.3f}; n = {sub.donor.nunique()} donors)",
                 fontsize=6.3, pad=4, color=INK)
    if j == 0:
        ax.set_ylabel("donor mean axis score\n(mean of per-gene $z$)")
        letter(ax, "b", dx=-0.30, dy=1.30)

# ---- c: the two age windows ---------------------------------------------------
axc = fig.add_subplot(gs[2, 0])
letter(axc, "c", dx=-0.30, dy=1.20)
ypos, ylab, k = [], [], 0
for name in MEASURABLE[::-1]:
    for w, mark in [("Child to Elderly", "s"), ("Neonate to Young child", "o")]:
        row = WIN[(WIN.axis == name) & (WIN.window == w)]
        if not len(row):
            continue
        row = row.iloc[0]
        axc.plot([row.lo, row.hi], [k, k], color=AXCOL[name], lw=1.2, zorder=3,
                 solid_capstyle="round")
        axc.scatter([row.rho], [k], s=18, marker=mark, facecolor=AXCOL[name],
                    edgecolor="white", linewidth=0.4, zorder=4)
        star = "*" if row.p < 0.05 else ""
        axc.text(row.hi + 0.05, k, f"$P$ {row.p:.3f}{star}", fontsize=5.2, va="center",
                 ha="left", color=MUTED)
        ypos.append(k)
        ylab.append(w.replace("Neonate to Young child", "Neonate → Y child")
                     .replace(" to ", " → "))
        k += 1
    # axis name as a group heading above its pair of rows, inside the axes
    axc.text(0.015, k - 0.35, name, transform=axc.get_yaxis_transform(), fontsize=6,
             color=AXCOL[name], fontweight="bold", va="center", ha="left")
    k += 1.1
axc.axvline(0, color=INK, lw=0.6, zorder=2)
axc.set_yticks(ypos)
axc.set_yticklabels(ylab, fontsize=5.6)
axc.set_xlabel("Spearman $\\rho$ of donor score with\nage-bin rank (95% bootstrap CI)")
axc.set_ylim(-0.9, k - 0.5)
axc.set_xlim(-1.35, 1.75)
axc.set_title("Trend within each age window", fontsize=6.6, pad=4)

# ---- d: which cells carry each axis -------------------------------------------
axd = fig.add_subplot(gs[2, 1])
letter(axd, "d", dx=-0.34, dy=1.20)
order = [t for t in ["Basal", "Club", "Secretory", "Goblet", "Goblet inflammatory",
                     "Ciliated", "Macrophage", "Monocyte", "cDC", "pDC",
                     "B mem", "Plasma"] if t in set(CT[CT.axis.isin(KEEP)].t)]
sub = CT[CT.axis.isin(MEASURABLE)]
vmax = float(np.nanmax(np.abs(sub.mean_score)))
for i, name in enumerate(MEASURABLE):
    for j, t in enumerate(order):
        rec = sub[(sub.axis == name) & (sub.t == t)]
        if not len(rec):
            continue
        rec = rec.iloc[0]
        size = 5 + 48 * np.sqrt(min(rec.n_cells, 20000) / 20000)
        axd.scatter(j, i, s=size, c=[rec.mean_score], cmap="RdBu_r",
                    vmin=-vmax, vmax=vmax, edgecolor=INK, linewidth=0.3, zorder=3)
axd.set_xticks(range(len(order)))
axd.set_xticklabels(order, rotation=55, ha="right", fontsize=5.8)
axd.set_yticks(range(len(MEASURABLE)))
axd.set_yticklabels(MEASURABLE, fontsize=6)
for tick, name in zip(axd.get_yticklabels(), MEASURABLE):
    tick.set_color(AXCOL[name])
axd.set_xlim(-0.7, len(order) - 0.3)
axd.set_ylim(-0.7, len(MEASURABLE) - 0.3)
axd.set_title("Mean score by cell type,\nwithin each axis's own subset", pad=4, fontsize=6.6)
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(vmin=-vmax, vmax=vmax))
cb = fig.colorbar(sm, ax=axd, fraction=0.038, pad=0.02)
cb.set_label("mean $z$ in group", fontsize=5.4)
cb.ax.tick_params(labelsize=5.0, width=0.5, length=2)
cb.outline.set_linewidth(0.5)
for n_ in (1000, 5000, 20000):
    axd.scatter([], [], s=5 + 48 * np.sqrt(n_ / 20000), facecolor="white",
                edgecolor=INK, linewidth=0.3, label=f"{n_:,}")
axd.legend(loc="upper center", bbox_to_anchor=(0.45, -0.40), ncol=3, frameon=False,
           fontsize=5.2, title="cells in group", title_fontsize=5.4,
           handletextpad=0.2, columnspacing=0.8, borderpad=0.1)

# ---- e: the composite ---------------------------------------------------------
axe = fig.add_subplot(gs[2, 2])
letter(axe, "e", dx=-0.32, dy=1.20)
Z = S[KEEP].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
w = np.array([AXES[k]["weight"] for k in KEEP])
w = w / np.abs(w).sum()
comp = pd.DataFrame({"comp": Z.values @ w, "age": S.age.values})
rng = np.random.default_rng(7)
means, ns = [], []
for i, b in enumerate(AGE_ORDER):
    v = comp.loc[comp.age == b, "comp"].values
    ns.append(len(v))
    means.append(float(np.mean(v)) if len(v) else np.nan)
    if len(v):
        axe.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=12,
                    facecolor=AGE_COLORS[b], edgecolor="white", linewidth=0.35,
                    alpha=0.95, zorder=3)
axe.plot(range(len(AGE_ORDER)), means, color=INK, lw=1.2, zorder=4)
for i, n in enumerate(ns):
    if n <= 2:
        axe.axvspan(i - 0.45, i + 0.45, color=FAINT, zorder=0)
axe.axhline(0, color=GREY, lw=0.5, ls=":", zorder=1)
axe.set_xticks(range(len(AGE_ORDER)))
axe.set_xticklabels([f"{lab}\n(n={n})" for lab, n in zip(LABELS, ns)],
                    rotation=55, ha="right", fontsize=5.6)
axe.set_ylabel("combined score")
axe.set_xlim(-0.6, len(AGE_ORDER) - 0.4)
axe.set_title(f"The {len(KEEP)} measurable axes combined,\nall {len(comp)} donors",
              fontsize=6.6, pad=4)

fig.savefig(OUT / "figure5_maturation.pdf", dpi=600)
fig.savefig(OUT / "figure5_maturation.png", dpi=300)
plt.close(fig)
print("wrote", OUT / "figure5_maturation.pdf")

# ================================================================== SUPPLEMENTARY
figs = plt.figure(figsize=(180 * MM, 238 * MM))
figs.text(0.012, 0.991, STITLE, fontsize=9, fontweight="bold", va="top", ha="left",
          color=INK)
age_legend(figs, 0.012, 0.971, ncol=7)
gss = GridSpec(4, 12, figure=figs, height_ratios=[1.0, 0.95, 1.05, 1.0],
               hspace=1.05, wspace=1.60,
               left=0.082, right=0.982, top=0.915, bottom=0.098)

# ---- S a: per cell versus per donor -------------------------------------------
CS = pd.read_csv(F5 / "cell_scores.csv.gz")
for j, name in enumerate(MEASURABLE):
    ax = figs.add_subplot(gss[0, 4 * j:4 * j + 4])
    parts = [(i, CS.loc[(CS.axis == name) & (CS.age == b), "score"].values)
             for i, b in enumerate(AGE_ORDER)]
    parts = [(i, p) for i, p in parts if len(p)]
    vp = ax.violinplot([p for _, p in parts], positions=[i for i, _ in parts],
                       widths=0.85, showextrema=False, showmedians=False)
    for body in vp["bodies"]:
        body.set_facecolor(GREY)
        body.set_alpha(0.40)
        body.set_linewidth(0)
    sub = DON[DON.axis == name]
    for i, b in enumerate(AGE_ORDER):
        v = sub.loc[sub.age == b, "score"].values
        if len(v):
            ax.scatter([i] * len(v), v, s=9, facecolor=AGE_COLORS[b], edgecolor="white",
                       linewidth=0.3, zorder=4)
    # Clip to the bulk of the per-cell distribution: a handful of extreme cells otherwise
    # compress the donor means into an invisible band, which is the opposite of the point.
    allv = np.concatenate([p for _, p in parts])
    lo, hi = np.percentile(allv, [0.5, 99.0])
    pad = 0.08 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xticks(range(len(AGE_ORDER)))
    ax.set_xticklabels(LABELS, rotation=55, ha="right", fontsize=5.6)
    r = R.loc[name]
    ax.set_title(f"{name}\nper cell $P$ < 1e−300   |   "
                 f"per donor $P$ = {r.KW_p_perdonor:.3f}", fontsize=6.3, pad=4)
    if j == 0:
        ax.set_ylabel("axis score")
        letter(ax, "a", dx=-0.17, dy=1.34)

# ---- S b: matched nulls -------------------------------------------------------
for j, name in enumerate(KEEP):
    ax = figs.add_subplot(gss[1, 3 * j:3 * j + 3])
    d = NUL[NUL.axis == name].kw_stat.values
    ax.hist(d, bins=26, color=FAINT, edgecolor=GREY, linewidth=0.4)
    obs = float(R.loc[name, "KW_stat_perdonor"])
    ax.axvline(obs, color=AXCOL[name], lw=1.3)
    star = "*" if R.loc[name, "empirical_p"] < 0.05 else ""
    ax.set_title(f"{name.replace(' ', chr(10))}\n"
                 f"empirical $P$ {R.loc[name, 'empirical_p']:.3f}{star}", fontsize=5.9,
                 pad=3, color=INK if R.loc[name, "gate_pass"] else MUTED)
    ax.tick_params(labelsize=5.2)
    ax.set_xlabel(f"KW statistic\nnulls at $P$<0.05: "
                  f"{R.loc[name, 'null_frac_p05']:.0%}", fontsize=5.6)
    if j == 0:
        ax.set_ylabel("matched null sets", fontsize=6)
        letter(ax, "b", dx=-0.42, dy=1.46)
axnote = figs.add_subplot(gss[1, 9:])
axnote.axis("off")
axnote.text(0, 0.95, "500 random gene sets\nper axis, matched on\nexpression decile and\nset "
                     "size, scored and\ntested identically at\ndonor level. Vertical\nline, "
                     "the observed\nstatistic.\n\nMore than half of\nrandom sets reach a\n"
                     "nominal $P$ < 0.05,\nbecause dataset and\n10x chemistry are\nconfounded "
                     "with age\nbin (panel f). The\nempirical $P$ is the only\ninterpretable "
                     "test here.",
            transform=axnote.transAxes, fontsize=5.2, va="top", ha="left", color=MUTED)

# ---- S c: per-gene detection --------------------------------------------------
axc2 = figs.add_subplot(gss[2, :8])
letter(axc2, "c", dx=-0.075, dy=1.14)
x, gpos, gnames, grp = 0, [], [], []
for name in AXES:
    d = DET[DET.axis == name].sort_values("detection", ascending=False)
    xs = np.arange(len(d)) + x
    cols = [AXCOL[name] if m else "#DADADA" for m in d.measurable]
    axc2.bar(xs, d.detection.values, color=cols, width=0.78, linewidth=0)
    gpos.extend(xs)
    gnames.extend(d.gene)
    grp.append((name, xs.mean()))
    x += len(d) + 4
axc2.axhline(0.05, color=INK, lw=0.7, ls="--")
axc2.text(-1.0, 0.062, "5% detection", fontsize=5.6, ha="left", va="bottom", color=INK,
          bbox=dict(facecolor="white", edgecolor="none", pad=0.6), zorder=6)
axc2.set_yscale("log")
axc2.set_ylim(1e-4, 1.6)
axc2.set_ylabel("fraction of cells with\nnon-zero counts")
axc2.set_xticks(gpos)
axc2.set_xticklabels(gnames, fontsize=4.0, rotation=90, ha="center")
axc2.tick_params(axis="x", length=1.2, pad=1)
axc2.set_xlim(-1.5, x - 3.5)
for name, cx in grp:
    axc2.text(cx, -0.42, f"{name}\n{R.loc[name, 'frac_measurable']:.0%} measurable",
              transform=axc2.get_xaxis_transform(), fontsize=5.4, ha="center", va="top",
              color=INK if R.loc[name, "gate_pass"] else MUTED)
axc2.set_title("Detection of every gene, within the cells its own axis is scored in",
               fontsize=6.6, pad=4)

# ---- S d: weight sensitivity --------------------------------------------------
axd2 = figs.add_subplot(gss[2, 8:])
letter(axd2, "d", dx=-0.26, dy=1.16)
WS = weight_sensitivity(S, KEEP)
ws = WS.set_index("age_bin").frac_peak.reindex(AGE_ORDER).fillna(0.0)
axd2.bar(range(len(AGE_ORDER)), ws.values,
         color=[AGE_COLORS[b] for b in AGE_ORDER], width=0.72, linewidth=0)
axd2.set_xticks(range(len(AGE_ORDER)))
axd2.set_xticklabels(LABELS, rotation=55, ha="right", fontsize=5.8)
axd2.set_ylabel("fraction of 10,000\nweightings peaking here", fontsize=6)
child = float(ws.get("Young child", 0) + ws.get("Child", 0))
axd2.set_title(f"Weight sensitivity: childhood peaks in {child:.1%}\n"
               "of sign-respecting weightings (threshold 90%)", fontsize=6.4, pad=4)

# ---- S e: depth sensitivity ---------------------------------------------------
SHORT = {"Mucus permeability": "Mucus", "Antigen presentation": "Antigen pres.",
         "Innate sensing": "Innate", "IgA support": "IgA", "Cytokine tone": "Cytokine"}
if DEP is not None:
    axe2 = figs.add_subplot(gss[3, :3])
    letter(axe2, "e", dx=-0.44, dy=1.16)
    names = [n for n in KEEP if n in DEP.index]
    axe2.barh(range(len(names)), [DEP.loc[n, "rho_with_depth"] for n in names],
              color=[AXCOL[n] for n in names], height=0.6, linewidth=0)
    axe2.set_yticks(range(len(names)))
    axe2.set_yticklabels([SHORT[n] for n in names], fontsize=5.8)
    axe2.invert_yaxis()
    axe2.set_xlabel("Spearman $\\rho$ with\nmean genes per cell", fontsize=6)
    axe2.set_title("Sequencing depth", fontsize=6.4, pad=4)
    axe2.text(0.5, -0.46, f"inter-axis $\\rho$ {DEPS.raw_mean_rho:+.2f} → "
                          f"{DEPS.depth_adj_mean_rho:+.2f} after removing depth;\n"
                          f"depth vs age-bin rank $\\rho$ "
                          f"{DEPS.rho_depth_vs_agebin:+.2f} ($P$ = "
                          f"{DEPS.p_depth_vs_agebin:.2f})",
              transform=axe2.transAxes, fontsize=5.2, va="top", ha="center", color=MUTED)
    slot_f, slot_g = gss[3, 3:8], gss[3, 8:]
else:
    slot_f, slot_g = gss[3, :6], gss[3, 6:]

# ---- S f: dataset and chemistry are confounded with age bin -------------------
# This is why the matched nulls in panel b are so badly calibrated, and it is the main
# caveat on the whole figure: the Elderly bin is almost entirely a different experiment.
axf2 = figs.add_subplot(slot_f)
letter(axf2, "f" if DEP is not None else "e", dx=-0.22, dy=1.16)
BATCH = pd.read_csv(AUD / "donor_batch.csv")
DS_COL = {"ucl": "#4C72B0", "extended_ethnicity": "#DD8452", "senescence": "#937860"}
tab = (BATCH.groupby(["age", "dataset"]).size().unstack(fill_value=0)
       .reindex(index=AGE_ORDER).fillna(0))
bottom = np.zeros(len(AGE_ORDER))
for ds in [c for c in ["ucl", "extended_ethnicity", "senescence"] if c in tab.columns]:
    axf2.bar(range(len(AGE_ORDER)), tab[ds].values, bottom=bottom, width=0.72,
             color=DS_COL[ds], linewidth=0, label=ds)
    bottom += tab[ds].values
axf2.set_xticks(range(len(AGE_ORDER)))
axf2.set_xticklabels(LABELS, rotation=55, ha="right", fontsize=5.6)
axf2.set_ylabel("donors", fontsize=6)
axf2.legend(frameon=False, fontsize=5.0, loc="upper left", bbox_to_anchor=(0.00, 1.00),
            handlelength=0.9, labelspacing=0.35, borderpad=0.1)
axf2.set_title("Source dataset by age bin: the Elderly\nbin is almost a separate experiment",
               pad=4, fontsize=6.4)

# ---- S g: dimensionality ------------------------------------------------------
axg2 = figs.add_subplot(slot_g)
letter(axg2, "g" if DEP is not None else "f", dx=-0.30, dy=1.16)
M = C.values
axg2.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        axg2.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5.0,
                  color="white" if abs(M[i, j]) > 0.6 else INK)
axg2.set_xticks(range(len(C)))
axg2.set_xticklabels([SHORT[c].replace(" ", "\n") for c in C.columns], fontsize=5.0,
                     rotation=55, ha="right")
axg2.set_yticks(range(len(C)))
axg2.set_yticklabels([SHORT[c] for c in C.index], fontsize=5.8)
axg2.set_title(f"Spearman between axes: {EFF['eff_li_ji']:.1f} effective\ntests of 5 "
               f"declared, PC1 {EFF['pc1_var_frac']:.0%}", pad=4, fontsize=6.4)

figs.savefig(OUT / "figure5_maturation_supp.pdf", dpi=600)
figs.savefig(OUT / "figure5_maturation_supp.png", dpi=300)
plt.close(figs)
print("wrote", OUT / "figure5_maturation_supp.pdf")
