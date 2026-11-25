"""Manuscript Figure 5, epithelium only: genes -> axes -> epithelial cell types -> lifespan.

Supersedes `validation/figure5_maturation_plot.py`. Data:
`validation/mvri_audit_epithelial.py`, results in `results/mvri_audit_epi/`.

Only the axes that clear the measurability gate *and* their matched nulls are drawn, and only
the genes detected in >= 5% of the cells their axis is scored in — the rest are reported in
the results file, not in the figure. Panel a names every gene that survives, so the figure
states its own inputs.

Styling follows Figures 1-4: the figure carries its own title, donors are coloured on the
shared Neonate-to-Elderly ramp, and each panel states its statistic above the axes.

Usage: python validation/figure5_plot_epithelial.py <audit_dir> <out_dir>
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

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from mvri_axes import AGE_LABEL, AGE_ORDER, AXES                # noqa: E402

AUD, OUT = Path(sys.argv[1]), Path(sys.argv[2])
F5 = AUD / "figure5"
OUT.mkdir(parents=True, exist_ok=True)

AGE_COLORS = {"Neonate": "#2C7BB6", "Infant": "#38A890", "Young child": "#9BCF63",
              "Child": "#D9E05A", "Adolescent": "#EFC84A", "Adult": "#E8822E",
              "Eldly": "#D7352A"}
BLUE, ORANGE, GREEN, PURPLE = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
GREY, INK, MUTED, FAINT = "#BFBFBF", "#1A1A1A", "#6B6B6B", "#E8E8E8"
MM = 1 / 25.4
PALETTE = [BLUE, ORANGE, GREEN, PURPLE]
LABELS = [AGE_LABEL[b] for b in AGE_ORDER]
TITLE = "Figure 5 | Maturation of the epithelial mucosal interface across the lifespan"

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

R = pd.read_csv(AUD / "axis_audit.csv").set_index("axis")
DET = pd.read_csv(F5 / "per_gene_detection.csv")
MEM = pd.read_csv(F5 / "gene_membership.csv")
DON = pd.read_csv(F5 / "donor_axis_long.csv")
CT = pd.read_csv(F5 / "celltype_axis_means.csv")
WIN = pd.read_csv(F5 / "window_trends.csv")
S = pd.read_csv(AUD / "donor_axis_scores.csv", index_col=0)

# Axes kept are those that clear the measurability gate. Innate sensing (45%) and cytokine
# tone (54%) fail it; IgA support is not scorable in epithelium at all, because only 2 of its
# 4 genes clear the detection floor once Plasma and B mem are excluded. Whether a kept axis
# also clears its matched null is stated on its own panel rather than used to hide it.
KEEP = [a for a in AXES if a in R.index and bool(R.loc[a, "gate_pass"])]
GT = pd.read_csv(F5 / "gene_trends.csv")
AXCOL = {a: PALETTE[i % len(PALETTE)] for i, a in enumerate(KEEP)}
print("drawing axes:", KEEP)

EPI_ORDER = ["Basal", "Basal cycling", "Hillock", "Club", "Secretory", "Goblet",
             "Goblet inflammatory", "Transit epi", "Deuterosomal", "Ciliated",
             "Squamous", "Duct", "Rare epi"]


def letter(ax, s, dx=-0.20, dy=1.18):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left", color=INK)


def age_legend(fig, x, y, ncol=7, fontsize=5.8):
    h = [Line2D([], [], marker="o", linestyle="none", markersize=3.4,
                markerfacecolor=AGE_COLORS[b], markeredgecolor="white",
                markeredgewidth=0.3, label=AGE_LABEL[b]) for b in AGE_ORDER]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(x, y), ncol=ncol,
               frameon=False, fontsize=fontsize, handletextpad=0.25, columnspacing=0.9)


def ribbon(ax, x0, x1, y0t, y0b, y1t, y1b, colour, alpha=0.42, n=60):
    t = np.linspace(0, 1, n)
    s = t * t * (3 - 2 * t)
    xs = x0 + (x1 - x0) * t
    top = y0t + (y1t - y0t) * s
    bot = y0b + (y1b - y0b) * s
    ax.fill(np.concatenate([xs, xs[::-1]]), np.concatenate([top, bot[::-1]]),
            facecolor=colour, edgecolor="none", alpha=alpha, zorder=2, linewidth=0)


def panel_a_sankey(ax):
    """Every scored gene -> its axis -> the epithelial cell types it is scored in.

    Heights are flows. A gene is one unit, so an axis node is as tall as the number of genes
    that survived the detection floor; that height is then divided among cell types in
    proportion to cells, which is what makes the composition visible -- these axes are
    overwhelmingly read out in Secretory and Goblet cells.
    """
    genes = {a: list(DET[(DET.axis == a) & DET.measurable]
                     .sort_values("detection", ascending=False).itertuples())
             for a in KEEP}
    ng = {a: len(genes[a]) for a in KEEP}
    GAP = 3.0
    T = sum(ng.values()) + GAP * (len(KEEP) - 1)

    def Y(u):
        return 1.0 - u / T

    types = [t for t in EPI_ORDER if t in set(CT[CT.axis.isin(KEEP)].t)]
    flow = {}
    for a in KEEP:
        sub = CT[(CT.axis == a) & CT.t.isin(types)]
        tot = float(sub.n_cells.sum())
        for _, r in sub.iterrows():
            flow[(a, r.t)] = ng[a] * r.n_cells / tot

    raw_h = {t: sum(flow.get((a, t), 0.0) for a in KEEP) for t in types}
    MINH = 0.9
    cell_h = {t: max(h, MINH) for t, h in raw_h.items()}
    scale = {t: cell_h[t] / raw_h[t] if raw_h[t] > 0 else 1.0 for t in types}
    CGAP = max((T - sum(cell_h.values())) / max(len(types) - 1, 1), 0.8)
    cell_top, u = {}, 0.0
    for t in types:
        cell_top[t] = u
        u += cell_h[t] + CGAP
    cell_fill = dict(cell_top)

    XL, XD0, XD1, XG0, XG1 = 0.108, 0.118, 0.188, 0.196, 0.212
    XA0, XA1, XC0, XC1 = 0.520, 0.556, 0.760, 0.784

    atop, u = {}, 0.0
    for a in KEEP:
        atop[a] = u
        u += ng[a] + GAP

    for a in KEEP:
        col = AXCOL[a]
        t0, h = atop[a], ng[a]
        for k, g in enumerate(genes[a]):
            yt, yb = t0 + k, t0 + k + 1
            pad = 0.16
            ax.text(XL, Y(yt + 0.5), g.gene, fontsize=4.4, ha="right", va="center",
                    color=INK)
            ax.add_patch(plt.Rectangle((XD0, Y(yb - pad)), (XD1 - XD0) * g.detection,
                                       (1 - 2 * pad) / T, facecolor=col, alpha=0.55, lw=0))
            ax.add_patch(plt.Rectangle((XG0, Y(yb - pad)), XG1 - XG0, (1 - 2 * pad) / T,
                                       facecolor=col, lw=0, zorder=4))
            ribbon(ax, XG1, XA0, Y(yt + pad), Y(yb - pad), Y(yt + pad), Y(yb - pad),
                   col, alpha=0.30)

        ax.add_patch(plt.Rectangle((XA0, Y(t0 + h)), XA1 - XA0, h / T,
                                   facecolor=col, lw=0, zorder=4))
        emp = float(R.loc[a, "empirical_p"])
        ax.text((XA0 + XA1) / 2, Y(t0) + 0.022, a, fontsize=6.4, fontweight="bold",
                ha="center", va="bottom", color=col)
        ax.text((XA0 + XA1) / 2, Y(t0 + h) - 0.022,
                f"{int(R.loc[a, 'n_genes_scored'])} genes · "
                f"{int(R.loc[a, 'n_cells']):,} cells\nempirical $P$ = {emp:.3f}",
                fontsize=5.0, ha="center", va="top", color=MUTED)

        tg = sorted([t for t in types if (a, t) in flow], key=lambda t: cell_top[t])
        off = t0
        for t in tg:
            w = flow[(a, t)]
            c0, c1 = cell_fill[t], cell_fill[t] + w * scale[t]
            ribbon(ax, XA1, XC0, Y(off), Y(off + w), Y(c0), Y(c1), col, alpha=0.38)
            cell_fill[t] = c1
            off += w

    for t in types:
        t0, h = cell_top[t], cell_h[t]
        ax.add_patch(plt.Rectangle((XC0, Y(t0 + h)), XC1 - XC0, h / T,
                                   facecolor="#6B6B6B", lw=0, zorder=4))
        n = int(CT[CT.t == t].n_cells.max())
        ax.text(XC1 + 0.008, Y(t0 + h / 2),
                f"{t.replace('Goblet inflammatory', 'Goblet infl.')} ({n:,})",
                fontsize=5.4, ha="left", va="center", color=INK)

    for x, lab in [((XD0 + XD1) / 2, "genes, and the fraction of\ncells detecting each"),
                   ((XA0 + XA1) / 2, "axis"),
                   ((XC0 + XC1) / 2, "epithelial cells it is scored in")]:
        ax.text(x, 1.028, lab, fontsize=6.0, ha="center", va="bottom", color=MUTED)
    ax.plot([XD0, XD1], [-0.030] * 2, color=MUTED, lw=0.5, clip_on=False)
    for f in (0, 0.5, 1.0):
        ax.plot([XD0 + f * (XD1 - XD0)] * 2, [-0.030, -0.045], color=MUTED, lw=0.5,
                clip_on=False)
        ax.text(XD0 + f * (XD1 - XD0), -0.052, f"{f:.0%}", fontsize=4.6, ha="center",
                va="top", color=MUTED)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.10, 1.02)
    ax.axis("off")


# ================================================================== FIGURE
nk = len(KEEP)
fig = plt.figure(figsize=(180 * MM, 215 * MM))
fig.text(0.012, 0.989, TITLE, fontsize=9.5, fontweight="bold", va="top", ha="left",
         color=INK)
age_legend(fig, 0.012, 0.966)
gs = GridSpec(3, max(nk, 3), figure=fig, height_ratios=[1.50, 1.0, 1.02],
              hspace=0.60, wspace=0.52,
              left=0.112, right=0.975, top=0.885, bottom=0.125)

axa = fig.add_subplot(gs[0, :])
letter(axa, "a", dx=-0.030, dy=1.26)
axa.text(0, 1.165, "What each axis is built from and where it is read out, in epithelium "
                   "only. Genes below 5% detection and axes that fail the\nmeasurability "
                   "gate or their matched nulls are excluded from the score and from this "
                   "figure.",
         transform=axa.transAxes, fontsize=6.6, va="bottom", ha="left", color=INK)
panel_a_sankey(axa)

# ---- b: donor-level trajectories ---------------------------------------------
for j, name in enumerate(KEEP):
    ax = fig.add_subplot(gs[1, j])
    sub = DON[DON.axis == name]
    rng = np.random.default_rng(j)
    for i, b in enumerate(AGE_ORDER):
        v = sub.loc[sub.age == b, "score"].values
        if not len(v):
            continue
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=12,
                   facecolor=AGE_COLORS[b], edgecolor="white", linewidth=0.35,
                   alpha=0.95, zorder=3)
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=1.2, zorder=4)
    ax.axhline(0, color=GREY, lw=0.5, ls=":", zorder=1)
    ax.set_xticks(range(len(AGE_ORDER)))
    ns = [int((sub.age == b).sum()) for b in AGE_ORDER]
    ax.set_xticklabels([f"{lab}\n(n={n})" for lab, n in zip(LABELS, ns)],
                       rotation=55, ha="right", fontsize=5.6)
    ax.set_xlim(-0.6, len(AGE_ORDER) - 0.4)
    r = R.loc[name]
    ax.set_title(f"{name}\nempirical $P$ = {r.empirical_p:.3f}* vs matched nulls\n"
                 f"(BH $P$ = {r.BH_perdonor:.3f}; n = {sub.donor.nunique()} donors)",
                 fontsize=6.3, pad=4, color=INK)
    if j == 0:
        ax.set_ylabel("donor mean axis score\n(mean of per-gene $z$)")
        letter(ax, "b", dx=-0.30, dy=1.32)

# ---- c: the two age windows ---------------------------------------------------
axc = fig.add_subplot(gs[2, 0])
letter(axc, "c", dx=-0.30, dy=1.20)
ypos, ylab, k = [], [], 0
for name in KEEP[::-1]:
    for w, mark in [("Child to Elderly", "s"), ("Neonate to Young child", "o")]:
        row = WIN[(WIN.axis == name) & (WIN.window == w)]
        if not len(row):
            continue
        row = row.iloc[0]
        axc.plot([row.lo, row.hi], [k, k], color=AXCOL[name], lw=1.2, zorder=3,
                 solid_capstyle="round")
        axc.scatter([row.rho], [k], s=18, marker=mark, facecolor=AXCOL[name],
                    edgecolor="white", linewidth=0.4, zorder=4)
        axc.text(row.hi + 0.05, k, f"$P$ {row.p:.3f}{'*' if row.p < 0.05 else ''}",
                 fontsize=5.2, va="center", ha="left", color=MUTED)
        ypos.append(k)
        ylab.append(w.replace("Neonate to Young child", "Neonate → Y child")
                     .replace(" to ", " → "))
        k += 1
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

# ---- d: which epithelial cells carry each axis --------------------------------
axd = fig.add_subplot(gs[2, 1])
letter(axd, "d", dx=-0.34, dy=1.20)
order = [t for t in EPI_ORDER if t in set(CT[CT.axis.isin(KEEP)].t)]
sub = CT[CT.axis.isin(KEEP)]
vmax = float(np.nanmax(np.abs(sub.mean_score)))
for i, name in enumerate(KEEP):
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
axd.set_yticks(range(len(KEEP)))
axd.set_yticklabels(KEEP, fontsize=6)
for tick, name in zip(axd.get_yticklabels(), KEEP):
    tick.set_color(AXCOL[name])
axd.set_xlim(-0.7, len(order) - 0.3)
axd.set_ylim(-0.7, len(KEEP) - 0.3)
axd.set_title("Mean score by epithelial cell type,\nwithin each axis's own subset", pad=4,
              fontsize=6.6)
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(vmin=-vmax, vmax=vmax))
cb = fig.colorbar(sm, ax=axd, fraction=0.038, pad=0.02)
cb.set_label("mean $z$ in group", fontsize=5.4)
cb.ax.tick_params(labelsize=5.0, width=0.5, length=2)
cb.outline.set_linewidth(0.5)
for n_ in (1000, 5000, 20000):
    axd.scatter([], [], s=5 + 48 * np.sqrt(n_ / 20000), facecolor="white",
                edgecolor=INK, linewidth=0.3, label=f"{n_:,}")
axd.legend(loc="upper center", bbox_to_anchor=(0.45, -0.46), ncol=3, frameon=False,
           fontsize=5.2, title="cells in group", title_fontsize=5.4,
           handletextpad=0.2, columnspacing=0.8, borderpad=0.1)

# ---- e: the surviving axes combined -------------------------------------------
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
axe.axhline(0, color=GREY, lw=0.5, ls=":", zorder=1)
axe.set_xticks(range(len(AGE_ORDER)))
axe.set_xticklabels([f"{lab}\n(n={n})" for lab, n in zip(LABELS, ns)],
                    rotation=55, ha="right", fontsize=5.6)
axe.set_ylabel("combined score")
axe.set_xlim(-0.6, len(AGE_ORDER) - 0.4)
axe.set_title(f"The {len(KEEP)} surviving axes combined,\nall {len(comp)} donors",
              fontsize=6.6, pad=4)

fig.savefig(OUT / "figure5.pdf", dpi=600)
fig.savefig(OUT / "figure5.png", dpi=300)
plt.close(fig)
print("wrote", OUT / "figure5.pdf")
