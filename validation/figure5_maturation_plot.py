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
WS = pd.read_csv(AUD / "weight_sensitivity.csv")
EFF = effective_tests(C)
N_TOTAL = DON.donor.nunique()      # donors carrying at least one axis subset

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


# ================================================================== MAIN FIGURE
fig = plt.figure(figsize=(180 * MM, 193 * MM))
fig.text(0.012, 0.988, TITLE, fontsize=10, fontweight="bold", va="top", ha="left",
         color=INK)
age_legend(fig, 0.012, 0.963, ncol=7)
gs = GridSpec(3, 3, figure=fig, height_ratios=[0.60, 1.0, 1.02],
              hspace=0.58, wspace=0.52,
              left=0.112, right=0.975, top=0.880, bottom=0.135)

# ---- a: the five axes, and what is measurable ---------------------------------
axa = fig.add_subplot(gs[0, :])
axa.set_xlim(0, 5)
axa.set_ylim(0, 1)
axa.axis("off")
letter(axa, "a", dx=-0.032, dy=1.42)
axa.text(0, 1.30, "The five axes of the index, and the fraction of each one's genes detected "
                  "in ≥ 5% of the cells it is scored in",
         transform=axa.transAxes, fontsize=7, va="bottom", ha="left", color=INK)

for i, name in enumerate(AXES):
    r = R.loc[name]
    ok = bool(r.gate_pass)
    edge = AXCOL[name] if ok else GREY
    axa.add_patch(FancyBboxPatch((i + 0.04, 0.08), 0.92, 0.82,
                                 boxstyle="round,pad=0.0,rounding_size=0.04",
                                 linewidth=0.9, edgecolor=edge,
                                 facecolor="#F6F8FA" if ok else "#F2F2F2",
                                 transform=axa.transData, zorder=2))
    tcol = INK if ok else MUTED
    axa.text(i + 0.50, 0.785, name, ha="center", va="center", fontsize=7,
             fontweight="bold", color=tcol, zorder=3)
    axa.text(i + 0.50, 0.645, f"{int(r.n_genes_found)} genes", ha="center", va="center",
             fontsize=6, color=MUTED, zorder=3)
    axa.text(i + 0.50, 0.545, SUBSET_LABEL[name], ha="center", va="center",
             fontsize=5.3, color=MUTED, zorder=3, style="italic")
    x0, w = i + 0.14, 0.72
    axa.add_patch(plt.Rectangle((x0, 0.30), w, 0.085, facecolor=FAINT, lw=0, zorder=3))
    axa.add_patch(plt.Rectangle((x0, 0.30), w * float(r.frac_measurable), 0.085,
                                facecolor=edge, lw=0, zorder=4))
    axa.plot([x0 + w * GATE] * 2, [0.275, 0.41], color=INK, lw=0.8, zorder=5)
    if i == 0:
        axa.text(x0 + w * GATE, 0.445, "60% gate", ha="center", va="center",
                 fontsize=5.2, color=INK, zorder=5)
    axa.text(i + 0.50, 0.165, f"{r.frac_measurable:.0%} detected"
             + ("" if ok else "  —  below floor"),
             ha="center", va="center", fontsize=6, color=tcol,
             fontweight="normal" if ok else "bold", zorder=3)

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
                     "B mem", "Plasma"] if t in set(CT.t)]
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
Z = S[list(AXES)].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
w = np.array([AXES[k]["weight"] for k in AXES])
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
axe.set_ylabel("weighted composite index")
axe.set_xlim(-0.6, len(AGE_ORDER) - 0.4)
axe.set_title(f"Weighted composite of all five axes\ncomputable for "
              f"{len(comp)} of {N_TOTAL} donors", fontsize=6.6, pad=4)

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
for j, name in enumerate(AXES):
    ax = figs.add_subplot(gss[1, 2 * j:2 * j + 2])
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
axnote = figs.add_subplot(gss[1, 10:])
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
    names = [n for n in AXES if n in DEP.index]
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
