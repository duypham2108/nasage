"""Compose Figure 5 (main) and Supplementary Figure S5.

This is the cilium-ageing replication figure and it keeps the Figure 5 slot. The
mucosal-vaccine-readiness figure is Figure 6 (validation/figure6_mvri_plot.py).

Every panel is drawn from a file in D:\\research\\nasage_validation\\results\\, produced by
a plan that was committed before the data were accessed. Nothing here recomputes a model
statistic; the script reads results and arranges them. The one exception is the Spearman
correlation annotated on the scatters, computed from the plotted points themselves.

Main Figure 5  — the results that replicated.
Supp. Fig. S5  — the nulls, the measurement caveats, and the pre-registration table.

Usage: python validation/figure5_plot.py <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import stats

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"
GREY, INK, MUTED = "#BFBFBF", "#1A1A1A", "#6B6B6B"
MM = 1 / 25.4

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


def letter(ax, s, dx=-0.22, dy=1.12):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left", color=INK)


def scatter_vs_age(ax, x, y, colour, fit_colour=INK):
    ax.scatter(x, y, s=9, facecolor=colour, edgecolor="white", linewidth=0.3,
               alpha=0.85, zorder=3)
    b, a = np.polyfit(x, y, 1)
    xx = np.linspace(x.min(), x.max(), 100)
    n = len(x)
    se = np.sqrt(np.sum((y - (a + b * x)) ** 2) / (n - 2))
    sx = np.sqrt(np.sum((x - x.mean()) ** 2))
    ci = 1.96 * se * np.sqrt(1 / n + (xx - x.mean()) ** 2 / sx ** 2)
    ax.fill_between(xx, a + b * xx - ci, a + b * xx + ci, color=fit_colour, alpha=0.10,
                    lw=0, zorder=2)
    ax.plot(xx, a + b * xx, color=fit_colour, lw=1.2, zorder=4)


def dot_contrast(ax, ds, aw, module, col, cts):
    rng = np.random.default_rng(1)
    ticklabs = []
    for i, ct in enumerate(cts):
        sub = ds[(ds.cell_type == ct) & ds.eligible & ds[col].notna()]
        for j, (grp, c) in enumerate((("child", ORANGE), ("adult", BLUE))):
            v = sub.loc[sub.age_bin == grp, col].values
            xpos = i + (-0.19 if j == 0 else 0.19)
            ax.scatter(xpos + rng.uniform(-0.065, 0.065, len(v)), v, s=9, facecolor=c,
                       edgecolor="white", linewidth=0.3, alpha=0.9, zorder=3)
            ax.plot([xpos - 0.14, xpos + 0.14], [np.median(v)] * 2, color=INK, lw=1.3, zorder=4)
        rec = aw[(aw.module_id == module) & (aw.cell_type == ct)].iloc[0]
        star = "" if pd.isna(rec.empirical_p) or rec.empirical_p >= 0.05 else "*"
        ticklabs.append(f"{ct}\nBH {rec.p_adj_BH:.0e}\nemp. $p$ {rec.empirical_p:.3f}{star}")
    ax.set_xticks(range(len(cts)))
    ax.set_xticklabels(ticklabs, fontsize=5.9, linespacing=1.5)
    ax.set_xlim(-0.5, len(cts) - 0.5)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.18 * (hi - lo))


# ================================================================ shared data
ms = pd.read_csv(R / "gse336399" / "module_scores.csv", index_col=0)
pa = pd.read_csv(R / "gse336399" / "primary_age_effects.csv").set_index("score")
nullc = pd.read_csv(R / "gse336399" / "matched_null_control.csv").set_index("module")
draws = pd.read_csv(R / "figure5" / "gse336399_null_draws.csv.gz")
ds = pd.read_csv(R / "loske2022_v2" / "donor_celltype_scores.csv")
aw = pd.read_csv(R / "loske2022_v2" / "age_effects_with_null.csv")
CTS = ["basal", "secretory", "ciliated"]

# ================================================================ MAIN FIGURE 5
fig = plt.figure(figsize=(180 * MM, 132 * MM))
gs = GridSpec(2, 6, figure=fig, height_ratios=[1.0, 1.0], hspace=0.60, wspace=1.25,
              left=0.075, right=0.985, top=0.930, bottom=0.105)

# -- a: M02 vs age
ax = fig.add_subplot(gs[0, 0:2])
x = ms["age"].values
scatter_vs_age(ax, x, ms["M02_mucociliary_axoneme"].values, BLUE)
r = pa.loc["M02_mucociliary_axoneme"]
lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.26 * (hi - lo))
ax.set_xlabel("age (years)")
ax.set_ylabel("M02 mucociliary axoneme\n(mean z, 18 genes)")
ax.text(0.03, 0.98, f"$\\rho$ = {r.rho:.2f}, $p$ = {r.p_rho:.3f}\n"
                    f"adj. $\\beta$ = {r.beta_age:+.4f}/yr, $p$ = {r.p_adj:.3f}\n"
                    f"BH = {r.BH_adj:.2f} in the family of 13",
        transform=ax.transAxes, va="top", ha="left", fontsize=6, color=MUTED, linespacing=1.5)
letter(ax, "a", dx=-0.30)
ax.set_title("Curated module, unbiased", pad=6, loc="left")

# -- b: CILIUM28 vs age
ax = fig.add_subplot(gs[0, 2:4])
yc = ms["CILIUM28_manuscript"].values
scatter_vs_age(ax, x, yc, GREEN, fit_colour=INK)
rc = stats.spearmanr(x, yc)
bc = nullc.loc["CILIUM28_manuscript"]
lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.26 * (hi - lo))
ax.set_xlabel("age (years)")
ax.set_ylabel("CILIUM28, this study\n(mean z, 27 genes)")
ax.text(0.03, 0.98, f"$\\rho$ = {rc.statistic:.2f}, $p$ = {rc.pvalue:.3f}\n"
                    f"adj. $\\beta$ = {bc.beta_age:+.4f}/yr, $p$ = {bc.p_adj:.3f}\n"
                    f"derived in the primary cohort",
        transform=ax.transAxes, va="top", ha="left", fontsize=6, color=MUTED, linespacing=1.5)
letter(ax, "b", dx=-0.30)
ax.set_title("Manuscript gene set", pad=6, loc="left")

# -- c: matched null
ax = fig.add_subplot(gs[0, 4:6])
nb = draws.loc[draws.module == "M02_mucociliary_axoneme", "beta"].values
ax.hist(nb, bins=34, color=GREY, edgecolor="white", linewidth=0.3, zorder=2)
obs2 = pa.loc["M02_mucociliary_axoneme", "beta_age"]
obs28 = nullc.loc["CILIUM28_manuscript", "beta_age"]
for v, c in ((obs28, GREEN), (obs2, BLUE)):
    ax.axvline(v, color=c, lw=1.4, zorder=4)
ymax = ax.get_ylim()[1] * 1.26
ax.set_ylim(0, ymax)
ax.text(obs2 - 0.0004, ymax * 0.99, "M02\n+0.0234", ha="right", va="top", fontsize=6.2,
        color=BLUE, linespacing=1.3)
ax.text(obs28 - 0.0004, ymax * 0.70, "CILIUM28\n+0.0180", ha="right", va="top", fontsize=6.2,
        color=GREEN, linespacing=1.3)
xlo, xhi = ax.get_xlim(); ax.set_xlim(xlo, xhi + 0.05 * (xhi - xlo))
ax.set_xlabel("age coefficient $\\beta$ per year")
ax.set_ylabel("500 expression-matched\nrandom gene sets")
ax.text(0.03, 0.98, "0 of 500 random sets\nreached $p$ < 0.05\n(nominal 5%)",
        transform=ax.transAxes, va="top", ha="left", fontsize=6.2, color=INK, linespacing=1.5)
letter(ax, "c", dx=-0.30)
ax.set_title("Both sets clear the matched null", pad=6, loc="left")

# -- d, e: childhood contrast, two modules
for k, (mod, col, lab, ltr, span) in enumerate((
        ("M05_interferon_stimulated_genes", "M05_interferon_stimulated_genes_UCell",
         "M05 interferon-stimulated genes\n(UCell, donor mean)", "d", slice(0, 3)),
        ("M06_mhc_I_processing", "M06_mhc_I_processing_UCell",
         "M06 MHC-I processing\n(UCell, donor mean)", "e", slice(3, 6)))):
    ax = fig.add_subplot(gs[1, span])
    dot_contrast(ax, ds, aw, mod, col, CTS)
    ax.set_ylabel(lab)
    if k == 0:
        h = [plt.Line2D([], [], marker="o", ls="", color=ORANGE, ms=4, label="children (n = 18)"),
             plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=4, label="adults (n = 23)")]
        ax.legend(handles=h, loc="upper right", frameon=False, handletextpad=0.3,
                  borderpad=0.1, ncol=2, columnspacing=0.9)
    letter(ax, ltr, dx=-0.175)
    ax.set_title("Childhood contrast replicates (Loske 2022)" if k == 0
                 else "A second module, same direction", pad=6, loc="left")

# Panels f and g (programme x cell-type specificity) were REMOVED 2026-11-25: they
# duplicated Figure 6a/b exactly. The specificity analysis belongs to the mucosal
# programme story, so it lives in Figure 6 (validation/figure6_mvri_plot.py) only.
# Figure 5 is the cilium-ageing replication and ends at panel e.

fig.savefig(OUT / "figure5.pdf", dpi=600)
fig.savefig(OUT / "figure5.png", dpi=400)
print("wrote", OUT / "figure5.pdf")

# ================================================================ SUPPLEMENTARY S5
fig = plt.figure(figsize=(180 * MM, 150 * MM))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.85], hspace=0.38, wspace=0.34,
              left=0.085, right=0.985, top=0.940, bottom=0.060)

# -- a: no gradient within older adults
inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 0], hspace=0.62)
p2 = pd.read_csv(R / "pneumo_gse261277" / "module_scores_all_samples.csv", index_col=0)
p2 = p2[p2.timepoint == "baseline"]
for k, (lab, xv, yv) in enumerate((
        ("GSE336399, n = 120", ms["age"].values, ms["M05_interferon_stimulated_genes"].values),
        ("GSE261277, n = 39", p2["age"].values, p2["M05_interferon_stimulated_genes"].values))):
    axk = fig.add_subplot(inner[k])
    axk.scatter(xv, yv, s=7, facecolor=BLUE, edgecolor="white", linewidth=0.25,
                alpha=0.8, zorder=3)
    bb, aa = np.polyfit(xv, yv, 1)
    xs = np.linspace(xv.min(), xv.max(), 50)
    axk.plot(xs, aa + bb * xs, color=MUTED, lw=1.1, ls=(0, (3, 2)), zorder=4)
    rr = stats.spearmanr(xv, yv)
    axk.text(0.97, 0.95, f"{lab}\n$\\rho$ = {rr.statistic:+.2f}, $p$ = {rr.pvalue:.2f}",
             transform=axk.transAxes, ha="right", va="top", fontsize=6, color=MUTED,
             linespacing=1.4)
    axk.set_ylabel("M05 ISG", labelpad=2)
    if k == 1:
        axk.set_xlabel("age (years)")
    if k == 0:
        letter(axk, "a", dx=-0.19, dy=1.22)
        axk.set_title("No interferon gradient within older adults", pad=6, loc="left")

# -- b: depth confound
cd = pd.read_csv(R / "crownlands_figrep" / "ciliated_cdkn2a.csv.gz")
ax = fig.add_subplot(gs[0, 1])
bins = np.linspace(0, 9000, 46)
for lab, c, sel in (("CDKN2A$^{-}$", BLUE, cd.cdkn2a == 0),
                    ("CDKN2A$^{+}$", ORANGE, cd.cdkn2a > 0)):
    v = cd.loc[sel, "n_genes"].values
    ax.hist(v, bins=bins, density=True, histtype="step", lw=1.2, color=c, label=lab, zorder=3)
    ax.axvline(np.median(v), color=c, lw=0.8, ls=(0, (2, 2)), zorder=2)
ax.set_xlabel("detected genes per cell")
ax.set_ylabel("density")
lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.30 * (hi - lo))
ax.legend(loc="upper right", frameon=False, handlelength=1.2, borderpad=0.1)
ax.text(0.03, 0.97, "median 2,888 vs 4,265\n+606 genes in 181/187 donors",
        transform=ax.transAxes, va="top", ha="left", fontsize=6, color=INK, linespacing=1.5)
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax.yaxis.get_offset_text().set_fontsize(6)
letter(ax, "b", dx=-0.19)
ax.set_title("p16$^{+}$ status is largely a depth event", pad=6, loc="left")

# -- c: SASP detection floor
MODS = ["M02_mucociliary_axoneme", "M05_interferon_stimulated_genes",
        "M11_senescence_arrest", "M12_inflammatory_sasp"]
SHORT = ["M02\naxoneme", "M05\nISG", "M11\narrest", "M12\nSASP"]
own = pd.read_csv(R / "own_cohort" / "module_measurability.csv")
own = own[own.cell_type == "Ciliated"]
cover = {"primary (10x 3')": own.groupby("module").frac_detected.apply(
    lambda s: (s >= 0.05).mean())}
for key, lab in (("ziegler2021_v2", "Ziegler 2021 (10x 3')"),
                 ("loske2022_v2", "Loske 2022 (10x 3')")):
    d = pd.read_csv(R / key / "module_coverage.csv")
    cover[lab] = d[d.cell_type == "ciliated"].set_index("module_id")["frac_genes_detected"]
d = pd.read_csv(R / "crownlands" / "measurability_gate.csv")
cover["olfactory (10x Flex)"] = d[d.cell_type == "ciliated"].set_index("module")["frac_detected"]

ax = fig.add_subplot(gs[1, 0])
xpos = np.arange(len(MODS))
mk = ["o", "s", "^", "D"]
for j, lab in enumerate(cover):
    vals = [cover[lab].get(m, np.nan) for m in MODS]
    ax.scatter(xpos + (j - 1.5) * 0.17, vals, s=16, marker=mk[j], facecolor="white",
               edgecolor=BLUE if j < 3 else MUTED, linewidth=1.0, label=lab,
               zorder=4, clip_on=False)
ax.axhline(0.6, color=ORANGE, lw=0.9, ls=(0, (3, 2)), zorder=2)
ax.text(3.42, 0.63, "measurability gate (60%)", fontsize=5.8, color=ORANGE, ha="right",
        va="bottom")
ax.set_xticks(xpos); ax.set_xticklabels(SHORT, fontsize=6.2, linespacing=1.4)
ax.set_xlim(-0.45, 3.45); ax.set_ylim(-0.04, 1.14)
ax.set_ylabel("fraction of module genes\nat $\\geq$ 5% detection")
ax.legend(loc="lower left", frameon=False, handletextpad=0.1, borderpad=0.1,
          labelspacing=0.25, fontsize=5.8, bbox_to_anchor=(-0.02, -0.03))
letter(ax, "c", dx=-0.19)
ax.set_title("SASP below the droplet detection floor", pad=6, loc="left")

# -- d: prereg table
ax = fig.add_subplot(gs[1, 1])
ax.axis("off")
letter(ax, "d", dx=-0.10, dy=1.12)
ax.text(-0.02, 1.12, "Outcome-anchored tests, as pre-registered",
        transform=ax.transAxes, fontsize=7.5, va="top", ha="left")
blocks = [
    ("GSE261277 · pneumococcal challenge", "NULL", MUTED,
     "17 POS / 22 NEG, baseline nasal\n"
     "locked direction: M05 +, M12 +, M01 −\n"
     "all $p$ ≥ 0.20, stable across scoring methods;\n"
     "underpowered below $d$ ≈ 0.8"),
    ("LAIV Gambia · live intranasal vaccine", "DESCRIPTIVE", BLUE,
     "open outcome table, participant link severed\n"
     "locked direction: innate / defence axis −\n"
     "OR 0.15 and 0.48 per SD, $p$ = 0.007, BH 0.064;\n"
     "pathways are the published model's own features"),
    ("GSE253686 · NOVA influenza challenge", "NOT SPECIFIC", ORANGE,
     "53 shedders / 21 non-shedders, baseline nasal\n"
     "locked direction: P5 mucin barrier −, one-sided\n"
     "OR 0.47 (0.22–1.02), $p$ = 0.028 — but 11.4% of\n"
     "matched random sets as extreme, so no claim"),
]
y = 0.98
for head, status, scol, body in blocks:
    ax.text(0.00, y, head, transform=ax.transAxes, fontsize=6.0, va="top", ha="left",
            color=INK, fontweight="bold")
    ax.text(1.00, y, status, transform=ax.transAxes, fontsize=6.0, va="top", ha="right",
            color=scol, fontweight="bold")
    ax.text(0.00, y - 0.075, body, transform=ax.transAxes, fontsize=5.5, va="top",
            ha="left", linespacing=1.45, color=MUTED)
    y -= 0.295
    ax.plot([0, 1], [y + 0.030] * 2, transform=ax.transAxes, color="#EEEEEE", lw=0.5,
            clip_on=False)
ax.text(0.00, y - 0.055, "Predictions, direction and multiplicity were fixed by a committed\n"
                        "plan before the data were accessed.",
        transform=ax.transAxes, fontsize=5.5, va="top", ha="left", linespacing=1.6,
        color=MUTED, style="italic")

# Panels e and f (the P5 three-cohort forest and the NOVA matched null) were REMOVED
# 2026-11-25: they duplicated Figure 6d and Supplementary S6f. P5 is a mucosal
# programme result and belongs to the Figure 6 story, not to this ageing supplement.

fig.savefig(OUT / "figure5_supp.pdf", dpi=600)
fig.savefig(OUT / "figure5_supp.png", dpi=400)
print("wrote", OUT / "figure5_supp.pdf")
