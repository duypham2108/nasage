"""Compose the cilium-ageing replication figure and its supplement.

**Renamed 2026-11-25.** This script used to write figure5.pdf / figure5_supp.pdf. Figure 5 of
the manuscript is now the mucosal-vaccine-readiness figure (validation/figure5_mvri_plot.py),
so this one writes figure5_cilium_replication.pdf / _supp.pdf instead. The content is
unchanged and nothing was deleted; only the output filenames moved, so the two scripts cannot
overwrite each other.


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
fig = plt.figure(figsize=(180 * MM, 196 * MM))
gs = GridSpec(3, 6, figure=fig, height_ratios=[1.0, 1.0, 1.15], hspace=0.52, wspace=1.25,
              left=0.075, right=0.985, top=0.952, bottom=0.060)

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

# -- f: programme x cell type, the specificity matrix
cs = pd.read_csv(R / "celltype_spec" / "donor_celltype_scores.csv")
cm = pd.read_csv(R / "celltype_spec" / "meta.csv").set_index("program")
PROGS = ["P1_basal_progenitor_repair", "P3_secretory_antimicrobial", "P4_iga_transport",
         "P5_mucin_barrier", "P7_ciliary_beat_machinery", "P6_entry_factor_availability"]
PLAB = ["P1 basal\nprogenitor", "P3 secretory\nantimicrobial", "P4 IgA\ntransport",
        "P5 mucin\nbarrier", "P7 ciliary\nbeat", "P6 entry\nfactors"]
OWN = ["Basal", "Secretory", "Secretory", "Secretory", "Ciliated", "Ciliated"]
FAMS = ["Basal", "Secretory", "Ciliated"]

ax = fig.add_subplot(gs[2, 0:3])
Mx = np.array([[cs.loc[cs.__fam == f, p].mean() for f in FAMS] for p in PROGS])
im = ax.imshow(Mx, cmap="RdBu_r", vmin=-np.abs(Mx).max(), vmax=np.abs(Mx).max(),
               aspect="auto")
ax.set_xticks(range(3)); ax.set_xticklabels(FAMS, fontsize=6.5)
ax.set_yticks(range(len(PROGS))); ax.set_yticklabels(PLAB, fontsize=5.8, linespacing=1.3)
for i, p in enumerate(PROGS):
    j = FAMS.index(OWN[i])
    ok = cm.loc[p, "verdict"].startswith("SPECIFIC")
    ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                               edgecolor=INK if ok else ORANGE,
                               lw=1.6 if ok else 1.8,
                               ls="-" if ok else (0, (2, 1.4)), zorder=5))
ax.text(0.0, -0.17, "box = the cell type the programme was curated for; "
                    "dashed orange = misassigned",
        transform=ax.transAxes, fontsize=5.4, color=MUTED, ha="left", va="top")
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cb.set_label("mean programme score (z)", fontsize=5.6)
cb.ax.tick_params(labelsize=5.4)
letter(ax, "f", dx=-0.13, dy=1.14)
ax.set_title("Programmes score highest in their own cell type", pad=6, loc="left")

# -- g: specificity contrast, own cell type vs the other two
ax = fig.add_subplot(gs[2, 3:6])
ypos = np.arange(len(PROGS))[::-1]
for yy, p in zip(ypos, PROGS):
    r = cm.loc[p]
    spec = r.verdict.startswith("SPECIFIC")
    col = INK if spec else ORANGE
    ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=1.2, zorder=3,
            solid_capstyle="butt")
    ax.scatter([r.smd], [yy], s=22, facecolor=col, edgecolor="white", linewidth=0.5,
               zorder=4)
ax.axvline(0, color=GREY, lw=0.8, ls=(0, (2, 2)), zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels([])          # rows are in the same order as panel f
ax.set_ylim(-0.9, len(PROGS) - 0.2)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("own cell type − other two (SMD)")
ax.set_xlim(-3.2, 7.6)
for _yy, _p in zip(ypos, PROGS):
    _r = cm.loc[_p]
    ax.text(7.4, _yy, f"{_r.smd:+.2f}", fontsize=5.8, ha="right", va="center",
            color=INK if _r.verdict.startswith("SPECIFIC") else ORANGE)
ax.text(7.4, ypos[-1] - 0.58, "misassigned", fontsize=5.5, color=ORANGE, ha="right",
        va="center")
ax.text(0.01, 0.99, "rows as in f · 3 cohorts, 224 donor × cell-type units\nsign unanimous in "
                    "every programme and cohort",
        transform=ax.transAxes, fontsize=5.4, color=MUTED, va="top", ha="left",
        linespacing=1.5)
letter(ax, "g", dx=-0.10, dy=1.14)
ax.set_title("Five specific, one points the wrong way", pad=6, loc="left")

# Renamed 2026-11-25. Figure 5 of the manuscript is now the mucosal-vaccine-readiness figure
# (validation/figure5_mvri_plot.py), which writes figure5.pdf. This cilium-replication figure
# keeps its own name so the two scripts cannot silently overwrite each other.
fig.savefig(OUT / "figure5_cilium_replication.pdf", dpi=600)
fig.savefig(OUT / "figure5_cilium_replication.png", dpi=400)
print("wrote", OUT / "figure5_cilium_replication.pdf")

# ================================================================ SUPPLEMENTARY S5
fig = plt.figure(figsize=(180 * MM, 210 * MM))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.0, 0.85, 0.80], hspace=0.52, wspace=0.34,
              left=0.085, right=0.985, top=0.952, bottom=0.042)

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

# -- e: three-cohort forest for P5, with the NOVA caveat on its face
fo = pd.read_csv(R / "nova" / "forest_three_cohort.csv")
pr = pd.read_csv(R / "nova" / "primary.csv").iloc[0]
ax = fig.add_subplot(gs[2, 0])
ypos = np.arange(len(fo))[::-1]
for yy, (_, r) in zip(ypos, fo.iterrows()):
    pooled_row = r.cohort.startswith("pooled")
    col = INK if pooled_row else (ORANGE if r.cohort.startswith("NOVA") else BLUE)
    ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=1.1, zorder=3,
            solid_capstyle="butt")
    ax.scatter([r.OR], [yy], s=30 if pooled_row else 20,
               marker="D" if pooled_row else "o", facecolor=col, edgecolor="white",
               linewidth=0.5, zorder=4)
ax.axvline(1.0, color=GREY, lw=0.8, ls=(0, (2, 2)), zorder=1)
ax.set_xscale("log")
ax.set_xlim(0.02, 30.0)
ax.set_xticks([0.05, 0.2, 1.0])
ax.set_xticklabels(["0.05", "0.2", "1"])
SHORTC = ["RSV", "pneumococcal", "NOVA", "pooled"]
ax.set_yticks(ypos)
ax.set_yticklabels(SHORTC, fontsize=6.2)
ax.set_ylim(-0.8, len(fo) - 0.2)
ax.set_xlabel("odds ratio per SD of baseline P5 (log scale)")
ax.tick_params(axis="y", length=0)
ax.spines["bottom"].set_bounds(0.02, 1.6)
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
for yy, (_, r) in zip(ypos, fo.iterrows()):
    col = ORANGE if r.cohort.startswith("NOVA") else INK
    ax.text(28, yy, f"{r.OR:.2f} ({r.ci_lo:.2f}–{r.ci_hi:.2f})", fontsize=5.8,
            color=col, ha="right", va="center")
ax.text(28, ypos[2] - 0.42, "fails matched null (panel f)", fontsize=5.5, color=ORANGE,
        ha="right", va="center")
ax.text(28, ypos[3] - 0.46, f"$I^2$ = {pr.meta_I2:.0f}%", fontsize=5.8, color=MUTED,
        ha="right", va="center")
letter(ax, "e", dx=-0.19, dy=1.16)
ax.set_title("Protective direction in three cohorts", pad=6, loc="left")

# -- f: why the NOVA cohort cannot carry it
nd = pd.read_csv(R / "nova" / "matched_null_draws.csv")
ax = fig.add_subplot(gs[2, 1])
obs = float(np.log(pr.OR))
ax.hist(nd.logOR.values, bins=34, color=GREY, edgecolor="white", linewidth=0.3, zorder=2)
ax.axvline(obs, color=ORANGE, lw=1.3, zorder=4)
lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi * 1.34)
ax.annotate("P5 observed", xy=(obs, hi * 0.30), xytext=(obs - 0.30, hi * 0.52),
            fontsize=6, color=ORANGE, ha="right", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.7, color=ORANGE))
ax.set_xlabel("log odds ratio, expression-matched random gene sets")
ax.set_ylabel("number of gene sets (of 500)")
ax.text(0.02, 0.97,
        f"empirical $p$ = {pr.empirical_p:.3f}  (FAIL)\n"
        f"{pr.null_frac_p05 * 100:.1f}% of random sets reach $p$ < 0.05\n"
        f"against a nominal 5%; 0/500 in both\n"
        f"earlier cohorts",
        transform=ax.transAxes, va="top", ha="left", fontsize=6, color=INK,
        linespacing=1.6)
letter(ax, "f", dx=-0.19, dy=1.16)
ax.set_title("Specificity fails in the largest cohort", pad=6, loc="left")

fig.savefig(OUT / "figure5_cilium_replication_supp.pdf", dpi=600)
fig.savefig(OUT / "figure5_cilium_replication_supp.png", dpi=400)
print("wrote", OUT / "figure5_cilium_replication_supp.pdf")
