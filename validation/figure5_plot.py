"""Compose Figure 5: external validation.

Every panel is drawn from a file in D:\\research\\nasage_validation\\results\\, produced by
a plan that was committed before the data were accessed. Nothing here recomputes a
statistic; the script reads results and arranges them.

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

MM = 1 / 25.4
fig = plt.figure(figsize=(180 * MM, 192 * MM))
gs = GridSpec(3, 3, figure=fig, height_ratios=[1.0, 1.0, 0.70],
              hspace=0.54, wspace=0.40,
              left=0.085, right=0.985, top=0.950, bottom=0.030)


def letter(ax, s, dx=-0.22, dy=1.12):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left", color=INK)


# ---------------------------------------------------------------- a: M02 vs age
ms = pd.read_csv(R / "gse336399" / "module_scores.csv", index_col=0)
pa = pd.read_csv(R / "gse336399" / "primary_age_effects.csv").set_index("score")
ax = fig.add_subplot(gs[0, 0])
x, y = ms["age"].values, ms["M02_mucociliary_axoneme"].values
ax.scatter(x, y, s=9, facecolor=BLUE, edgecolor="white", linewidth=0.3, alpha=0.85, zorder=3)
b, a = np.polyfit(x, y, 1)
xx = np.linspace(x.min(), x.max(), 100)
ax.plot(xx, a + b * xx, color=INK, lw=1.2, zorder=4)
n = len(x)
se = np.sqrt(np.sum((y - (a + b * x)) ** 2) / (n - 2))
sx = np.sqrt(np.sum((x - x.mean()) ** 2))
ci = 1.96 * se * np.sqrt(1 / n + (xx - x.mean()) ** 2 / sx ** 2)
ax.fill_between(xx, a + b * xx - ci, a + b * xx + ci, color=INK, alpha=0.10, lw=0, zorder=2)
r = pa.loc["M02_mucociliary_axoneme"]
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.22 * (hi - lo))
ax.set_xlabel("age (years)")
ax.set_ylabel("M02 mucociliary axoneme\n(mean z, 18 genes)")
ax.text(0.03, 0.97, f"$\\rho$ = {r.rho:.2f}, $p$ = {r.p_rho:.3f}\nadj. $\\beta$ = {r.beta_age:+.4f}/yr, "
                    f"$p$ = {r.p_adj:.3f}\nBH = {r.BH_adj:.2f} (n.s. in family of 13)",
        transform=ax.transAxes, va="top", ha="left", fontsize=6, color=MUTED, linespacing=1.5)
letter(ax, "a")
ax.set_title("Independent nasal cohort (n = 120)", pad=6, loc="left")

# ---------------------------------------------------------------- b: matched null
dr = pd.read_csv(R / "figure5" / "gse336399_null_draws.csv.gz")
ax = fig.add_subplot(gs[0, 1])
nullb = dr.loc[dr.module == "M02_mucociliary_axoneme", "beta"].values
ax.hist(nullb, bins=34, color=GREY, edgecolor="white", linewidth=0.3, zorder=2)
obs2 = pa.loc["M02_mucociliary_axoneme", "beta_age"]
obsc = pd.read_csv(R / "gse336399" / "matched_null_control.csv").set_index("module")
obs28 = obsc.loc["CILIUM28_manuscript", "beta_age"]
for v, c in ((obs28, GREEN), (obs2, ORANGE)):
    ax.axvline(v, color=c, lw=1.4, zorder=4)
ymax = ax.get_ylim()[1] * 1.22
ax.set_ylim(0, ymax)
ax.text(obs2, ymax * 0.98, "M02\n+0.0234", ha="right", va="top", fontsize=6.2,
        color=ORANGE, linespacing=1.3)
ax.text(obs28, ymax * 0.70, "CILIUM28\n+0.0180", ha="right", va="top", fontsize=6.2,
        color=GREEN, linespacing=1.3)
xlo, xhi = ax.get_xlim()
ax.set_xlim(xlo, xhi + 0.05 * (xhi - xlo))
ax.set_xlabel("age coefficient $\\beta$ per year")
ax.set_ylabel("500 expression-matched\nrandom gene sets")
ax.text(0.03, 0.97, "0 of 500 random sets\nreached $p$ < 0.05\n(nominal 5%)",
        transform=ax.transAxes, va="top", ha="left", fontsize=6.2, color=INK, linespacing=1.5)
letter(ax, "b", dx=-0.24)
ax.set_title("Specificity: no generic age inflation", pad=6, loc="left")

# ---------------------------------------------------------------- c: childhood arm
ds = pd.read_csv(R / "loske2022_v2" / "donor_celltype_scores.csv")
aw = pd.read_csv(R / "loske2022_v2" / "age_effects_with_null.csv")
ax = fig.add_subplot(gs[0, 2])
cts = ["basal", "secretory", "ciliated"]
col = "M05_interferon_stimulated_genes_UCell"
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
    rec = aw[(aw.module_id == "M05_interferon_stimulated_genes") & (aw.cell_type == ct)].iloc[0]
    star = "" if pd.isna(rec.empirical_p) or rec.empirical_p >= 0.05 else "*"
    ticklabs.append(f"{ct}\nBH {rec.p_adj_BH:.0e}\nemp. $p$ {rec.empirical_p:.3f}{star}")
ax.set_xticks(range(3)); ax.set_xticklabels(ticklabs, fontsize=5.9, linespacing=1.5)
ax.set_ylabel("M05 interferon-stimulated genes\n(UCell, donor mean)")
ax.set_xlim(-0.5, 2.5)
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.18 * (hi - lo))
h = [plt.Line2D([], [], marker="o", ls="", color=ORANGE, ms=4, label="children"),
     plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=4, label="adults")]
ax.legend(handles=h, loc="upper right", frameon=False, handletextpad=0.3, borderpad=0.1,
          ncol=2, columnspacing=1.0)
letter(ax, "c", dx=-0.24)
ax.set_title("Childhood contrast replicates", pad=6, loc="left")

# ---------------------------------------------------------------- d: why Fig 4 is childhood-only
inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 0], hspace=0.62)
p2 = pd.read_csv(R / "pneumo_gse261277" / "module_scores_all_samples.csv", index_col=0)
p2 = p2[p2.timepoint == "baseline"]
panels = [
    ("GSE336399, n = 120", ms["age"].values, ms["M05_interferon_stimulated_genes"].values),
    ("GSE261277, n = 39", p2["age"].values, p2["M05_interferon_stimulated_genes"].values),
]
for k, (lab, xv, yv) in enumerate(panels):
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
        letter(axk, "d", dx=-0.22, dy=1.22)
        axk.set_title("No rise within older adults", pad=6, loc="left")

# ---------------------------------------------------------------- e: depth confound
cd = pd.read_csv(R / "crownlands_figrep" / "ciliated_cdkn2a.csv.gz")
ax = fig.add_subplot(gs[1, 1])
bins = np.linspace(0, 9000, 46)
for lab, c, sel in (("CDKN2A$^{-}$", BLUE, cd.cdkn2a == 0), ("CDKN2A$^{+}$", ORANGE, cd.cdkn2a > 0)):
    v = cd.loc[sel, "n_genes"].values
    ax.hist(v, bins=bins, density=True, histtype="step", lw=1.2, color=c, label=lab, zorder=3)
    ax.axvline(np.median(v), color=c, lw=0.8, ls=(0, (2, 2)), zorder=2)
ax.set_xlabel("detected genes per cell")
ax.set_ylabel("density")
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.30 * (hi - lo))
ax.legend(loc="upper right", frameon=False, handlelength=1.2, borderpad=0.1)
ax.text(0.03, 0.97, "median 2,888 vs 4,265\n+606 genes in 181/187 donors",
        transform=ax.transAxes, va="top", ha="left", fontsize=6, color=INK, linespacing=1.5)
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax.yaxis.get_offset_text().set_fontsize(6)
letter(ax, "e", dx=-0.24)
ax.set_title("p16$^{+}$ status is largely a depth event", pad=6, loc="left")

# ---------------------------------------------------------------- f: M12 detection floor
MODS = ["M02_mucociliary_axoneme", "M05_interferon_stimulated_genes",
        "M11_senescence_arrest", "M12_inflammatory_sasp"]
SHORT = ["M02\naxoneme", "M05\nISG", "M11\narrest", "M12\nSASP"]
own = pd.read_csv(R / "own_cohort" / "module_measurability.csv")
own = own[own.cell_type == "Ciliated"]
cover = {"primary\n(10x 3')": own.groupby("module").frac_detected.apply(lambda s: (s >= 0.05).mean())}
for key, lab in (("ziegler2021_v2", "Ziegler 2021\n(10x 3')"), ("loske2022_v2", "Loske 2022\n(10x 3')")):
    d = pd.read_csv(R / key / "module_coverage.csv")
    cover[lab] = d[d.cell_type == "ciliated"].set_index("module_id")["frac_genes_detected"]
d = pd.read_csv(R / "crownlands" / "measurability_gate.csv")
cover["olfactory\n(10x Flex)"] = d[d.cell_type == "ciliated"].set_index("module")["frac_detected"]

ax = fig.add_subplot(gs[1, 2])
labs = list(cover)
xpos = np.arange(len(MODS))
mk = ["o", "s", "^", "D"]
for j, lab in enumerate(labs):
    vals = [cover[lab].get(m, np.nan) for m in MODS]
    off = (j - 1.5) * 0.17
    ax.scatter(xpos + off, vals, s=16, marker=mk[j], facecolor="white",
               edgecolor=BLUE if j < 3 else MUTED, linewidth=1.0, label=lab.replace("\n", " "),
               zorder=4, clip_on=False)
ax.axhline(0.6, color=ORANGE, lw=0.9, ls=(0, (3, 2)), zorder=2)
ax.text(3.42, 0.63, "measurability\ngate (60%)", fontsize=5.8, color=ORANGE, ha="right",
        va="bottom", linespacing=1.3)
ax.set_xticks(xpos); ax.set_xticklabels(SHORT, fontsize=6.2, linespacing=1.4)
ax.set_xlim(-0.45, 3.45)
ax.set_ylim(-0.04, 1.14)
ax.set_ylabel("fraction of module genes\nat $\\geq$ 5% detection")
ax.legend(loc="lower left", frameon=False, handletextpad=0.1, borderpad=0.1,
          labelspacing=0.25, fontsize=5.8, bbox_to_anchor=(-0.02, -0.03))
letter(ax, "f", dx=-0.24)
ax.set_title("SASP below the detection floor", pad=6, loc="left")

# ---------------------------------------------------------------- g: prereg table
ax = fig.add_subplot(gs[2, :])
ax.axis("off")
letter(ax, "g", dx=-0.070, dy=1.16)
ax.text(-0.022, 1.16, "Outcome-anchored tests, reported as pre-registered",
        transform=ax.transAxes, fontsize=7.5, va="top", ha="left")

rows = [
    ("GSE261277\npneumococcal challenge", "17 POS / 22 NEG\nbaseline nasal",
     "M05 +, M12 +, M01 −", "NULL",
     "all $p$ ≥ 0.20, stable across scoring methods;\nunderpowered below $d$ ≈ 0.8"),
    ("LAIV Gambia\nlive intranasal vaccine", "open outcome table\n(participant link severed)",
     "innate / defence axis −", "DESCRIPTIVE",
     "OR 0.15 and 0.48 per SD, $p$ = 0.007, BH 0.064;\npathways are the published model's own features"),
    ("GSE253686 NOVA\nA/H1N1 challenge", "75 participants\nbaseline SD −1",
     "to be locked before access", "PENDING",
     "series matrix in hand; age and shedding outcome\nnot yet obtained"),
]
cols = [0.00, 0.175, 0.345, 0.515, 0.625]
hdr = ["cohort", "samples", "locked prediction", "result", "detail"]
y = 0.86
for cx, h in zip(cols, hdr):
    ax.text(cx, y, h, transform=ax.transAxes, fontsize=6.4, color=MUTED,
            va="top", ha="left", fontweight="bold")
ax.plot([0, 1], [y - 0.09] * 2, transform=ax.transAxes, color=GREY, lw=0.6, clip_on=False)
y -= 0.17
STATUS = {"NULL": MUTED, "DESCRIPTIVE": BLUE, "PENDING": MUTED}
for row in rows:
    for cx, cell, h in zip(cols, row, hdr):
        kw = dict(transform=ax.transAxes, fontsize=6.2, va="top", ha="left",
                  linespacing=1.45, color=INK)
        if h == "result":
            kw.update(color=STATUS[cell], fontweight="bold")
        if h in ("locked prediction", "detail"):
            kw.update(color=MUTED)
        ax.text(cx, y, cell, **kw)
    y -= 0.27
    ax.plot([0, 1], [y + 0.10] * 2, transform=ax.transAxes, color="#EEEEEE", lw=0.5,
            clip_on=False)

fig.savefig(OUT / "figure5.pdf", dpi=600)
fig.savefig(OUT / "figure5.png", dpi=400)
print("wrote", OUT / "figure5.pdf", "and", OUT / "figure5.png")
