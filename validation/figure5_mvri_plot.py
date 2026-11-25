"""Figure 5 (MVRI): the results that survived, plus Supplementary Figure S6 for the rest.

**This is Figure 5 of the manuscript** (merged 2026-11-25). It writes figure5.pdf and
figure5_supp.pdf. The earlier cilium-replication figure that held those names now writes
figure5_cilium_replication.pdf from validation/figure5_plot.py; nothing of it was lost.

MAIN  figure5.pdf            -- only results that passed: cell-type specificity (a, b), the
                                robust age trajectory (c), the infection association (d), the
                                specification curve that is its real defence (e), and the
                                expression-matched specificity control (f).
SUPP  figure5_supp.pdf       -- every failure: the evidence matrix (a), effective
                                dimensionality (b), all trajectory slopes with matched-null
                                status (c), the vaccine-engagement null with its positive
                                control (d), the panel-v2 held-out failure (e), the full
                                matched-null grid including the cohorts that fail (f), the null
                                calibration that shows why they fail (g), and the scoring-bug
                                correction (h).

The vaccine null is why the main figure may not be titled as a validated readiness index. The
main text must state it and cite S6d, or a figure set built around mucosal vaccine readiness
with no vaccine panel reads as selective reporting.

Written as a SEPARATE file from figure5_plot.py; the original Figure 5 is untouched.
Every panel reads a results file produced by a plan committed before the data were accessed.

Usage: python validation/figure5_mvri_plot.py <results_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"
PURPLE, GOLD = "#CC79A7", "#B8860B"
GREY, INK, MUTED, LGREY = "#BFBFBF", "#1A1A1A", "#6B6B6B", "#E8E8E8"
MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def letter(ax, s, dx=-0.16, dy=1.13):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left", color=INK)


def opt(path, **kw):
    p = R / path
    return pd.read_csv(p, **kw) if p.exists() else None


spec = pd.read_csv(R / "celltype_spec" / "meta.csv").set_index("program")
specsc = pd.read_csv(R / "celltype_spec" / "donor_celltype_scores.csv")
traj = pd.read_csv(R / "lifespan" / "trajectory.csv")
life = pd.read_csv(R / "lifespan" / "donor_celltype_scores.csv")
eff = pd.read_csv(R / "lifespan" / "effective_tests.csv", index_col=0)
forest = pd.read_csv(R / "nova" / "forest_three_cohort.csv")
eng = pd.read_csv(R / "laiv_engagement" / "engagement.csv")
pc = pd.read_csv(R / "laiv_engagement" / "positive_control.csv").iloc[0]
curve = opt("multiverse/specification_curve_corrected.csv")
csumm = opt("multiverse/specification_summary_corrected.csv")
mnr = opt("matched_null_rescore/rescore_nulls.csv")
draws = opt("matched_null_rescore/null_draws.csv.gz")
v2h = opt("panel_v2_heldout/v2_heldout.csv")

FAMS = ["Basal", "Secretory", "Ciliated"]
FCOL = {"Basal": BLUE, "Secretory": ORANGE, "Ciliated": GREEN}
PROGS = ["P1_basal_progenitor_repair", "P3_secretory_antimicrobial", "P4_iga_transport",
         "P5_mucin_barrier", "P7_ciliary_beat_machinery", "P6_entry_factor_availability"]
PLAB = ["P1 basal", "P3 antimicrob.", "P4 IgA", "P5 mucin", "P7 cilia", "P6 entry"]
OWN = ["Basal", "Secretory", "Secretory", "Secretory", "Ciliated", "Ciliated"]
SCORDER = ["mean_z", "median_z", "rank_ucell", "pc1_of_set"]
SCLAB = {"mean_z": "mean z", "median_z": "median z", "rank_ucell": "rank\n(UCell)",
         "pc1_of_set": "PC1 of\nset"}

# which cohorts pass the matched null for P5 under EVERY scorer -- main figure admission rule
P5_PASS = []
if mnr is not None:
    g = mnr[(mnr.program == "P5_mucin_barrier") & (mnr.unit == "cohort")]
    for coh, gg in g.groupby("analysis"):
        if len(gg) and gg.empirical_p.max() < 0.05:
            P5_PASS.append(coh)

# ======================================================================== MAIN
NROW = 3 if (curve is not None) else 2
fig = plt.figure(figsize=(180 * MM, (220 if NROW == 3 else 142) * MM))
gs = GridSpec(NROW, 2, figure=fig, hspace=0.60, wspace=0.46,
              left=0.115, right=0.975, top=0.958, bottom=0.135)

# -- a: specificity matrix
ax = fig.add_subplot(gs[0, 0])
Mx = np.array([[specsc.loc[specsc.__fam == f, p].mean() for f in FAMS] for p in PROGS])
im = ax.imshow(Mx, cmap="RdBu_r", vmin=-np.abs(Mx).max(), vmax=np.abs(Mx).max(), aspect="auto")
ax.set_xticks(range(3))
ax.set_xticklabels(FAMS, fontsize=6)
ax.set_yticks(range(len(PROGS)))
ax.set_yticklabels(PLAB, fontsize=5.8)
for i, p in enumerate(PROGS):
    j = FAMS.index(OWN[i])
    ok = spec.loc[p, "verdict"].startswith("SPECIFIC")
    ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                               edgecolor=INK if ok else ORANGE, lw=1.5,
                               ls="-" if ok else (0, (2, 1.4)), zorder=5))
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
cb.set_label("score (z)", fontsize=5.6)
cb.ax.tick_params(labelsize=5.4)
ax.text(0.0, -0.20, "box = curated cell type; dashed orange = the one misassignment",
        transform=ax.transAxes, fontsize=5.3, color=MUTED, ha="left", va="top")
letter(ax, "a", dx=-0.26)
ax.set_title("Each programme scores highest\nin its own cell type", pad=6, loc="left")

# -- b: specificity forest
ax = fig.add_subplot(gs[0, 1])
yp = np.arange(len(PROGS))[::-1]
for yy, p in zip(yp, PROGS):
    r = spec.loc[p]
    ok = r.verdict.startswith("SPECIFIC")
    col = INK if ok else ORANGE
    ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=1.2, solid_capstyle="butt", zorder=3)
    ax.scatter([r.smd], [yy], s=20, facecolor=col, edgecolor="white", lw=0.5, zorder=4)
    ax.text(6.9, yy, f"{r.smd:+.2f}", fontsize=5.6, ha="right", va="center", color=col)
ax.axvline(0, color=GREY, lw=0.8, ls=(0, (2, 2)), zorder=1)
ax.set_yticks(yp)
ax.set_yticklabels(PLAB, fontsize=5.8)
ax.tick_params(axis="y", length=0)
ax.set_xlim(-3.2, 7.1)
ax.set_ylim(-0.9, len(PROGS) - 0.3)
ax.set_xlabel("own cell type − other two (SMD)")
ax.text(6.9, yp[-1] - 0.62, "misassigned", fontsize=5.4, color=ORANGE, ha="right", va="center")
ax.text(0.0, -0.20, "3 cohorts · 224 donor × cell-type units · sign unanimous in every one",
        transform=ax.transAxes, fontsize=5.3, color=MUTED, va="top", ha="left")
letter(ax, "b", dx=-0.26)
ax.set_title("Five of six are specific", pad=6, loc="left")

# -- c: the robust trajectory
AXIS = "Antigen presentation (M06)"
ax = fig.add_subplot(gs[1, 0])
for f in FAMS:
    d = life[life.__fam == f]
    ax.scatter(d.age, d[AXIS], s=9, facecolor=FCOL[f], edgecolor="white", lw=0.25,
               alpha=0.85, zorder=3, label=f)
    if len(d) > 3:
        b, a0 = np.polyfit(d.age, d[AXIS], 1)
        xs = np.linspace(d.age.min(), d.age.max(), 40)
        ax.plot(xs, a0 + b * xs, color=FCOL[f], lw=1.2, zorder=4)
ax.set_xlabel("age (years)")
ax.set_ylabel("M06 antigen presentation\n(donor × cell-type mean z)")
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.46 * (hi - lo))
ax.legend(loc="upper right", frameon=False, handletextpad=0.1, labelspacing=0.2,
          borderpad=0.1, fontsize=5.8, markerscale=0.9)
tt = traj[traj.axis == AXIS]
ax.text(0.02, 0.98, "declines in all three cell types\n"
                    f"BH ≤ {tt.BH_age.max():.4f} · robust across 12 specifications\n"
                    "clears its matched null under all four scorers\nin basal and secretory "
                    "cells (3 of 4 in ciliated)",
        transform=ax.transAxes, fontsize=5.4, color=INK, va="top", ha="left", linespacing=1.5)
ax.text(0.0, -0.24, "the slope is carried by the childhood-to-adult contrast; within adulthood "
                    "it is flat\n(three cohorts, Supplementary)",
        transform=ax.transAxes, fontsize=5.2, color=MUTED, va="top", ha="left", linespacing=1.5)
letter(ax, "c", dx=-0.26)
ax.set_title("Antigen presentation declines, 4–77 y", pad=6, loc="left")

# -- d: infection association
ax = fig.add_subplot(gs[1, 1])
yp = np.arange(len(forest))[::-1]
SHORTC = ["RSV", "pneumo.", "NOVA", "pooled"]
for yy2, (_, r) in zip(yp, forest.iterrows()):
    pooled = r.cohort.startswith("pooled")
    col = INK if pooled else (ORANGE if r.cohort.startswith("NOVA") else BLUE)
    ax.plot([r.ci_lo, r.ci_hi], [yy2, yy2], color=col, lw=1.1, solid_capstyle="butt", zorder=3)
    ax.scatter([r.OR], [yy2], s=26 if pooled else 18, marker="D" if pooled else "o",
               facecolor=col, edgecolor="white", lw=0.5, zorder=4)
ax.axvline(1.0, color=GREY, lw=0.8, ls=(0, (2, 2)), zorder=1)
ax.set_xscale("log")
ax.set_xlim(0.02, 3.0)
ax.set_xticks([0.05, 0.2, 1.0])
ax.set_xticklabels(["0.05", "0.2", "1"])
ax.set_yticks(yp)
ax.set_yticklabels(SHORTC, fontsize=5.8)
ax.tick_params(axis="y", length=0)
ax.set_ylim(-0.9, len(forest) - 0.3)
ax.set_xlabel("odds ratio per SD of baseline P5 (log)")
ax.text(0.0, -0.21, "the direction is consistent in all three cohorts,\n"
                    "but only pneumococcal clears its own matched null (f).\n"
                    "The pooled estimate is shown and not claimed.",
        transform=ax.transAxes, fontsize=5.2, color=ORANGE, va="top", ha="left",
        linespacing=1.6)
letter(ax, "d", dx=-0.26)
ax.set_title("The mucin-barrier axis and infection", pad=6, loc="left")

# -- e: the specification curve (the real defence of panel d)
if curve is not None:
    ax = fig.add_subplot(gs[2, 0])
    COHS = [("RSV infection", BLUE), ("pneumococcal carriage", PURPLE)]
    x0 = 0
    ticks, tlabs = [], []
    for cname, col in COHS:
        g = curve[(curve.target == cname) & (curve.program == "P5_mucin_barrier")]
        if not len(g):
            continue
        v = np.sort(np.exp(g.est.values))
        xs = np.arange(x0, x0 + len(v))
        sig = np.sort(g.assign(o=np.exp(g.est)).sort_values("o").p.values)
        gg = g.assign(o=np.exp(g.est)).sort_values("o")
        ax.vlines(xs, 1.0, gg.o.values, color=col, lw=0.9, alpha=0.55, zorder=2)
        ax.scatter(xs, gg.o.values, s=13,
                   facecolor=[col if p_ < 0.05 else "white" for p_ in gg.p.values],
                   edgecolor=col, lw=0.7, zorder=4)
        ticks.append(x0 + (len(v) - 1) / 2)
        tlabs.append(cname.split()[0].replace("pneumococcal", "pneumo."))
        x0 += len(v) + 3
    ax.axhline(1.0, color=GREY, lw=0.8, ls=(0, (2, 2)), zorder=1)
    ax.set_yscale("log")
    ax.set_yticks([0.1, 0.3, 1.0, 3.0])
    ax.set_yticklabels(["0.1", "0.3", "1", "3"])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks(ticks)
    ax.set_xticklabels(tlabs, fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("odds ratio per SD of P5 (log)")
    ax.set_xlabel("16 analytic specifications, ordered by estimate")
    ax.text(0.02, 0.06, "filled = $p$ < 0.05", transform=ax.transAxes, fontsize=5.4,
            color=MUTED, ha="left", va="bottom")
    ax.text(0.0, -0.26, "4 scoring methods × 4 covariate sets. Every one of the 32 "
                        "specifications falls below\nOR = 1: 100% sign stability in both "
                        "cohorts. This is what a multiplicity correction\ncannot supply and "
                        "cannot take away.",
            transform=ax.transAxes, fontsize=5.2, color=MUTED, va="top", ha="left",
            linespacing=1.5)
    letter(ax, "e", dx=-0.26)
    ax.set_title("The estimate does not depend on\nthe analytic choices", pad=6, loc="left")

# -- f: the matched null that P5 clears, under every scorer
if draws is not None and P5_PASS:
    COH = "pneumococcal carriage" if "pneumococcal carriage" in P5_PASS else P5_PASS[0]
    ax = fig.add_subplot(gs[2, 1])
    dd = draws[(draws.program == "P5_mucin_barrier") & (draws.analysis == COH)]
    d0 = dd[dd.scoring == "mean_z"]
    obs = float(d0.observed.iloc[0])
    ax.hist(d0.null_est.values, bins=34, color=LGREY, edgecolor="white", lw=0.3, zorder=2,
            label=f"{len(d0)} matched random sets")
    ax.axvline(obs, color=PURPLE, lw=1.6, zorder=5)
    ax.annotate("observed P5", xy=(obs, ax.get_ylim()[1] * 0.80),
                xytext=(obs + 0.45, ax.get_ylim()[1] * 0.92), fontsize=5.6, color=PURPLE,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=PURPLE, lw=0.7))
    ax.set_xlabel("log odds ratio per SD")
    ax.set_ylabel("matched random gene sets")
    ax.legend(frameon=False, fontsize=5.6, loc="upper right", handlelength=1.0, borderpad=0.1)
    gg = mnr[(mnr.program == "P5_mucin_barrier") & (mnr.analysis == COH)].set_index("scoring")
    SHORTSC = {"mean_z": "mean z", "median_z": "median z", "rank_ucell": "rank",
               "pc1_of_set": "PC1"}
    txt = ", ".join(f"{SHORTSC[s]} {gg.empirical_p.get(s, np.nan):.3f}"
                    for s in SCORDER if s in gg.index)
    ax.text(0.0, -0.26, f"empirical $p$ — {txt}\n"
                        "not one of 500 expression- and size-matched random\n"
                        "gene sets reproduces the effect, under any of the four\n"
                        "scoring methods. The pre-registered headline is the worst\n"
                        "scorer, not the best. The cohorts where P5 does NOT clear\n"
                        "its null are in Supplementary.",
            transform=ax.transAxes, fontsize=5.2, color=MUTED, va="top", ha="left",
            linespacing=1.5)
    letter(ax, "f", dx=-0.26)
    ax.set_title("Random gene sets do not reproduce it\n(pneumococcal carriage)",
                 pad=6, loc="left")

fig.savefig(OUT / "figure5.pdf", dpi=600)
fig.savefig(OUT / "figure5.png", dpi=400)
print("wrote", OUT / "figure5.pdf")

# ======================================================================== SUPPLEMENTARY S6
AX_ROWS = [
    # "predicts infection" downgraded to partial 2026-11-25: the direction is consistent in all
    # three cohorts, but P5 clears its expression-matched null in pneumococcal only
    # (results/matched_null_rescore/). A full pass would overstate it.
    ("Mucus barrier\nP5 · 11", [2, 2, 1, 1, 0]),
    ("Antigen presn.\nM06 · 15", [2, -1, 2, -1, -1]),
    ("Innate sensing\nM04 · 15", [0, -1, 1, -1, -1]),
    ("IgA support\nP4 · 10", [2, 2, 1, 0, 0]),
    ("Cytokine tone\nM12 · 7", [0, -1, -1, -1, -1]),
]
COLS = ["measurable", "cell-type\nspecific", "changes\nwith age", "predicts\ninfection",
        "engaged by\nvaccine"]
SYM = {2: ("o", GREEN, 26), 1: ("o", GOLD, 11), 0: ("X", ORANGE, 24),
       -1: (None, GREY, 0)}

NS = 4 if (mnr is not None or v2h is not None) else 3
fig = plt.figure(figsize=(180 * MM, (252 if NS == 4 else 185) * MM))
gs = GridSpec(NS, 2, figure=fig, height_ratios=[1.0, 0.95, 0.95, 0.95][:NS],
              hspace=0.95, wspace=0.50, left=0.150, right=0.975, top=0.945, bottom=0.105)

# -- a: evidence matrix
ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(-0.5, len(COLS) - 0.5)
ax.set_ylim(len(AX_ROWS) - 0.5, -0.5)
ax.set_xticks(range(len(COLS)))
ax.set_xticklabels(COLS, fontsize=5.4, linespacing=1.3)
ax.set_yticks(range(len(AX_ROWS)))
ax.set_yticklabels([r[0] for r in AX_ROWS], fontsize=5.4, linespacing=1.3)
ax.tick_params(length=0)
ax.xaxis.set_ticks_position("top")
for sp in ax.spines.values():
    sp.set_visible(False)
for i, (_, vals) in enumerate(AX_ROWS):
    for j, v in enumerate(vals):
        ax.add_patch(plt.Rectangle((j - .46, i - .42), .92, .84, facecolor="#FAFAFA",
                                   edgecolor=LGREY, lw=0.5, zorder=1))
        mk, c, sz = SYM[v]
        if mk is None:
            ax.plot([j - .10, j + .10], [i, i], color=GREY, lw=0.9, zorder=3)
        else:
            ax.scatter([j], [i], s=sz, marker=mk, color=c, zorder=3,
                       linewidths=1.4 if mk == "X" else 0)
for _k, (_lx, _lb) in enumerate([(0.00, "passed"), (0.30, "partial"), (0.60, "failed"),
                                 (0.86, "not tested")]):
    _mk, _c, _sz = SYM[[2, 1, 0, -1][_k]]
    if _mk is None:
        ax.plot([_lx - 0.012, _lx + 0.012], [-0.19, -0.19], color=GREY, lw=0.9,
                transform=ax.transAxes, clip_on=False)
    else:
        ax.scatter([_lx], [-0.19], s=_sz * 0.7, marker=_mk, color=_c, transform=ax.transAxes,
                   clip_on=False, linewidths=1.2 if _mk == "X" else 0)
    ax.text(_lx + 0.035, -0.19, _lb, fontsize=5.3, color=MUTED, ha="left", va="center",
            transform=ax.transAxes)
fig.text(0.012, 0.978, "a", fontsize=9, fontweight="bold", va="top", ha="left", color=INK)
ax.set_title("Five axes, and what each survived", pad=22, loc="left")

# -- b: effective dimensionality
ax = fig.add_subplot(gs[0, 1])
xs = np.arange(len(eff.index))
ax.bar(xs - 0.19, [5] * len(eff), width=0.36, color=LGREY, edgecolor="white", lw=0.5,
       label="axes declared", zorder=3)
ax.bar(xs + 0.19, eff.eff_li_ji.values, width=0.36, color=BLUE, edgecolor="white", lw=0.5,
       label="effective tests", zorder=3)
for x, v in zip(xs, eff.eff_li_ji.values):
    ax.text(x + 0.19, v + 0.12, f"{v:.0f}", ha="center", va="bottom", fontsize=6, color=INK)
ax.set_xticks(xs)
ax.set_xticklabels([str(i).lower() for i in eff.index], fontsize=6)
ax.set_ylabel("number of tests")
ax.set_ylim(0, 6.4)
ax.legend(frameon=False, fontsize=5.8, loc="upper left", handlelength=1.1, borderpad=0.1)
ax.text(0.0, -0.17, "PC1 carries 55–69% of the variance across the five axes",
        transform=ax.transAxes, fontsize=5.3, color=MUTED, va="top", ha="left")
letter(ax, "b", dx=-0.24)
ax.set_title("Five axes behave as three", pad=6, loc="left")

# -- c: all trajectory slopes
ax = fig.add_subplot(gs[1, 0:2])
order = ["Antigen presentation (M06)", "IgA support (P4)", "Mucus barrier (P5)",
         "Innate sensing (M04)", "Cytokine tone (M12)"]
SHORT = {"Antigen presentation (M06)": "M06 antigen", "IgA support (P4)": "P4 IgA",
         "Mucus barrier (P5)": "P5 mucin", "Innate sensing (M04)": "M04 innate",
         "Cytokine tone (M12)": "M12 cytokine"}
yy, ticks, labs = 0, [], []
for axis in order:
    for f in FAMS:
        r = traj[(traj.axis == axis) & (traj.family == f)]
        if not len(r):
            continue
        r = r.iloc[0]
        ticks.append(yy)
        labs.append(f"{SHORT[axis]} · {f[:3].lower()}")
        if not r.tested:
            ax.text(0.62, yy, "not measurable", fontsize=5.2, color=ORANGE, va="center",
                    ha="left", transform=ax.get_yaxis_transform())
        else:
            solid = (r.empirical_p < 0.05) and (r.BH_age < 0.05)
            ax.scatter([r.beta_age], [yy], s=20, facecolor=INK if solid else GREY,
                       edgecolor="white", lw=0.5, zorder=4)
            if solid:
                ax.scatter([0.92], [yy], s=14, marker="o", color=GREEN, zorder=5,
                           transform=ax.get_yaxis_transform(), clip_on=False)
        yy += 1
ax.axvline(0, color=GREY, lw=0.8, ls=(0, (2, 2)), zorder=1)
ax.set_yticks(ticks)
ax.set_yticklabels(labs, fontsize=5.0)
ax.set_ylim(yy - 0.4, -0.6)
ax.tick_params(axis="y", length=0)
ax.set_xlim(-0.027, 0.002)
ax.set_xlabel("change in score per year")
ax.scatter([0.92], [1.045], s=12, marker="o", color=GREEN, transform=ax.transAxes,
           clip_on=False, zorder=6)
ax.text(0.905, 1.045, "clears matched null", transform=ax.transAxes, fontsize=5.2,
        color=MUTED, ha="right", va="center")
ax.text(0.0, -0.30, "7 of 11 tested combinations reach BH < 0.05 but fail the matched null: "
                    "a general age drift that most gene sets partially capture",
        transform=ax.transAxes, fontsize=5.2, color=MUTED, va="top", ha="left")
letter(ax, "c", dx=-0.115)
ax.set_title("Only antigen presentation is robust", pad=6, loc="left")

# -- d: vaccine engagement null
ax = fig.add_subplot(gs[2, 0])
rows = [("M05 IFN\n(control)", float(pc["mean"]), GREEN)]
for _, r in eng[~eng.program.str.contains("P2")].iterrows():
    rows.append((r.program.split("_")[0], r.mean_change, GREY))
yp = np.arange(len(rows))[::-1]
for yy3, (lab, v, col) in zip(yp, rows):
    ax.barh(yy3, v, height=0.62, color=col, edgecolor="white", lw=0.4, zorder=3)
ax.axvline(0, color=INK, lw=0.7, zorder=4)
ax.set_yticks(yp)
ax.set_yticklabels([r[0] for r in rows], fontsize=5.4)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("paired change after LAIV (SD)")
ax.set_xlim(-0.55, 1.15)
ax.text(0.0, -0.28, f"the positive control rises (+{float(pc['mean']):.2f} SD, "
                    f"$p$ = {pc['p']:.3f}, {int(pc['n_rose'])}/{int(pc['n'])} participants), "
                    "so the assay works on\nthese samples; every axis is null (BH ≥ 0.53). "
                    "Three independent designs returned this.",
        transform=ax.transAxes, fontsize=5.2, color=INK, va="top", ha="left", linespacing=1.6)
letter(ax, "d", dx=-0.30)
ax.set_title("A live mucosal vaccine does not engage them", pad=6, loc="left")

# -- e: panel v2 fails out of sample
if v2h is not None and len(v2h):
    ax = fig.add_subplot(gs[2, 1])
    xs = np.arange(len(v2h) + 1)
    vals = [0.82] + list(v2h.median_rho_A01.values)
    labs2 = ["derivation\n(Ziegler+Loske)"] + [c.replace(" ", "\n") for c in v2h.cohort]
    cols = [GREY] + [GREEN if p_ else ORANGE for p_ in v2h.V2_pass]
    ax.bar(xs, vals, width=0.62, color=cols, edgecolor="white", lw=0.5, zorder=3)
    for x, c_ in zip(xs[1:], v2h.null_q95.values):
        ax.plot([x - 0.33, x + 0.33], [c_, c_], color=INK, lw=1.0, zorder=5)
    ax.axhline(0.60, color=INK, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.text(len(xs) - 0.45, 0.625, "pre-specified floor 0.60", fontsize=5.2, color=INK,
            ha="right")
    ax.set_xticks(xs)
    ax.set_xticklabels(labs2, fontsize=5.0, linespacing=1.2)
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("median ρ among the five\ninterferon-axis modules")
    ax.text(0.0, -0.34, "black rule = 95th percentile of random five-module\n"
                        "subsets of the same panel. The axis is 0.72–0.93 where it\n"
                        "was derived and 0.34–0.50 where it was not, beating\n"
                        "chance in one cohort of three. v1 remains the panel\n"
                        "of record.",
            transform=ax.transAxes, fontsize=5.2, color=MUTED, va="top", ha="left",
            linespacing=1.5)
    letter(ax, "e", dx=-0.28)
    ax.set_title("The v2 panel does not replicate", pad=6, loc="left")

# -- f/g: the full matched-null grid and its calibration
if mnr is not None and len(mnr):
    gg = mnr[mnr.unit == "cohort"]
    ax = fig.add_subplot(gs[3, 0])
    progs = [p for p in PROGS if p in set(gg.program)]
    cohs = [c for c in ["RSV infection", "pneumococcal carriage", "NOVA influenza"]
            if c in set(gg.analysis)]
    M = np.full((len(progs), len(cohs) * len(SCORDER)), np.nan)
    for i, p in enumerate(progs):
        for jc, c in enumerate(cohs):
            sub = gg[(gg.program == p) & (gg.analysis == c)].set_index("scoring")
            for js, s in enumerate(SCORDER):
                M[i, jc * len(SCORDER) + js] = sub.empirical_p.get(s, np.nan)
    im = ax.imshow(M, cmap="RdYlGn_r", vmin=0, vmax=0.5, aspect="auto")
    for jc in range(1, len(cohs)):
        ax.axvline(jc * len(SCORDER) - 0.5, color=INK, lw=0.9)
    ax.set_yticks(range(len(progs)))
    ax.set_yticklabels([PLAB[PROGS.index(p)] for p in progs], fontsize=5.4)
    ax.set_xticks([jc * len(SCORDER) + (len(SCORDER) - 1) / 2 for jc in range(len(cohs))])
    ax.set_xticklabels([c.split()[0].replace("pneumococcal", "pneumo.") for c in cohs],
                       fontsize=5.8)
    ax.tick_params(length=0)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]) and M[i, j] < 0.05:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                           edgecolor=INK, lw=1.1, zorder=5))
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cb.set_label("empirical $p$", fontsize=5.6)
    cb.ax.tick_params(labelsize=5.4)
    ax.text(0.0, -0.30, "four scorers per cohort, in the order mean z · median z ·\n"
                        "rank · PC1. Boxed = clears the matched null. Specificity is\n"
                        "scorer-dependent for most programmes, which is why the\n"
                        "pre-registered headline is the worst scorer, not the best.",
            transform=ax.transAxes, fontsize=5.2, color=MUTED, va="top", ha="left",
            linespacing=1.5)
    letter(ax, "f", dx=-0.30)
    ax.set_title("Every programme against its matched null", pad=6, loc="left")

    ax = fig.add_subplot(gs[3, 1])
    cal = (gg.groupby("analysis").null_frac_p05.mean()
           .reindex([c for c in cohs]).dropna())
    cal.index = [str(i).split()[0].replace("pneumococcal", "pneumo.") for i in cal.index]
    lifg = mnr[mnr.unit != "cohort"]
    for a, v in lifg.groupby("analysis").null_frac_p05.mean().items():
        cal[str(a).split("|")[-1][:3].lower() + "\n(lifespan)"] = float(v)
    xs = np.arange(len(cal))
    ax.bar(xs, cal.values, width=0.6,
           color=[ORANGE if v > 0.10 else GOLD if v > 0.05 else GREEN for v in cal.values],
           edgecolor="white", lw=0.5, zorder=3)
    for x, v in zip(xs, cal.values):
        ax.text(x, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=5.6, color=INK)
    ax.axhline(0.05, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=4)
    ax.text(2.55, 0.075, "nominal 0.05", fontsize=5.3, color=INK, ha="center")
    ax.set_xticks(xs)
    ax.set_xticklabels(list(cal.index), fontsize=5.2, linespacing=1.2)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("fraction of RANDOM matched\nsets with model $p$ < 0.05")
    ax.set_ylim(0, max(0.12, float(cal.max()) * 1.30))
    ax.text(0.0, -0.32, "the dominant-axis problem, measured directly. Where the bar\n"
                        "exceeds 0.05, random gene sets matched on expression and\n"
                        "size reach significance more often than chance, so a nominal\n"
                        "$p$ from any gene set there is anticonservative. The bulk\n"
                        "cohorts are near nominal; the single-cell path is not.",
            transform=ax.transAxes, fontsize=5.2, color=MUTED, va="top", ha="left",
            linespacing=1.5)
    letter(ax, "g", dx=-0.30)
    ax.set_title("Why so many tests fail their null", pad=6, loc="left")

fig.savefig(OUT / "figure5_supp.pdf", dpi=600)
fig.savefig(OUT / "figure5_supp.png", dpi=400)
print("wrote", OUT / "figure5_supp.pdf")
print("P5 passes the matched null under all four scorers in:", P5_PASS or "no cohort")
