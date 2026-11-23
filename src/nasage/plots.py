"""Figures. Each returns a matplotlib Figure; nothing is written implicitly."""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .scoring import score_column


def _bins(summary):
    c = pd.Categorical(summary.age_bin)
    return [str(x) for x in c.categories]


def composition_figure(prop: pd.DataFrame, donor_meta: pd.DataFrame, cell_types=None,
                       highlight=("ciliated", "secretory", "Tcell", "dendritic")):
    """Stacked per-donor composition plus per-population proportions by age bin."""
    cell_types = list(cell_types or prop.columns)
    bins = list(pd.Categorical(donor_meta.age_bin).categories)
    order = sorted(prop.index, key=lambda d: (bins.index(str(donor_meta.loc[d, "age_bin"])), str(d)))
    pr = prop.loc[order, cell_types]
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, 20))[:2 * len(cell_types):2]
    fig = plt.figure(figsize=(7.4, 6.0))
    gsA = fig.add_gridspec(1, 1, left=0.085, right=0.985, top=0.93, bottom=0.57)
    gsB = fig.add_gridspec(1, len(highlight), left=0.085, right=0.985, top=0.36,
                           bottom=0.10, wspace=0.42)
    ax = fig.add_subplot(gsA[0])
    bottom = np.zeros(len(pr))
    for i, ct in enumerate(cell_types):
        ax.bar(range(len(pr)), pr[ct].values, bottom=bottom, color=colors[i], width=0.82,
               label=ct, linewidth=0)
        bottom += pr[ct].values
    ax.set_xlim(-0.6, len(pr) - 0.4); ax.set_ylim(0, 1)
    ax.set_xticks(range(len(pr))); ax.set_xticklabels(order, rotation=90, fontsize=5.5)
    ax.set_ylabel("fraction of recovered cells")
    ax.set_title("Cell-type composition per donor", pad=18)
    for ab in bins:
        xs = [k for k, d in enumerate(order) if str(donor_meta.loc[d, "age_bin"]) == ab]
        if not xs:
            continue
        ax.plot([min(xs) - 0.4, max(xs) + 0.4], [1.03, 1.03], color="0.25", lw=1.4, clip_on=False)
        ax.text(np.mean(xs), 1.05, ab, ha="center", va="bottom", fontsize=7)
    ax.legend(ncol=len(cell_types), frameon=False, fontsize=6.2, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), columnspacing=1.1, handlelength=1.0)
    rng = np.random.default_rng(0)
    for j, ct in enumerate(highlight):
        axj = fig.add_subplot(gsB[j])
        for k, ab in enumerate(bins):
            v = pr.loc[[d for d in order if str(donor_meta.loc[d, "age_bin"]) == ab], ct].values
            axj.scatter(np.full(len(v), k) + rng.uniform(-.13, .13, len(v)), v, s=11,
                        color=colors[cell_types.index(ct)], alpha=.85, linewidth=0)
            if len(v):
                axj.plot([k - .26, k + .26], [np.median(v)] * 2, color="0.15", lw=1.6)
        axj.set_xticks(range(len(bins))); axj.set_xticklabels(bins, rotation=45, ha="right", fontsize=6)
        axj.set_title(ct, fontsize=7); axj.set_ylim(bottom=0); axj.margins(y=0.12)
        if j == 0:
            axj.set_ylabel("proportion of donor's cells")
    fig.text(0.085, 0.005, "Proportions are relative: a rise in one population forces "
                           "a fall in the rest", ha="left", fontsize=7.5)
    return fig


def donor_score_figure(summary: pd.DataFrame, combos, ncol=3):
    """Donor-level score by age bin for selected (module, cell type) pairs."""
    bins = _bins(summary)
    nrow = int(np.ceil(len(combos) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.5 * ncol, 2.3 * nrow), squeeze=False)
    rng = np.random.default_rng(3)
    for ax, (mod, ct) in zip(axes.ravel(), combos):
        col = score_column(mod)
        d = summary[summary.cell_type == ct]
        ok, miss = d[d.eligible], d[~d.eligible]
        for k, ab in enumerate(bins):
            v = ok.loc[ok.age_bin.astype(str) == ab, col].dropna().values
            ax.scatter(np.full(len(v), k) + rng.uniform(-.14, .14, len(v)), v, s=14,
                       color="#2b6cb0", alpha=.85, linewidth=0)
            if len(v):
                ax.plot([k - .27, k + .27], [np.median(v)] * 2, color="0.1", lw=1.7)
            nm = int((miss.age_bin.astype(str) == ab).sum())
            if nm:
                ax.annotate(f"{nm} n.d.", (k, 0.02), xycoords=("data", "axes fraction"),
                            ha="center", fontsize=5.8, color="0.45")
        ax.set_xticks(range(len(bins))); ax.set_xticklabels(bins, rotation=45, ha="right", fontsize=6)
        ax.set_xlim(-0.5, len(bins) - 0.5); ax.margins(y=0.18)
        ax.set_title(f"{mod.split('_', 1)[1]}\nin {ct}", fontsize=7)
        ax.set_ylabel({"_UCell": "UCell score", "_logit": "logit IgA fraction"}.get(col[col.rfind("_"):], "log1p(mean CP10K)"), fontsize=7)
    for ax in axes.ravel()[len(combos):]:
        ax.axis("off")
    fig.text(0.5, 0.0, "n.d. = donor below the minimum cell threshold; recorded as missing, not zero",
             ha="center", fontsize=6.5, color="0.35")
    fig.tight_layout(rect=[0, 0.03, 1, 1]); fig.subplots_adjust(hspace=0.8, wspace=0.4)
    return fig


def effect_heatmap(omnibus: pd.DataFrame, contrasts: pd.DataFrame, ref_bin: str):
    """Standardised age-bin differences within cell types, with BH marks."""
    ctr = [c for c in contrasts.contrast.unique()]
    rows = sorted(omnibus.itertuples(index=False), key=lambda r: (r.module_id, r.cell_type))
    M = np.full((len(rows), len(ctr)), np.nan)
    for i, r in enumerate(rows):
        for j, c in enumerate(ctr):
            q = contrasts[(contrasts.module_id == r.module_id) &
                          (contrasts.cell_type == r.cell_type) & (contrasts.contrast == c)]
            if len(q):
                M[i, j] = q.std_effect.values[0]
    vmax = np.ceil(np.nanmax(np.abs(M)) * 10) / 10 if np.isfinite(M).any() else 1.0
    fig, ax = plt.subplots(figsize=(1.1 + 1.5 * len(ctr), 0.19 * len(rows) + 1.6))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    labels = []
    for i, r in enumerate(rows):
        sig = pd.notna(getattr(r, "p_adj_BH", np.nan)) and r.p_adj_BH < 0.05
        labels.append(f"{r.module_id.split('_', 1)[1][:26]} · {r.cell_type}"
                      + ("  (expl.)" if r.tier != "primary" else "") + ("  *" if sig else ""))
        if not r.estimable:
            ax.add_patch(mpl.patches.Rectangle((-0.5, i - 0.5), len(ctr), 1, facecolor="0.88",
                                               edgecolor="white", hatch="///", lw=0))
            ax.text((len(ctr) - 1) / 2, i, "not assessed (too few donors)", ha="center",
                    va="center", fontsize=5.6, color="0.35")
            continue
        for j in range(len(ctr)):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=5.8,
                        color="white" if abs(M[i, j]) > 0.62 * vmax else "0.15")
    ax.set_xticks(range(len(ctr))); ax.set_xticklabels([c.replace(" vs ", "\nvs ") for c in ctr], fontsize=6.5)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels, fontsize=6)
    for i, r in enumerate(rows):
        if r.tier != "primary":
            ax.get_yticklabels()[i].set_color("0.45")
    ax.set_title("Age-bin differences in module score within each cell type", fontsize=8.5, pad=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label(f"standardised difference vs {ref_bin} (estimate / residual SD)", fontsize=6.5)
    return fig


def forest_plot(contrasts: pd.DataFrame, combos, ref_bin: str, unit="UCell score"):
    """Effect estimates with 95% CIs for selected combinations."""
    sub = contrasts[contrasts.apply(lambda r: (r.module_id, r.cell_type) in set(combos), axis=1)].copy()
    sub["key"] = sub.module_id.str.split("_", n=1).str[1] + " · " + sub.cell_type
    ctr = list(sub.contrast.unique())
    cmap = dict(zip(ctr, plt.get_cmap("coolwarm")(np.linspace(0.05, 0.95, len(ctr)))))
    keys = sorted(sub.key.unique(), reverse=True)
    fig, ax = plt.subplots(figsize=(5.6, 0.55 * len(keys) + 1.6))
    offs = np.linspace(0.26, -0.26, len(ctr))
    for y, k in enumerate(keys):
        d = sub[sub.key == k]
        for off, c in zip(offs, ctr):
            q = d[d.contrast == c]
            if not len(q):
                continue
            ax.errorbar(q.estimate.values[0], y + off,
                        xerr=[[q.estimate.values[0] - q.ci_low.values[0]],
                              [q.ci_high.values[0] - q.estimate.values[0]]],
                        fmt="o", ms=4.2, lw=1.4, color=cmap[c], capsize=2.2,
                        label=c if y == 0 else None)
    ax.axvline(0, color="0.5", lw=0.9, zorder=0)
    ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys, fontsize=7)
    ax.set_xlabel(f"difference vs {ref_bin} (95% CI), {unit}")
    ax.margins(x=0.12, y=0.08)
    ax.legend(frameon=False, fontsize=6.5, loc="best")
    return fig


def pooled_vs_within_figure(pooled_res: pd.DataFrame, within_res: pd.DataFrame):
    """Side-by-side omnibus -log10 p for pooled vs within-cell-type analyses."""
    m = pooled_res.merge(within_res, on="module_id", suffixes=("_pooled", "_within"))
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    x = np.arange(len(m))
    ax.bar(x - 0.2, -np.log10(m.p_omnibus_pooled), width=0.4, label="pooled (ignores cell identity)",
           color="#c23b22")
    ax.bar(x + 0.2, -np.log10(m.p_omnibus_within), width=0.4, label="within cell type", color="#2b6cb0")
    ax.axhline(-np.log10(0.05), color="0.4", ls="--", lw=0.9)
    ax.set_xticks(x); ax.set_xticklabels([s.split("_", 1)[1] for s in m.module_id], rotation=30,
                                         ha="right", fontsize=6.5)
    ax.set_ylabel("$-\\log_{10}$ omnibus p")
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title("Pooling across epithelial subtypes can manufacture an age effect", fontsize=8.5)
    return fig
