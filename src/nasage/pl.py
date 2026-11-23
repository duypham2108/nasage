"""Figures from the results stored by ``nas.tl`` (each returns a matplotlib Figure)."""
from __future__ import annotations

from . import plots as _plots
from .tl import _need


def composition(adata, **kw):
    """Per-donor composition, and proportions of key populations by age group."""
    prop = _need(adata, "proportions", "composition")
    meta = _need(adata, "design", "donor_summary").set_index("donor")
    meta.index = meta.index.astype(str)
    prop = prop.copy()
    prop.index = prop.index.astype(str)
    hl = [c for c in kw.pop("highlight", ("ciliated", "secretory", "Tcell", "dendritic"))
          if c in prop.columns]
    return _plots.composition_figure(prop, meta, highlight=tuple(hl), **kw)


def effect_heatmap(adata):
    """Standardised age-group differences within cell types, BH-significant marked."""
    st = _need(adata, "params", "donor_summary")
    return _plots.effect_heatmap(_need(adata, "omnibus", "age_effects"),
                                 _need(adata, "contrasts", "age_effects"), st["reference"])


def donor_scores(adata, combos, ncol=3):
    """Donor-level values by age group for [(module_id, cell_type), ...]."""
    return _plots.donor_score_figure(_need(adata, "donor_summary", "donor_summary"), combos,
                                     ncol=ncol)


def forest(adata, combos, unit="UCell score"):
    """Estimates with 95% CIs vs the reference group for [(module_id, cell_type), ...]."""
    st = _need(adata, "params", "donor_summary")
    return _plots.forest_plot(_need(adata, "contrasts", "age_effects"), combos, st["reference"],
                              unit=unit)
