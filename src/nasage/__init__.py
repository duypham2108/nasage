"""
nasage: age-associated nasal mucosal immune-state profiling from
single-cell RNA-seq, scanpy-style.

    import scanpy as sc
    import nasage as nas

    adata = sc.read_h5ad("nasal.h5ad")
    nas.pp.filter_cells(adata, celltype_key="population", layer="counts")
    nas.tl.score_modules(adata, layer="counts")
    nas.tl.donor_summary(adata, donor_key="donor", celltype_key="population",
                         age_key="age_group", batch_key="batch", sex_key="sex")
    nas.tl.age_effects(adata, reference="18-39")

Results live in adata.uns["nasage"]. The package does NOT produce a
composite "vaccine readiness" score: see docs/INTERPRETATION.md.
"""
__version__ = "0.2.0"

from . import pl, pp, tl
from .modules import drop_hla, load_applicability, load_modules
from .report import write_report

__all__ = ["pp", "tl", "pl", "load_modules", "load_applicability", "drop_hla",
           "write_report", "__version__"]
