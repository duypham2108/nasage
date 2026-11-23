"""Gene modules, applicability matrix and symbol harmonisation."""
from __future__ import annotations

import json
from importlib import resources
from pathlib import Path

import pandas as pd

#: Legacy -> current HGNC symbols, each verified against mygene.info.
#: Extend with care; every entry should be checked, not guessed.
VERIFIED_ALIASES = {
    "DDX58": "RIGI",      # RNA sensor RIG-I
    "PIFO": "CIMAP3",     # ciliary microtubule associated protein 3
    "MB21D1": "CGAS",
    "TMEM173": "STING1",
}


def _resource(name: str) -> Path:
    # Resources ship inside the package (see [tool.setuptools.package-data]), so
    # this works for both editable and regular installs.
    p = Path(str(resources.files("nasage") / "resources" / name))
    if not p.exists():
        raise FileNotFoundError(f"packaged resource not found: {p}")
    return p


def load_modules(path: str | Path | None = None) -> dict:
    """Return {module_id: {genes, name, definition, interpretation, ...}}."""
    with open(path or _resource("gene_modules.json")) as fh:
        return json.load(fh)


def load_applicability(path: str | Path | None = None) -> dict:
    """Return {module_id: {cell_type: 1.0 primary | 0.5 exploratory | 0.0 n/a}}."""
    with open(path or _resource("module_celltype_applicability.json")) as fh:
        return json.load(fh)


def load_provenance(path: str | Path | None = None) -> pd.DataFrame:
    return pd.read_csv(path or _resource("gene_set_provenance_modules.csv"))


def harmonise_symbols(var_names, modules: dict, online: bool = False,
                      email: str | None = None) -> tuple[dict, pd.DataFrame]:
    """Map module genes onto the symbols present in a dataset.

    Exact matches first, then the verified alias table. With ``online=True`` the
    remaining symbols are resolved through mygene.info (one batched request);
    this requires network access and is off by default so that a run is
    reproducible offline.

    Returns (modules_mapped, report_dataframe).
    """
    present = set(map(str, var_names))
    rev = {v: k for k, v in VERIFIED_ALIASES.items()}
    out, rows = {}, []
    unresolved = []
    for mod, meta in modules.items():
        genes = meta["genes"] if isinstance(meta, dict) else meta
        mapped = []
        for g in genes:
            if g in present:
                mapped.append(g); status = "exact"
            elif VERIFIED_ALIASES.get(g) in present:
                mapped.append(VERIFIED_ALIASES[g]); status = f"alias->{VERIFIED_ALIASES[g]}"
            elif rev.get(g) in present:
                mapped.append(rev[g]); status = f"alias->{rev[g]}"
            else:
                status = "absent"; unresolved.append(g)
            rows.append(dict(module_id=mod, gene=g, status=status))
        out[mod] = mapped
    if online and unresolved:
        extra = _resolve_online(sorted(set(unresolved)), present, email)
        for mod, meta in modules.items():
            genes = meta["genes"] if isinstance(meta, dict) else meta
            for g in genes:
                if g in extra and extra[g] not in out[mod]:
                    out[mod].append(extra[g])
                    for r in rows:
                        if r["module_id"] == mod and r["gene"] == g:
                            r["status"] = f"mygene->{extra[g]}"
    return out, pd.DataFrame(rows)


def _resolve_online(symbols, present, email=None):
    import json as _json
    import urllib.parse
    import urllib.request
    body = urllib.parse.urlencode({"q": ",".join(symbols), "scopes": "symbol,alias",
                                   "fields": "symbol", "species": "human"}).encode()
    req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
    hits = _json.load(urllib.request.urlopen(req, timeout=60))
    out = {}
    for h in hits:
        s = h.get("symbol")
        if s in present:
            out[h.get("query")] = s
    return out


def drop_hla(modules: dict) -> dict:
    """Copy of `modules` without polymorphic HLA genes (HLA-*), for sensitivity.

    HLA alleles align to the reference genome with allele-dependent efficiency,
    so HLA-based scores partly track a donor's HLA type and ancestry. Rerun the
    analysis with this and check that age effects in M06/M07 survive.
    B2M is kept: it is not polymorphic.
    """
    out = {}
    for mod, meta in modules.items():
        if isinstance(meta, dict):
            out[mod] = {**meta, "genes": [g for g in meta["genes"] if not g.startswith("HLA-")]}
        else:
            out[mod] = [g for g in meta if not g.startswith("HLA-")]
    return out


def module_report(modules_mapped: dict, modules: dict) -> pd.DataFrame:
    """How much of each module survived mapping onto the dataset."""
    rows = []
    for mod, genes in modules_mapped.items():
        defined = modules[mod]["genes"] if isinstance(modules[mod], dict) else modules[mod]
        rows.append(dict(module_id=mod, n_defined=len(defined), n_found=len(genes),
                         frac_found=len(genes) / max(len(defined), 1),
                         missing=";".join(sorted(set(defined) - set(genes)))))
    df = pd.DataFrame(rows)
    df["usable_as_module"] = df.n_found >= 3
    return df
