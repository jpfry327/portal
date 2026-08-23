"""Data-access layer for the chordoma portal.

All queries hit Parquet files lazily through DuckDB, so nothing large is
resident in memory: a single-gene lookup reads only the row groups that
contain that gene (~ms even on the 756k-probe methylation matrix).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import duckdb
import pandas as pd

DATA = Path(__file__).resolve().parent / "data"

RNA_COHORTS = ["discovery", "zhang", "inhouse", "bai"]
COHORT_LABELS = {
    "discovery": "Discovery (n=135, multi-omic)",
    "zhang": "Zhang proteogenomic (validation)",
    "inhouse": "In-house (validation)",
    "bai": "Bai (validation)",
}
# viridis(4), matching the manuscript main figures (C1 purple ... C4 yellow)
SUBTYPE_COLORS = {
    "C1": "#440154",   # hypermetabolic
    "C2": "#31688E",   # neural proliferative
    "C3": "#35B779",   # immunogenic
    "C4": "#FDE725",   # LYST-activated
}
SUBTYPE_NAMES = {
    "C1": "C1 hypermetabolic",
    "C2": "C2 neural proliferative",
    "C3": "C3 immunogenic",
    "C4": "C4 LYST-activated",
}


def _con() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("SET enable_progress_bar = false")
    return con


def pq(name: str) -> str:
    return str(DATA / f"{name}.parquet")


@lru_cache(maxsize=1)
def samples() -> pd.DataFrame:
    return _con().execute(f"SELECT * FROM read_parquet('{pq('samples')}')").df()


@lru_cache(maxsize=1)
def gene_list() -> list[str]:
    q = f"SELECT DISTINCT gene FROM read_parquet('{pq('gene_index')}') ORDER BY gene"
    return _con().execute(q).df()["gene"].tolist()


@lru_cache(maxsize=1)
def gene_availability() -> pd.DataFrame:
    """gene x (cohort, modality) presence."""
    return _con().execute(f"SELECT * FROM read_parquet('{pq('gene_index')}')").df()


def _melt_row(df: pd.DataFrame, id_col: str, value_name: str) -> pd.DataFrame:
    long = df.melt(id_vars=[id_col], var_name="sample_id", value_name=value_name)
    return long.drop(columns=[id_col]).dropna(subset=[value_name])


@lru_cache(maxsize=256)
def expression(gene: str, cohort: str) -> pd.DataFrame:
    """Long-form expression for one gene in one cohort, joined to subtype."""
    f = pq(f"expr_{cohort}")
    df = _con().execute(f"SELECT * FROM read_parquet('{f}') WHERE gene = ?", [gene]).df()
    if df.empty:
        return pd.DataFrame(columns=["sample_id", "value", "subtype", "cohort"])
    long = _melt_row(df, "gene", "value")
    s = samples()[samples().cohort == cohort][["sample_id", "subtype", "subtype_prob"]]
    out = long.merge(s, on="sample_id", how="left")
    out["cohort"] = cohort
    return out


@lru_cache(maxsize=256)
def protein(gene: str) -> pd.DataFrame:
    df = _con().execute(
        f"SELECT * FROM read_parquet('{pq('protein_zhang')}') WHERE gene_symbol = ?", [gene]
    ).df()
    if df.empty:
        return pd.DataFrame(columns=["sample_id", "value", "subtype"])
    long = _melt_row(df, "gene_symbol", "value")
    s = samples()[samples().cohort == "zhang"][["sample_id", "subtype"]]
    out = long.merge(s, on="sample_id", how="left")
    out["cohort"] = "zhang"
    return out


@lru_cache(maxsize=64)
def methylation(probe: str) -> pd.DataFrame:
    df = _con().execute(
        f"SELECT * FROM read_parquet('{pq('meth_discovery')}') WHERE meth_site = ?", [probe]
    ).df()
    if df.empty:
        return pd.DataFrame(columns=["sample_id", "value", "subtype"])
    long = _melt_row(df, "meth_site", "value")
    s = samples()[samples().cohort == "discovery"][["sample_id", "subtype"]]
    return long.merge(s, on="sample_id", how="left")


def probe_exists(probe: str) -> bool:
    n = _con().execute(
        f"SELECT COUNT(*) FROM read_parquet('{pq('meth_discovery')}') WHERE meth_site = ?",
        [probe],
    ).fetchone()[0]
    return n > 0


@lru_cache(maxsize=4)
def clinical(cohort: str = "discovery") -> pd.DataFrame:
    """Clinical + subtype for 'discovery' (corrected OS, RFS) or 'zhang' (OS, PFS)."""
    df = _con().execute(f"SELECT * FROM read_parquet('{pq(f'clinical_{cohort}')}')").df()
    s = samples()[samples().cohort == cohort][["sample_id", "subtype"]]
    return df.merge(s, on="sample_id", how="inner")


@lru_cache(maxsize=1)
def cell_proportions() -> pd.DataFrame:
    df = _con().execute(
        f"SELECT * FROM read_parquet('{pq('cell_proportions_discovery')}')"
    ).df()
    s = samples()[samples().cohort == "discovery"][["sample_id", "subtype"]]
    return df.merge(s, on="sample_id", how="left")


def pair_expression(gene_a: str, gene_b: str, cohort: str) -> pd.DataFrame:
    """Two genes side by side for correlation, same cohort."""
    a = expression(gene_a, cohort).rename(columns={"value": "x"})[["sample_id", "x", "subtype"]]
    b = expression(gene_b, cohort).rename(columns={"value": "y"})[["sample_id", "y"]]
    return a.merge(b, on="sample_id", how="inner")


@lru_cache(maxsize=1)
def cellline_meta() -> pd.DataFrame:
    return _con().execute(f"SELECT * FROM read_parquet('{pq('cellline_meta')}')").df()


@lru_cache(maxsize=1)
def cellline_genes() -> list[str]:
    q = f"SELECT gene_name FROM read_parquet('{pq('tpm_cellline')}') ORDER BY gene_name"
    return _con().execute(q).df()["gene_name"].tolist()


@lru_cache(maxsize=128)
def cellline_tpm(gene: str) -> pd.DataFrame:
    """Long-form TPM for one gene across cell line / PDX replicates."""
    df = _con().execute(
        f"SELECT * FROM read_parquet('{pq('tpm_cellline')}') WHERE gene_name = ?", [gene]
    ).df()
    if df.empty:
        return pd.DataFrame(columns=["replicate", "tpm", "sample_id", "sample_type"])
    long = df.melt(id_vars=["gene_name"], var_name="replicate", value_name="tpm")
    long = long.drop(columns=["gene_name"]).dropna(subset=["tpm"])
    return long.merge(cellline_meta(), on="replicate", how="inner")


def rna_protein_pairs(gene: str) -> pd.DataFrame:
    """RNA vs protein for one gene (Zhang cohort)."""
    r = expression(gene, "zhang").rename(columns={"value": "rna"})[["sample_id", "rna", "subtype"]]
    p = protein(gene).rename(columns={"value": "prot"})[["sample_id", "prot"]]
    return r.merge(p, on="sample_id", how="inner")
