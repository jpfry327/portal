"""ETL: convert chordoma cohort TSVs into Parquet + DuckDB for the portal.

Outputs (portal/data/):
  samples.parquet                      harmonized sample metadata (all cohorts)
  expr_<cohort>.parquet                wide RNA matrix, gene-sorted, float32
  protein_zhang.parquet                wide protein matrix, gene-sorted, float32
  meth_discovery.parquet               wide methylation M-values, probe-sorted, float32
  clinical_discovery.parquet           cleaned discovery clinical (OS/RFS)
  cell_proportions_discovery.parquet   BayesPrism cell-type fractions
  embeddings.parquet                   16-dim GNN embeddings (zhang + bai spaces)
  gene_index.parquet                   gene -> which cohorts/modalities carry it
  chordoma.duckdb                      views over all parquet files

Run:  python portal/etl.py
"""
from __future__ import annotations

import re
from pathlib import Path

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent   # ChordomaPortal/
RAW = ROOT / "data"
OUT = Path(__file__).resolve().parent / "data"
OUT.mkdir(exist_ok=True)

COHORTS_RNA = {
    "discovery": RAW / "discovery/rna_seq_logrpkm.tsv",
    "zhang": RAW / "validation_zhang/proteogenomics_rna.tsv",
    "inhouse": RAW / "validation_inhouse/inhouse_rna.tsv",
    "bai": RAW / "validation_bai/bai_rna.tsv",
}

con = duckdb.connect()


def wide_tsv_to_parquet(src: Path, dest: Path, id_col: str, rename: dict[str, str] | None = None):
    """Read a wide numeric TSV, drop empty cols, cast float32, sort by id, write parquet."""
    rel = con.read_csv(str(src), sep="\t", header=True)
    cols = rel.columns
    keep = [c for c in cols if c.strip() != ""]
    sample_cols = [c for c in keep if c != id_col]
    sel = [f'"{id_col}"']
    for c in sample_cols:
        new = (rename or {}).get(c, c)
        sel.append(f'CAST("{c}" AS FLOAT) AS "{new}"')
    q = f'SELECT {", ".join(sel)} FROM rel WHERE "{id_col}" IS NOT NULL ORDER BY "{id_col}"'
    con.execute(
        f"COPY ({q}) TO '{dest}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 20000)"
    )
    n = con.execute(f"SELECT COUNT(*) FROM read_parquet('{dest}')").fetchone()[0]
    print(f"  {dest.name}: {n} rows x {len(sample_cols)} samples "
          f"({dest.stat().st_size/1e6:.1f} MB)")
    return [str((rename or {}).get(c, c)) for c in sample_cols]


# ---------------------------------------------------------------- sample map
print("== sample metadata ==")

# Zhang: map HRR accession (RNA cols) -> canonical sample id (protein col naming)
zmeta = pd.read_csv(RAW / "validation_zhang/proteogenomics_metadata.csv", encoding="utf-8-sig")
acc2sample = {
    a: re.sub(r"_RNA$", "", t)
    for a, t in zip(zmeta["Accession"], zmeta["Run title"])
    if isinstance(t, str) and t.endswith("_RNA")
}

frames = []

disc_meta = pd.read_csv(RAW / "discovery/discovery_samples_meta.tsv", sep="\t")
frames.append(pd.DataFrame({
    "sample_id": disc_meta["chordoma_id"],
    "cohort": "discovery",
    "subtype": disc_meta["subtype"],
    "subtype_source": "assigned (multi-omic clustering)",
    "subtype_prob": pd.NA,
}))

# Classifier calls below MIN_PRED_PROB are low-confidence and excluded: the
# model only assigns a subtype at > 0.6 probability (paper methods). Drops 6
# Zhang and 1 Bai sample.
MIN_PRED_PROB = 0.6

for cohort, pred_path, id_col in [
    ("zhang", RAW / "validation_zhang/proteogenomics_rna_predictions.tsv", "sample"),
    ("inhouse", RAW / "validation_inhouse/inhouse_rna_methylation_predictions.tsv", "chordoma_id"),
    ("bai", RAW / "validation_bai/bai_rna_predictions.tsv", "sample"),
]:
    pred = pd.read_csv(pred_path, sep="\t")
    top = pred[["C1", "C2", "C3", "C4"]].max(axis=1)
    n_low = (top <= MIN_PRED_PROB).sum()
    if n_low:
        print(f"  {cohort}: excluding {n_low} low-confidence prediction(s) "
              f"(top prob <= {MIN_PRED_PROB})")
    pred = pred[top > MIN_PRED_PROB]
    sid = pred[id_col].map(acc2sample).fillna(pred[id_col]) if cohort == "zhang" else pred[id_col]
    frames.append(pd.DataFrame({
        "sample_id": sid,
        "cohort": cohort,
        "subtype": pred["prediction"],
        "subtype_source": "predicted (classifier)",
        "subtype_prob": pred[["C1", "C2", "C3", "C4"]].max(axis=1),
    }))

samples = pd.concat(frames, ignore_index=True)
samples["subtype_prob"] = pd.to_numeric(samples["subtype_prob"], errors="coerce")
samples.to_parquet(OUT / "samples.parquet", index=False)
print(f"  samples.parquet: {len(samples)} samples")
print(samples.groupby(["cohort", "subtype"]).size().unstack(fill_value=0))

# ---------------------------------------------------------------- expression
print("== RNA expression ==")
gene_rows = []
for cohort, src in COHORTS_RNA.items():
    rename = acc2sample if cohort == "zhang" else None
    wide_tsv_to_parquet(src, OUT / f"expr_{cohort}.parquet", "gene", rename)
    g = con.execute(
        f"SELECT DISTINCT gene FROM read_parquet('{OUT}/expr_{cohort}.parquet')"
    ).df()["gene"]
    gene_rows.append(pd.DataFrame({"gene": g, "cohort": cohort, "modality": "rna"}))

# ---------------------------------------------------------------- proteomics
print("== proteomics (zhang) ==")
wide_tsv_to_parquet(
    RAW / "validation_zhang/proteogenomics_prot_corrected.tsv",
    OUT / "protein_zhang.parquet", "gene_symbol",
)
g = con.execute(
    f"SELECT DISTINCT gene_symbol AS gene FROM read_parquet('{OUT}/protein_zhang.parquet')"
).df()["gene"]
gene_rows.append(pd.DataFrame({"gene": g, "cohort": "zhang", "modality": "protein"}))

pd.concat(gene_rows, ignore_index=True).to_parquet(OUT / "gene_index.parquet", index=False)

# ---------------------------------------------------------------- methylation
print("== methylation (discovery, ~756k probes — takes a minute) ==")
wide_tsv_to_parquet(
    RAW / "discovery/methylation_mvalues_normalized_batchcorrected.tsv",
    OUT / "meth_discovery.parquet", "meth_site",
)

# ---------------------------------------------------------------- clinical
print("== clinical ==")
# Discovery: time_os_nr / status_os_fr are the CORRECTED overall-survival columns.
clin = pd.read_csv(RAW / "discovery/clinical_data_os.tsv", sep="\t")
keep = [
    "chordoma_id", "age_at_initial_diagnosis_resection_years", "gender_0_female_1_male",
    "histology", "location_1_clival_2_petroclival_3_c_spine_4_seeding_5_sacral",
    "time_os_nr", "status_os_fr", "time_rfs", "status_rfs",
]
clin = clin[[c for c in keep if c in clin.columns]].rename(columns={
    "chordoma_id": "sample_id",
    "age_at_initial_diagnosis_resection_years": "age",
    "gender_0_female_1_male": "sex_male",
    "location_1_clival_2_petroclival_3_c_spine_4_seeding_5_sacral": "location_code",
    "time_os_nr": "time_os",       # corrected OS
    "status_os_fr": "status_os",   # corrected OS event indicator
})
for c in clin.columns:
    if c not in ("sample_id", "histology"):
        clin[c] = pd.to_numeric(clin[c], errors="coerce")
clin.to_parquet(OUT / "clinical_discovery.parquet", index=False)
print(f"  clinical_discovery.parquet: {len(clin)} rows "
      f"(OS n={clin.time_os.notna().sum()}, events={int(clin.status_os.sum(skipna=True))})")

# Zhang: OS and PFS from the processed proteogenomics clinical table.
# tumor_sample_type: 1=primary, 2=recurrent w/o RT, 3=recurrent after RT.
# NOTE: recurrent samples inherit the matched primary's time-to-death (time
# zero = diagnosis, not this resection), so OS analyses must be restricted to
# tumor_sample_type == 1 (as in the paper's methods).
zclin = pd.read_csv(RAW / "validation_zhang/proteogenomics_clinical_data_processed.txt", sep="\t")
zclin = zclin.rename(columns={
    "Sample.ID": "sample_id",
    "OS..day.": "time_os", "OS": "status_os",
    "PFS..day.": "time_pfs", "PFS": "status_pfs",
    "Tumor.sample.type.based.on.clinical.therapy.and.recurrent.status..1..Primary.tumor..2..Recurrent.tumor.without.radiotherapy..3..Recurrent.tumor.after.radiotherapy.":
        "tumor_sample_type",
})[["sample_id", "time_os", "status_os", "time_pfs", "status_pfs", "tumor_sample_type"]]
for c in zclin.columns:
    if c != "sample_id":
        zclin[c] = pd.to_numeric(zclin[c], errors="coerce")
zclin.to_parquet(OUT / "clinical_zhang.parquet", index=False)
print(f"  clinical_zhang.parquet: {len(zclin)} rows "
      f"(OS n={zclin.time_os.notna().sum()}, PFS n={zclin.time_pfs.notna().sum()}, "
      f"primary n={(zclin.tumor_sample_type == 1).sum()})")

# ------------------------------------------------------- cell proportions
props = pd.read_csv(RAW / "discovery/bulk_proportions_from_bayesprism.tsv", sep="\t")
props = props.rename(columns={"chordoma_id": "sample_id"})
props.to_parquet(OUT / "cell_proportions_discovery.parquet", index=False)
print(f"  cell_proportions_discovery.parquet: {len(props)} rows")

# ---------------------------------------------------------------- embeddings
emb_frames = []
for space, path in [("zhang", RAW / "validation_zhang/proteogenomics_rna_embeddings.csv"),
                    ("bai", RAW / "validation_bai/bai_rna_imbeddings.csv")]:
    e = pd.read_csv(path)
    e = e.rename(columns={e.columns[0]: "sample_id"})
    e["sample_id"] = e["sample_id"].map(lambda s: acc2sample.get(s, s))
    e["space"] = space
    emb_frames.append(e)
emb = pd.concat(emb_frames, ignore_index=True)
emb.columns = [f"d{c}" if str(c).isdigit() else c for c in emb.columns]
emb.to_parquet(OUT / "embeddings.parquet", index=False)
print(f"  embeddings.parquet: {len(emb)} rows")

# ---------------------------------------------------------------- cell lines
print("== cell lines (CCL Tracker) ==")
wide_tsv_to_parquet(RAW / "cell_line/tpm.tsv", OUT / "tpm_cellline.parquet", "gene_name")
clmeta = pd.read_csv(RAW / "cell_line/meta.tsv", sep="\t")
clmeta.to_parquet(OUT / "cellline_meta.parquet", index=False)
print(f"  cellline_meta.parquet: {len(clmeta)} replicates, "
      f"{clmeta.sample_id.nunique()} samples "
      f"({(clmeta.sample_type == 'Cell Line').sum()} CL / {(clmeta.sample_type == 'PDX').sum()} PDX reps)")

# ---------------------------------------------------------------- duckdb views
print("== chordoma.duckdb ==")
dbpath = OUT / "chordoma.duckdb"
dbpath.unlink(missing_ok=True)
db = duckdb.connect(str(dbpath))
views = {
    "samples": "samples.parquet",
    "expr_discovery": "expr_discovery.parquet",
    "expr_zhang": "expr_zhang.parquet",
    "expr_inhouse": "expr_inhouse.parquet",
    "expr_bai": "expr_bai.parquet",
    "protein_zhang": "protein_zhang.parquet",
    "meth_discovery": "meth_discovery.parquet",
    "clinical_discovery": "clinical_discovery.parquet",
    "clinical_zhang": "clinical_zhang.parquet",
    "cell_proportions_discovery": "cell_proportions_discovery.parquet",
    "embeddings": "embeddings.parquet",
    "gene_index": "gene_index.parquet",
    "tpm_cellline": "tpm_cellline.parquet",
    "cellline_meta": "cellline_meta.parquet",
}
if (OUT / "meth_probe_index.parquet").exists():
    views["meth_probe_index"] = "meth_probe_index.parquet"
for _c in range(1, 23):
    if (OUT / f"meth_chr{_c}.parquet").exists():
        views[f"meth_chr{_c}"] = f"meth_chr{_c}.parquet"

for name, f in views.items():
    db.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{OUT / f}')")
db.close()
print("done. total size:",
      f"{sum(p.stat().st_size for p in OUT.glob('*.parquet'))/1e6:.0f} MB parquet")
