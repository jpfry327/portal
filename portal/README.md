# Chordoma Multi-omic Atlas — companion portal

Interactive companion to the chordoma multi-omic subtyping manuscript. Users can
query gene expression (4 cohorts), proteomics (Zhang), and DNA methylation
(discovery, 756k probes) by molecular subtype (C1–C4), run gene–gene and
RNA–protein correlations, and draw Kaplan–Meier curves stratified by subtype or
by gene expression.

## Layout

```
portal/
  etl.py           TSV -> Parquet + DuckDB conversion (run once, rerun when data changes)
  portal_data.py   data-access layer (DuckDB lazy queries over Parquet)
  app.py           Streamlit app (5 pages)
  data/            generated Parquet files + chordoma.duckdb (518 MB total; NOT in git)
```

## Quick start

```bash
# 1. environment
pip install duckdb pyarrow pandas streamlit plotly scipy lifelines statsmodels

# 2. build the columnar store from ../data/*.tsv  (~1-2 min, methylation dominates)
python etl.py

# 3. run
streamlit run app.py
```

## Why Parquet + DuckDB

* Raw TSVs: **2.5 GB** -> Parquet (zstd, float32): **518 MB** (methylation 2.0 GB -> 453 MB).
* Nothing is loaded into RAM up-front. Every query is a DuckDB predicate-pushdown
  scan over Parquet row groups; a single-gene or single-probe lookup returns in
  **~30 ms**, including on the 756,448-probe methylation matrix.
* `data/chordoma.duckdb` exposes all tables as SQL views, so power users (or a
  future API) can query the same files directly:
  `duckdb data/chordoma.duckdb "SELECT * FROM samples LIMIT 5"`.

## Data model

| table | contents |
|---|---|
| `samples` | sample_id, cohort, subtype (assigned or classifier-predicted), prediction probability. Classifier calls with top probability <= 0.6 are excluded as low-confidence (6 Zhang, 1 Bai) — the model only assigns a subtype above 0.6 (paper methods) |
| `expr_discovery/zhang/inhouse/bai` | wide gene x sample RNA matrices (log scale) |
| `protein_zhang` | wide protein matrix (corrected log abundance) |
| `meth_discovery` | 756k probes x 161 arrays, M-values (normalized, batch-corrected) |
| `clinical_discovery` | corrected OS (time_os_nr/status_os_fr) + RFS, age, sex, histology, location |
| `clinical_zhang` | OS + PFS from the Zhang proteogenomics clinical table |
| `cell_proportions_discovery` | BayesPrism cell-type fractions |
| `embeddings` | 16-dim GNN embeddings (zhang & bai spaces) |
| `gene_index` | gene -> cohort/modality availability |

Zhang RNA columns (HRR accessions) are mapped to the proteomics sample IDs via
the HRA run-title metadata, so RNA, protein, and subtype join on one ID.

## Deployment options (next step)

1. **Streamlit Community Cloud / Posit Connect** — simplest; push repo, mount
   `portal/data` (518 MB fits comfortably).
2. **Docker on a lab VM** — `python:3.13-slim` + pip deps, ~1 GB image with data
   baked in; put nginx/caddy in front for TLS + a citable domain name.
3. **Static-first upgrade path** — because everything is Parquet, the same files
   can back a serverless portal later (DuckDB-WASM in the browser reading
   range-requested Parquet from object storage) with zero backend to maintain.

## Survival analysis conventions

- Cohorts: discovery (corrected OS = `time_os_nr`/`status_os_fr`, plus RFS) and
  Zhang (OS, PFS). In-house/Bai clinical data are too sparse and are excluded.
- Per-gene analysis: the KM plot is drawn at the **optimal cutpoint** (maximally
  selected log-rank statistic, 10th–90th percentile scan, >=10 per arm), but the
  **reported** statistics are from Cox PH regression on continuous expression —
  likelihood-ratio p and HR per 1 SD (comparable across genes). The cutpoint
  log-rank p is biased by cutpoint selection and is shown only with a warning.

## Style

Subtype colors are viridis(4) matching the manuscript: C1 #440154 (purple),
C2 #31688E, C3 #35B779, C4 #FDE725 (yellow). Expression/methylation panels are
violins with black outlines and an inner grey box, Kruskal-Wallis p annotated
top-left, `simple_white` template.

## Caveats / TODO

- Discovery methylation currently includes replicate arrays (161 columns for 135
  tumors); replicates are shown as-is.
- Protein rows without an RNA-mapped sample (27 Zhang samples) show protein-only.
- Consider adding: probe->gene annotation (450k/EPIC manifest) so users can query
  methylation by gene; batch download endpoints; subtype-vs-subtype differential
  tables precomputed.
