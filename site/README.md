# Chordoma Atlas — static site (GitHub Pages stack)

Fully client-side portal: DuckDB-WASM range-reads Parquet over HTTP, all
statistics (Kruskal-Wallis, Spearman/Pearson, Kaplan-Meier, k-group log-rank,
Cox PH with Efron ties, maximally-selected cutpoint) run in the browser via
`stats.js` — validated to <=1e-4 relative error against scipy/lifelines on the
real cohort data.

## Local development

```bash
cd site
ln -s ../portal/data data     # already done; parquet files from portal/etl.py
python serve.py 8000          # dev server WITH HTTP Range support
# open http://localhost:8000
```

`python -m http.server` will NOT work — DuckDB-WASM needs Range requests.
`serve.py` implements them so local dev behaves exactly like production.

## Files

```
index.html   shell, theme (viridis header), subtype chips, nav
app.js       DuckDB-WASM boot + 4 pages (expression, survival, correlation, methylation)
stats.js     dependency-free statistics library (validated vs scipy/lifelines)
config.js    DATA_BASE: "./data" for dev; swap for production URL
serve.py     local dev server with byte-range support
```

## Production hosting — verified facts

The site (HTML+JS, <100 KB) goes on **GitHub Pages** for free. The Parquet
files (~520 MB total) need a host that supports HTTP Range + CORS:

| Option | Range? | CORS? | Limits / cost | Verdict |
|---|---|---|---|---|
| **Zenodo** | YES — verified 206 Partial Content on per-file URLs | YES — `access-control-allow-origin: *` (verified) | 50 GB/record, free, DOI included. Files do NOT need to be tarred — each file in a record has its own URL: `zenodo.org/api/records/<ID>/files/<name>/content` | **Recommended**: free + DOI + range-reads work |
| **GitHub Releases** | YES (verified 206) | Effectively usable via redirect to S3-backed CDN | 2 GB per file (fine), unlimited total, free | Good backup; no DOI |
| **Cloudflare R2** | YES | configurable | Free tier: 10 GB storage, **zero egress fees**; ~$0.015/GB/mo after | Best latency; small cost |
| **AWS S3** | YES | configurable | ~$0.023/GB/mo storage + **$0.09/GB egress** — egress is the trap: 1000 heavy users x 50 MB = $4.50/mo, but a viral month could be $50+ | Works; R2 strictly cheaper |
| GitHub Pages itself | — | — | 1 GB site limit, soft 100 GB/mo bandwidth — too small for the methylation file + bandwidth risk | Site only, not data |

Per-query transfer: row groups are 2048 rows (~1.2 MB compressed), so a
single-gene expression query moves ~1-2 MB; a methylation probe query on the
453 MB file moves the 4 MB footer once (cached) + ~1.2 MB per probe.

### Zenodo notes
- Upload the parquet files as individual files in one record (not a tarball).
  Batch-download convenience is preserved (Zenodo offers record-level zip).
- New record version = new file URLs; pin `config.js` to a record ID.
- Zenodo rate-limits (~60-100 req/min/IP) are fine for this access pattern
  (a handful of range requests per user action).

## Deploy checklist

1. `python portal/etl.py` -> parquet files
2. Create Zenodo record, upload all `portal/data/*.parquet`, publish -> DOI
3. Edit `config.js`: `DATA_BASE = "https://zenodo.org/api/records/<ID>/files"`,
   `ZENODO_MODE = true` (appends `/content` to file URLs)
4. Push `site/` to the repo, enable GitHub Pages
5. Put the DOI + citation in `index.html` (`#citation`)
