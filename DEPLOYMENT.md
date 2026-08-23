# Deploying the Chordoma Atlas

## Pre-publication access controls (REMOVE AT PUBLICATION)

The manuscript is under review, so the deployed site is deliberately
uninformative to anyone who has not been given the passphrase.

**This is obscurity, not security.** The repo is public (GitHub Pages will not
serve a private repo on the free tier) and the site is static, so anyone who
opens developer tools or browses the repo can reach the data. These measures
stop *accidental* discovery — a search engine, a forwarded link, a guessed
URL. They would not stop someone who wants the data.

What is in place:

| Measure | File | Effect |
| --- | --- | --- |
| Passphrase gate | `site/gate.js` | Hides the page and blocks app boot until "chordoma" is entered; SHA-256 compared, verified to accept the passphrase and reject variants. Unlock persists for the browser session only. |
| Randomised data filenames | `portal/obfuscate_names.py`, `site/datamap.json` | Published files are opaque tokens (`d6114d3912a63.parquet`), so a URL or directory listing does not describe its contents. Tokens are a salted hash of the logical name, so they are stable across rebuilds and do not churn git history. |
| No identifying HTML | `site/index.html` | Title is "Research Data Portal"; the heading, description and citation are injected by JS after unlock, so the served HTML names neither the disease nor the study. |
| Crawler directives | `site/robots.txt`, `<meta robots>` | `Disallow: /` plus `noindex, nofollow, noarchive, nosnippet`. |
| CI leak checks | `.github/workflows/deploy-pages.yml` | The deploy **fails** if identifying terms appear in rendered HTML or if any data file is left untokenized. |

To restore the public version at publication: delete `site/gate.js` and its
`<script>` tag, delete `site/robots.txt` and the robots `<meta>` tags, inline
the heading/citation back into `index.html`, and drop the rename step from the
workflow (the app falls back to readable filenames when `datamap.json` is
absent).

If you need genuine access control before then, host on Cloudflare Pages
(which *does* serve private repos) behind Cloudflare Access — free for up to
50 users, and unauthenticated requests never reach the content. Note that
Cloudflare Pages caps single assets at 25 MiB, so six of the current data
files would have to move to R2.

---


The portal is a static site: DuckDB-WASM queries Parquet directly from the
browser. Data ships **in this repo** alongside the site, so updating the
atlas is a normal `git push` — no separate upload step.

## One-time setup

1. Create the GitHub repo and add it as a remote:
   ```bash
   git remote add origin https://github.com/<user>/<repo>.git
   ```
2. In the repo: **Settings → Pages → Build and deployment → Source:
   GitHub Actions**.
3. Push `main`. The workflow in `.github/workflows/deploy-pages.yml`
   assembles and publishes the site.

The site appears at `https://<user>.github.io/<repo>/`.

## Everyday workflow

Edit code in `site/`, or regenerate data into `portal/data/`, then:

```bash
git add -A && git commit -m "update" && git push
```

The workflow republishes automatically. It fails the build if any file
exceeds 100 MB or if a required data file is missing, so a broken deploy is
caught in CI rather than in the browser.

## Why the data is laid out this way

**GitHub Pages honours HTTP Range requests** (returns `HTTP 206` with a
correct `content-range`, including the suffix ranges DuckDB-WASM issues to
read a Parquet footer). A probe lookup therefore transfers ~1–2 MB, not a
whole file. Because data and page share an origin, no CORS setup is needed.

**Git's hard limit is 100 MB per file.** The combined methylation matrix
(`meth_discovery.parquet`, ~460 MB) exceeds it and cannot be pushed, so it is
gitignored and split by chromosome:

- `portal/split_methylation.py` writes `meth_chr1..22.parquet`
  (largest: chr1, ~45 MB), each sorted by genomic position with small row
  groups so range queries prune to the row groups they need.
- `meth_probe_index.parquet` maps `meth_site -> (chr, pos)` so the app knows
  which chromosome file a probe lives in.

Splitting by chromosome matches how the data is used: no view ever needs two
chromosomes at once.

Regenerate after re-running the ETL:

```bash
python portal/etl.py               # writes meth_discovery.parquet (local only)
python portal/split_methylation.py # writes the per-chromosome files
```

## Limits to watch

| Limit | Value | Current |
| --- | --- | --- |
| Per file (git) | 100 MB hard | 45 MB (chr1) |
| Published site | 1 GB | ~520 MB |
| Bandwidth | 100 GB/month soft | far below |

Note that every committed revision of a data file is stored permanently in
git history, and Parquet does not delta-compress — re-committing all 22
chromosome files repeatedly will grow the repo. If history becomes unwieldy,
move `portal/data/` to a separate data repo or to Hugging Face and point
`DATA_BASE` in `site/config.js` at it (Range + CORS verified there).

## Alternative hosts

Recorded so the evaluation does not have to be repeated:

| Host | Range | CORS | Mutable | Verdict |
| --- | --- | --- | --- | --- |
| GitHub Pages (this repo) | yes | same-origin | yes | **in use** |
| Hugging Face | yes | yes | yes | best fallback |
| Cloudflare R2 | yes | configurable | yes | if >1 GB |
| Zenodo | yes | yes | **no** (DOI-frozen) | final archival |
| GitHub Releases | yes | **no header** | yes | unusable in-browser |
| Git LFS on Pages | — | — | — | serves pointer text, not data |

At publication, deposit the frozen dataset in Zenodo for a DOI and cite that
alongside the live portal.
