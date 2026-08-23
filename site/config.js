/* Data location.
 *
 * "./data" works for BOTH local dev and the GitHub Pages deployment: the
 * deploy workflow copies portal/data/*.parquet into the published site, so
 * data is same-origin with the page. Same-origin means no CORS to configure,
 * and GitHub Pages honours HTTP Range requests (verified: HTTP 206 with a
 * correct content-range, including the suffix ranges DuckDB-WASM uses to read
 * a Parquet footer), so the browser fetches only the row groups it needs.
 *
 * Alternative hosts, if the data ever outgrows the 1 GB Pages site limit:
 *   Hugging Face (mutable, git-based; Range + CORS verified on both hops):
 *     "https://huggingface.co/datasets/<user>/<repo>/resolve/main"
 *     Pin `main` to a commit SHA for a reproducible citation.
 *   Cloudflare R2 / S3 / GCS: "https://<bucket>.r2.dev/chordoma-v1"
 *   Zenodo (archival, DOI-frozen — not for iterative updates):
 *     "https://zenodo.org/api/records/<RECORD_ID>/files"
 *     -> set ZENODO_MODE = true to append the required /content suffix.
 *
 * NOT usable: GitHub Release assets serve Range requests but send no
 * access-control-allow-origin header, so browser fetches are blocked. Git LFS
 * is also out — Pages serves the LFS pointer text, not the file contents.
 */
export const DATA_BASE = "./data";
export const ZENODO_MODE = false;
