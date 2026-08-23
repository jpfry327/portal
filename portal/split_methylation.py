"""Split the wide methylation matrix into one Parquet file per chromosome.

Why: GitHub enforces a hard 100 MB per-file limit, and the combined
meth_discovery.parquet is ~460 MB — it cannot be pushed. Splitting by
chromosome yields 22 files whose largest (chr1) is ~45 MB, so the whole
dataset can live in the repo and be served by GitHub Pages, which honours
HTTP Range requests (verified: returns HTTP 206 with a correct
content-range, including suffix ranges).

Splitting by chromosome is also the natural access pattern: no view in the
portal ever needs two chromosomes at once.

Each output is sorted by genomic position with small row groups, so a
genomic-range query prunes to the row groups it needs via Parquet
statistics rather than scanning the file.

Inputs  (portal/data/):
  meth_discovery.parquet     wide M-values, one row per probe
  meth_probe_index.parquet   meth_site -> (chr, pos), sorted by meth_site

Outputs (portal/data/):
  meth_chr1.parquet ... meth_chr22.parquet
      columns: meth_site, pos, <162 sample columns>; sorted by pos

Probe coordinates come from the Bioconductor annotation package
IlluminaHumanMethylationEPICanno.ilm10b4.hg19 (EPIC v1, build hg19).
All 756,447 probes in the matrix are annotated by it. The matrix contains
autosomes only — chrX/chrY probes were filtered upstream.

Run:  python portal/split_methylation.py
"""
from __future__ import annotations

from pathlib import Path

import duckdb

OUT = Path(__file__).resolve().parent / "data"
SRC = OUT / "meth_discovery.parquet"
IDX = OUT / "meth_probe_index.parquet"

# 2048 rows/group keeps each group near ~1.2 MB, so a range query transfers
# roughly that much per group touched instead of the whole file.
ROW_GROUP = 2048


def main() -> None:
    if not SRC.exists():
        raise SystemExit(f"missing {SRC} — run etl.py first")
    if not IDX.exists():
        raise SystemExit(f"missing {IDX} — build the probe index first")

    con = duckdb.connect()
    chroms = [r[0] for r in con.execute(
        f"SELECT DISTINCT chr FROM read_parquet('{IDX}') ORDER BY "
        f"TRY_CAST(REPLACE(chr, 'chr', '') AS INTEGER), chr"
    ).fetchall()]
    print(f"== splitting {SRC.name} into {len(chroms)} chromosome files ==")

    total = 0
    for c in chroms:
        dest = OUT / f"meth_{c}.parquet"
        # Join on the index to attach pos, then sort by position so that
        # Parquet row-group min/max statistics make range queries prunable.
        con.execute(f"""
            COPY (
              SELECT i.pos AS pos, m.*
              FROM read_parquet('{SRC}') m
              JOIN read_parquet('{IDX}') i USING (meth_site)
              WHERE i.chr = '{c}'
              ORDER BY i.pos
            ) TO '{dest}' (FORMAT PARQUET, COMPRESSION ZSTD,
                           ROW_GROUP_SIZE {ROW_GROUP})
        """)
        n = con.execute(f"SELECT COUNT(*) FROM read_parquet('{dest}')").fetchone()[0]
        mb = dest.stat().st_size / 1e6
        total += mb
        flag = "" if mb < 100 else "  <-- OVER GITHUB 100 MB LIMIT"
        print(f"  {dest.name:20s} {n:>7d} probes  {mb:6.1f} MB{flag}")

    print(f"total {total:.0f} MB across {len(chroms)} files")


if __name__ == "__main__":
    main()
