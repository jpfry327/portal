"""Map readable Parquet filenames to opaque tokens for the published site.

Pre-publication the portal is deployed from a PUBLIC repo (GitHub Pages will
not serve a private repo on the free tier). Randomising the data filenames
keeps a directory listing or a stray link from advertising what the files
contain — "expr_zhang.parquet" tells a reader a lot; "d6114d3912a63.parquet"
does not.

This is obscurity, not security. Anyone who loads the site can read the
resolved URLs out of the JS. It exists to prevent accidental discovery
(search engines, guessed URLs), not to withstand someone who wants the data.

Tokens are a salted hash of the logical name, so they are STABLE: re-running
the ETL regenerates identical filenames, which keeps git from storing a fresh
copy of every file on each deploy.

Writes datamap.json (logical name -> token), consumed by app.js at runtime.

IMPORTANT: the map belongs to the PUBLISHED build only. Renaming happens in
the deploy workflow, so the files on disk during local development keep their
readable names — a datamap.json sitting in site/ would make the local app
request tokens that do not exist there and 404. app.js treats a missing map as
"use readable names", which is exactly right for local dev.

The deploy workflow calls this with an explicit output path inside its build
directory. Running it with no argument writes to ./datamap.json so it cannot
accidentally poison site/.

Run (deploy workflow):  python portal/obfuscate_names.py _site/datamap.json
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "portal" / "data"
# Never defaults into site/ — see the note above.
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("datamap.json")

# Changing SALT re-randomises every filename (forcing a full re-upload).
SALT = "chordoma-portal-v1"

# Kept under its real name: it is the ETL input, is gitignored, and is never
# published (>100 MB, over GitHub's per-file limit).
EXCLUDE = {"meth_discovery"}


def token(name: str) -> str:
    return "d" + hashlib.sha256((SALT + name).encode()).hexdigest()[:12]


def build_map() -> dict[str, str]:
    names = sorted(p.stem for p in DATA.glob("*.parquet") if p.stem not in EXCLUDE)
    if not names:
        raise SystemExit(f"no parquet files in {DATA} — run etl.py first")
    mapping = {n: token(n) for n in names}
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit("token collision — lengthen the hash prefix")
    return mapping


def main() -> None:
    mapping = build_map()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(mapping, indent=1, sort_keys=True) + "\n")
    print(f"wrote {OUT} ({len(mapping)} files)")
    for n in list(mapping)[:3]:
        print(f"  {n} -> {mapping[n]}")


if __name__ == "__main__":
    main()
