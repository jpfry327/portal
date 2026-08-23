"""Local dev server for the static atlas with HTTP Range support.

Python's http.server does NOT support Range requests, which DuckDB-WASM needs.
This wrapper adds single-range support so local dev behaves like GitHub
Releases / Zenodo / S3.

Run:  python serve.py [port]     (default 8000)
Then open http://localhost:<port>
"""
import os
import re
import sys
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class RangeHandler(SimpleHTTPRequestHandler):
    def send_head(self):
        path = self.translate_path(self.path)
        if os.path.isdir(path):
            return super().send_head()
        rng = self.headers.get("Range")
        if not rng:
            return super().send_head()
        m = re.match(r"bytes=(\d*)-(\d*)$", rng.strip())
        if not m or not os.path.isfile(path):
            return super().send_head()
        size = os.path.getsize(path)
        start = int(m.group(1)) if m.group(1) else None
        end = int(m.group(2)) if m.group(2) else None
        if start is None:          # suffix range: last N bytes
            start = max(0, size - (end or 0))
            end = size - 1
        else:
            end = min(end if end is not None else size - 1, size - 1)
        if start > end or start >= size:
            self.send_response(416)
            self.send_header("Content-Range", f"bytes */{size}")
            self.end_headers()
            return None
        f = open(path, "rb")
        f.seek(start)
        self.send_response(206)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(end - start + 1))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self._range_span = end - start + 1
        return f

    def copyfile(self, source, outputfile):
        if hasattr(self, "_range_span"):
            remaining = self._range_span
            del self._range_span
            while remaining > 0:
                chunk = source.read(min(65536, remaining))
                if not chunk:
                    break
                outputfile.write(chunk)
                remaining -= len(chunk)
        else:
            super().copyfile(source, outputfile)

    # Dev ergonomics: never let the browser cache the app code. ES modules
    # (<script type="module">) are cached aggressively enough that a plain
    # reload keeps serving a stale app.js after an edit. Parquet is left
    # cacheable so DuckDB-WASM range reads stay fast.
    NO_CACHE_EXT = (".js", ".html", ".htm", ".css", ".json", ".mjs")

    def end_headers(self):
        path = self.translate_path(self.path)
        if path.lower().endswith(self.NO_CACHE_EXT) or os.path.isdir(path):
            self.send_header("Cache-Control", "no-store, must-revalidate")
            self.send_header("Pragma", "no-cache")
            self.send_header("Expires", "0")
        super().end_headers()


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    server = ThreadingHTTPServer(("127.0.0.1", port), partial(RangeHandler))
    print(f"Serving atlas at http://localhost:{port} (Range requests enabled)")
    server.serve_forever()
