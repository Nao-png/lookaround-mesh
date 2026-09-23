"""out/web の中身を手元で配るだけの簡単なサーバー（export_web.py の出力を見る用）。

    python scripts/serve_web.py out/web [ポート]

ポートは引数、なければ環境変数 PORT、それもなければ 8000。
"""
import functools
import http.server
import os
import sys

d = sys.argv[1] if len(sys.argv) > 1 else "out/web"
port = int(sys.argv[2]) if len(sys.argv) > 2 else int(os.environ.get("PORT", 8000))
handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=d)
print(f"http://localhost:{port}/  ({d})", flush=True)
http.server.ThreadingHTTPServer(("127.0.0.1", port), handler).serve_forever()
