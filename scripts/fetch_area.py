"""ある地点のまわりの Look Around パノラマの mt/7 をまとめて取る。

    python scripts/fetch_area.py 36.5722 136.6680 --radius 200 --out out/kanazawa

被覆タイル（3×3）を引いて、中心から --radius メートル以内のパノラマを順に取る。
出力先に `<panoid>_mt7.bin` と、一覧の `index.json`（id, build_id, lat, lon, date,
dist_m, bytes）を書く。すでにあるファイルは飛ばすので、途中で止めても再開できる。

取得したデータは Apple の著作物（samples/NOTICE.md）。リポジトリには入れないこと。
"""
import argparse
import json
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from lookaround_mesh.fetch import fetch_mt  # noqa: E402


def nearby(lat, lon, radius):
    """中心から radius メートル以内のパノラマを、近い順に。"""
    from streetlevel.lookaround import lookaround as la

    tile = la.get_coverage_tile_by_latlon(lat, lon)
    x, y = tile.x, tile.y
    seen, out = set(), []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            t = la.get_coverage_tile(x + dx, y + dy)
            for p in (t.panos if hasattr(t, "panos") else t):
                if p.id in seen:
                    continue
                seen.add(p.id)
                d = math.hypot((p.lat - lat) * 111320,
                               (p.lon - lon) * 111320 * math.cos(math.radians(lat)))
                if d <= radius:
                    out.append((d, p))
    return [(d, p) for d, p in sorted(out, key=lambda t: t[0])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lat", type=float)
    ap.add_argument("lon", type=float)
    ap.add_argument("--radius", type=float, default=200.0, help="中心からの距離 [m]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=None, help="取る枚数の上限（近い順）")
    ap.add_argument("--delay", type=float, default=0.3, help="1 枚ごとに待つ秒数")
    a = ap.parse_args()

    found = nearby(a.lat, a.lon, a.radius)[:a.limit]
    os.makedirs(a.out, exist_ok=True)
    print(f"{len(found)} 枚")
    meta, done = [], 0
    for i, (d, p) in enumerate(found):
        path = os.path.join(a.out, f"{p.id}_mt7.bin")
        if os.path.exists(path):
            n = os.path.getsize(path)
        else:
            try:
                b = fetch_mt(p.id, p.build_id)
            except Exception as e:  # 取れないパノラマは飛ばす
                print(f"  {p.id}: {e}")
                continue
            open(path, "wb").write(b)
            n = len(b)
            done += 1
            time.sleep(a.delay)
        meta.append({"id": p.id, "build_id": p.build_id, "lat": p.lat, "lon": p.lon,
                     "date": str(p.date), "dist_m": round(d, 1), "bytes": n})
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(found)}  {sum(m['bytes'] for m in meta) / 1e6:.1f} MB")
    json.dump(meta, open(os.path.join(a.out, "index.json"), "w"), indent=1)
    print(f"{len(meta)} 枚（うち新しく取得 {done}）, "
          f"{sum(m['bytes'] for m in meta) / 1e6:.1f} MB -> {a.out}")


if __name__ == "__main__":
    main()
