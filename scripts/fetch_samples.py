"""渋谷スクランブル交差点の Look Around パノラマ 8 枚の mt/7 を取り直す。

    python scripts/fetch_samples.py [出力先]

同じ走行で 4.5 m 間隔の 6 枚（中心＋5 枚）と、翌日の別走行で中心から 1.7 m と
3.3 m の 2 枚。一覧は samples/neighbours.json。取得したデータは Apple の著作物
（samples/NOTICE.md）。
"""
import json
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from lookaround_mesh.fetch import fetch_mt  # noqa: E402

CENTER = 14755389547370281509   # 2026-04-21 の車載撮影、交差点の最寄り
LAT, LON = 35.6595, 139.7005


def pick():
    from streetlevel.lookaround import lookaround as la
    tile = la.get_coverage_tile_by_latlon(LAT, LON)
    panos = tile.panos if hasattr(tile, "panos") else tile
    a = next(p for p in panos if p.id == CENTER)

    def dist(p):
        return math.hypot((p.lat - a.lat) * 111320, (p.lon - a.lon) * 111320 * math.cos(math.radians(a.lat)))
    same = sorted((p for p in panos if p.id != CENTER and abs((p.date - a.date).total_seconds()) < 300), key=dist)[:5]
    other = sorted((p for p in panos if abs((p.date - a.date).total_seconds()) > 3600), key=dist)[:2]
    return [(p, dist(p), i <= len(same)) for i, p in enumerate([a] + same + other)]


def main(out):
    os.makedirs(out, exist_ok=True)
    meta = []
    for k, (p, d, same) in enumerate(pick()):
        b = fetch_mt(p.id, p.build_id)
        open(os.path.join(out, f"pano{k}_mt7.bin"), "wb").write(b)
        meta.append({"k": k, "id": p.id, "build_id": p.build_id, "lat": p.lat, "lon": p.lon,
                     "date": str(p.date), "dist_m": round(d, 1), "same_drive": same, "bytes": len(b)})
        print(f"pano{k}: {'同じ走行' if same else '別の走行'} {d:5.1f} m  {p.date}  {len(b)} B")
        time.sleep(0.5)
    json.dump(meta, open(os.path.join(out, "neighbours.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "samples"))
