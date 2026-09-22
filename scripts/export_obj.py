"""mt/7 の 6 面のメッシュを 1 つの OBJ に書き出す。

    python scripts/export_obj.py samples/pano0_mt7.bin pano0.obj [--max-edge 3]

座標はパノラマの位置を原点にした東・北・上（メートル）。空（w = 0）の頂点を
含む三角形は除く。--max-edge を付けると、それより長い辺を持つ三角形（物陰を
またいで手前と奥をつなぐもの）も除く。面ごとに OBJ のグループ（face0〜5）。
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lookaround_mesh as lm  # noqa: E402


def enu_basis(p):
    """ECEF の点 p での東・北・上の単位ベクトル（行）。"""
    up = p / np.linalg.norm(p)
    east = np.cross([0.0, 0.0, 1.0], up)
    east /= np.linalg.norm(east)
    return np.array([east, np.cross(up, east), up])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mt7")
    ap.add_argument("obj")
    ap.add_argument("--max-edge", type=float, default=None)
    a = ap.parse_args()

    meshes = lm.meshes(open(a.mt7, "rb").read())
    origin = lm.header(meshes[0])[0][:3, 3]
    R = enu_basis(origin)
    lines, base = [], 1
    for f, raw in enumerate(meshes):
        tris, P = lm.face_mesh(raw)
        ok = P[:, 2] > 0
        X = (lm.to_ecef(raw, np.where(ok[:, None], P, 1)) - origin) @ R.T
        keep = [t for t in tris if ok[list(t)].all()]
        if a.max_edge:
            keep = [t for t in keep
                    if max(np.linalg.norm(X[t[i]] - X[t[i - 1]]) for i in range(3)) <= a.max_edge]
        lines.append(f"g face{f}")
        lines += [f"v {x:.3f} {y:.3f} {z:.3f}" for x, y, z in X]
        lines += [f"f {t[0] + base} {t[1] + base} {t[2] + base}" for t in keep]
        base += len(X)
    open(a.obj, "w").write("\n".join(lines) + "\n")
    print(f"{a.obj}: {base - 1} vertices")


if __name__ == "__main__":
    main()
