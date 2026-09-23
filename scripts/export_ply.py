"""mt/7 のメッシュから色付きの点群を作り、PLY で書き出す。

    python scripts/export_ply.py samples/*.bin out/cloud.ply [--spacing 0.1] [--max-range 50]

各三角形の面に点をまき、面の低解像度画像（MCP4 の種類 3）から色を引く。点の
間隔は近くで --spacing（m）、遠くでは角度が一定（--angular ラジアン）になる。
面の画像の 1 画素が約 5.44 mrad なので、それより細かくしても情報は増えない。複数のファイルを渡すと 1 つの点群にまとめ、最初のファイルの
パノラマの位置を原点にした東・北・上（メートル）で書き出す。重なった部分は
--spacing の格子でまとめる（位置と色は平均）。

除くもの: 空（w = 0）を含む三角形、視線とほぼ平行な三角形（物陰をまたいで手前と
奥をつなぐもの。法線と視線のなす角の余弦が --min-facing 未満）、--max-range より
遠い点。--max-edge を付けると、それより長い辺を持つ三角形も除く。

側面の面には、画像では空が写っているのに遠くの面が作られている部分がある（空の
帽子の面 4 は w = 0 だが、側面の空はそうならない）。既定では、色が空らしく
（--sky-min より明るいか青が赤より 20 以上強く、青 ≥ 赤）、かつ --sky-range より遠い点を除く。
--sky-range 0 で無効。

画像の対応: 面の画像は、その面のメッシュの角度の範囲をちょうど覆う（1 画素
およそ 5.44 mrad）。docs/FORMAT.md「種類 9」の範囲の表を参照。
"""
import argparse
import io
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lookaround_mesh as lm  # noqa: E402

# 面ごとの (u 最小, u 最大, v 最小, v 最大)。サンプル 8 枚で共通（docs/FORMAT.md）
FACE_BOUNDS = {0: (-1078, 1078, -1220, 609), 1: (-588, 588, -1220, 609),
               2: (-1078, 1078, -1220, 609), 3: (-588, 588, -1220, 609),
               4: (-522, 522, -522, 522), 5: (-1062, 1062, -1132, 1132)}


def enu_basis(p):
    """ECEF の点 p での東・北・上の単位ベクトル（行）。"""
    up = p / np.linalg.norm(p)
    east = np.cross([0.0, 0.0, 1.0], up)
    east /= np.linalg.norm(east)
    return np.array([east, np.cross(up, east), up])


def face_image(blob, f):
    from pillow_heif import register_heif_opener
    from PIL import Image

    register_heif_opener()
    return np.asarray(Image.open(io.BytesIO(lm.face_images(blob)[f])).convert("RGB"))


def sample_triangles(edge, spacing):
    """三角形ごとに、最長辺を spacing（三角形ごとの配列）で割った数だけ分割した
    バリ中心座標の格子。

    Returns (三角形の番号, バリ中心座標 (n, 3))。"""
    k = np.clip(np.ceil(edge / spacing).astype(int), 1, 200)
    idx, bary = [], []
    for kk in np.unique(k):
        g = np.array([(i / kk, j / kk) for i in range(kk + 1) for j in range(kk + 1 - i)])
        w = np.c_[1 - g.sum(1), g]
        where = np.flatnonzero(k == kk)
        idx.append(np.repeat(where, len(w)))
        bary.append(np.tile(w, (len(where), 1)))
    return np.concatenate(idx), np.concatenate(bary)


def cloud(path, origin, R, a):
    """1 つの mt/7 から (点 (n, 3), 色 (n, 3))。点は origin, R の東・北・上。"""
    blob = open(path, "rb").read()
    pts, cols = [], []
    for f, raw in enumerate(lm.meshes(blob)):
        tris, P = lm.face_mesh(raw)
        ok = P[:, 2] > 0
        keep = np.array([t for t in tris if ok[list(t)].all()])
        if not len(keep):
            continue
        X = lm.to_local(np.where(ok[:, None], P, 1))
        A, B, C = X[keep[:, 0]], X[keep[:, 1]], X[keep[:, 2]]
        edge = np.maximum.reduce([np.linalg.norm(A - B, axis=1), np.linalg.norm(B - C, axis=1),
                                  np.linalg.norm(C - A, axis=1)])
        n = np.cross(B - A, C - A)
        mid = (A + B + C) / 3
        facing = np.abs((n * mid).sum(1)) / (np.linalg.norm(n, axis=1) * np.linalg.norm(mid, axis=1) + 1e-9)
        sel = facing >= a.min_facing
        if a.max_edge:
            sel &= edge <= a.max_edge
        keep, edge = keep[sel], edge[sel]
        if not len(keep):
            continue
        r = np.linalg.norm(mid[sel], axis=1)
        idx, bary = sample_triangles(edge, np.maximum(a.spacing, a.angular * r))
        uvw = (P[keep[idx]] * bary[:, :, None]).sum(1)

        img = face_image(blob, f)
        H, W = img.shape[:2]
        u0, u1, v0, v1 = FACE_BOUNDS[f]
        x = np.clip(((uvw[:, 0] - u0) / (u1 - u0) * W).astype(int), 0, W - 1)
        y = np.clip(((uvw[:, 1] - v0) / (v1 - v0) * H).astype(int), 0, H - 1)
        col = img[y, x]
        mat, _, _ = lm.header(raw)
        X3 = lm.to_local(uvw)
        if a.sky_range:  # 画像では空なのに面が作られている部分を除く
            c = col.astype(np.int16)
            sky = ((c.min(1) > a.sky_min) | (c[:, 2] - c[:, 0] > 20)) & (c[:, 2] >= c[:, 0])                 & (np.linalg.norm(X3, axis=1) > a.sky_range)
            X3, col = X3[~sky], col[~sky]
        cols.append(col)
        pts.append((X3 @ mat[:3, :3].T + mat[:3, 3] - origin) @ R.T)
    p, c = np.vstack(pts), np.vstack(cols)
    r = np.linalg.norm(p - (lm.header(lm.meshes(blob)[0])[0][:3, 3] - origin) @ R.T, axis=1)
    return p[r <= a.max_range], c[r <= a.max_range]


def voxel_average(p, c, size):
    """size の格子ごとに 1 点にまとめる（位置と色は平均）。"""
    key = np.round(p / size).astype(np.int64)
    _, inv, n = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    out = np.zeros((len(n), 3)), np.zeros((len(n), 3))
    for k in range(3):
        out[0][:, k] = np.bincount(inv, p[:, k]) / n
        out[1][:, k] = np.bincount(inv, c[:, k]) / n
    return out[0], out[1].round().astype(np.uint8)


def write_ply(path, p, c):
    head = ("ply\nformat binary_little_endian 1.0\n"
            f"element vertex {len(p)}\n"
            "property float x\nproperty float y\nproperty float z\n"
            "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            "end_header\n")
    rec = np.zeros(len(p), dtype=[("x", "<f4"), ("y", "<f4"), ("z", "<f4"),
                                  ("r", "u1"), ("g", "u1"), ("b", "u1")])
    for i, k in enumerate("xyz"):
        rec[k] = p[:, i]
    for i, k in enumerate("rgb"):
        rec[k] = c[:, i]
    with open(path, "wb") as fh:
        fh.write(head.encode("ascii"))
        fh.write(rec.tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mt7", nargs="+")
    ap.add_argument("ply")
    ap.add_argument("--spacing", type=float, default=0.1, help="点の間隔 [m]")
    ap.add_argument("--angular", type=float, default=0.0027,
                    help="遠くでの点の間隔（ラジアン）。面の画像の 1 画素は約 0.00544")
    ap.add_argument("--max-range", type=float, default=100.0, help="パノラマからの距離の上限 [m]")
    ap.add_argument("--max-edge", type=float, default=None, help="これより長い辺を持つ三角形を除く [m]")
    ap.add_argument("--sky-range", type=float, default=50.0,
                    help="この距離より遠く、色が空らしい点を除く [m]。0 で無効")
    ap.add_argument("--sky-min", type=int, default=165, help="空とみなす色の暗い方の成分の下限")
    ap.add_argument("--min-facing", type=float, default=0.08,
                    help="視線に対する面の向き（余弦）の下限。物陰をまたぐ三角形を除く")
    a = ap.parse_args()

    origin = lm.header(lm.meshes(open(a.mt7[0], "rb").read())[0])[0][:3, 3]
    R = enu_basis(origin)
    pts, cols = [], []
    for path in a.mt7:
        p, c = cloud(path, origin, R, a)
        pts.append(p)
        cols.append(c)
        print(f"{os.path.basename(path)}: {len(p)} points")
    p, c = voxel_average(np.vstack(pts), np.vstack(cols), a.spacing)
    write_ply(a.ply, p, c)
    print(f"{a.ply}: {len(p)} points, "
          f"{p[:, 0].min():.0f}..{p[:, 0].max():.0f} E, {p[:, 1].min():.0f}..{p[:, 1].max():.0f} N, "
          f"{p[:, 2].min():.0f}..{p[:, 2].max():.0f} U [m]")


if __name__ == "__main__":
    main()
