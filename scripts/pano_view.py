"""Look Around のパノラマから、ふつうのカメラで撮ったような画像を切り出す。

    python scripts/pano_view.py out/kanazawa/<id>_mt7.bin view.jpg --yaw 105 --pitch -3 --fov 60
    python scripts/pano_view.py out/kanazawa/<id>_mt7.bin view.jpg --yaw 105 --zoom 2   # 高解像度

--yaw は方位（北 0、東 90）、--pitch は上向きが正、--fov は縦の画角（three.js の
PerspectiveCamera と同じ）、いずれも度。向きはパノラマの位置での東・北・上で決める。

面の画像は、その面のメッシュの角度の範囲をちょうど覆う（docs/FORMAT.md）。面ごとに
向きの行列が違うので、視線を面の局所座標に直してから (u, v) = (方位, 下向きの仰角) を
求め、画像を引く。使うのは経度の帯 4 面（上下の帽子の面は使わない。帯は下 30°〜上 62° を覆う）。

既定では mt/7 に入っている低解像度の画像（1 画素約 5.4 mrad）を使う。--zoom を付けると
Apple のサーバーから高解像度の面画像を取る（0 が最高、7 が最低。3 で帯 1 枚 1568×1168）。
取った画像は mt/7 と同じフォルダに <id>_face<面>_z<zoom>.heic として保存し、次からはそれを使う。
高解像度の画像も mt/7 の画像と同じ範囲をずれなく覆う。
"""
import argparse
import io
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lookaround_mesh as lm  # noqa: E402
from lookaround_mesh.vertices import ANGLE  # noqa: E402

STRIPS = range(4)


def enu_basis(p):
    up = p / np.linalg.norm(p)
    east = np.cross([0.0, 0.0, 1.0], up)
    east /= np.linalg.norm(east)
    return np.array([east, np.cross(up, east), up])


def load_image(data):
    from PIL import Image
    from pillow_heif import register_heif_opener

    register_heif_opener()
    return np.asarray(Image.open(io.BytesIO(data)).convert("RGB"), np.float32)


def face_images(path, blob, zoom):
    """帯 4 面の画像。zoom が None なら mt/7 の中のもの。"""
    if zoom is None:
        return [load_image(lm.face_images(blob)[f]) for f in STRIPS]
    pid = int(os.path.basename(path).split("_")[0])
    index = json.load(open(os.path.join(os.path.dirname(path), "index.json")))
    bid = next(d["build_id"] for d in index if d["id"] == pid)
    out = []
    for f in STRIPS:
        cache = os.path.join(os.path.dirname(path), f"{pid}_face{f}_z{zoom}.heic")
        if not os.path.exists(cache):
            from streetlevel.lookaround import get_panorama_face
            from streetlevel.lookaround.auth import Authenticator

            open(cache, "wb").write(get_panorama_face((pid, bid), f, zoom, Authenticator()))
        out.append(load_image(open(cache, "rb").read()))
    return out


def bilinear(img, x, y):
    H, W = img.shape[:2]
    x = np.clip(x, 0, W - 1.001)
    y = np.clip(y, 0, H - 1.001)
    x0, y0 = x.astype(int), y.astype(int)
    fx, fy = (x - x0)[:, None], (y - y0)[:, None]
    return ((img[y0, x0] * (1 - fx) + img[y0, x0 + 1] * fx) * (1 - fy)
            + (img[y0 + 1, x0] * (1 - fx) + img[y0 + 1, x0 + 1] * fx) * fy)


def camera_rays(yaw, pitch, fov, w, h):
    """東・北・上での視線 (h*w, 3)。"""
    y, p = np.radians(yaw), np.radians(pitch)
    fwd = np.array([np.sin(y) * np.cos(p), np.cos(y) * np.cos(p), np.sin(p)])
    right = np.array([np.cos(y), -np.sin(y), 0.0])
    up = np.cross(right, fwd)
    tv = np.tan(np.radians(fov) / 2)
    th = tv * w / h
    xs = ((np.arange(w) + 0.5) / w * 2 - 1) * th
    ys = (1 - (np.arange(h) + 0.5) / h * 2) * tv
    X, Y = np.meshgrid(xs, ys)
    return (fwd + X[..., None] * right + Y[..., None] * up).reshape(-1, 3)


def render(path, yaw, pitch, fov, w, h, zoom=None):
    blob = open(path, "rb").read()
    meshes = lm.meshes(blob)
    imgs = face_images(path, blob, zoom)
    pos = lm.header(meshes[0])[0][:3, 3]
    d = camera_rays(yaw, pitch, fov, w, h) @ enu_basis(pos)      # ECEF の向き
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    best = np.full(len(d), np.inf)
    out = np.zeros((len(d), 3), np.float32)
    for f in STRIPS:
        M = lm.header(meshes[f])[0][:3, :3]
        _, P = lm.face_mesh(meshes[f])
        u0, u1, v0, v1 = P[:, 0].min(), P[:, 0].max(), P[:, 1].min(), P[:, 1].max()
        loc = d @ M                                               # 面の局所座標（x 前、y 横、z 下）
        u = np.arctan2(loc[:, 1], loc[:, 0]) / ANGLE
        v = np.arcsin(np.clip(loc[:, 2], -1, 1)) / ANGLE
        score = np.abs(u - (u0 + u1) / 2) / ((u1 - u0) / 2)      # 面の中心に近いほどよい
        ok = (u >= u0) & (u <= u1) & (v >= v0) & (v <= v1) & (score < best)
        img = imgs[f]
        H, W = img.shape[:2]
        x = (u[ok] - u0) / (u1 - u0) * W - 0.5
        y = (v[ok] - v0) / (v1 - v0) * H - 0.5
        out[ok] = bilinear(img, x, y)
        best[ok] = score[ok]
    return out.reshape(h, w, 3).clip(0, 255).astype(np.uint8), pos


def main():
    from PIL import Image

    ap = argparse.ArgumentParser()
    ap.add_argument("mt7")
    ap.add_argument("out")
    ap.add_argument("--yaw", type=float, default=0.0, help="方位 [度]（北 0、東 90）")
    ap.add_argument("--pitch", type=float, default=0.0, help="上向きが正 [度]")
    ap.add_argument("--fov", type=float, default=60.0, help="縦の画角 [度]")
    ap.add_argument("--size", default="1600x900")
    ap.add_argument("--zoom", type=int, default=None, help="高解像度の面画像のズーム（0〜7）")
    a = ap.parse_args()
    w, h = map(int, a.size.lower().split("x"))
    img, pos = render(a.mt7, a.yaw, a.pitch, a.fov, w, h, a.zoom)
    Image.fromarray(img).save(a.out, quality=90)
    print(f"{a.out}: yaw {a.yaw} pitch {a.pitch} fov {a.fov}, カメラ ECEF {pos.round(2).tolist()}")


if __name__ == "__main__":
    main()
