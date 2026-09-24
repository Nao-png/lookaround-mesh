"""色付きの点群（export_ply.py の出力）を、学習なしで 3D Gaussian Splatting の形にする。

    python scripts/export_splat.py out/street_crop.ply out/splat_street [--spacing 0.05]
    python scripts/serve_web.py out/splat_street 8000

点を --spacing の格子でまとめ、1 点を 1 つのガウスにする。点群は面の上にまいた点
なので、ガウスは周りの点（--neighbors 個）の主成分分析で求めた面に沿って平たくする:
面の中の 2 方向の大きさは近くの点までの平均距離から、法線方向はその面からの
ばらつきから決める（葉など面にならないところは丸くなる）。色は球面調和の 0 次
だけ（向きによる色の変化はない）、不透明度は --opacity。

写真から学習した 3DGS ではないので、見た目は点群を滑らかにした程度。
学習の初期値にも使える。

出力: scene.ply（Inria の 3DGS と同じ形式。SuperSplat などでもそのまま開ける）、
index.html（Spark のビューア）、meta.json。座標は範囲の中心を原点にした
東・北・上（メートル）で、元の点群の座標は meta.json の offset を足すと戻る。
データは手元のサーバーから読むだけで、どこにも送らない。
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from export_ply import voxel_average  # noqa: E402
from export_web import read_ply  # noqa: E402

SH_C0 = 0.28209479177387814   # 球面調和の 0 次の係数

HTML = """<!DOCTYPE html>
<meta charset="utf-8">
<title>gaussian splats</title>
<style>
  html, body { margin: 0; height: 100%; background: #0b0d10; color: #cfd6dd;
               font: 13px/1.5 system-ui, sans-serif; overflow: hidden; }
  #ui { position: fixed; left: 12px; top: 12px; padding: 10px 12px; border-radius: 8px;
        background: rgba(20,24,28,.82); }
  #ui b { color: #fff; font-weight: 600; }
</style>
<div id="ui">
  <b id="n">読み込み中…</b>
  <div style="margin-top:6px;opacity:.75">左ドラッグ: 回転 / 右ドラッグ: 平行移動 / ホイール: 拡大</div>
</div>
<script type="importmap">
{ "imports": { "three": "https://cdn.jsdelivr.net/npm/three@0.180.0/build/three.module.js",
               "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.180.0/examples/jsm/",
               "@sparkjsdev/spark": "https://cdn.jsdelivr.net/npm/@sparkjsdev/spark@2.2.0/dist/spark.module.js" } }
</script>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { SparkRenderer, SplatMesh } from '@sparkjsdev/spark';

const meta = await (await fetch('meta.json')).json();
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0b0d10);
const cam = new THREE.PerspectiveCamera(60, innerWidth / innerHeight, 0.05, 2000);
cam.up.set(0, 0, 1);
const r = meta.radius;
cam.position.set(0, -r * 0.5, r * 0.35);

const ren = new THREE.WebGLRenderer({ antialias: false });
ren.setPixelRatio(Math.min(devicePixelRatio, 2));
ren.setSize(innerWidth, innerHeight);
document.body.appendChild(ren.domElement);
scene.add(new SparkRenderer({ renderer: ren }));

const splats = new SplatMesh({ url: 'scene.ply' });
scene.add(splats);
await splats.initialized;
document.getElementById('n').textContent = meta.count.toLocaleString() + ' ガウス（学習なし）';

const ctl = new OrbitControls(cam, ren.domElement);
ctl.target.set(0, 0, 0);
ctl.enableDamping = true;
ctl.update();
window.view = { cam, ctl, meta };   // コンソールから視点を動かす用

addEventListener('resize', () => {
  cam.aspect = innerWidth / innerHeight; cam.updateProjectionMatrix();
  ren.setSize(innerWidth, innerHeight);
});
ren.setAnimationLoop(() => { ctl.update(); ren.render(scene, cam); });
</script>
"""


def shapes(p, k, spacing):
    """点ごとの (回転行列 (n, 3, 3)、標準偏差 (n, 3))。列は面の中の 2 方向と法線。"""
    from scipy.spatial import cKDTree

    tree = cKDTree(p)
    R = np.empty((len(p), 3, 3))
    S = np.empty((len(p), 3))
    for s in range(0, len(p), 500_000):
        q = p[s:s + 500_000]
        d, j = tree.query(q, k=k + 1, workers=-1)
        nb = p[j[:, 1:]] - p[j[:, 1:]].mean(1, keepdims=True)
        lam, vec = np.linalg.eigh(np.einsum("nki,nkj->nij", nb, nb) / k)   # 小さい順
        lam = np.maximum(lam, 0)
        n, t1 = vec[:, :, 0], vec[:, :, 2]
        t2 = np.cross(n, t1)                     # 右手系にそろえる（行列式 +1）
        R[s:s + len(q)] = np.stack([t1, t2, n], axis=2)
        side = np.clip(d[:, 1:].mean(1) * 0.6, spacing * 0.3, spacing * 4)
        S[s:s + len(q), 0] = side
        S[s:s + len(q), 1] = side
        S[s:s + len(q), 2] = np.clip(np.sqrt(lam[:, 0]), side * 0.1, side)
    return R, S


def write_splat_ply(path, p, c, R, S, opacity):
    from scipy.spatial.transform import Rotation

    q = Rotation.from_matrix(R).as_quat()        # (x, y, z, w)
    names = (["x", "y", "z", "nx", "ny", "nz", "f_dc_0", "f_dc_1", "f_dc_2", "opacity",
              "scale_0", "scale_1", "scale_2", "rot_0", "rot_1", "rot_2", "rot_3"])
    rec = np.zeros(len(p), dtype=[(k, "<f4") for k in names])
    for i, k in enumerate("xyz"):
        rec[k] = p[:, i]
    for i in range(3):
        rec[f"f_dc_{i}"] = (c[:, i] / 255.0 - 0.5) / SH_C0
        rec[f"scale_{i}"] = np.log(S[:, i])
    rec["opacity"] = np.log(opacity / (1 - opacity))   # 保存するのはシグモイドの前の値
    rec["rot_0"], rec["rot_1"], rec["rot_2"], rec["rot_3"] = q[:, 3], q[:, 0], q[:, 1], q[:, 2]
    head = ("ply\nformat binary_little_endian 1.0\n"
            f"element vertex {len(p)}\n"
            + "".join(f"property float {k}\n" for k in names)
            + "end_header\n")
    with open(path, "wb") as fh:
        fh.write(head.encode("ascii"))
        fh.write(rec.tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ply")
    ap.add_argument("out")
    ap.add_argument("--spacing", type=float, default=0.05, help="まとめる格子の大きさ [m]")
    ap.add_argument("--neighbors", type=int, default=12, help="面の向きを求める近傍の点の数")
    ap.add_argument("--opacity", type=float, default=0.95)
    a = ap.parse_args()

    p, c = read_ply(a.ply)
    print(f"{len(p)} 点を読み込み", flush=True)
    p, c = voxel_average(p, c, a.spacing)
    p = p.astype(np.float64)
    offset = (p.min(0) + p.max(0)) / 2
    p -= offset
    print(f"{len(p)} 点にまとめた。面の向きを求める…", flush=True)
    R, S = shapes(p, a.neighbors, a.spacing)
    os.makedirs(a.out, exist_ok=True)
    write_splat_ply(os.path.join(a.out, "scene.ply"), p, c, R, S, a.opacity)
    json.dump({"count": len(p), "spacing": a.spacing, "offset": [float(v) for v in offset],
               "radius": float(np.abs(p[:, :2]).max())},
              open(os.path.join(a.out, "meta.json"), "w"))
    open(os.path.join(a.out, "index.html"), "w", encoding="utf-8").write(HTML)
    size = os.path.getsize(os.path.join(a.out, "scene.ply"))
    print(f"{len(p)} ガウス -> {a.out}（scene.ply {size / 1e6:.0f} MB）")
    print(f"見るには: python scripts/serve_web.py {a.out} 8000")


if __name__ == "__main__":
    main()
