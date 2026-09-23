"""PLY の点群を、ブラウザで見るための一式に変換する。

    python scripts/export_web.py out/kanazawa.ply out/web [--spacing 0.25]

出力先に `index.html`（three.js のビューア）、`cloud.bin`（点の位置と色）、
`meta.json` を書く。見るときはその場で簡単なサーバーを立てる:

    python -m http.server 8000 --directory out/web

点は --spacing の格子でまとめてから、位置を int16 に丸める（1 点 9 バイト）。
ブラウザで扱える量にするため、既定では元より粗くする。

データは手元のサーバーから読むだけで、どこにも送らない。
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from export_ply import voxel_average  # noqa: E402

HTML = """<!DOCTYPE html>
<meta charset="utf-8">
<title>point cloud</title>
<style>
  html, body { margin: 0; height: 100%; background: #0b0d10; color: #cfd6dd;
               font: 13px/1.5 system-ui, sans-serif; overflow: hidden; }
  #ui { position: fixed; left: 12px; top: 12px; padding: 10px 12px; border-radius: 8px;
        background: rgba(20,24,28,.82); }
  #ui b { color: #fff; font-weight: 600; }
  label { display: block; margin-top: 6px; }
  input[type=range] { width: 150px; vertical-align: middle; }
</style>
<div id="ui">
  <b id="n">読み込み中…</b>
  <label>点の大きさ <input id="size" type="range" min="1" max="12" step="1" value="3"></label>
  <label>明るさ <input id="gain" type="range" min="50" max="250" step="5" value="100"></label>
  <div style="margin-top:6px;opacity:.75">左ドラッグ: 回転 / 右ドラッグ: 平行移動 / ホイール: 拡大</div>
</div>
<script type="importmap">
{ "imports": { "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
               "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/" } }
</script>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const meta = await (await fetch('meta.json')).json();
const buf = await (await fetch('cloud.bin')).arrayBuffer();
const n = meta.count;
const qi = new Int16Array(buf, 0, n * 3);
const ci = new Uint8Array(buf, n * 6, n * 3);

const pos = new Float32Array(n * 3);
for (let i = 0; i < n * 3; i++) pos[i] = qi[i] * meta.scale;
const col = new Float32Array(n * 3);
for (let i = 0; i < n * 3; i++) col[i] = ci[i] / 255;

const geo = new THREE.BufferGeometry();
geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
geo.setAttribute('color', new THREE.BufferAttribute(col, 3));
geo.computeBoundingSphere();

const mat = new THREE.PointsMaterial({ size: 3, sizeAttenuation: false, vertexColors: true });
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0b0d10);
scene.add(new THREE.Points(geo, mat));

const cam = new THREE.PerspectiveCamera(60, innerWidth / innerHeight, 0.1, 5000);
const r = geo.boundingSphere.radius;
cam.position.set(0, -r * 0.9, r * 0.55);
cam.up.set(0, 0, 1);

const ren = new THREE.WebGLRenderer({ antialias: false });
ren.setPixelRatio(devicePixelRatio);
ren.setSize(innerWidth, innerHeight);
document.body.appendChild(ren.domElement);

const ctl = new OrbitControls(cam, ren.domElement);
ctl.target.set(0, 0, 0);
ctl.enableDamping = true;
ctl.update();

document.getElementById('n').textContent =
  n.toLocaleString() + ' 点 · ' + meta.spacing.toFixed(2) + ' m 間隔';
document.getElementById('size').oninput = e => { mat.size = +e.target.value; };
document.getElementById('gain').oninput = e => {
  const g = +e.target.value / 100;
  for (let i = 0; i < n * 3; i++) col[i] = Math.min(1, ci[i] / 255 * g);
  geo.attributes.color.needsUpdate = true;
};
addEventListener('resize', () => {
  cam.aspect = innerWidth / innerHeight; cam.updateProjectionMatrix();
  ren.setSize(innerWidth, innerHeight);
});
(function loop() { requestAnimationFrame(loop); ctl.update(); ren.render(scene, cam); })();
</script>
"""


def read_ply(path):
    with open(path, "rb") as fh:
        head = b""
        while b"end_header\n" not in head:
            head += fh.read(1024)
        cut = head.index(b"end_header\n") + len(b"end_header\n")
        n = int([ln for ln in head[:cut].split(b"\n") if ln.startswith(b"element vertex")][0].split()[2])
        raw = head[cut:] + fh.read()
    rec = np.frombuffer(raw, dtype=[("x", "<f4"), ("y", "<f4"), ("z", "<f4"),
                                    ("r", "u1"), ("g", "u1"), ("b", "u1")], count=n)
    return (np.c_[rec["x"], rec["y"], rec["z"]].astype(np.float64),
            np.c_[rec["r"], rec["g"], rec["b"]].astype(np.float64))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ply")
    ap.add_argument("out")
    ap.add_argument("--spacing", type=float, default=0.25, help="まとめる格子の大きさ [m]")
    a = ap.parse_args()

    p, c = read_ply(a.ply)
    print(f"{len(p)} 点を読み込み")
    p, c = voxel_average(p, c, a.spacing)
    p = p - (p.min(0) + p.max(0)) / 2      # 中心を原点に
    scale = float(np.abs(p).max()) / 32767.0
    q = np.round(p / scale).astype("<i2")

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "cloud.bin"), "wb") as fh:
        fh.write(q.tobytes())
        fh.write(np.asarray(c, np.uint8).tobytes())
    json.dump({"count": len(p), "scale": scale, "spacing": a.spacing},
              open(os.path.join(a.out, "meta.json"), "w"))
    open(os.path.join(a.out, "index.html"), "w", encoding="utf-8").write(HTML)
    size = os.path.getsize(os.path.join(a.out, "cloud.bin"))
    print(f"{len(p)} 点 -> {a.out}  cloud.bin {size / 1e6:.1f} MB")
    print(f"見るには: python -m http.server 8000 --directory {a.out}")


if __name__ == "__main__":
    main()
