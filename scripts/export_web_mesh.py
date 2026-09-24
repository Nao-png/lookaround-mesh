"""複数のパノラマのメッシュを、面の画像を貼った 1 つのメッシュにして、ブラウザで見る一式を書き出す。

    python scripts/export_web_mesh.py out/kanazawa out/web_mesh [--max-range 25] [--jobs 12]
    python scripts/serve_web.py out/web_mesh 8000

パノラマのメッシュはそれぞれ自分の位置から見た奥行きなので、全部をそのまま重ねると
同じ壁が何十枚も重なる。そこで三角形ごとに、それをいちばん正確に写している
パノラマのものだけを残す。正確さは、画像の 1 画素が面の上で占める大きさ
（距離 × 1 画素の角度 ÷ 視線と法線のなす角の余弦）と、奥行きの刻み（距離の逆数 w が
整数なので、距離 r では約 r² / 1000 m。20 m で 0.4 m、50 m で 2.5 m）を合わせた
誤差で比べる。候補は水平に近い --candidates 枚で、
面を裏から見る位置にあるものは除く（壁の向こうの通りのパノラマなど）。遮るものは
考えない。路面だけなら、ほぼ水平距離でいちばん近いパノラマを選ぶことになるが、高い
ビルの壁は足元の歩道からより通りの向こうからの方が細かく写るので、そちらが選ばれる。
隣どうしの三角形の切れ目はかみ合わないので、判定は重心でなく 3 つの頂点で行い、
どれか 1 つでも自分が最良なら残す。境目では両側の三角形が 1 列ずつ重なり、すき間が
できない（重心で決めると、どちらにも属さない細いすき間が残る）。誤差が最良の
(1 + --overlap) 倍以内なら最良とみなす。

面の画像は、その面のメッシュの角度の範囲をちょうど覆う（export_ply.py と同じ対応）
ので、頂点の (u, v) がそのまま画像の座標になる。残った三角形が使う部分だけを
64 px の区画ごとに切り出し、幅 4096 px のアトラスに詰める。

除くもの: 空（w = 0）を含む三角形、視線とほぼ平行な三角形（--min-facing。どの面も
よく写すパノラマから見るので、点群のときより厳しくしてよい）、
--max-edge より長い辺を持つ三角形、パノラマから --max-range より遠い三角形、
--sky-range より遠くて色が空らしい三角形。

出力: index.html（three.js のビューア）、mesh.bin、atlas*.jpg、meta.json。
座標は全パノラマの中心を原点にした東・北・上（メートル）。
データは手元のサーバーから読むだけで、どこにも送らない。
"""
import argparse
import glob
import json
import multiprocessing as mp
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lookaround_mesh as lm  # noqa: E402
from export_ply import enu_basis, face_image  # noqa: E402

ATLAS = 4096
PIXEL = 0.0054   # 面の画像の 1 画素 [rad]
PAD = 2      # 切り出しの余白 [px]（隣の切り出しがにじまないように）
TILE = 64    # 切り出しの区画 [px]

HTML = """<!DOCTYPE html>
<meta charset="utf-8">
<title>mesh</title>
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
  <label>明るさ <input id="gain" type="range" min="50" max="250" step="5" value="100"></label>
  <label><input id="wire" type="checkbox"> 三角形の辺</label>
  <label><input id="panos" type="checkbox"> パノラマの位置</label>
  <div style="margin-top:6px;opacity:.75">左ドラッグ: 回転 / 右ドラッグ: 平行移動 / ホイール: 拡大<br>
  ダブルクリック: そこを回転の中心に</div>
</div>
<script type="importmap">
{ "imports": { "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
               "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/" } }
</script>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const meta = await (await fetch('meta.json')).json();
const buf = await (await fetch('mesh.bin')).arrayBuffer();
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0b0d10);

const loader = new THREE.TextureLoader();
const mats = [], meshes = [];
let tris = 0;
const texs = await Promise.all(meta.groups.map(g => loader.loadAsync(g.file)));   // 並行して読む
for (const [k, g] of meta.groups.entries()) {
  const tex = texs[k];
  tex.flipY = false;
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.generateMipmaps = false;               // 隣の切り出しがにじまないように
  tex.minFilter = THREE.LinearFilter;
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(buf, g.pos, g.nv * 3), 3));
  geo.setAttribute('uv', new THREE.BufferAttribute(new Float32Array(buf, g.uv, g.nv * 2), 2));
  geo.setIndex(new THREE.BufferAttribute(new Uint32Array(buf, g.idx, g.nt * 3), 1));
  geo.computeBoundingSphere();
  const mat = new THREE.MeshBasicMaterial({ map: tex, side: THREE.DoubleSide });
  const m = new THREE.Mesh(geo, mat);
  scene.add(m); mats.push(mat); meshes.push(m);
  tris += g.nt;
}

const pp = new THREE.BufferGeometry();
pp.setAttribute('position', new THREE.BufferAttribute(new Float32Array(meta.panos), 3));
const panos = new THREE.Points(pp, new THREE.PointsMaterial({ color: 0xff5040, size: 6, sizeAttenuation: false }));
panos.visible = false;
scene.add(panos);

const cam = new THREE.PerspectiveCamera(60, innerWidth / innerHeight, 0.1, 5000);
const r = meta.radius;
cam.position.set(0, -r * 0.6, r * 0.45);
cam.up.set(0, 0, 1);

const ren = new THREE.WebGLRenderer({ antialias: true });
ren.setPixelRatio(devicePixelRatio);
ren.setSize(innerWidth, innerHeight);
document.body.appendChild(ren.domElement);

const ctl = new OrbitControls(cam, ren.domElement);
ctl.target.set(0, 0, 0);
ctl.enableDamping = true;
ctl.update();
window.view = { cam, ctl, meta };   // コンソールから視点を動かす用

document.getElementById('n').textContent =
  tris.toLocaleString() + ' 三角形 · ' + (meta.label || 'パノラマ ' + meta.count + ' 枚');
document.getElementById('gain').oninput = e => {
  for (const m of mats) m.color.setScalar(+e.target.value / 100);
};
document.getElementById('wire').onchange = e => {
  for (const m of mats) m.wireframe = e.target.checked;
};
document.getElementById('panos').onchange = e => { panos.visible = e.target.checked; };

const ray = new THREE.Raycaster();
ren.domElement.addEventListener('dblclick', e => {
  ray.setFromCamera(new THREE.Vector2(e.clientX / innerWidth * 2 - 1, -e.clientY / innerHeight * 2 + 1), cam);
  const hit = ray.intersectObjects(meshes, false)[0];
  if (!hit) return;
  cam.position.add(hit.point.clone().sub(ctl.target));
  ctl.target.copy(hit.point);
});
addEventListener('resize', () => {
  cam.aspect = innerWidth / innerHeight; cam.updateProjectionMatrix();
  ren.setSize(innerWidth, innerHeight);
});
(function loop() { requestAnimationFrame(loop); ctl.update(); ren.render(scene, cam); })();
</script>
"""


def pano_position(path):
    return lm.header(lm.meshes(open(path, "rb").read())[0])[0][:3, 3]


_JOB = {}


def _init(pos, origin, R, a):
    from scipy.spatial import cKDTree

    _JOB.update(pos=pos, origin=origin, R=R, a=a, tree=cKDTree(pos[:, :2]))


def best_view(i, cen, nrm, a):
    """点 cen（法線 nrm の面の上）のうち、パノラマ i がほかの候補と比べて
    (1 + overlap) 倍以内の誤差で写しているもの（真偽の配列）。"""
    pos = _JOB["pos"]
    k = min(a.candidates, len(pos))
    _, J = _JOB["tree"].query(cen[:, :2], k=k)
    J = J.reshape(len(cen), k)
    nn = np.linalg.norm(nrm, axis=1) + 1e-12

    def score(V):   # V: 視点から点へ (..., 3)。面の上での誤差 [m]
        r = np.linalg.norm(V, axis=-1) + 1e-9
        cos = (V * nrm[..., None, :] if V.ndim == 3 else V * nrm).sum(-1)
        cos = cos / (r * (nn[:, None] if V.ndim == 3 else nn))
        side = r * PIXEL / np.maximum(np.abs(cos), 1e-6)   # 1 画素が面の上で占める大きさ
        depth = r * r / lm.vertices.DEPTH                  # w が整数なので奥行きの刻み
        return np.hypot(side, depth), cos, r

    own, oc, _ = score(cen - pos[i])
    s, c, r = score(cen[:, None, :] - pos[J])
    bad = (np.sign(c) != np.sign(oc)[:, None]) | (np.abs(c) < a.min_facing) | (r > a.max_range)
    s[bad] = np.inf
    return own <= np.minimum(s.min(1), own) * (1 + a.overlap)


def trimmed(i, path):
    """パノラマ i の、自分がよく写している部分の三角形。面ごとに
    (頂点 (m, 3) 東・北・上, 切り出した画像での画素座標 (m, 2), 三角形 (k, 3), 切り出した画像)。"""
    a, origin, R = _JOB["a"], _JOB["origin"], _JOB["R"]
    blob = open(path, "rb").read()
    out = []
    for f, raw in enumerate(lm.meshes(blob)):
        tris, P = lm.face_mesh(raw)
        ok = P[:, 2] > 0
        T = np.asarray(tris)
        if not len(T):
            continue
        T = T[ok[T].all(1)]
        if not len(T):
            continue
        X = lm.to_local(np.where(ok[:, None], P, 1))
        A, B, C = X[T[:, 0]], X[T[:, 1]], X[T[:, 2]]
        edge = np.maximum.reduce([np.linalg.norm(A - B, axis=1), np.linalg.norm(B - C, axis=1),
                                  np.linalg.norm(C - A, axis=1)])
        n = np.cross(B - A, C - A)
        mid = (A + B + C) / 3
        facing = np.abs((n * mid).sum(1)) / (np.linalg.norm(n, axis=1) * np.linalg.norm(mid, axis=1) + 1e-9)
        sel = (facing >= a.min_facing) & (np.linalg.norm(mid, axis=1) <= a.max_range)
        if a.max_edge:
            sel &= edge <= a.max_edge
        mat, _, _ = lm.header(raw)
        E = (X @ mat[:3, :3].T + mat[:3, 3] - origin) @ R.T
        cen = (E[T[:, 0]] + E[T[:, 1]] + E[T[:, 2]]) / 3
        nrm = np.cross(E[T[:, 1]] - E[T[:, 0]], E[T[:, 2]] - E[T[:, 0]])
        best = best_view(i, E[T.T.ravel()], np.tile(nrm, (3, 1)), a).reshape(3, -1)
        sel &= best.any(0)
        img = face_image(blob, f)
        H, W = img.shape[:2]
        u0, u1 = P[:, 0].min(), P[:, 0].max()   # 画像はこの範囲をちょうど覆う
        v0, v1 = P[:, 1].min(), P[:, 1].max()
        px = np.c_[(P[:, 0] - u0) / (u1 - u0) * W, (P[:, 1] - v0) / (v1 - v0) * H]
        if a.sky_range:  # 画像では空なのに面が作られている部分を除く
            c = (px[T[:, 0]] + px[T[:, 1]] + px[T[:, 2]]) / 3
            col = img[np.clip(c[:, 1].astype(int), 0, H - 1),
                      np.clip(c[:, 0].astype(int), 0, W - 1)].astype(np.int16)
            sky = ((col.min(1) > a.sky_min) | (col[:, 2] - col[:, 0] > 20)) & (col[:, 2] >= col[:, 0])
            sel &= ~(sky & (np.linalg.norm(mid, axis=1) > a.sky_range))
        T = T[sel]
        if not len(T):
            continue
        # 残る三角形は画像の中で散らばるので、重心の区画ごとに切り出す（面ごとの
        # 外接矩形で切ると空白まで詰めることになる）
        c = (px[T[:, 0]] + px[T[:, 1]] + px[T[:, 2]]) / 3
        cell = (c[:, 1] // TILE).astype(int) * 1000 + (c[:, 0] // TILE).astype(int)
        for k in np.unique(cell):
            t = T[cell == k]
            used, t = np.unique(t, return_inverse=True)
            q = px[used]
            x0, y0 = np.maximum(np.floor(q.min(0)).astype(int) - PAD, 0)
            x1 = min(int(np.ceil(q[:, 0].max())) + PAD, W)
            y1 = min(int(np.ceil(q[:, 1].max())) + PAD, H)
            out.append((E[used].astype(np.float32), (q - [x0, y0]).astype(np.float32),
                        t.reshape(-1, 3).astype(np.uint32), np.ascontiguousarray(img[y0:y1, x0:x1])))
    return out


def _one(job):
    try:
        return trimmed(*job)
    except Exception as e:  # 失敗したパノラマは飛ばす
        return f"{os.path.basename(job[1])}: {e}"


class Shelves:
    """切り出した画像を ATLAS 四方の画像に、棚を積むように詰める。"""

    def __init__(self):
        self.pages = []

    def put(self, h, w):
        if not self.pages or not self._fit(h, w):
            self.pages.append({"shelves": [], "top": 0})
            if not self._fit(h, w):
                raise ValueError(f"切り出しが大きすぎる: {w}x{h}")
        return self.last

    def _fit(self, h, w):
        page = self.pages[-1]
        for s in page["shelves"]:
            if h <= s["h"] and s["x"] + w <= ATLAS:
                self.last = (len(self.pages) - 1, s["x"], s["y"])
                s["x"] += w
                return True
        if page["top"] + h <= ATLAS and w <= ATLAS:
            page["shelves"].append({"y": page["top"], "h": h, "x": w})
            self.last = (len(self.pages) - 1, 0, page["top"])
            page["top"] += h
            return True
        return False


def main():
    from PIL import Image

    ap = argparse.ArgumentParser()
    ap.add_argument("mt7", nargs="*", help="mt/7 ファイル、またはそれが入ったフォルダ")
    ap.add_argument("out")
    ap.add_argument("--list", help="mt/7 のパスを 1 行 1 つ書いた一覧")
    ap.add_argument("--max-range", type=float, default=25.0, help="パノラマからの距離の上限 [m]")
    ap.add_argument("--max-edge", type=float, default=None, help="これより長い辺を持つ三角形を除く [m]")
    ap.add_argument("--min-facing", type=float, default=0.3,
                    help="視線に対する面の向き（余弦）の下限。物陰をまたぐ三角形を除く")
    ap.add_argument("--sky-range", type=float, default=10.0,
                    help="この距離より遠く、色が空らしい三角形を除く [m]。0 で無効")
    ap.add_argument("--sky-min", type=int, default=165, help="空とみなす色の暗い方の成分の下限")
    ap.add_argument("--overlap", type=float, default=0.1,
                    help="誤差が最良のこの割合増しまでなら残す（境目のすき間を埋める）")
    ap.add_argument("--candidates", type=int, default=16,
                    help="三角形ごとに比べる、水平に近いパノラマの数")
    ap.add_argument("--quality", type=int, default=90, help="アトラスの JPEG の画質")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    a = ap.parse_args()
    if a.list:
        a.mt7 = [ln.strip() for ln in open(a.list) if ln.strip()] + a.mt7
    a.mt7 = [f for p in a.mt7
             for f in (sorted(glob.glob(os.path.join(p, "*.bin"))) if os.path.isdir(p) else [p])]
    if not a.mt7:
        sys.exit("mt/7 ファイルがない")

    with mp.Pool(a.jobs) as pool:
        ecef = np.array(pool.map(pano_position, a.mt7, chunksize=16))
    mid = ecef.mean(0)
    R = enu_basis(mid)
    pos = (ecef - mid) @ R.T
    origin = mid + R.T @ [*((pos.min(0) + pos.max(0)) / 2)[:2], 0]   # 水平の中心、高さはパノラマの平均
    pos = (ecef - origin) @ R.T
    print(f"{len(a.mt7)} 枚のパノラマ、{np.ptp(pos[:, 0]):.0f} x {np.ptp(pos[:, 1]):.0f} m", flush=True)

    faces, done = [], 0
    with mp.Pool(a.jobs, initializer=_init, initargs=(pos, origin, R, a)) as pool:
        for k, res in enumerate(pool.imap_unordered(_one, enumerate(a.mt7), chunksize=4)):
            if isinstance(res, str):
                print(res)
                continue
            faces += res
            if (k + 1) % 100 == 0:
                print(f"  {k + 1}/{len(a.mt7)}  {sum(len(f[2]) for f in faces) / 1e6:.2f}M 三角形",
                      flush=True)
    if not faces:
        sys.exit("三角形が 1 つも残らなかった")

    # 高い順に詰めると棚のすき間が小さい
    faces.sort(key=lambda f: -f[3].shape[0])
    shelves = Shelves()
    placed = [(shelves.put(*f[3].shape[:2]), f) for f in faces]

    os.makedirs(a.out, exist_ok=True)
    groups, blobs, off = [], [], 0
    for page in range(len(shelves.pages)):
        mine = [(x, y, f) for (pg, x, y), f in placed if pg == page]
        h = shelves.pages[page]["top"]   # 使った高さまでで切る
        atlas = np.zeros((h, ATLAS, 3), np.uint8)
        for x, y, (_, _, _, img) in mine:
            atlas[y:y + img.shape[0], x:x + img.shape[1]] = img
        name = f"atlas{page}.jpg"
        Image.fromarray(atlas).save(os.path.join(a.out, name), quality=a.quality)
        base, idx = 0, []
        for _, _, (V, _, T, _) in mine:
            idx.append(T + base)
            base += len(V)
        V = np.vstack([f[0] for _, _, f in mine]).astype("<f4")
        uv = np.vstack([(f[1] + [x, y]) / [ATLAS, h] for x, y, f in mine]).astype("<f4")
        idx = np.vstack(idx).astype("<u4")
        g = {"file": name, "nv": len(V), "nt": len(idx)}
        for key, arr in (("pos", V), ("uv", uv), ("idx", idx)):
            g[key] = off
            blobs.append(arr.tobytes())
            off += len(blobs[-1])
        groups.append(g)
    with open(os.path.join(a.out, "mesh.bin"), "wb") as fh:
        for b in blobs:
            fh.write(b)
    radius = float(np.abs(pos[:, :2]).max()) + a.max_range
    json.dump({"count": len(a.mt7), "groups": groups, "radius": radius,
               "panos": [round(float(v), 3) for v in pos.ravel()]},
              open(os.path.join(a.out, "meta.json"), "w"))
    open(os.path.join(a.out, "index.html"), "w", encoding="utf-8").write(HTML)
    nt = sum(g["nt"] for g in groups)
    size = sum(os.path.getsize(os.path.join(a.out, f)) for f in os.listdir(a.out))
    print(f"{nt / 1e6:.2f}M 三角形、アトラス {len(groups)} 枚 -> {a.out}（{size / 1e6:.0f} MB）")
    print(f"見るには: python scripts/serve_web.py {a.out} 8000")


if __name__ == "__main__":
    main()
