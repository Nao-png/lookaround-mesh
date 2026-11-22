"""samples/ の mt/7 ファイルで、解けている層の主張をすべて検証する。

    python scripts/check.py [samples ディレクトリ]

各項目を「満たしたメッシュ（またはパノラマ）の数 / 全体」で表示する。
"""
import collections
import json
import os
import struct
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import lookaround_mesh as lm  # noqa: E402
from lookaround_mesh.container import (  # noqa: E402
    CHNK_CLERS, CHNK_FACE_ID, CHNK_PARAMS, CHNK_PROJECTION, CHNK_VERTICES, FACE_INDEX, INDEX_TABLE)
from lookaround_mesh.vertices import BODY, ESCAPE  # noqa: E402


def nearest(A, B, chunk=256):
    """A の各点から B の最も近い点までの距離（scipy なしの総当たり）。"""
    out = []
    for i in range(0, len(A), chunk):
        d = ((A[i:i + chunk, None, :] - B[None, :, :]) ** 2).sum(-1)
        out.append(np.sqrt(d.min(1)))
    return np.concatenate(out)


def surface_samples(pano, centre, radius=30.0, longest=3.0, n=6):
    """半径内の三角形の面上に点をまく。物陰をまたぐ長い三角形（最長辺 > longest m）は除く。"""
    bary = np.array([(i / n, j / n) for i in range(n + 1) for j in range(n + 1 - i)])
    out = []
    for raw, tris, P in pano:
        ok = P[:, 2] > 0
        E = lm.to_ecef(raw, np.where(ok[:, None], P, 1))
        T = np.array([t for t in tris if ok[list(t)].all()])
        if not len(T):  # 空の面は全頂点が無限遠
            continue
        a, b, c = E[T[:, 0]], E[T[:, 1]], E[T[:, 2]]
        edge = np.maximum.reduce([np.linalg.norm(a - b, axis=1), np.linalg.norm(b - c, axis=1),
                                  np.linalg.norm(c - a, axis=1)])
        keep = (edge < longest) & (np.linalg.norm(a - centre, axis=1) < radius)
        a, b, c = a[keep], b[keep], c[keep]
        out.append((a[:, None] + bary[None, :, :1] * (b - a)[:, None] + bary[None, :, 1:] * (c - a)[:, None]).reshape(-1, 3))
    return np.vstack(out)


def geometry(faces, meta, check):
    """座標を地心座標に置いたときの主張: 面ごとの範囲、地面の平面、別走行との重なり。"""
    ext = collections.defaultdict(list)
    for k, pano in faces.items():
        for f, (raw, tris, P) in enumerate(pano):
            ext[f].append(np.r_[P[:, :2].min(0), P[:, :2].max(0)])
    for f, e in ext.items():
        check("座標: 面ごとの u, v の範囲が 8 枚で共通（差 2 以内）", np.ptp(np.array(e), axis=0).max() <= 2)

    centre = {}
    for k, pano in faces.items():
        raw, tris, P = pano[5]  # 地面の面
        E = lm.to_ecef(raw, P[P[:, 2] > 0])
        c = lm.header(raw)[0][:3, 3]
        up = c / np.linalg.norm(c)
        G = E[np.linalg.norm(E - c, axis=1) < 3]
        h = (G - c) @ up
        check("座標: カメラ高（半径 3 m の地面の高さの中央値）が 2.2〜2.4 m", 2.2 < -np.median(h) < 2.4)
        F = G[np.abs(h - np.median(h)) < 0.1]  # 人や車など地面から浮いた点を除く
        _, S, Vt = np.linalg.svd(F - F.mean(0))
        n = Vt[2]
        check("座標: 半径 3 m の地面が平面（RMS 4 cm 未満）で法線が鉛直から 3° 以内",
              S[2] / np.sqrt(len(F)) < 0.04 and np.degrees(np.arccos(abs(n @ up))) < 3)
        centre[k] = (c, up)

    c0, up = centre[0]
    A = np.vstack([lm.to_ecef(raw, P[P[:, 2] > 0]) for raw, tris, P in faces[0]])
    A = A[(np.linalg.norm(A - c0, axis=1) < 20) & ((A - c0) @ up > -1.3)]  # 20 m 以内、地面より 1 m 以上上
    for m in meta:
        if m["same_drive"]:
            continue
        cB, upB = centre[m["k"]]
        S = surface_samples(faces[m["k"]], cB)
        K = np.array([[0, -upB[2], upB[1]], [upB[2], 0, -upB[0]], [-upB[1], upB[0], 0]])
        R = np.eye(3) + K + K @ K  # 鉛直軸まわりに 90° 回した対照
        d, dc = nearest(A, S), nearest(A, (S - cB) @ R.T + cB)
        check("座標: 別走行のパノラマと重なる（建物の点の中央値 0.5 m 未満、90° 回した対照は 2 m 超）",
              np.median(d) < 0.5 and np.median(dc) > 2)


def main(samples):
    meta = json.load(open(os.path.join(samples, "neighbours.json")))
    tally = collections.Counter()
    total = collections.Counter()

    def check(name, ok):
        total[name] += 1
        tally[name] += bool(ok)

    params = set()
    origin = None
    faces = {}
    for m in meta:
        b = open(os.path.join(samples, f"pano{m['k']}_mt7.bin"), "rb").read()
        recs = lm.mcp4_records(b)  # 表の直後から順にすき間なく続くことは中で照合済み
        check("MCP4: 0x08 からの 20 記録がファイルをすき間なく覆う",
              len(recs) == 20 and recs[0][1] == 0x08 + 10 * 20 and recs[-1][1] + recs[-1][2] == len(b))
        check("MCP4: 記録の並び（種類 5、メッシュ 6、種類 2・3 の組 6、種類 0）",
              "".join(str(t) for t, _, _ in recs) == "5" + "1" * 6 + "23" * 6 + "0")
        faces_of = list(range(6)) + [f for f in range(6) for _ in (2, 3)]
        idx = next(b[o:o + s] for t, o, s in recs if t == INDEX_TABLE)
        check("記録 0: 117 バイト、9 バイトの後の 6 バイト項目 18 個の先頭が各記録の面番号",
              len(idx) == 117 and idx[:9] == bytes(9)
              and [idx[9 + 6 * i] for i in range(18)] == faces_of
              and all(idx[9 + 6 * i + 1:9 + 6 * i + 6] == bytes(5) for i in range(18)))
        check("種類 2: 1 バイト目が面番号",
              [b[o + 1] for t, o, s in recs if t == FACE_INDEX] == list(range(6)))
        info = lm.pano_info(b)
        check("種類 5: パノラマ id が一覧の id と一致", info["panoid"] == m["id"])
        check("種類 5: 定数 0.01・180・127、値の数 180",
              abs(info["scale"] - 0.01) < 1e-7 and len(info["values"]) == 180)
        check("種類 5: 180 個の値の累積和が 0〜100 に収まる（サンプルはすべて 0）",
              ((info["profile"] >= 0) & (info["profile"] <= 100)).all())
        meshes = lm.meshes(b)  # LZFSE の展開後の大きさは中で照合済み
        check("LZFSE: 6 面とも宣言どおりの大きさに展開", len(meshes) == 6)
        mats = []
        faces[m["k"]] = []
        for f, raw in enumerate(meshes):
            mat, nv, nf = lm.header(raw)
            mats.append(mat)
            rot = mat[:3, :3]
            check("行列: 回転が正規直交（誤差 1e-9 未満）で行列式 +1",
                  np.abs(rot @ rot.T - np.eye(3)).max() < 1e-9 and abs(np.linalg.det(rot) - 1) < 1e-9)
            c = lm.chnks(raw)
            params.add(c[CHNK_PARAMS])
            check("種類 0: 行列の直前の u64 が面番号", struct.unpack_from("<Q", c[0], 22)[0] == f)
            check("種類 10: 先頭 u32 1・u16 1・u16 面番号、続く 134 バイトが種類 9 の見出しの写し",
                  struct.unpack_from("<IHH", c[CHNK_PROJECTION], 0) == (1, 1, f)
                  and c[CHNK_PROJECTION][8:142] == c[CHNK_VERTICES][16:150])
            check("種類 9: 投影の位置は恒等 (1, 1, 0, 0)",
                  struct.unpack_from("<4d", c[CHNK_VERTICES], 150) == (1.0, 1.0, 0.0, 0.0))
            check("種類 4: 末尾 12 バイトが (F, V, 1)", struct.unpack_from("<3I", c[CHNK_CLERS], 24 + nf) == (nf, nv, 1))
            check("種類 5: u32 F、u32 1、F バイトの (面番号 + 1)",
                  struct.unpack_from("<2I", c[CHNK_FACE_ID], 0) == (nf, 1) and len(c[CHNK_FACE_ID]) == 8 + nf
                  and set(c[CHNK_FACE_ID][8:]) == {f + 1})
            s, bnd = lm.clers(raw)
            cnt = collections.Counter(s)
            check("CLERS: C/L/E/R/S の 5 文字だけで長さ F", len(s) == nf and set(s) <= set("CLERS"))
            check("CLERS: E の数 = S の数 + 1", cnt["E"] == cnt["S"] + 1)
            need, end, _ = lm.split_lengths(s)
            check("最初のループの長さ（文字列から計算）= V - Cの数", need == bnd and end == len(s))
            tris, _, n = lm.decode(s, bnd)
            und, dire = collections.Counter(), collections.Counter()
            for t in tris:
                for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                    dire[e] += 1
                    und[frozenset(e)] += 1
            boundary = [e for e, k in und.items() if k == 1]
            check("復元: 三角形数と頂点数が一致", len(tris) == nf and n == nv)
            check("復元: 潰れた三角形がない", all(len(set(t)) == 3 for t in tris))
            check("復元: どの辺も 2 枚以下、向きが一貫", max(und.values()) <= 2 and max(dire.values()) == 1)
            check("復元: 境界の辺が b 本、境界の頂点が 0..b-1",
                  len(boundary) == bnd and {v for e in boundary for v in e} == set(range(bnd)))
            body = np.frombuffer(c[CHNK_VERTICES][BODY:BODY + 3 * nv], np.uint8)
            n_exc = (len(c[CHNK_VERTICES]) - BODY - 3 * nv) / 4
            check("頂点: 例外表の長さ = 本体の 0xff の数", n_exc == int((body == ESCAPE).sum()))
            P = lm.decode_vertices(raw, tris, bnd)
            val = lm.vertex_values(raw).reshape(3, nv).T
            ids = lm.apple_ids(tris, bnd)
            check("座標: 最初の三角形の 3 頂点は値そのもの（番号 0 と b-1 を含む）",
                  all((P[j] == val[ids[j]]).all() for j in tris[0]) and {0, bnd - 1} <= {ids[j] for j in tris[0]})
            check("座標: 奥行きの逆数 w がすべて 0 以上", (P[:, 2] >= 0).all())
            su, sv, ou, ov = lm.projection(raw)
            ndc = P[:, :2] * 1e-3 * (su, sv) + (ou, ov)
            check("種類 10 の投影: 頂点の正規化座標が [-1, 1] に収まり、両端まで 1% 以内",
                  (np.abs(ndc) <= 1.0).all() and (np.abs(ndc.min(0) + 1) < 0.01).all() and (np.abs(ndc.max(0) - 1) < 0.01).all())
            faces[m["k"]].append((raw, tris, P))
        T = np.array([x[:3, 3] for x in mats])
        check("位置: 6 面が同じ平行移動を共有（1 mm 未満）", np.abs(T - T[0]).max() < 1e-3)
        check("種類 5: ECEF 位置が種類 0 の平行移動と一致（1 mm 未満）", np.abs(info["ecef"] - T[0]).max() < 1e-3)
        origin = T[0] if origin is None else origin
        check("位置: パノラマ間の距離が GPS と 0.2 m 以内", abs(np.linalg.norm(T[0] - origin) - m["dist_m"]) < 0.2)
    check("種類 11: 全メッシュで同一", len(params) == 1)
    geometry(faces, meta, check)

    width = max(len(k) for k in total)
    for k in total:
        mark = "OK " if tally[k] == total[k] else "NG "
        print(f"{mark}{k:<{width}}  {tally[k]}/{total[k]}")
    return all(tally[k] == total[k] for k in total)


if __name__ == "__main__":
    ok = main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "samples"))
    sys.exit(0 if ok else 1)
