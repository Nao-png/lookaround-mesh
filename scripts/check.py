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
from lookaround_mesh.container import CHNK_HEADER_COPY, CHNK_PARAMS, CHNK_VERTICES  # noqa: E402
from lookaround_mesh.vertices import BODY, ESCAPE  # noqa: E402


def main(samples):
    meta = json.load(open(os.path.join(samples, "neighbours.json")))
    tally = collections.Counter()
    total = collections.Counter()

    def check(name, ok):
        total[name] += 1
        tally[name] += bool(ok)

    params = set()
    origin = None
    for m in meta:
        b = open(os.path.join(samples, f"pano{m['k']}_mt7.bin"), "rb").read()
        recs = lm.mcp4_records(b)
        check("MCP4: 記録がファイルをすき間なく覆う", recs[0][1] == 0x1b0 and recs[-1][1] + recs[-1][2] == len(b))
        meshes = lm.meshes(b)  # LZFSE の展開後の大きさは中で照合済み
        check("LZFSE: 6 面とも宣言どおりの大きさに展開", len(meshes) == 6)
        mats = []
        for raw in meshes:
            mat, nv, nf = lm.header(raw)
            mats.append(mat)
            rot = mat[:3, :3]
            check("行列: 回転が正規直交（誤差 1e-9 未満）で行列式 +1",
                  np.abs(rot @ rot.T - np.eye(3)).max() < 1e-9 and abs(np.linalg.det(rot) - 1) < 1e-9)
            c = lm.chnks(raw)
            params.add(c[CHNK_PARAMS])
            check("種類 10: 種類 9 の見出し（134 バイト）の写し", c[CHNK_HEADER_COPY][8:142] == c[CHNK_VERTICES][16:150])
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
        T = np.array([x[:3, 3] for x in mats])
        check("位置: 6 面が同じ平行移動を共有（1 mm 未満）", np.abs(T - T[0]).max() < 1e-3)
        origin = T[0] if origin is None else origin
        check("位置: パノラマ間の距離が GPS と 0.2 m 以内", abs(np.linalg.norm(T[0] - origin) - m["dist_m"]) < 0.2)
    check("種類 11: 全メッシュで同一", len(params) == 1)

    width = max(len(k) for k in total)
    for k in total:
        mark = "OK " if tally[k] == total[k] else "NG "
        print(f"{mark}{k:<{width}}  {tally[k]}/{total[k]}")
    return all(tally[k] == total[k] for k in total)


if __name__ == "__main__":
    ok = main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "samples"))
    sys.exit(0 if ok else 1)
