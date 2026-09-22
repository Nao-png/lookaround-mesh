"""頂点ストリーム（CHNK 種類 9）。

解けているのはエントロピー層だけ: 見出し 146 バイト（u32 V、f32 0.001、
f32 0.001、f32 1000、u16 0、地心→局所の 4x4 行列）、全メッシュ共通の 40
バイト、そのあと V×3 バイトの本体（値 = バイト - 127、0xff は「例外表から
次の値を取る」印）と、本体の 0xff と同じ数の int32 の例外表。

値から座標を組み立てる規則（予測の方式と座標の意味）は未解明。
`reconstruct_experimental` は現時点で最も筋の良い試み（初めて触れた順＋
平行四辺形予測）で、正しい座標ではない（docs/HYPOTHESES.md）。
"""
import numpy as np

from .container import CHNK_VERTICES, chnks, header

BODY = 146 + 40  # 種類 9 の本体の中で、1 バイト値の列が始まる位置
ESCAPE = 0xFF
BIAS = 127


def vertex_values(m):
    """V×3 個の整数（予測との差とみられる）。48/48 のメッシュで、例外表の長さ
    が本体の 0xff の数と一致する。"""
    _, nv, _ = header(m)
    p = np.frombuffer(chnks(m)[CHNK_VERTICES][BODY:], np.uint8)
    body = p[:3 * nv].astype(np.int64) - BIAS
    esc = p[:3 * nv] == ESCAPE
    exc = np.frombuffer(p[3 * nv:].tobytes(), "<i4")
    if len(exc) != int(esc.sum()):
        raise ValueError(f"{len(exc)} exception values for {int(esc.sum())} escapes")
    body[esc] = exc
    return body


def reconstruct_experimental(res, tris, order, gate_of):
    """未完成。頂点ごとに格納 (V, 3) の残差を、初めて触れた順に、ゲート (a, b)
    とその向こう側の頂点 o から a + b - o で予測して積み上げる。境界は範囲内に
    収まるが、内側は早い段階から境界の 5〜20 倍に広がるので、この規則は正しくない。"""
    P, known, ti = {}, {}, 0
    for i, v in enumerate(order):
        g = gate_of[v]
        if i == 0:
            P[v] = res[0].copy()
            continue
        if g is None:
            P[v] = P[order[i - 1]] + res[i]
            continue
        while ti < len(tris) and all(u in P for u in tris[ti]):
            t = tris[ti]
            for k in range(3):
                known[(t[k], t[(k + 1) % 3])] = t[(k + 2) % 3]
            ti += 1
        a, b = g
        o = known.get((b, a))
        pred = (P[a] + P[b]) // 2 if o is None else P[a] + P[b] - P[o]
        P[v] = pred + res[i]
    return np.array([P[i] for i in range(len(order))])
