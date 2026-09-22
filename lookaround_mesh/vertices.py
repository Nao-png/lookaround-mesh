"""頂点ストリーム（CHNK 種類 9）と、そこから座標を組み立てる規則。

エントロピー層: 見出し 146 バイト（u32 V、f32 0.001、f32 0.001、f32 1000、
u16 0、地心→局所の 4x4 行列）、全メッシュ共通の 40 バイト、そのあと V×3
バイトの本体（値 = バイト - 127、0xff は「例外表から次の値を取る」印）と、
本体の 0xff と同じ数の int32 の例外表。

本体は**軸ごと**に並ぶ（V 個の u、V 個の v、V 個の w）。各軸の中は Apple の
頂点番号の順。座標は面ごとの球面座標で、u = 経度、v = 緯度（下向きが正）、
単位はミリラジアン。w は距離の逆数で、距離 = 1000 / w メートル（空は w = 0）。

組み立て（docs/FORMAT.md）:
- 頂点番号: 内側の頂点は C の順に b, b+1, …（edgebreaker.decode と同じ）。
  境界の頂点は、最初の三角形の角 1 から境界を歩いて b-1 から数え下げる
  （`apple_ids`）。
- 最初の三角形の 3 頂点は値そのもの。
- 残りは面の走査で決まる。各三角形の角を (先端, a, b) = 0, 1, 2 とし、面を
  取り出すたびに角 0, 1, 2 の順に、その角の向かいの辺の先にある未訪問の面を
  スタックに積む（積んだ時点で訪問済み）。積んだ面を取り出したとき、その先端の
  頂点がまだなら、積んだ側の面から平行四辺形で予測し、値を足す。

Apple Flyover の C3M 形式の Edgebreaker 変種（retroplasma/flyover-reverse-
engineering）と同じ系統だが、値を引かずに足す点と、頂点を積む時ではなく
取り出す時に決める点が違う。
"""
import numpy as np

from .container import CHNK_VERTICES, chnks, header

BODY = 146 + 40  # 種類 9 の本体の中で、1 バイト値の列が始まる位置
ESCAPE = 0xFF
BIAS = 127
ANGLE = 1e-3     # u, v の単位（ラジアン）
DEPTH = 1000.0   # 距離 = DEPTH / w（メートル）


def vertex_values(m):
    """V×3 個の整数（軸ごと: u×V, v×V, w×V）。48/48 のメッシュで、例外表の長さ
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


def _next(c):
    return 3 * (c // 3) + (c + 1) % 3


def _prev(c):
    return 3 * (c // 3) + (c + 2) % 3


def corner_table(tris):
    """(角 → 頂点, 角 → 向かいの角 or -1)。三角形 (a, b, x) の角 0, 1, 2 は x, a, b。"""
    vert = [0] * (3 * len(tris))
    for w, (a, b, x) in enumerate(tris):
        vert[3 * w:3 * w + 3] = x, a, b
    edge = {(vert[_next(c)], vert[_prev(c)]): c for c in range(len(vert))}
    opp = [edge.get((vert[_prev(c)], vert[_next(c)]), -1) for c in range(len(vert))]
    return vert, opp


def apple_ids(tris, b):
    """{edgebreaker.decode の頂点番号: 値の列での番号}。内側はそのまま。境界は、
    最初の三角形の角 1 の頂点の扇を回って境界の辺に出てから、境界を歩いて
    b-1, b-2, … と振る。"""
    vert, opp = corner_table(tris)
    vt = [v if v >= b else -1 for v in vert]
    c = _next(1)
    while True:
        r = _next(c)
        c = opp[r]
        if c < 0:
            break
    n = b - 1
    while True:
        vt[_next(r)] = n
        r = _prev(r)
        o = opp[r]
        while o >= 0:  # 同じ頂点の扇を回る
            vt[_next(o)] = n
            r = _prev(o)
            o = opp[r]
        n -= 1
        if vt[_next(r)] != -1:
            break
    ids = {vert[c]: vt[c] for c in range(len(vert))}
    if n != -1 or sorted(ids.values()) != list(range(len(ids))):
        raise ValueError("boundary walk did not number every boundary vertex once")
    return ids


def decode_vertices(m, tris, b):
    """(V, 3) の整数座標 (u, v, w)。行は edgebreaker.decode の頂点番号の順。"""
    _, nv, _ = header(m)
    val = vertex_values(m).reshape(3, nv).T
    ids = apple_ids(tris, b)
    vert, opp = corner_table(tris)
    P = np.zeros((nv, 3), np.int64)
    done = np.zeros(nv, bool)
    for c in (0, 1, 2):
        P[vert[c]] = val[ids[vert[c]]]
        done[vert[c]] = True
    seen = [False] * len(tris)
    seen[0] = True
    stack = []

    def push(f):
        for c in range(3 * f, 3 * f + 3):
            o = opp[c]
            if o >= 0 and not seen[o // 3]:
                seen[o // 3] = True
                stack.append(o)

    push(0)
    while stack:
        c = stack.pop()
        k = vert[c]
        if not done[k]:
            s = opp[c]  # 積んだ側の面で、共有する辺の向かいの角
            P[k] = P[vert[_next(s)]] + P[vert[_prev(s)]] - P[vert[s]] + val[ids[k]]
            done[k] = True
        push(c // 3)
    if not done.all():
        raise ValueError("traversal did not reach every vertex")
    return P


def to_local(P):
    """(u, v, w) → 面の局所座標 (x 前, y 横, z 下) のメートル。w = 0（空）は inf。"""
    lon, lat = P[:, 0] * ANGLE, P[:, 1] * ANGLE
    with np.errstate(divide="ignore"):
        r = DEPTH / P[:, 2].astype(float)
    d = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], 1)
    return d * r[:, None]


def to_ecef(m, P):
    """(u, v, w) → 地心座標 ECEF（種類 0 の行列）。w = 0 の頂点は inf になる。"""
    mat, _, _ = header(m)
    with np.errstate(invalid="ignore"):
        return to_local(P) @ mat[:3, :3].T + mat[:3, 3]
