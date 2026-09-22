"""Edgebreaker の接続情報（CLERS 文字列）の復元。

各面のメッシュは位相的に円盤で、その境界（b = V - Cの数 個の頂点）が最初の
アクティブループになる。境界の頂点を 0..b-1（ループ順）、C が生むたびに
内側の頂点 b, b+1, ... と番号を振る。ループはリストで、先頭の辺
(loop[0], loop[1]) がゲート。S のあとは右側を先に処理する。

この約束で 48/48 のメッシュについて、三角形数・頂点数が一致し、潰れた三角形
がなく、どの辺も 2 枚以下の三角形に属し、向きが一貫し、境界の辺がちょうど
b 本になる（L と R を入れ替えるとこれらが全部崩れる）。
"""
import sys

from .container import CHNK_CLERS, chnks, header

sys.setrecursionlimit(100000)


def clers(m):
    """(CLERS 文字列, 境界の頂点数 b)。種類 4 の本体 = 24 バイトの見出し
    (u32 F, u32 V, u32 F, 0, u32 12, 0)、F 文字、12 バイトの末尾。"""
    _, nv, nf = header(m)
    s = chnks(m)[CHNK_CLERS][24:24 + nf].decode("ascii")
    return s, nv - s.count("C")


def split_lengths(s):
    """各 S が先に切り出すループの長さと、文字列全体が必要とする最初のループの
    長さ。文字列だけから決まる: 分岐の中で C は +1、L と R は -1、入れ子の S は
    1 - (切り出す長さ)、E は長さ 3 のループを閉じる。

    Returns (必要な最初のループの長さ, 読み終えた位置, {S の位置: 長さ})。"""
    lens = {}
    need, end = _walk(s, lens)
    return need, end, lens


def _walk(s, lens):
    """split_lengths の本体（明示的なスタックで、深い入れ子でも落ちない）。"""
    # 各フレーム: [この分岐の delta, この分岐を開いた S の位置 or None]
    stack = [[0, None]]
    i = 0
    while i < len(s):
        op = s[i]
        top = stack[-1]
        if op == "C":
            top[0] += 1
        elif op in "LR":
            top[0] -= 1
        elif op == "S":
            stack.append([0, i])
        elif op == "E":
            need = 3 - top[0]
            stack.pop()
            if not stack:
                return need, i + 1
            if top[1] is not None:
                lens[top[1]] = need
                stack[-1][0] += 1 - need
        else:
            raise ValueError(f"bad op {op!r} at {i}")
        i += 1
    raise ValueError("CLERS string ended before its last E")


def decode(s, b):
    """(三角形の頂点番号の組のリスト, 内側の頂点ごとのゲート (a, b), 頂点数)。"""
    need, end, lens = split_lengths(s)
    if need != b or end != len(s):
        raise ValueError(f"string needs an initial loop of {need}, the mesh has {b}")
    loop = list(range(b))
    stack, tris, gates, nxt = [], [], {}, b
    for i, op in enumerate(s):
        v0, v1 = loop[0], loop[1]
        if op == "C":
            x = nxt
            nxt += 1
            gates[x] = (v0, v1)
            tris.append((v0, v1, x))
            loop = [x] + loop[1:] + [v0]
        elif op == "L":
            x = loop[-1]
            tris.append((v0, v1, x))
            loop = [x] + loop[1:-1]
        elif op == "R":
            x = loop[2]
            tris.append((v0, v1, x))
            loop = [v0] + loop[2:]
        elif op == "S":
            k = lens[i]
            x = loop[k]
            tris.append((v0, v1, x))
            stack.append([v0] + loop[k:])
            loop = [x] + loop[1:k]
        elif op == "E":
            if len(loop) != 3:
                raise ValueError(f"E at {i} on a loop of {len(loop)}")
            tris.append((v0, v1, loop[2]))
            loop = stack.pop() if stack else []
    return tris, gates, nxt


def first_touch(tris):
    """走査で初めて触れた順の頂点と、各頂点を初めて触れた三角形のゲート。
    最初のゲートの 2 頂点のあとは、新しく触れる頂点は必ず三角形の 3 つ目。"""
    order, seen, gate_of = [], set(), {}
    a0, b0, _ = tris[0]
    for v in (a0, b0):
        order.append(v)
        seen.add(v)
        gate_of[v] = None
    for a, b, x in tris:
        if x not in seen:
            order.append(x)
            seen.add(x)
            gate_of[x] = (a, b)
    return order, gate_of
