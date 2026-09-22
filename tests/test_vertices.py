import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from lookaround_mesh import vertices as vx  # noqa: E402
import lookaround_mesh as lm  # noqa: E402


def encode(tris, b, P):
    """decode_vertices の逆: 座標 P から値の列（軸ごと）を作る。走査は同じ。"""
    ids = vx.apple_ids(tris, b)
    vert, opp = vx.corner_table(tris)
    val = np.zeros_like(P)
    for c in (0, 1, 2):
        val[ids[vert[c]]] = P[vert[c]]
    seen, stack, done = [True] + [False] * (len(tris) - 1), [], set(vert[:3])

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
        if k not in done:
            s = opp[c]
            val[ids[k]] = P[k] - (P[vert[vx._next(s)]] + P[vert[vx._prev(s)]] - P[vert[s]])
            done.add(k)
        push(c // 3)
    return val


def test_boundary_numbering_is_a_permutation_that_starts_on_the_first_triangle():
    for s, b in (("CRRE", 4), ("SERE", 6), ("CCRCRCSERCRCRRCSERCRSECRRRSERE", 14)):
        tris, _, n = lm.decode(s, b)
        ids = vx.apple_ids(tris, b)
        assert sorted(ids.values()) == list(range(n))
        assert all(ids[v] == v for v in range(b, n))  # 内側は C の順のまま
        assert {0, b - 1} <= {ids[v] for v in tris[0]}


def test_values_round_trip_through_the_traversal(monkeypatch):
    s, b = "CCRCRCSERCRCRRCSERCRSECRRRSERE", 14
    tris, _, n = lm.decode(s, b)
    P = np.random.default_rng(0).integers(-500, 500, (n, 3))
    val = encode(tris, b, P)
    monkeypatch.setattr(vx, "header", lambda m: (None, n, len(tris)))
    monkeypatch.setattr(vx, "vertex_values", lambda m: val.T.reshape(-1))
    assert (vx.decode_vertices(None, tris, b) == P).all()


def test_depth_and_angles_become_metres():
    # 経度 0・緯度 0・w = 100 → 前方 10 m。緯度 +90°（下）→ z = +10 m
    X = vx.to_local(np.array([[0, 0, 100], [0, 1571, 100]]))
    assert np.allclose(X[0], [10, 0, 0]) and np.allclose(X[1], [0, 0, 10], atol=1e-2)
