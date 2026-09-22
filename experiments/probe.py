"""Swap the prediction rule at one stream position and see how far the clean
decode then reaches."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SAMPLES = os.path.join(ROOT, "samples")

import lookaround_mesh as lamesh  # noqa: E402
from lookaround_mesh import edgebreaker as eb  # noqa: E402
import collections  # noqa: E402

import numpy as np  # noqa: E402



def decode(m, override=None):
    """override: (stream index, rule name) for one vertex; rules use the
    gate (a, b), the vertex o across it, and the previous stream vertex."""
    nv, tris = m["nv"], m["tris"]
    order, gate_of = eb.first_touch(tris)
    res = m["vals"].reshape(nv, 3)
    P, known, ti = {}, {}, 0
    for i, v in enumerate(order):
        g = gate_of[v]
        if i == 0:
            P[v] = res[0].copy(); continue
        if g is None:
            P[v] = P[order[i - 1]] + res[i]; continue
        while ti < len(tris) and all(u in P for u in tris[ti]):
            t = tris[ti]
            for k in range(3):
                known[(t[k], t[(k + 1) % 3])] = t[(k + 2) % 3]
            ti += 1
        a, b = g
        o = known.get((b, a))
        rule = override[1] if (override and override[0] == i) else "para"
        A, B, prev = P[a], P[b], P[order[i - 1]]
        O = P[o] if o is not None else (A + B) // 2
        pred = {"para": A + B - O if o is not None else (A + B) // 2,
                "mid": (A + B) // 2, "a": A, "b": B, "o": O, "prev": prev,
                "2a-b": 2 * A - B, "2b-a": 2 * B - A, "zero": np.zeros(3, np.int64),
                "a+o-b": A + O - B, "b+o-a": B + O - A}[rule]
        P[v] = pred + res[i]
    return P, order


def first_bad(P, order, tris, nv, start=0, k=6.0):
    X = np.array([P[u] for u in range(nv)], float)
    nb = collections.defaultdict(set)
    for t in tris:
        for x in t:
            nb[x] |= set(t) - {x}
    L = [np.linalg.norm(X[u] - X[w]) for t in tris for u, w in ((t[0], t[1]), (t[1], t[2]))]
    med = np.median(L)
    for i in range(start, len(order)):
        u = order[i]
        if min(np.linalg.norm(X[u] - X[w]) for w in nb[u]) > k * med:
            return i
    return len(order)
