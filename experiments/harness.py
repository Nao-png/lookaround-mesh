"""Score candidate vertex decodings: where in the stream the coordinates
start to run away (a decoding mistake accumulates; real geometry stays within
the extent the clean part already showed)."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SAMPLES = os.path.join(ROOT, "samples")

import lookaround_mesh as lamesh  # noqa: E402
from lookaround_mesh import edgebreaker as eb  # noqa: E402
import json  # noqa: E402

import numpy as np  # noqa: E402


def meshes(min_v=100):
    out = []
    for m in json.load(open(os.path.join(SAMPLES, "neighbours.json"))):
        for raw in lamesh.meshes(open(os.path.join(SAMPLES, f"pano{m['k']}_mt7.bin"), "rb").read()):
            c = lamesh.chnks(raw)
            _, nv, nf = lamesh.header(raw)
            if nv < min_v:
                continue
            s = c[4][24:24 + nf].decode()
            b = nv - s.count("C")
            tris, gates, _ = eb.decode(s, b)
            out.append({"nv": nv, "b": b, "s": s, "tris": tris, "gates": gates,
                        "vals": lamesh.vertex_values(raw)})
    return out


def runaway(P, order, lead=0.25, factor=4.0):
    """Fraction of the stream decoded before any coordinate leaves `factor`
    times the extent of the first `lead` of it (per component)."""
    Q = np.array([P[v] for v in order], float)
    n0 = max(10, int(lead * len(Q)))
    lo, hi = Q[:n0].min(0), Q[:n0].max(0)
    ext = np.maximum(hi - lo, 1.0)
    c = (lo + hi) / 2
    out = np.abs(Q - c) > factor * ext
    bad = np.where(out.any(1))[0]
    first = bad[0] / len(Q) if len(bad) else 1.0
    per_axis = [(np.where(out[:, k])[0][0] / len(Q) if out[:, k].any() else 1.0) for k in range(3)]
    return first, per_axis
