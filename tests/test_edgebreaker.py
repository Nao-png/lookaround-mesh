import collections
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import lookaround_mesh as lm  # noqa: E402


def edges(tris):
    und, dire = collections.Counter(), collections.Counter()
    for t in tris:
        for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            dire[e] += 1
            und[frozenset(e)] += 1
    return und, dire


def test_a_fan_around_one_interior_vertex():
    # 境界 0-1-2-3 の四角形の中に内側の頂点 4 が 1 つ
    need, end, _ = lm.split_lengths("CRRE")
    assert (need, end) == (4, 4)
    tris, gates, n = lm.decode("CRRE", 4)
    assert tris == [(0, 1, 4), (4, 1, 2), (4, 2, 3), (4, 3, 0)]
    assert gates == {4: (0, 1)} and n == 5
    und, dire = edges(tris)
    assert max(dire.values()) == 1
    assert sorted(tuple(sorted(e)) for e, k in und.items() if k == 1) == [(0, 1), (0, 3), (1, 2), (2, 3)]


def test_a_split_goes_right_first():
    # 境界 6 頂点、S が 3 頂点のループを右に切り出して E で閉じ、残りを R E で閉じる
    need, end, lens = lm.split_lengths("SERE")
    assert need == 6 and lens == {0: 3}
    tris, _, n = lm.decode("SERE", 6)
    assert tris == [(0, 1, 3), (3, 1, 2), (0, 3, 4), (0, 4, 5)]
    und, dire = edges(tris)
    assert max(dire.values()) == 1 and sum(1 for k in und.values() if k == 1) == 6


def test_a_string_that_does_not_fit_its_boundary_is_refused():
    with pytest.raises(ValueError):
        lm.decode("CRRE", 5)


def test_first_touch_starts_with_the_first_gate():
    tris, _, _ = lm.decode("SERE", 6)
    order, gate_of = lm.first_touch(tris)
    assert order[:3] == [0, 1, 3] and gate_of[0] is None and gate_of[3] == (0, 1)
    assert sorted(order) == list(range(6))


SAMPLES = os.path.join(ROOT, "samples")


@pytest.mark.skipif(not os.path.exists(os.path.join(SAMPLES, "neighbours.json")), reason="no samples")
def test_every_claim_holds_on_the_samples():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import check
    assert check.main(SAMPLES)
