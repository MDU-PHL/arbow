"""
Test count bases
"""

from arbow.__main__ import count_base


def test_count_base(seq_gen):
    a = 5
    c = 10
    g = 8
    t = 2
    n = 0
    seq = seq_gen(a, c, g, t, n)
    for base, count in zip(list('acgtn'), [a, c, g, t, n]):
        assert count_base(seq, base) == count
