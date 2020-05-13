"""
Test seq2series function from __main__
"""

from arbow.__main__ import seq2series


def test_seq2series(seq_gen):
    a = 2
    c = 5
    g = 3
    t = 8
    n = 3
    seq = seq_gen(a, c, g, t, n)
    seq_series = seq2series(seq, 'seq1')
    assert seq_series.name == 'seq1'
    assert len(seq_series) == len(seq)
    for i, b in enumerate(seq):
        assert seq_series[i] == b
