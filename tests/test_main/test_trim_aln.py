"""
Test alignment trimming
"""

from arbow.__main__ import trim_aln


def test_trim_aln(gen_aln):
    aln, *_ = gen_aln
    trim = trim_aln(aln, "ref", 3, 16)
    print(trim)
    assert 0