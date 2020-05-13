"""
Test fasta2df function
"""

from arbow.__main__ import fasta2df


def test_fasta2df_cols(fasta_file):
    """Check the number of columns is correct"""
    df = fasta2df(fasta_file)
    assert df.shape[1] == 10


def test_fasta2df_rows(fasta_file):
    """Check the number of rows is correct"""
    df = fasta2df(fasta_file)
    assert df.shape[0] == 10

