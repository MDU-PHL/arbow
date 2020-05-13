"""
Configure tests and fixtures
"""

import random
import pytest
import pandas as pd


@pytest.fixture
def fasta_file(tmp_path, n_seqs=10, seq_len=10, proportion_missing=0, seed=42):
    random.seed(seed)
    fasta_file = tmp_path / 'seqs.fa'
    seqs = ""
    if proportion_missing == 0:
        cum_weights = [25, 50, 75, 100, 100]
    else:
        missing_weight = 100 * proportion_missing
        base_weights = (100 - missing_weight)/4
        cum_weights = [i * base_weights for i in range(1, 5)] + [100]
    for i in range(0, n_seqs):
        seq = "".join(random.choices(list('acgtn'), cum_weights=cum_weights, k=seq_len))
        seqs += f'>seq{i}\n{seq}\n'
    fasta_file.write_text(seqs)
    return fasta_file


@pytest.fixture
def seq_gen():
    def fasta_seq(count_a, count_c, count_g, count_t, count_n, seed=42):
        random.seed(seed)
        a = ['a'] * count_a
        c = ['c'] * count_c
        g = ['g'] * count_g
        t = ['t'] * count_t
        n = ['n'] * count_n
        b = a + c + g + t + n
        random.shuffle(b)
        seq = "".join(b)
        return seq
    return fasta_seq


def mutate_seq(seq):
    base = random.choice(list('acgt'))
    while True:
        pos = random.choice(list(range(0, len(seq))))
        if base != seq[pos]:
            new_seq = seq[0:pos] + base + seq[(pos + 1):]
            break
    return new_seq, base, pos


@pytest.fixture
def gen_aln(request, seq_gen):
    ref = seq_gen(5, 5, 5, 5, 0)
    mut, base, pos = mutate_seq(ref)
    aln = pd.DataFrame.from_dict({'ref': list(ref), 'alt': list(mut)})
    return aln.transpose(), base, pos
