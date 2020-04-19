import logging
import subprocess
import datetime
import multiprocessing
import click
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from .version import __version__ as version

TOTAL_CORES = multiprocessing.cpu_count()
MAX_CORES = 16 if TOTAL_CORES > 16 else TOTAL_CORES

logging.basicConfig(level=logging.INFO)
logger = logging

n_seqs = 0
aln_length = 0


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"arbow {version}")
    ctx.exit()


def count_base(seq, base):
    return str(seq.seq).lower().count(base)


def seq2series(seq):
    seq_str = str(seq.seq).replace("-", "n").replace(" ", "n")
    return pd.Series(list(seq_str), name=seq.id)


def fasta2df(fn, labels=[0.5, 1.0, 5.0]):
    logger.info("Loading FASTA alignment...")
    global n_seqs
    global aln_length
    recs = []
    print("[SEQ]\tID\tA\tC\tG\tT\tN\tPROP_MISS\tLENGTH\tSTATUS")
    with open(fn) as fasta:
        for rec in SeqIO.parse(fasta, format="fasta"):
            seq_len = len(rec)
            if n_seqs == 0:
                aln_length = seq_len
            tA, tC, tG, tT = [count_base(rec, base) for base in ["a", "c", "g", "t"]]
            total_valid_bases = sum([tA, tC, tG, tT])
            tN = seq_len - total_valid_bases
            prop_missing = 100 * tN / seq_len
            lab = "*" * sum(np.array(labels) < prop_missing)
            print(
                f"[SEQ]\t{rec.id}\t{tA}\t{tC}\t{tG}\t{tT}\t{tN}\t{100*tN/seq_len:.03}\t{seq_len}\t{lab}"
            )
            recs.append(seq2series(rec))
            n_seqs += 1
    logger.info(f"Loaded {n_seqs} sequences")
    return pd.DataFrame(recs)


def summary(row, letters, stream=True, fmt_string=None):
    bases = dict(zip(letters, [0] * len(letters)))
    for index, value in row.items():
        try:
            bases[value] += 1
        except Exception:
            bases["n"] += 1
    if stream and fmt_string:
        print(fmt_string.format(pos=row.name, **bases))
    return bases


def get_per_column_summary(tab, all_iupac=True, stream=True):
    logger.info("Getting column stats...")
    alphabet = IUPACAmbiguousDNA()
    nucs = sorted(list(alphabet.letters[0:4]))
    if all_iupac:
        nucs += sorted(list(alphabet.letters[5:-1]))
    nucs += ["N"]
    letters = [l.lower() for l in nucs]
    if stream:
        fmt_string = "[ALN]\t{pos}\t" + "\t".join([f"{{{l}}}" for l in letters])
        header = "[ALN]\tpos\t" + "\t".join(letters)
        print(header)
    else:
        fmt_string = None
    return pd.DataFrame(
        tab.apply(lambda x: summary(x, letters, stream, fmt_string), axis=0).to_list()
    )


# def apply_hard_filter(aln, col_data, max_missing):
#     logger.info("Applying hard filter...")
#     ix = col_data.query(f"n<{max_missing}").index.to_list()
#     aln = aln.loc[:,ix]
#     logger.info(f"Total sites removed {len(ix)}...")
#     logger.info(f"Total sites remaining {aln.shape[1]}...")
#     return aln


def test_major_allele_freq(aln_col, major_freq):
    mask = aln_col.index.isin(["n"])
    bases = aln_col[~mask]
    total = sum(bases)
    freq = bases / total
    return freq.max() >= major_freq


def is_const(col_data, major_allele_freq, allow_missing_data=True):
    logger.info("Finding constant sites...")
    if allow_missing_data and major_allele_freq == 1:
        const = col_data.apply(lambda x: any(x == (n_seqs - x.n)), axis=1)
    elif allow_missing_data and major_allele_freq < 1.0:
        const = col_data.apply(
            lambda x: test_major_allele_freq(x, major_allele_freq), axis=1
        )
    elif major_allele_freq < 1.0:
        const = col_data.apply(
            lambda x: x.n == 0 and test_major_allele_freq(x, major_allele_freq), axis=1
        )
    else:
        const = col_data.apply(lambda x: any(x == (n_seqs)) and x.n == 0, axis=1)
    logger.info(f"Total constant sites: {sum(const)}")
    logger.info(f"Total variable sites: {sum(~const)}")
    return const


def include_sites(col_data, max_missing=20):
    included_ix = col_data.eval(f"n<{max_missing}")
    logger.info(f"Total included sites: {sum(included_ix)}")
    logger.info(f"Total removed sites: {sum(~included_ix)}")
    return included_ix


def get_base_counts(col_stats, const_sites):
    const_sites = col_stats[const_sites]
    base_counts = (
        const_sites.apply(lambda x: x.idxmax(), axis=1).value_counts().to_dict()
    )
    logger.info(
        "Constant site counts - A:{a} | C:{c} | G:{g} | T:{g}".format(**base_counts)
    )
    return base_counts


def output_variable_aln(aln, pos, outfile="var.aln"):
    aln = aln.loc[:, pos]
    fasta_str = "\n".join(
        aln.apply(lambda x: f">{x.name}\n{''.join(x.to_list())}", axis=1)
    )
    with open(outfile, "w") as fasta:
        fasta.write(fasta_str)
    logger.info(f"Wrote variable sitest to {outfile}")


def run_iqtree(
    aln_fn,
    base_counts,
    pre,
    mset="HKY,TIM2,GTR",
    mfreq="F",
    mrate="G,R",
    cmax=5,
    bb=1000,
    alrt=1000,
    threads=8,
):
    fconst = "{a},{c},{g},{t}".format(**base_counts)
    cmd = f"iqtree -s {aln_fn} -pre {pre} -mset {mset} -mfreq {mfreq} -mrate {mrate} -bb {bb} -alrt {alrt} -fconst {fconst} -nt {threads}"
    logging.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, capture_output=True)


def default_prefix(file_type):
    return f"{file_type}-" + datetime.datetime.strftime(
        datetime.datetime.now(), "%Y-%m-%d-%H%M%S"
    )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("aln")
@click.option(
    "--version", is_flag=True, callback=print_version, expose_value=False, is_eager=True
)
@click.option(
    "-i",
    "--all-iupac",
    is_flag=True,
    help="Print count of all IUPAC code for column stats?",
)
@click.option("-s", "--no-stream", is_flag=True, help="Stop streaming stats to console")
@click.option(
    "-mm",
    "--max-missing",
    default=20,
    help="Remove sites with 'mm' missing sites or more",
    show_default=True,
)
@click.option(
    "-x",
    "--major-allele-freq",
    help="If major allele frequency is equal or larger than consider the site constant.",
    default=0.99,
    show_default=True,
)
@click.option(
    "-o",
    "--out-var-aln",
    default=default_prefix("aln") + ".aln",
    help="Filename for alignment of variable sites.",
    show_default=True,
)
@click.option(
    "-p",
    "--prefix",
    default=default_prefix("tree"),
    show_default=True,
    help="Prefix to append to IQTree output files.",
)
@click.option(
    "-t",
    "--iqtree-threads",
    default=1,
    show_default=True,
    help="Number of cores to run IQtree",
)
@click.option(
    "-m",
    "--iqtree-models",
    default="HKY,TIM2,GTR",
    show_default=True,
    help="Substitution models to test.",
)
@click.option(
    "-f",
    "--iqtree-freq",
    default="F",
    show_default=True,
    help="Base frequency models to test.",
)
@click.option(
    "-r",
    "--iqtree-rates",
    default="G,R",
    show_default=True,
    help="Rate category models to test.",
)
@click.option(
    "-b",
    "--iqtree-bb",
    default=1000,
    show_default=True,
    help="Maximum number of UltraFast Bootstrap iterations to attempt.",
)
@click.option(
    "-a",
    "--iqtree-alrt",
    default=1000,
    show_default=True,
    help="Number of replicates to perform SH-aLRT.",
)
@click.option(
    "-c",
    "--iqtree-cmax",
    default=5,
    show_default=True,
    help="Maximum number of rate categories to test.",
)
def main(
    aln,
    all_iupac,
    no_stream,
    max_missing,
    major_allele_freq,
    out_var_aln,
    prefix,
    iqtree_threads,
    iqtree_models,
    iqtree_freq,
    iqtree_rates,
    iqtree_cmax,
    iqtree_bb,
    iqtree_alrt,
):
    aln = fasta2df(aln)
    col_stats = get_per_column_summary(aln, all_iupac, not no_stream)
    included_sites_ix = include_sites(col_stats, max_missing=max_missing)
    included_col_stats = col_stats[included_sites_ix]
    const_sites_ix = is_const(included_col_stats, major_allele_freq)
    var_sites_ix = const_sites_ix.index[~const_sites_ix].to_list()
    output_variable_aln(aln, var_sites_ix, outfile=out_var_aln)
    base_counts = get_base_counts(included_col_stats, const_sites_ix)
    run_iqtree(
        out_var_aln,
        base_counts,
        prefix,
        mset=iqtree_models,
        mfreq=iqtree_freq,
        mrate=iqtree_rates,
        bb=iqtree_bb,
        alrt=iqtree_alrt,
        threads=iqtree_threads,
        cmax=iqtree_cmax,
    )
    logger.info(
        "Support values on tree are SH-aLRT (good support is >= 80%) and UFboot (good support is >= 95%)"
    )


if __name__ == "__main__":
    main()
