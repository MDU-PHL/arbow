import logging
import subprocess
import datetime
import multiprocessing
import re
import pathlib
import click
import pandas as pd
import numpy as np
from scipy import stats
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
    seq_str = str(seq.seq).lower().replace("-", "n").replace(" ", "n")
    return pd.Series(list(seq_str), name=seq.id)


def clean_seqs(fasta, outfasta):
    logger.info("Cleaning sequences...")
    invalid_bases = re.compile("[^actg]", re.IGNORECASE)
    find_n_strings = re.compile("n+", re.IGNORECASE)
    find_n5_strings = re.compile("n{5,}", re.IGNORECASE)
    with open(fasta) as fa, open(outfasta, "w") as ofa:
        for rec in SeqIO.parse(fa, "fasta"):
            rec_id = rec.id
            seq = str(rec.seq)
            len_seq = len(seq)
            # make sure we only have `actg` and `n`
            seq = re.sub(invalid_bases, "n", seq)
            # find all the gaps
            gaps = find_n_strings.findall(seq)
            n_gaps = len(gaps)
            # need to take care of the three cases:
            #  1. if we have more than 1 gap
            #  2. if we have only one gap
            #  3. if we have no gaps
            if n_gaps > 1:
                len_gaps = [len(g) for g in gaps]
                stat_summary = stats.describe(len_gaps)
                print(
                    f"[RSQ]{rec_id}\t{len_seq}\t{n_gaps}\t{stat_summary.mean:0.3}\t{stat_summary.variance:0.3}\t{stat_summary.minmax[0]}\t{stat_summary.minmax[1]}"
                )
            elif n_gaps == 1:
                len_gaps = [len(g) for g in gaps]
                print(
                    f"[RSQ]{rec_id}\t{len_seq}\t{n_gaps}\t{len_gaps[0]}\t\t{len_gaps[0]}\t{len_gaps[0]}"
                )
            else:
                print(f"[RSQ]{rec_id}\t{len_seq}\t0\t\t\t\t")
            clean_seq = re.sub(find_n5_strings, "", seq)
            ofa.write(f">{rec_id}\n{clean_seq}\n")
    logger.info("Data clean...")


def run_maftt(fasta):
    logger.info("Running maftt...")
    logger.info("Alignment ready...")


def fasta2df(fn, labels=[0.5, 1.0, 5.0]):
    logger.info("Loading FASTA alignment...")
    global n_seqs
    global aln_length
    recs = []
    print("[SEQ]ID\tA\tC\tG\tT\tN\tPROP_MISS\tLENGTH\tSTATUS")
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
                f"[SEQ]{rec.id}\t{tA}\t{tC}\t{tG}\t{tT}\t{tN}\t{100*tN/seq_len:.03}\t{seq_len}\t{lab}"
            )
            recs.append(seq2series(rec))
            n_seqs += 1
    logger.info(f"Loaded {n_seqs} sequences")
    return pd.DataFrame(recs)


def trim_aln(aln_df, ref_seq="MN908947.3", five_prime_end=265, three_prime_start=29675):
    logger.info("Trimming according to ref sequence...")
    missing_in_ref = aln_df.loc[ref_seq, ].apply(lambda nuc: nuc == "n")
    logger.info(f"Found {sum(missing_in_ref)} introduced gaps into the ref...")
    df_ref = aln_df.loc[:, ~missing_in_ref]
    logger.info(f"New alignment length: {df_ref.shape[1]}...")
    logger.info("Trimming 5' and 3' UTR regions...")
    df_ref.columns = list(range(1, df_ref.shape[1] + 1))
    df_ref = df_ref.loc[:, (five_prime_end + 1): (three_prime_start - 1)]
    logger.info(f"Clean alignment length: {df_ref.shape[1]}")
    return df_ref


def scrub_seq(seq, ref, window, threshold, substitution="*"):
    new_seq = seq
    for i in range(0, len(seq)):
        q = seq[i:(i + window)]
        if len(q) == window:
            r = ref[i:(i + window)]
            total_var = 0
            total_n = 0
            var_pos = []
            for j, (x, y) in enumerate(zip(q, r)):
                if x == '-' or x == 'n':
                    total_n += 1
                    continue
                if x != y:
                    total_var += 1
                    var_pos += [j]
                    continue
            if (total_n + total_var) / window > threshold:
                subs = "".join([b if k not in var_pos else substitution for k, b in enumerate(q)])
                new_seq = new_seq[:i] + subs + new_seq[(i + 8):]
    return new_seq


def remove_gaps(seq, ix=None):
    if ix is None:
        ix = [i.start() for i in re.finditer("-", seq)]
    return "".join([b for i, b in enumerate(seq) if i not in ix]), ix


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
        fmt_string = "[ALN]{pos}\t" + "\t".join([f"{{{l}}}" for l in letters])
        header = "[ALN]pos\t" + "\t".join(letters)
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
    """
    Identify the sites to be included in the final alignment based on the 
    number of missing data points. 

    Input: a table of counts of bases per column in the alignment, including missing or 
           gappy sites
    Output: a boolean Pandas Series indicating the position of the sites to keep
    """
    included_ix = col_data.eval(f"n<={max_missing}")
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


def output_aln(aln, pos, outfile="var.aln", filter_const=True):
    if filter_const:
        aln = aln.iloc[:, pos]
    fasta_str = "\n".join(
        aln.apply(lambda x: f">{x.name}\n{''.join(x.to_list())}", axis=1)
    )
    with open(outfile, "w") as fasta:
        fasta.write(fasta_str)
    logger.info(f"Wrote variable sitest to {outfile}")


def output_separate_trees(prefix):
    logger.info("Writting out separate trees for BB and aLRT....")
    p = re.compile("\\)(\\d+\\.?\\d?\\d?)\\/(\\d+):")
    treefile = prefix + ".treefile"
    bbtree = prefix + "_bb.treefile"
    alrttree = prefix + "_alrt.treefile"
    with open(treefile) as tr, open(bbtree, "w") as bb, open(alrttree, "w") as alrt:
        tr = tr.read().strip()
        bb.write(re.sub(p, ")\\2:", tr) + "\n")
        alrt.write(re.sub(p, ")\\1:", tr) + "\n")


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
    include_const=False
):
    cmd = f"iqtree -s {aln_fn} -pre {pre} -mset {mset} -mfreq {mfreq} -mrate {mrate} -bb {bb} -alrt {alrt} -nt {threads}"
    if not include_const:
        fconst = "{a},{c},{g},{t}".format(**base_counts)
        cmd += f" -fconst {fconst}"
    logging.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)


def default_prefix(file_type, outdir=None):
    if outdir is not None:
        file_type = pathlib.Path(outdir).joinpath(file_type)
    return f"{file_type}-" + datetime.datetime.strftime(
        datetime.datetime.now(), "%Y-%m-%d-%H%M%S"
    )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("fasta")
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
    default=0.999,
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
@click.option(
    "-r",
    "--ref-id",
    default="MN908947.3",
    help="Sequence ID of the reference",
    show_default=True,
)
@click.option(
    "--five-prime-end",
    default=265,
    help="Last base of the 5' UTR region in 1-index in the ref sequence",
    show_default=True,
)
@click.option(
    "--three-prime-start",
    default=29675,
    help="First base of the 3' UTR region in 1-index in the ref sequence",
    show_default=True,
)
@click.option("--include-const", is_flag=True, help="When outputting the clean alignment, leave constant sites in the alignment. [default is to remove]")
def main(
    fasta,
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
    ref_id,
    five_prime_end,
    three_prime_start,
    include_const
):
    # outfa = default_prefix("arbow-clean-seqs") + ".fa"
    # fasta = clean_seqs(fasta, outfa)
    aln = fasta2df(fasta)
    aln = trim_aln(
        aln,
        ref_seq=ref_id,
        five_prime_end=five_prime_end,
        three_prime_start=three_prime_start,
    )
    col_stats = get_per_column_summary(aln, all_iupac, not no_stream)
    included_sites_ix = include_sites(col_stats, max_missing=max_missing)
    included_col_stats = col_stats[included_sites_ix]
    const_sites_ix = is_const(included_col_stats, major_allele_freq)
    var_sites_ix = const_sites_ix.index[~const_sites_ix].to_list()
    aln = aln.loc[:, included_sites_ix.to_list()]
    output_aln(aln, var_sites_ix, outfile=out_var_aln, filter_const=not include_const)
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
        include_const=include_const
    )
    output_separate_trees(prefix)
    logger.info(
        "Support values on tree are SH-aLRT (good support is >= 80%) and UFboot (good support is >= 95%)"
    )


if __name__ == "__main__":
    main()
