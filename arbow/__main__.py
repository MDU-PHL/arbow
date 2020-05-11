import logging
import subprocess
import datetime
import re
import pathlib
import warnings
import sys
import click
import pandas as pd
import numpy as np
from scipy import stats # noqa
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from .utils import IQTree
from .version import __version__ as version


# logging.captureWarnings(capture=True)
logger = logging.getLogger(__name__)
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='arbow:%(levelname)s:%(asctime)s:%(message)s', datefmt="%Y-%m-%d %H:%M:%S")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

n_seqs = 0
aln_length = 0


# noinspection PyUnusedLocal
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


def clean_seqs(fasta, outfasta, log):
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
                stat_summary: pd.Series = stats.describe(len_gaps)
                print(
                    f'[RSQ]{rec_id}\t{len_seq}\t{n_gaps}\t{stat_summary.mean:0.3}\t{stat_summary.variance:0.3}\t'
                    f'{stat_summary.minmax[0]}\t{stat_summary.minmax[1]}',
                    file=log,
                )
            elif n_gaps == 1:
                len_gaps = [len(g) for g in gaps]
                print(
                    f"[RSQ]{rec_id}\t{len_seq}\t{n_gaps}\t{len_gaps[0]}\t\t{len_gaps[0]}\t{len_gaps[0]}",
                    file=log,
                )
            else:
                print(f"[RSQ]{rec_id}\t{len_seq}\t0\t\t\t\t", file=log)
            clean_seq = re.sub(find_n5_strings, "", seq)
            ofa.write(f">{rec_id}\n{clean_seq}\n")
    logger.info("Data clean...")


def run_maftt(fasta):
    """
    TODO --- still need to be implemented
    """
    logger.info("Running maftt...")
    logger.info("Alignment ready...")


# noinspection PyPep8Naming
def fasta2df(fn, labels=(0.5, 1.0, 5.0), log=sys.stdout):
    logger.info("Loading FASTA alignment...")
    global n_seqs
    global aln_length
    recs = []
    print("[SEQ]ID\tA\tC\tG\tT\tN\tPROP_MISS\tLENGTH\tSTATUS", file=log)
    with open(fn) as fasta:
        for rec in SeqIO.parse(fasta, format="fasta"):
            seq_len = len(rec)
            if n_seqs == 0:
                aln_length = seq_len
            else:
                if seq_len != aln_length:
                    logger.error(f"Sequences do not appear to be aligned. "
                                  f"Please run an alignment program, such as MAFTT "
                                  f"before using arbow")
                    sys.exit()
            tA, tC, tG, tT = [count_base(rec, base) for base in ["a", "c", "g", "t"]]
            total_valid_bases = sum([tA, tC, tG, tT])
            tN = seq_len - total_valid_bases
            prop_missing = 100 * tN / seq_len
            lab = "*" * sum(np.array(labels) < prop_missing)
            print(
                f"[SEQ]{rec.id}\t{tA}\t{tC}\t{tG}\t{tT}\t{tN}\t{100*tN/seq_len:.03}\t{seq_len}\t{lab}",
                file=log,
            )
            recs.append(seq2series(rec))
            n_seqs += 1
    logger.info(f"Loaded {n_seqs} sequences")
    return pd.DataFrame(recs)


def trim_aln(aln_df, ref_seq="MN908947.3", five_prime_end=265, three_prime_start=29675):
    try:
        logger.info("Trimming according to ref sequence...")
        missing_in_ref = aln_df.loc[ref_seq, ].apply(lambda nuc: nuc == "n")
        logger.info(f"Found {sum(missing_in_ref)} introduced gaps into the ref...")
        df_ref = aln_df.loc[:, ~missing_in_ref]
        logger.info(f"New alignment length: {df_ref.shape[1]}...")
        logger.info("Trimming 5' and 3' UTR regions...")
        df_ref.columns = list(range(1, df_ref.shape[1] + 1))
        df_ref = df_ref.loc[:, (five_prime_end + 1): (three_prime_start - 1)]
        logger.info(f"Clean alignment length: {df_ref.shape[1]}")
    except KeyError:
        logger.critical(f"Did not find any sequence with id {ref_seq} in alignment!")
        sys.exit(1)
    except Exception as e:
        logger.critical(e)
        sys.exit(1)
    return df_ref


def scrub_seq(seq, ref, window, threshold, substitution="*"):
    """
    TODO --- implementation pending
    The idea here is to implement a rolling window to apply some filtering
    options.
    """
    new_seq = seq
    for i in range(0, len(seq)):
        q = seq[i: (i + window)]
        if len(q) == window:
            r = ref[i: (i + window)]
            total_var = 0
            total_n = 0
            var_pos = []
            for j, (x, y) in enumerate(zip(q, r)):
                if x == "-" or x == "n":
                    total_n += 1
                    continue
                if x != y:
                    total_var += 1
                    var_pos += [j]
                    continue
            if (total_n + total_var) / window > threshold:
                subs = "".join(
                    [b if k not in var_pos else substitution for k, b in enumerate(q)]
                )
                new_seq = new_seq[:i] + subs + new_seq[(i + 8):]
    return new_seq


def remove_gaps(seq, ix=None):
    if ix is None:
        ix = [i.start() for i in re.finditer('-', seq)]
    return "".join([b for i, b in enumerate(seq) if i not in ix]), ix


def summary(row, letters, stream=True, fmt_string=None, log=sys.stdout):
    bases = dict(zip(letters, [0] * len(letters)))
    for index, value in row.items():
        try:
            bases[value] += 1
        except KeyError:
            bases["n"] += 1
    bases['pos'] = row.name
    if stream and fmt_string:
        print(fmt_string.format(**bases), file=log)
    return bases


def get_per_column_summary(tab, all_iupac=True, stream=True, log=sys.stdout):
    logger.info("Getting column stats...")
    alphabet = IUPACAmbiguousDNA()
    nucs = sorted(list(alphabet.letters[0:4]))
    if all_iupac:
        nucs += sorted(list(alphabet.letters[5:-1]))
    nucs += ["N"]
    letters = [nuc.lower() for nuc in nucs]
    if stream:
        fmt_string = "[ALN]{pos}\t" + "\t".join([f"{{{letter}}}" for letter in letters])
        header = "[ALN]pos\t" + "\t".join(letters)
        print(header, file=log)
    else:
        fmt_string = None
    return pd.DataFrame(
        tab.apply(
            lambda x: summary(x, letters, stream, fmt_string, log=log), axis=0
        ).to_list()
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


def is_const_by_count(col_data, max_alt_count, allow_missing=True):
    """
    If the user supplies a maximum count for the alternate allele in
    order to consider the site a constant site, run this function
    """
    miss = "allowing" if allow_missing else "not allowing"
    logger.info(f"Filtering constant sites by {max_alt_count} max minor allele count and {miss} "
                 f"missing sites.")

    try:
        max_alt_count = int(max_alt_count)
    except ValueError:
        logger.critical("Maximum alternate count must be an integer")
    if max_alt_count < 1:
        logger.critical(f"Maximum alternate allele cannot be less than 1.")
        sys.exit(1)
    if max_alt_count > n_seqs/2:
        logger.critical(f"Maximum alternate allele cannot be 1/2 of {n_seqs}")
        sys.exit(1)
    if allow_missing:
        # if allow missing, then consider constant the sites those for which one of the
        # possible bases (a,c,t,g) the count is equal to or larger than the
        # total sequences in the alignment (number of entries in a column) minus the
        # number of observed missing elements in the column minus the maximum allowed
        # alternate alleles for example:
        # max_alt_allele = 2
        # allow_missing = True
        # total_sequences = 100
        # base count is:
        # a = 88, c = 2, g = 0, t = 0, n = 10 ==> constant
        # a = 89, c = 1, g = 0, t = 0, n = 10 ==> constant
        # a = 88, c = 1, g = 1, t = 0, n = 10 ==> constant
        # a = 90, c = 0, g = 0, t = 0, n = 10 ==> constant
        # a = 87, c = 3, g = 0, t = 0, n = 10 ==> variable
        # a = 95, c = 0, g = 0, t = 0, n = 5  ==> constant
        # a = 95, c = 0, g = 3, t = 2, n = 0  ==> variable

        const = col_data.apply(lambda x: any(x >= (n_seqs - max_alt_count - x.n)), axis=1)
    else:
        const = col_data.apply(lambda x: x.n == 0 and any(x >= (n_seqs - max_alt_count)), axis=1)
    return const


def is_const_by_freq(col_data, min_major_allele_freq, allow_missing=True):
    """
    If the user supplies a minimum allele frequency for the maximum allele in
    order to consider the site a constant site, run this function
    """
    miss = "allowing" if allow_missing else "not allowing"
    logger.info(f"Filtering constant sites by {min_major_allele_freq} minimum major allele frequency"
                 f" and {miss} missing sites.")
    try:
        min_major_allele_freq = float(min_major_allele_freq)
    except ValueError:
        logger.critical("Minimum major allele frequency must be a float")
    if min_major_allele_freq <= 0.5:
        logger.critical("Major allele frequency cannot be less than or equal to 0.5")
        sys.exit(1)
    if min_major_allele_freq > 1.0:
        logger.critical("Major allele frequency cannot be larger than 1.0")
        sys.exit(1)
    if allow_missing and min_major_allele_freq == 1.0:
        const = col_data.apply(lambda x: any(x == (n_seqs - x.n)), axis=1)
    elif allow_missing and min_major_allele_freq < 1.0:
        const = col_data.apply(
            lambda x: test_major_allele_freq(x, min_major_allele_freq), axis=1
        )
    else:
        const = col_data.apply(
            lambda x: x.n == 0 and test_major_allele_freq(x, min_major_allele_freq), axis=1
        )
    return const


def is_const(col_data, max_alt_count=None, min_major_allele_freq=None, allow_missing_data=True):
    logger.info("Finding constant sites...")
    if max_alt_count is not None:
        const = is_const_by_count(col_data, max_alt_count, allow_missing_data)
    elif min_major_allele_freq is not None:
        const = is_const_by_freq(col_data, min_major_allele_freq, allow_missing_data)
    else:
        const = col_data.apply(lambda x: any(x == n_seqs) and x.n == 0, axis=1)
    logger.info(f"Total constant sites: {sum(const)}")
    logger.info(f"Total variable sites: {sum(~const)}")
    if sum(~const) == 0:
        logger.critical("Current filtering parameters resulted in no variable sites remaining."
                         " Try relaxing your search.")
        sys.exit(1)
    return col_data[const]['pos'], col_data[~const]['pos']


def include_sites(col_data, max_missing_count=None, max_missing_proportion=None):
    """
    Identify the sites to be included in the final alignment based on the 
    number of missing data points.

    Give proportion of missing sites precedence.

    Input: a table of counts of bases per column in the alignment, including missing or 
           gappy sites
    Output: a boolean Pandas Series indicating the position of the sites to keep
    """
    if max_missing_proportion:
        # In this case, we need to calculate the number of missing sites based on
        # the total number of sequences, and then filter them out.
        logger.info(f"Filtering out sites with > {max_missing_proportion} proportion of "
                     f"missing sites.")
        n_miss = max_missing_proportion * n_seqs
        included_ix = col_data.eval(f"n<={n_miss}")
        total_included_sites = sum(included_ix)
        total_removed_sites = sum(~included_ix)
        logger.info(f"Removing {total_removed_sites} due to missing data.")
    elif max_missing_count:
        # This case is easy, we just query the data.frame for columns that have at most
        # max_missing_count `n`s
        logger.info(f"Filtering out sites with > {max_missing_count} sites.")
        included_ix = col_data.eval(f"n<={max_missing_count}")
        total_included_sites = sum(included_ix)
        total_removed_sites = sum(~included_ix)
        logger.info(f"Removing {total_removed_sites} due to missing data.")
    else:
        logger.warning("Not filtering out any sites with missing data.")
        included_ix = col_data.eval(f"n>=0")
        total_included_sites = sum(included_ix)

    logger.info(f"Total included sites: {total_included_sites}")
    return included_ix


def get_base_counts(col_stats, const_sites):
    const_sites = col_stats.set_index('pos').loc[const_sites]
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
    logger.info(f"Wrote variable sites to {outfile}")


def output_separate_trees(prefix):
    logger.info("Writing out separate trees for BB and aLRT....")
    p = re.compile("\\)(\\d+\\.?\\d?\\d?)/(\\d+):")
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
    iqtree: IQTree,
    mset="HKY,TIM2,GTR",
    mfreq="F",
    mrate="G,R",
    cmax=5,
    bb=1000,
    alrt=1000,
    threads=8,
    include_const=False,
    use_ref_as_outgroup=None,
    iqtree_outgroup=None
):
    if iqtree.version == 1:
        cmd = f"{iqtree} -s {aln_fn} -pre {pre} " \
              f"-mset {mset} -mfreq {mfreq} -mrate {mrate} -cmax {cmax} " \
              f"-bb {bb} -alrt {alrt} -nt {threads}"
    else:
        cmd = f"{iqtree} -s {aln_fn} --prefix {pre} -T {threads} " \
              f"--ufboot {bb} --alrt {alrt} " \
              f"--mset {mset} -mfreq {mfreq} -mrate {mrate} --cmax {cmax}"
    # There are two cases to consider when adding an outgroup
    # 1. The user wants to simply use the reference as the outgroup
    # 2. The user wants to specify a particular taxon (or list of taxa)
    #    as the outgroup
    # I am assuming that one would be more likely to want to use the ref
    # as the outgroup, however, I do add a check here in case the user
    # unintentionally specifies both options in the command line.
    # So, first we have the assumed predominant use case
    if use_ref_as_outgroup:
        cmd += f" -o {use_ref_as_outgroup}"
    # this is the use case we want to consider as well, but we must be
    # careful to make sure that use_ref_as_outgroup is not specified as
    # well, otherwise we are unclear what the user would want
    if iqtree_outgroup and not use_ref_as_outgroup:
        cmd += f" -o {iqtree_outgroup}"
    # if the user's intentions are unclear, make a guess and issue
    # a warning to alert the user of what we have done
    if iqtree_outgroup and use_ref_as_outgroup:
        logger.warning("You specified both --use-ref-as-outgroup and --iqtree-outgroup-id "
                       "Using the reference as the outgroup!"
                       "If that is incorrect, please only specify --iqtree-outgroup-id")
    if not include_const:
        fconst = "{a},{c},{g},{t}".format(**base_counts)
        cmd += f" -fconst {fconst}"
    logger.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)


def default_prefix(file_type, outdir: str = None):
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
    "-mc",
    "--max-missing-count",
    default=None,
    help="Sites with more than max-missing-count will be removed from the alignment. "
         "This option is ignored if -mp is used.",
    show_default=True,
    type=int,
)
@click.option(
    "-mm",
    "--max-missing-count-deprecated",
    default=None,
    help="Sites with more than max-missing-count will be removed from the alignment. "
         "This option is ignored if -mp is used. "
         "THIS OPTION IS NOW DEPRECATED IN FAVOUR OF -mc",
    show_default=True,
    type=int,
)
@click.option(
    "-mp",
    "--max-missing-proportion",
    default=None,
    help="Sites with more than max-missing-proportion of the sites missing will be "
         "removed from the alignment. This option takes precedence over -mp.",
    show_default=True,
    type=click.FloatRange(min=0, max=1),
)
@click.option(
    "-x",
    "--min-major-allele-freq",
    help="The minimum major allele frequency required to consider a site constant. Needs to be a float from 0.5 "
         "to 1.0. (e.g., -x 0.99 means that in an alignment with 100 sequences, if a site has 99 As and 1 T, the site "
         "will be considered constant). Ignored if -c is used.",
    default=None,
    show_default=True,
    type=click.FloatRange(min=0.5, max=1),
)
@click.option(
    "-c",
    "--max-alt-allele-count",
    help="The maximum count of all alternate alleles for a site to be considered constant (the alternate alleles are "
         "defined as all the bases with count less than the base with the largest count). (e.g., -c 1 means that "
         "in an alignment with 100 sequences a site with 99 As and 1 T is considered constant."
         "Needs to be an integer from 1 to 1/2 total seqs",
    default=1,
    show_default=True,
    type=int,
)
@click.option(
    "-o",
    "--out-var-aln",
    default=default_prefix("arbow-aln") + ".aln",
    help="Filename for alignment of variable sites.",
    show_default=True,
)
@click.option(
    "-p",
    "--prefix",
    default=default_prefix("arbow-tree"),
    show_default=True,
    help="Prefix to append to IQTree output files.",
)
@click.option(
    "-iv",
    "--iqtree-version",
    default="2",
    show_default=True,
    help="Version of IQTree to use.",
    type=click.Choice(["1", "2"])
)
@click.option(
    "-ip",
    "--iqtree-path",
    show_default=True,
    help="Path to iqtree executable",
    envvar="PATH",
    type=click.Path(),
    multiple=True
)
@click.option(
    "-ie",
    "--iqtree-exec",
    show_default=True,
    help="The IQTree executable",
    default="iqtree2"
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
    "--iqtree-cmax",
    default=5,
    show_default=True,
    help="Maximum number of rate categories to test.",
)
@click.option(
    "-io",
    "--iqtree-outgroup",
    default=None,
    show_default=True,
    help="ID(s) of samples to be used as outgroup in IQTree. (e.g., single: sample_1, "
         "multiple: 'sample_2,sample_3'. This option is ignore if --use-ref-as-outgroup is "
         "selected.",
)
@click.option(
    "-r",
    "--ref-id",
    default="MN908947.3",
    help="Sequence ID of the reference",
    show_default=True,
)
@click.option(
    "-u",
    "--use-ref-as-outgroup",
    is_flag=True,
    help="Use the reference sequence as the outgroup in IQTree",
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
@click.option(
    "--include-const",
    is_flag=True,
    help="When outputting the clean alignment, leave constant sites in the alignment. [default is to remove]",
)
@click.option(
    "-l",
    "--log",
    type=click.File("w"),
    default=default_prefix("arbow") + ".log",
    help="Log file to store output. Use '-' to log to stdout",
    show_default=True,
)
def main(
    fasta,
    all_iupac,
    no_stream,
    max_missing_count,
    max_missing_count_deprecated,
    max_missing_proportion,
    min_major_allele_freq,
    max_alt_allele_count,
    out_var_aln,
    prefix,
    iqtree_version,
    iqtree_path,
    iqtree_exec,
    iqtree_threads,
    iqtree_models,
    iqtree_freq,
    iqtree_rates,
    iqtree_cmax,
    iqtree_bb,
    iqtree_alrt,
    iqtree_outgroup,
    ref_id,
    use_ref_as_outgroup,
    five_prime_end,
    three_prime_start,
    include_const,
    log,
):
    iqtree = IQTree(iqtree_path, iqtree_exec, iqtree_version)
    iqtree.check_version()
    # outfa = default_prefix("arbow-clean-seqs") + ".fa"
    # fasta = clean_seqs(fasta, outfa)
    if max_missing_count_deprecated is not None:
        warnings.warn("-mm has been deprecated in favour of -mc. "
                      "-mm will be removed in v0.8.0", FutureWarning)
        max_missing_count = max_missing_count_deprecated
    aln = fasta2df(fasta, log=log)
    aln = trim_aln(
        aln,
        ref_seq=ref_id,
        five_prime_end=five_prime_end,
        three_prime_start=three_prime_start,
    )
    col_stats = get_per_column_summary(aln, all_iupac, not no_stream, log=log)
    included_sites_ix = include_sites(col_stats,
                                      max_missing_count=max_missing_count,
                                      max_missing_proportion=max_missing_proportion)
    included_col_stats = col_stats[included_sites_ix]
    if max_missing_count is None and max_missing_proportion is None:
        allow_missing = False
    elif max_missing_count is not None and max_missing_count > 0:
        allow_missing = True
    elif max_missing_proportion is not None and max_missing_proportion > 0:
        allow_missing = True
    else:
        allow_missing = False
    const_sites_pos, var_sites_pos = is_const(included_col_stats,
                                              max_alt_count=max_alt_allele_count,
                                              min_major_allele_freq=min_major_allele_freq,
                                              allow_missing_data=allow_missing)
    output_aln(aln, var_sites_pos, outfile=out_var_aln, filter_const=not include_const)
    base_counts = get_base_counts(included_col_stats, const_sites_pos)
    run_iqtree(
        out_var_aln,
        base_counts,
        prefix,
        iqtree=iqtree,
        mset=iqtree_models,
        mfreq=iqtree_freq,
        mrate=iqtree_rates,
        bb=iqtree_bb,
        alrt=iqtree_alrt,
        threads=iqtree_threads,
        cmax=iqtree_cmax,
        include_const=include_const,
        use_ref_as_outgroup=ref_id if use_ref_as_outgroup else None,
        iqtree_outgroup=iqtree_outgroup
    )
    output_separate_trees(prefix)
    logger.info(
        'Support values on tree are SH-aLRT (good support is >= 80%) and UFboot (good support is >= 95%)'
    )


if __name__ == "__main__":
    main()
