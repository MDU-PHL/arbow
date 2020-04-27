# Arbow: cultivate your multiple sequence aligment to get better trees

## Name

We named this tool `arbow` as that would be the phonetic pronounciation of the short, endearing, 
term for an [arborist](https://en.wikipedia.org/wiki/Arborist) in Australia.

## What it does

The goal of `arbow` is to automate and simplify the production of trees from multiple sequence alignments. The tool 
has been developed in the context of viral phylogenomics.

In the current version (`0.5.*`) it:

1. Reads an alignment in `multiFASTA` format
2. Calculates stats for each sequence in the alignment
3. Trims 5/3 prime UTR regions --- defaults set to SARS-COV-2 (Genbank accession: `NC_045512.2`)
4. Calculates stats per column in the alignment
5. Allows the user to set a threshold of tolerable missing data in a column, and removes all non-conforming columns from the alignment
6. From the remaining columns, `arbow` finds all the `constant` columns according to two `user` defined criteria: `allow missing data` (i.e., a column with missing data can still count to towards `constant` sites if it meets other criteria), and the frequency of the major allow is equal to or larger than a trheshold (i.e., if the threshold is set to 0.99 and there are 100 samples, 99 of which are `A` and one is `G`, that column would be counted as a constant `A`). Filtering by frequency allows one to remove potential sequencing error.
7. It then filters out all the `variable` columns, and outputs the variable alignment as a `multiFASTA` alignment.
8. It runs `IQTree` with a few sensible `presets`

Currently, in step `4` above, columns that have a single observed `nucleotide` (e.g., `C`) but still have missing data that were not filtered out in step `3` are counted towards the overall frequency of that `base` in the alignment. In other words, if a `user` specifies a maximum number of 20 missing bases, and a column with 5 missing bases but with `A` in all other samples, that column will count towards the overall frequency of `A` in the alignment (i.e., majority consensus imputation). This assumptions is less risky the larger the number of samples in the alignment.

For step `5`, missing data (i.e., `-` and `N`) are all codes as `N`.

Tests are underway to figure out how these assumptions might affect the output.

## Dependencies

1. Python >=3.6
2. IQTree 1.6+ (not tested on IQTree 2 as it is not production ready yet)
3. BioPython
4. Pandas
5. NumPy

## Installation

### Brew

```
brew install iqtree
pip<3> install arbow
```

### Conda

```
conda install -c bioconda iqtree
pip<3> install arbow
```

## Running

1. Generate a mulitple sequence alignment with your favourite aligner (e.g., MAFTT). Output a `multiFASTA` file.
2. Run `arbow <aln.fa>`
3. Open `tree-YYYY-MM-DD_HHMMSS.treefile` in your favourite tree viewer (e.g, FigTree)
4. Open `tree-YYYY-MM-DD_HHMMSS_bb.treefile` or `tree-YYYY-MM-DD_HHMMSS_alrt.treefile` for branches with `ultra-fast bootstrap` support or `SH-aLRT` support only, respectively.

## Data stream

When running `arbow`, by default a stream is output to the console (`stdout`). 

Data about the each sequence in the alignment is prefixed with `[SEQ]`, and is followed by:

1. Count of each base (`A`, `C`, `G`, `T`, and `N` – `N` is any character other than `ACGT`)
2. Percent missing data
3. A status column that has 0, 1, 2, or 3 `*` depending on whether the percent missing data is `<0.5`, `>=0.5 and <1.0`, `>=1.0 and <5.0`, or `>=5`, respectively.

Data about each column in the alignment is prefixed with `[ALN]`, and is followed by:

1. Position in the alignment
2. Count of each base (bases counted will depend on whether all IUPAC codes are allowed or not - see below in usage)


## Command line

## Usage

```
Usage: arbow [OPTIONS] ALN

Options:
  --version
  -i, --all-iupac               Print count of all IUPAC code for column
                                stats?

  -s, --no-stream               Stop streaming stats to console
  -mm, --max-missing INTEGER    Remove sites with 'mm' missing sites or more
                                [default: 20]

  -x, --major-allele-freq FLOAT  If major allele frequency is equal or larger
                                 than consider the site constant.  [default:
                                 0.99]

  -o, --out-var-aln TEXT        Filename for alignment of variable sites.
                                [default: aln-2020-04-07-150443.aln]

  -p, --prefix TEXT             Prefix to append to IQTree output files.
                                [default: tree-2020-04-07-150443]

  -t, --iqtree-threads INTEGER  Number of cores to run IQtree  [default: 4]
  -m, --iqtree-models TEXT      Substitution models to test.  [default:
                                HKY,TIM2,GTR]

  -f, --iqtree-freq TEXT        Base frequency models to test.  [default: F]
  -r, --iqtree-rates TEXT       Rate category models to test.  [default: G,R]
  -b, --iqtree-bb INTEGER       Maximum number of UltraFast Bootstrap
                                iterations to attempt.  [default: 1000]

  -a, --iqtree-alrt INTEGER     Number of replicates to perform SH-aLRT.
                                [default: 1000]

  -c, --iqtree-cmax INTEGER     Maximum number of rate categories to test.
                                [default: 5]

  -r, --ref-id TEXT              Sequence ID of the reference  [default:
                                 MN908947.3]

  --five-prime-end INTEGER       Last base of the 5' UTR region in 1-index in
                                 the ref sequence  [default: 265]

  --three-prime-start INTEGER    First base of the 3' UTR region in 1-index in
                                 the ref sequence  [default: 29675]

  --include-const                When outputting the clean alignment, leave
                                 constant sites in the alignment. [default is
                                 to remove]

  --help                        Show this message and exit.
```

## Default behaviour explained

By default, `arbow` will remove any site in the alignment that has `20` missing data points or more, will trim the 5' and 3' UTR regions, and will consider as constant any site that has a major allele frequency larger or equal to 0.997.

### Remove sites with any gaps in the alignment 

Let us say that you wish to remove all sites in the alignment that have **any** missing data, and retain all complete columns:

```
arbow -x 1.0 -mm 0 <in.aln>
```

### Keep all sites in an alignment (i.e., skip any filtering)

Let us say that you wish to keep all sites in the alignment, and you have an alignment with 200 sequences:

```
arbow -x 1.0 -m 200 <in.aln>
```

### Keep constant sites in the clean alignment

```
arbow --include-const <in.aln>
```

### Get help

```
arbow <-h|--help>
```

### Get version
```
arbow --version
```





