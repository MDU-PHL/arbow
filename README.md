# Arbow: cultivate your multiple sequence aligment to get better trees

## Name

We called `arbow` as that would be the phonetic pronounciation of the short, endearing, 
terms for an [arborist](https://en.wikipedia.org/wiki/Arborist) in Australia.

## What it does

The goal of `arbow` is to automate and simplify the production of trees from multiple sequence alignments.

In this version it does:

1. Reads all the alignment
2. Calculates stats per column in the alignment
3. Allows the user to set a threshold of tolerable missing data in a column, and `arbow` removes all non-conforming columns from the alignment
4. From the remaining columns, `arbow` finds all the `constant` columns, and calculates the frequency of `A`, `C`, `G`, and `T`.
5. It then filters out all the `variable` columns, and outputs that as a `multiFASTA` alignment.
6. It runs `IQTree` with a few sensible `presets`

Currently, in step `4` above, columns that have a single `base` but still have missing data that did were not filtered out in step `3` are counted as being the `base`. In other words, if a `user` specifies a maximum number of 20 missing bases, and a column with 5 missing bases but with `A` in all samples, it will count towards `A` (i.e., majority consensus imputation).

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

## Data stream

When running `arbow`, by default a stream is output to the console (`stdout`). 

Data about the each sequence in the alignment is prefixed with `[SEQ]`, and is followed by:

1. Count of each base (`A`, `C`, `G`, `T`, and `N` – `N` is any character other than `ACGT`)
2. Percent missing data
3. A status column that has 0, 1, 2, or 3 `*` depending on whether the percent missing data is `<0.5`, `>=0.5 and <1.0`, `>=1.0 and <5.0`, or `>=5`.

Data about each column in the alignment is prefixed with `[ALN]`, and is followed by:

1. Position in the alignment
2. Count of each base (bases counted will depend on whether all IUPAC codes are allowed or not - see below in usage)



## Command line

## Usage

```
Usage: arbow [OPTIONS] ALN

Options:
  --version
  -a, --all-iupac               Print count of all IUPAC code for column
                                stats?

  -s, --no-stream               Stop streaming stats to console
  -mm, --max-missing INTEGER    Remove sites with 'mm' missing sites or more
                                [default: 20]

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

  --help                        Show this message and exit.
```

### Get help

```
arbow --help
```

### Get version
```
arbow --version
```





