# DNApi <sub>(version: 1.0)<sub>
*de novo* adapter prediction (iterative) algorithm for small RNA sequencing data.

## Introduction
DNApi is a software package that predicts 3′ adapter sequences *de
novo* to process any small RNA libraries. The package is composed of
the following two programs:
* [`dnap`](https://github.com/jnktsj/DNApi#dnap-3-adapter-prediction)
  predicts 3′ adapter sequences from an input FASTQ.
* [`dnai`](https://github.com/jnktsj/DNApi#dnai-iterative-3-adapter-search-and-quality-control)
  predicts 3′ adapter sequences iteratively and/or performs quality
  control for an input FASTQ.

If you want to integrate the adapter prediction algorithm into your
program, see: [API](https://github.com/jnktsj/DNApi#api)

For quick examples, see:
[Examples](https://github.com/jnktsj/DNApi/tree/master/examples#examples)

For other useful utilities, see:
[Utilities](https://github.com/jnktsj/DNApi#utilities)

Of course, (sadly) there are some limitations on 3′ adapter prediction
although DNApi gives better results. For the information, see:
[Limitations](https://github.com/jnktsj/DNApi#limitations)

## Requirement
DNApi requires Python >=2.5 under a Linux/Unix environment.

## Programs
To see the usage for each program, type:

    $ <dnap | dnai> -h

or

    $ <dnap | dnai> --help

`dnap` and `dnai` accept (un)compressed FASTQ files or redirected
standard input (`stdin`) as an input.

### `dnap`: 3′ adapter prediction
`dnap` predicts 3′ adapter sequences *de novo* from an input FASTQ.
#### Usage

    $ dnap [options] <fastq>

#### Options
###### -k BP
K-mer length to use to compute k-mer frequency in the input FASTQ.
The default value is 9 nucleotides (nt).
###### -r FLOAT
Cutoff ratio for filtering less frequent k-mers. For each k-mer, a
ratio of the frequency of the most abundant k-mer to the frequency of
a target k-mer will be computed. If a ratio is lower than the cutoff
specified with `-r`, the k-mer with the ratio will be discarded.
###### -a
This option shows other predicted 3′adapter candidates (if any).

### `dnai`: iterative 3′ adapter search and quality control
`dnai` searches 3′ adapter sequences iteratively and conducts quality
control for an input FASTQ by mapping reads after adapter removal. If
a 3′adapter sequence is specified with `-3`, the program only executes
quality control using a given genome mapping command. `dnai` also maps
reads without adapter removal to investigate whether the reads are
already processed or mappable to the genome.

`dnai` judges the input FASTQ quality is poor when the mapping rate is
below 20%.
#### Usage

    $ dnai [options] <mapping_cmd> <fastq>

`<mapping_cmd>` is the genome mapping command to be tested.
For this argument, any read mapping software package can be used. 
The requirements for this argument are:
* Specify FASTA as the input read format
* Specify the input read filename as `@in`
* Specify SAM as the output format for the mapping results
* Specify the output SAM filename as `@out`
* Pass `<mapping_cmd>` as a string in the command for `dnai`

For example, when you want to use
[Bowtie](http://bowtie-bio.sourceforge.net) as a mapping engine, the
entire command line for `dnai` will be:

    $ dnai "/path_to/bowtie /path_to/genome_index -p8 -v0 -k1 -S -f @in > @out" <fastq>

[Bowtie](http://bowtie-bio.sourceforge.net) options used:
* `-p <int>`: Number of `<int>` CPUs
* `-v <int>`: Number of `<int>` mismatches
* `-k <int>`: report up to `<int>` valid alignments
* `-S`: SAM output
* `-f`: FASTA input

The results will be printed in standard output (`stdout`). The length
of predicted 3′ adapter sequences will be the 3′ adapter prefix match
length specified by `-l` + 5nt.

#### Options
##### Adapter removal
###### -l BP
3′ adapter prefix match length. `dnai` only considers perfect adapter
matches. The default is 7nt. This option affects the length of
predicted 3′ adapter sequences in the final output.
###### -m BP
Minimum read length to keep for mapping. Extracted small RNA reads
will be discarded if the lengths are *shorter* than the specified
length with `-m`. The default is 16nt.
###### -x BP
Maximum read length to keep for mapping. Extracted small RNA reads
will be discarded if the lengths are *longer* than the specified
length with `-x`. The default is 36nt.
###### --trim-5p
Trim specified number of bases from 5′ ends after adapter removal.
This option will be combined with the adapter clipping process to
trim down specific number of bases additionally.
###### --trim-3p
Trim specified number of bases from 3′ ends after adapter removal.
This option can be combined with the adapter clipping process to
trim down specific number of bases additionally.
##### Evaluation of 3′ adapter candidates
###### -3 SEQ1,SEQ2,...
Comma-separated list of 3′ adapter(s) for quality control.  When the
option is specified, `dnai` maps the processed reads after clipping
each 3′ adapter in every run and checks the genome mapping rate.
##### Iterative 3′ adapter search
###### -p FLOAT
Subsampling fraction of reads in an input FASTQ.  In the default,
`dnai` uses all reads, i.e., `-p 1.0`.  Small read sets can make
`dnai` faster.
###### -k BEG:END:INT
K-mers to predict a 3′ adapter in the input FASTQ. `BEG` is the smallest
k-mer to start, `END` is the largest k-mer to end, and `INT` is an
interval of the k-mers. The default is `9:11:2`, i.e., from 9mer to
11mer in a 2nt interval (k = 9, 11).
###### -r BEG:END:INT
Cutoff ratios for filtering less abundant k-mers. As in option `-k`,
`BEG` is the smallest ratio to start, `END` is the largest ratio to
end, and `INT` is an interval of the ratios. The default is
`1.2:1.4:0.1`, i.e., from 1.2 to 1.4 in a 0.1 interval (r = 1.2, 1.3,
1.4).
###### --temp PATH
Path for the temporary directory. `dnai` creates a temporary directory
during a computation. In the default setting, the program makes the
directory in the current directory.

## API
You can access the adapter prediction algorithm once you import
`dnapilib.apred` in your python program. `adapterPrediction` is the
core function for adapter prediciton. It takes four arguments: FASTQ
file name, filtering ratio, k-mer size, and subsampling read count. As
the result, `adapterPrediction` returns the list of tuples containing
predicted 3′ adapters and the assembly scores. The returned list is
sorted by the scores.

```python
from dnapilib.apred import adapterPrediction

# adapterPrediction(FASTQ, ratio, k-mer, subsampleReadCount)
adapts = adapterPrediction("examples/good.fq", 1.4, 9, 50000)

# all predicted adapters
print adapts
# >> [('TGGAATTCTCGGGTGCCAAGGAACTCC', 913012)]

# predicted adapter with the highest score
print adapts[0][0]
# >> 'TGGAATTCTCGGGTGCCAAGGAACTCC'
```


## Utilities
In addition to `dnap` and `dnai`, there are potentially useful three
programs in the `utils` directory:

* `qual-offset` estimates ASCII-encoded quality score offsets of FASTQ
  files.
* `qual-trim` trims low quality bases in input FASTQ reads. The
  quality trimming algorithm in the program is the same as the one in
  BWA.
* `to-fasta` removes specified 5′ and/or 3′ adapter sequences, merges
  identical reads while retainig the counts, and writes the collapsed
  reads as FASTA in standard output (`stdout`).

To see the usage for each program, type:

    $ <program-name> -h

or

    $ <program-name> --help

## Limitations
DNApi has a few limitations on 3′ adapter prediction:
* Poly(A) or other low-complexity 3′ adapters can't be predicted due
  to the low-complexity k-mer filtering step.
* Prediction accuracy will drop if gel-extracted lengths of RNAs are
  long enough to be sequenced, i.e. if few 3′ adapters are in FASTQ.
* DNApi can't do de-multiplexing.
