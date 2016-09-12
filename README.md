# DNApi <sub>(version: 1.1)<sub>
*de novo* adapter prediction (iterative) algorithm for small RNA
sequencing data. DNApi requires Python (2 or 3) under a Linux/Unix
environment.

DNApi accept (un)compressed FASTQ files or redirected standard input
(`stdin`) as an input. You can simply run:

    $ python dnapi.py <fastq>

or

    $ <process-generates-fastq> | python dnapi.py -

To see the usage for each program, type:

    $ python dnapi.py [-h | --help]

For quick examples, see:
[Examples](https://github.com/jnktsj/DNApi/tree/master/examples#examples)

DNApi can predict most 3′ adapters correctly with the default
parameters. However, if you want to tweak the parameters or want to
run other prediction modes, see [prediction modes and parameters]
(https://github.com/jnktsj/DNApi#prediction-modes-and-parameters) for
mode defail.


## Miscellaneous
For other useful utilities, see:
[Utilities](https://github.com/jnktsj/DNApi#utilities)

If you want to integrate the adapter prediction algorithm into your
program, see: [API](https://github.com/jnktsj/DNApi#api)

Of course, (sadly) there are some limitations on 3′ adapter prediction
although DNApi gives near-perfect results. For the information, see:
[Limitations](https://github.com/jnktsj/DNApi#limitations)


## Prediction modes and parameters
The package covers three ways (hereafter modes) to predict adapters.
The prediction algorithm needs two main parameters
[`-k`](https://github.com/jnktsj/DNApi#-k-kmer_begkmer_endincrement--kmer_len)
(k-mer lengths) and
[`-r`](https://github.com/jnktsj/DNApi#-r-ratio_begratio_endintcrement--ratio)
(filtering ratio for less abundant kmers). The default is *iterative*
mode with `-k 9:11:2` and `-r 1.2:1.4:0.1`.  The default setting
already works well on any small RNA libraries, but you can tweak the
parametes with `-k` and `-r` (For more detail, see
[Options](https://github.com/jnktsj/DNApi#options)).

##### *Iterative* mode
*Iterative* mode runs the algorithm multiple times with different
combinations of *k* and *R* and refines the ranks of predicted adapter
candidates in subsequent iterations.

    $ python3 dnapi.py -k 9:11:2 -r 1.2:1.4:0.1 <fastq>

##### *Single* mode
*Single* mode runs a single adapter prediction algorithm with a
specific combination of *k* and *R*.

    $ python3 dnapi.py -k 9 -r 1.4 <fastq>

##### *Exaustive* mode
*Exhaustive* mode exhaustively searches an optimal 3´ adapter by
running the algorithm multiple times with different combinations of
*k* and *R* to obtain a non-redundant list of adapter candidates, and
incorporating adapter removal and read mapping. Only this mode can
judge whether input libraries are already clean (i.e. the 3′ adapter
sequences are already removed). To turn on this mode, you need to run
with
[`--map-command`](https://github.com/jnktsj/DNApi#--map-command-command)
(For more detail, see
[Options](https://github.com/jnktsj/DNApi#options)).

    $ python3 dnapi.py --map-command "<command>" <fastq>

You can also incorporate `-k` and `-r`. The default setting is
`-k 9:11:2` and `-r 1.2:1.4:0.1`.

This mode also outputs cleansed reads in FASTA format if the input
FASTQ is not processed. The reads in the output FASTA are
non-redundant, and the read counts are written in FASTA headers.

If a 3′adapter sequence is specified with `--adapter-seq`, DNApi
only executes quality control using a given genome mapping
command.

    $ python3 dnapi.py --map-command "<command>" --adapter-seq SEQ1 [SEQ2 SEQ3...] <fastq>

DNApi judges the input FASTQ quality is poor when the mapping rate is
below 20%.


#### Options

##### Adapter prediction parameters

###### -k [KMER_BEG:KMER_END:INCREMENT | KMER_LEN]
K-mer(s) to predict a 3′ adapter in the input FASTQ. When you specify
the longer argument with ":", DNApi performs *iterative* mode.  In the
longer argument, `KMER_BEG` is the smallest k-mer to start, `KNER_END`
is the largest k-mer to end, and `INCREMENT` is an interval of the
k-mers. The default is `9:11:2`, i.e., from 9mer to 11mer in a 2nt
interval (k = 9, 11). When you specify a single k-mer length
`KMER_LEN` with a single ratio (see `-r` below), DNApi runs *single*
mode.

###### -r [RATIO_BEG:RATIO_END:INTCREMENT | RATIO]
Cutoff ratio(s) for filtering less frequent k-mers. For each k-mer, a
ratio of the frequency of the most abundant k-mer to the frequency of
a target k-mer will be computed. If a ratio is lower than the cutoff
specified with `-r`, the k-mer with the ratio will be discarded. As in
option `-k` for *iterative* search, `RATIO_BEG` is the smallest ratio
to start, `RATIO_END` is the largest ratio to end, and `INCREMENT` is
an interval of the ratios. The default is `1.2:1.4:0.1`, i.e., from
1.2 to 1.4 in a 0.1 interval (r = 1.2, 1.3, 1.4). When you specify a
single ratio `RATIO` with a single k-mer length, DNApi runs *single*
mode.

###### --show-all
This option shows other predicted 3′adapter candidates (if any).

##### Exhaustive adapter search with mapping process

###### --map-command COMMAND
`COMMAND` is the genome mapping command to be tested.
For this argument, any read mapping software package can be used. 
The requirements for this argument are:
* Specify FASTA as the input read format
* Specify the input read filename as `@in`
* Specify SAM as the output format for the mapping results
* Specify the output SAM filename as `@out`
* Pass `COMMAND` as a string in the command for DNApi

For example, when you want to use
[Bowtie](http://bowtie-bio.sourceforge.net) as a mapping engine, the
entire command line for DNApi will be:

    $ python3 dnapi.py "/path_to/bowtie /path_to/genome_index -p8 -v0 -k1 -S -f @in > @out" <FASTQ>

[Bowtie](http://bowtie-bio.sourceforge.net) options used:
* `-p <int>`: Number of `<int>` CPUs
* `-v <int>`: Number of `<int>` mismatches
* `-k <int>`: report up to `<int>` valid alignments
* `-S`: SAM output
* `-f`: FASTA input

The results will be printed in standard output (`stdout`). The length
of predicted 3′ adapter sequences will be the 3′ adapter prefix match
length specified by `--prefix-match` + 5nt.

###### --subsample-rate FLOAT
Subsampling fraction of reads in an input FASTQ for *exhaustive* mode.
In the default, DNApi uses all reads (`--subsample-rate 1.0`).
Small read sets can make DNApi faster (For more detail, see [Tips for
making the exhaustive search mode faster]
(https://github.com/jnktsj/DNApi/tree/master/examples#tips-for-making-the-exhaustive-search-mode-faster)).

###### --output-dir DIRECTORY
Output directory for cleansed reads after a computation of
*exhaustive* mode. If the input FASTQ is not processed, DNApi removes
predicted 3′ adapters from the reads and generates a FASTA file
containing cleansed reads. In the default setting, DNApi creates the
output in the current directory as `./dnapi_out`.

###### --no-output-files
Suppress the output of the report and the cleansed reads, and only
display report on the screen.

###### --temp-dir DIRECTORY
Place for the temporary directory. DNApi creates a temporary directory
during a computation of *exhaustive* mode. In the default setting, the
program makes the directory in `/tmp`.

##### Evaluation of 3′ adapter candidates

###### --adapter-seq SEQ [SEQ ...]
A list of 3′ adapter(s) for quality control.  When the option is
specified, DNApi maps the processed reads after clipping each 3′
adapter in every run and checks the genome mapping rate.

##### Adapter removal parameters

###### --prefix-match LENGTH
3′ adapter prefix match length. DNApi only considers perfect adapter
matches. The default is 7nt. This option affects the length of
predicted 3′ adapter sequences in the final output.

###### --min-len LENGTH
Minimum read length to keep for mapping. Extracted small RNA reads
will be discarded if the lengths are *shorter* than the specified
length with `--min-len`. The default is 16nt.

###### --max-len LENGTH
Maximum read length to keep for mapping. Extracted small RNA reads
will be discarded if the lengths are *longer* than the specified
length with `--max-len`. The default is 36nt.

###### --trim-5p LENGTH
Trim specified number of bases from 5′ ends after adapter removal.
This option will be combined with the adapter clipping process to
trim down specific number of bases additionally.

###### --trim-3p LENGTH
Trim specified number of bases from 3′ ends after adapter removal.
This option can be combined with the adapter clipping process to
trim down specific number of bases additionally.


## API

You can access the adapter prediction algorithm once you import
`dnapilib.apred` in your python
program. `iterative_adapter_prediction` and `adapter_prediction` are
the core function for *iterative* and *single* adapter prediciton. It
takes four arguments: FASTQ file name, filtering ratio, k-mer size,
and subsampling read count. As the result, the two functions return
the list of tuples containing predicted 3′ adapters and the assembly
scores. The returned list is sorted by the scores.

```python
from dnapilib.apred import adapter_prediction
from dnapilib.apred import iterative_adapter_prediction

# [iterative mode]
# iterative_adapter_prediction(FASTQ, ratios, k-mers, subsample_read_count, length_to_print=12)
iterative_result = iterative_adapter_prediction("examples/good.fq", [1.2, 1.3, 1.4], [9, 11], 50000)

# [single mode]
# adapter_prediction(FASTQ, ratio, k-mer, subsample_read_count)
single_result = adapter_prediction("examples/good.fq", 1.4, 9, 50000)

# all predicted adapters
print(iterative_result)
# >> [('TGGAATTCTCGG', 200.0)]
print(single_result)
# >> [('TGGAATTCTCGGGTGCCAAGGAACTCC', 100.0)]

# predicted adapter with the highest score
print(iterative_result[0][0])
# >> 'TGGAATTCTCGGG'
print(single_result[0][0])
# >> 'TGGAATTCTCGGGTGCCAAGGAACTCC'
```


## Utilities
In addition to DNApi, there are potentially useful three programs in
the `utils` directory:
* `qual-offset.py` estimates ASCII-encoded quality score offsets of
  FASTQ files.
* `qual-trim.py` trims low quality bases in input FASTQ reads. The
  quality trimming algorithm in the program is the same as the one in
  BWA.
* `to-fasta.py` removes specified 5′ and/or 3′ adapter sequences,
  merges identical reads while retainig the counts, and writes the
  collapsed reads as FASTA in standard output (`stdout`).

To see the usage for each program, type:

    $ python3 <program-name> [-h | --help]


## Limitations
DNApi has a few limitations on 3′ adapter prediction:
* Poly(A) or other low-complexity 3′ adapters can't be predicted due
  to the low-complexity k-mer filtering step.
* Prediction accuracy will drop if gel-extracted lengths of RNAs are
  long enough to be sequenced, i.e. if few 3′ adapters are in FASTQ.
* DNApi can't do demultiplexing.
