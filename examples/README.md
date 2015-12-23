## Examples
This direcotry contains the following example FASTQ files:
* `good.fq`: good quality FASTQ (3′ adapter: `TGGAATTC`)
* `poor.fq`: poor quality FASTQ (3′ adapter: `CGCCTTGG`)
* `processed.fq`: reads already processed, i.e., FASTQ that do not
  contain any 3′ adapters

### `dnap`: Adapter prediction
To predict the 3′ adapter sequence in `good.fq`, type:

```shell
$ dnap good.fq
# you will get:
# Predicted_3'adapter_1=TGGAATTCTCGGGTGCCAAGGAACTCC
```

If you want to use different k-mer size and filtering rate,
type:

```shell
$ dnap -k 8 -r 1.1 good.fq
# you will get:
# Predicted_3'adapter_1=TGGAATTCTCGGGTGCCAAGGAACTC
```

Adapter search with lower k-mer sizes and lower filtering rates would
be lenient and would lead to (mostly) less accurate results.

### `dnapi`: Iterative adapter search and quality control
`dnapi` maps reads with a given mapping command after adapter removal.
For the mapping step, let's prepare the reference genome index first.
Although you can use any read mapping software packages for mapping,
we will use [Bowtie](http://bowtie-bio.sourceforge.net) in this example.

To generate the genome index, download human genome ([hg38]
(http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz)),
and then:

```shell
# concatenate all chromosomes into one file
cat *.fa > hg38.fa
bowtie-build hg38.fa <index_name>
```

##### Case 1: `good.fq`
With the following command:

    $ dnapi "/path_to/bowtie /path_to/genome_index -p <cpu_num> -v0 -k1 -S -f @in > @out" good.fq

You will get:

    Cleansed_reads=./good_TGGAATT.fa
    Optimal_3'adapter=TGGAATTCTCGG

    # Report: sampled_reads=100000 (total_reads * 1.00)
    # 3'adapter     reads_extracted  (reads_extracted/sampled_reads)%  reads_mapped  (reads_mapped/sampled_reads)%  params_k:r
      TGGAATTCTCGG   95845            95.84                            84409          84.41                         9:1.2;9:1.3;9:1.4;11:1.2;11:1.3;11:1.4
      RAW_INPUT     100000           100.00                               66           0.07                         NO_TREATMENT

##### Case 2: `processed.fq`
If the adapter was already clipped from the reads, you will get:

    Optimal_3'adapter=RAW_INPUT

    # Report: sampled_reads=100000 (total_reads * 1.00)
    # 3'adapter     reads_extracted  (reads_extracted/sampled_reads)%  reads_mapped  (reads_mapped/sampled_reads)%  params_k:r
      TGGCAGTGTCTT       0             0.00                                0           0.00                         9:1.3;11:1.2
      TAATACTGCCTG	 0             0.00                                0           0.00                         9:1.2;9:1.4;11:1.3;11:1.4
      RAW_INPUT     100000           100.00                            75031          75.03                         NO_TREATMENT
    # Input reads look already clean!

Note that the FASTA output will not be generated because the input
FASTQ were already processed.

##### Case 3: `poor.fq`
When the quality of reads in a FASTQ is poor, you will get:

    Cleansed_reads=./poor_CGCCTTG.fa
    Optimal_3'adapter=CGCCTTGGCCGT/POOR_QUALITY

    # Report: sampled_reads=100000 (total_reads * 1.00)
    # 3'adapter     reads_extracted  (reads_extracted/sampled_reads)%  reads_mapped  (reads_mapped/sampled_reads)%  params_k:r
      CGCCTTGGCCGT   15661            15.66                            4617           4.62                          9:1.2;9:1.3;9:1.4;11:1.2;11:1.3;11:1.4
      RAW_INPUT     100000           100.00                               3           0.00                          NO_TREATMENT

Even if the predicted 3' adapter sequence is correct, `/POOR_QUALITY`
will be attached when the mapping rate is below 20 for quality control.

#### Tips for making `dnapi` faster
In the default setting, `dnapi` processes all reads in an input FASTQ.
It is possible to make the computation faster if you subsample a
portion of reads with `-p`. The following shell script snippet shows
the example to use 1 million reads from an input FASTQ for `dnapi` run.

```shell
# !/bin/sh

# Bowtie command as a mapping engine
MAPCMD="/path_to/bowtie /path_to/genome_index -p${CPU} -v0 -k1 -S -f @in > @out"

SUBSAMPLE_RATE=1.0
SUBSAMPLE_READS=1000000
TOTAL_READS=`wc -l ${FASTQ} | awk '{print $1/4}'`
if [ ${TOTAL_READS} -gt ${SUBSAMPLE_READS} ]; then
    # Compute a subsample rate if the total reads are more than 1 million
    SUBSAMPLE_RATE=`echo ${SUBSAMPLE} ${TOTAL_READS} | awk '{print $1/$2}'`
fi

# Input the subsample rate with -p to dnapi
dnapi -p ${SUBSAMPLE_RATE} "${MAPCMD}" ${FASTQ}
````

One caveat is that the output FASTA will be made based on subsampled
reads. If you want to get all reads processed, you need to run `dnapi`
without subsampling reads or run other software packages with the
predicted adapter by `dnapi`.

Although this example used [Bowtie](http://bowtie-bio.sourceforge.net)
to map reads, you can try any command lines and read mapping software
packages you like.
