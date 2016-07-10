## Examples
This direcotry contains the following example FASTQ files:
* `good.fq`: good quality FASTQ (3′ adapter: `TGGAATTC`)
* `poor.fq`: poor quality FASTQ (3′ adapter: `CGCCTTGG`)
* `processed.fq`: reads already processed, i.e., FASTQ that do not
  contain any 3′ adapters


### Adapter prediction: *iterative* and *single* modes
To predict the 3′ adapter sequence in `good.fq`, type:

```shell
$ python3 dnapi.py good.fq
# you will get: TGGAATTCTCGG
```

The default setting is *iterative* search mode. If you want to use
different k-mer sizes and filtering ratios, type:

```shell
$ python3 dnapi.py -k 8:11:1 -r 1.1:1.4:0.1 good.fq
# you will get: TGGAATTCTCGG
```

If you want to see all predicted adapter candidates from the run, you
can add `--show-all` like below:

```shell
$ python3 dnapi.py --show-all poor.fq
# you will get:
# CGCCTTGGCCGT    score=200.00
# AGCAGAAGGGG     score=30.25
```

The option prints the predicted adapter sequences with the assembly
scores. The higher score indicates more reliable prediction results.

If you want to use only a single combination of k-mer size and
filtering ratio (i.e. *single* search mode), type:

```shell
$ python3 dnapi.py -k 9 -r 1.4 good.fq
# you will get: TGGAATTCTCGGGTGCCAAGGAACTCC
```

Adapter search with either lower k-mer sizes and lower filtering
ratios gives you mostly less accurate results.


### Adapter search with mapping process and quality control: *exhaustive* mode
DNApi executes a given mapping command for read mapping after adapter
removal. Adding `--map-command` to arguments turns on this
*exhaustive* search mode. For the mapping step, you need to prepare
the reference genome index first. Although you can use any read
mapping software packages for mapping, we will use
[Bowtie](http://bowtie-bio.sourceforge.net) in this example.

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

    $ python3 dnapi.py --map-command "/path_to/bowtie /path_to/genome_index -p <cpu_num> -v0 -k1 -S -f @in > @out" good.fq

You will get:

    optimal_3'adapter=TGGAATTCTCGG

    # sampled_reads=100000 (total_reads * 1.00)
    # 3'adapter     reads_extracted  (reads_extracted/sampled_reads)%  reads_mapped  (reads_mapped/sampled_reads)%  params_k:r
      TGGAATTCTCGG   95845            95.84                            84129          84.13                         9:1.2;9:1.3;9:1.4;11:1.2;11:1.3;11:1.4
      RAW_INPUT     100000           100.00                               65           0.07                         NO_TREATMENT

Each run generates the above report and the cleansed reads as files in
an output directory (default: `./dnapi_out`). You can control the path
and the name of the output directory by `--output-dir` or suppress the
output by `--no-output-files`.

##### Case 2: `processed.fq`
If the adapter was already clipped from the reads, you will get:

    optimal_3'adapter=RAW_INPUT

    # sampled_reads=100000 (total_reads * 1.00)
    # 3'adapter     reads_extracted  (reads_extracted/sampled_reads)%  reads_mapped  (reads_mapped/sampled_reads)%  params_k:r
      TAATACTGCCTG	 0             0.00                                0           0.00                         9:1.2;9:1.4;11:1.3;11:1.4
      TGGCAGTGTCTT       0             0.00                                0           0.00                         9:1.3;11:1.2
      RAW_INPUT     100000           100.00                            75031          75.03                         NO_TREATMENT
    # input reads look already clean!

Note that the FASTA output will not be generated because the input
FASTQ were already processed.

##### Case 3: `poor.fq`
When the quality of reads in a FASTQ is poor, you will get:

    optimal_3'adapter=CGCCTTGGCCGT/POOR_QUALITY

    # sampled_reads=100000 (total_reads * 1.00)
    # 3'adapter     reads_extracted  (reads_extracted/sampled_reads)%  reads_mapped  (reads_mapped/sampled_reads)%  params_k:r
      CGCCTTGGCCGT   15661            15.66                            4536           4.54                          9:1.2;9:1.3;9:1.4;11:1.2;11:1.3;11:1.4
      RAW_INPUT     100000           100.00                               3           0.00                          NO_TREATMENT

Even if the predicted 3' adapter sequence is correct, `/POOR_QUALITY`
will be attached when the mapping rate is below 20% for quality control.

##### Case 4: `user-input`
If you want to evaluate adapter sequences without prediction, you can
directly input the sequences with `--adapter-seq`. Using the option,
the following command shows the example of evaluating three adapter
sequences:

    $ python3 dnapi.py --adapter-seq TGGAATTCTCGG TAATACTGCCTG CGCCTTGGCCGT --map-command "/path_to/bowtie /path_to/genome_index -p <cpu_num> -v0 -k1 -S -f @in > @out" good.fq

From the above command, You will get:

    optimal_3'adapter=TGGAATTCTCGG

    # sampled_reads=100000 (total_reads * 1.00)
    # 3'adapter     reads_extracted  (reads_extracted/sampled_reads)%  reads_mapped  (reads_mapped/sampled_reads)%  params_k:r
      TAATACTGCCTG      2             0.00                                 2          0.00                          user-input
      CGCCTTGGCCGT   1206             1.21                              1062          1.06                          user-input
      TGGAATTCTCGG  95845            95.84                             84129         84.13                          user-input


#### Tips for making the exhaustive search mode faster
In the default setting, DNApi processes all reads in an input FASTQ.
It is possible to make the computation faster if you subsample a
portion of reads with `--subsample-rate`. The following shell script
snippet shows the example to use 1 million reads from an input FASTQ
for DNApi run.

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

# Input the subsample rate with -p to DNApi
python3 dnapi.py --subsample-rate ${SUBSAMPLE_RATE} --map-command "${MAPCMD}" ${FASTQ}
````

One caveat is that the output FASTA will be based on **only subsampled
reads**. If you want to get all reads processed, you need to run
DNApi without subsampling reads or run other software packages with
the predicted adapter by DNApi.

Although this example used [Bowtie](http://bowtie-bio.sourceforge.net)
to map reads, you can try any command lines and read mapping software
packages you like.
