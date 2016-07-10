"""Functions for exhaustive adapter search
   incorporating with read mapping process.

"""

import re
import os.path
import subprocess
import fileinput

from dnapilib.io_utils import get_file_obj
from dnapilib.io_utils import fastq_sequence
from dnapilib.io_utils import fastq_record


def rm_temp_dir(temp_dir):
    """Remove temporary directory.

    """
    if temp_dir:
        if os.path.exists(temp_dir):
            subprocess.call("rm -r {}".format(temp_dir).split())


def clip_adapter(fp, aseed, tm5, tm3, min_len, max_len):
    """Return adapter-clipped clean reads.

    """
    seed_len = len(aseed)
    pp = re.compile("(.*)"+aseed, re.IGNORECASE)
    for seq in fastq_sequence(fp):
        if len(seq) < tm5 or len(seq) < tm3:
            raise Exception("trimming length is too large")
        match = pp.search(seq)
        if not match:
            continue
        end = match.end() - seed_len
        clipped_seq = seq[tm5 : end-tm3]
        L = len(clipped_seq)
        if min_len <= L and L <= max_len:
            yield clipped_seq


def to_fasta(fastq, fasta, aseed, tm5, tm3, min_len, max_len):
    """Write FASTA containing clean reads, and return
       the number of the reads.

    """
    fq_obj = get_file_obj(fastq)
    if "RAW_INPUT".startswith(aseed):
        iterator = fastq_sequence(fq_obj)
    else:
        iterator = clip_adapter(fq_obj, aseed, tm5, tm3, min_len, max_len)
    fas = {}
    clean_read_count = 0
    for seq in iterator:
        fas[seq] = fas.get(seq, 0) + 1
    fa_obj = open(fasta, "w")
    for seq, cnt in fas.items():
        clean_read_count += cnt
        fa_obj.write(">{0}_{1}\n{0}\n".format(seq, cnt))
    fa_obj.close()
    fq_obj.close()
    return clean_read_count


def fastq_input_prep(fastq, ratio, temp_dir):
    """Write FASTQ in the temporary directory, and retrun
       (subsampled) FASTQ name, the total read count,
       standard deviation of read lengths.

    """
    num = int(1/ratio)
    read_count = 0.0
    stats = {}
    fq_out = "{}/input.fq".format(temp_dir)
    fq_obj = get_file_obj(fastq)
    fout = open(fq_out, "w")
    for i, rec in enumerate(fastq_record(fq_obj)):
        if i % num == 0:
            fout.write(rec)
            read_count += 1
            L = len(rec.split("\n")[1])
            stats[L] = stats.get(L,0) + 1
    fout.close()
    fq_obj.close()
    mean = sum([L*c for L,c in stats.items()]) / read_count
    sum_square = sum([(L-mean)**2 * c for L,c in stats.items()])
    sd = (sum_square / read_count)**0.5
    return fq_out, read_count, sd


def count_mapped_read_sam(samout):
    """Return the number of mapped reads to the genome.

    """
    if not os.path.exists(samout):
        raise Exception("can't open SAM")
    mapped = set()
    for x in fileinput.input(samout):
        if not x or x.startswith("@"):
            continue
        x = x.rstrip().split("\t")
        if x[2] != '*':
            mapped.add(x[0])
    cnt = sum([int(n.split('_')[1]) for n in mapped])
    return cnt


def map_clean_reads(fastq, adapter, tm5, tm3,
                    min_len, max_len, map_command, temp_dir):
    """Execute mapping command, and return the numbers
       of clean and mapped reads.

    """
    fasta = "{0}/insert_{1}.fa".format(temp_dir, adapter)
    samout = "{}/output.sam".format(temp_dir)
    clipped = to_fasta(fastq, fasta, adapter, tm5, tm3, min_len, max_len)
    map_command = map_command.replace("@in",fasta).replace("@out",samout)
    map_command += " 2> /dev/null"
    if subprocess.call(map_command, shell=True) != 0:
        raise Exception("mapping failed, check command line")
    mapped = count_mapped_read_sam(samout)
    return clipped, mapped


def make_stats_report(table, sampled_read, subsample_rate, prefix_match,
                      sd, fastq, output_dir, temp_dir, no_output_files):
    """Report read statistics with predicted adapters.

    """
    out = ["# sampled_reads={} (total_reads * {:.2f})".format(
              int(sampled_read), subsample_rate)]
    out.append("\t".join([
          "# 3'adapter",
          "reads_extracted",
          "(reads_extracted/sampled_reads)%",
          "reads_mapped",
          "(reads_mapped/sampled_reads)%",
          "params_k:r"]))
    max_mapped_read = -1
    max_index = -1
    for i, x in enumerate(table):
        if x[3] > max_mapped_read:
            max_mapped_read = x[3]
            max_index = i
        out.append("{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}".format(*x))
    optimal = [table[max_index][0]]

    fq_prefix = os.path.basename(fastq).split(".")[0]
    if table[max_index][4] < 20:
        optimal.append("/POOR_QUALITY")
    if optimal[0] == "RAW_INPUT":
        if sd:
            out.append("# input reads look already clean!")
        else:
            optimal.append("?")
    else:
        if no_output_files:
            pass
        else:
            if not os.path.exists(output_dir):
                subprocess.call("mkdir {}".format(output_dir).split())
            aseq = optimal[0][:prefix_match]
            fa_tmp = "{}/insert_{}.fa".format(temp_dir, aseq)
            fa_out = "{}/{}_{}.fa".format(output_dir, fq_prefix, aseq)
            subprocess.call(("mv {} {}".format(fa_tmp,fa_out)).split())

    out.insert(0, "optimal_3'adapter={}\n".format(''.join(optimal)))
    report = "\n".join(out)
    print(report)

    if not no_output_files:
        f = open("{}/{}_report.txt".format(output_dir, fq_prefix), "w")
        f.write(report + "\n")
        f.close()
