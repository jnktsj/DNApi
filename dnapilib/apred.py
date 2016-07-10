"""Functions for adapter prediction.

"""

from operator import itemgetter

from dnapilib.io_utils import get_file_obj, fastq_sequence
from dnapilib.kmer import count_kmers, filter_kmers, assemble_kmers


def adapter_prediction(fastq, ratio, kmer_len, sample_num):
    """Return a list of predicted adapters.

       Predict 3' adapter sequence with a combination of k and R.
    """
    fq_obj = get_file_obj(fastq)
    fq_seq = fastq_sequence(fq_obj)
    freq = count_kmers(fq_seq, kmer_len, sample_num)
    clean = filter_kmers(freq, kmer_len, ratio)
    assembl = sorted(assemble_kmers(clean, kmer_len//2),
                     key=itemgetter(1), reverse=True)
    fq_obj.close()
    return assembl

def iterative_adapter_prediction(fastq, ratios, kmer_lens,
                                 sample_num, keep_len=12):
    """Return a list of predicted adapters.

       Iteratively predict 3' adapter sequence with different
       combinations of k and R.
    """
    fq_seq = []
    fq_obj = get_file_obj(fastq)
    for i, s in enumerate(fastq_sequence(fq_obj)):
        if i == sample_num:
            break
        fq_seq.append(s)
    fq_obj.close()

    collection = {}
    for kmer_len in kmer_lens:
        curated = {}
        freq = count_kmers(fq_seq, kmer_len, sample_num)
        for ratio in ratios:
            clean = filter_kmers(freq, kmer_len, ratio)
            assembl = assemble_kmers(clean, kmer_len//2)
            for s, c in assembl:
                key = s[:keep_len]
                curated[key] = max(curated.get(key,0), c)
        for s, c in curated.items():
            collection[s] = round(collection.get(s,0)+c, 4)
    asmbl_min_len = min(map(len, collection.keys()))
    assembl = sorted(assemble_kmers(list(collection.items()), asmbl_min_len//2),
                     key=itemgetter(1), reverse=True)
    return assembl
