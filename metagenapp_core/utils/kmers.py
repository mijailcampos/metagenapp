def generate_kmers(seq, k):

    seq = seq.upper()

    n = len(seq) - k + 1

    for i in range(n):
        kmer = seq[i:i+k]

        if "N" in kmer:
            continue

        yield kmer
