def generate_minimizers(seq, k=8, m=4):

    seq = seq.upper()

    n = len(seq)

    for i in range(n - k + 1):

        kmer = seq[i:i+k]

        if "N" in kmer:
            continue

        minim = None

        for j in range(k - m + 1):

            sub = kmer[j:j+m]

            if minim is None or sub < minim:
                minim = sub

        if minim:
            yield minim
