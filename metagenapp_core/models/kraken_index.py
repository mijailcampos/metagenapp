def build_kraken_index(raw_model):

    k = raw_model["kmer_size"]
    taxa = raw_model["taxonomy_labels"]
    kmer_counts = raw_model["kmer_counts"]

    kmer_index = {}

    for taxid, (taxon, kmers) in enumerate(kmer_counts.items()):

        for kmer, weight in kmers.items():

            # ignorar kmers poco informativos
            if weight < 1e-5:
                continue

            if kmer not in kmer_index:
                kmer_index[kmer] = taxid

            else:
                prev_taxid = kmer_index[kmer]

                prev_weight = kmer_counts[taxa[prev_taxid]].get(kmer, 0)

                if weight > prev_weight:
                    kmer_index[kmer] = taxid

    return {
        "k": k,
        "taxonomy": taxa,
        "kmer_index": kmer_index
    }


from metagenapp_core.utils.minimizers import generate_minimizers


def build_minimizer_index(raw_model, k=8, m=4):

    taxa = raw_model["taxonomy_labels"]
    kmer_counts = raw_model["kmer_counts"]

    minimizer_index = {}

    for taxid, (taxon, kmers) in enumerate(kmer_counts.items()):

        for kmer in kmers.keys():

            minim = min(kmer[i:i+m] for i in range(k-m+1))

            if minim not in minimizer_index:
                minimizer_index[minim] = set()

            minimizer_index[minim].add(taxid)

    return {
        "k": k,
        "m": m,
        "taxonomy": taxa,
        "minimizer_index": minimizer_index
    }
