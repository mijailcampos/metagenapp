def build_kraken_index(raw_model, max_taxa_per_kmer=50):

    k = raw_model["kmer_size"]
    taxa = raw_model["taxonomy_labels"]
    kmer_counts = raw_model["kmer_counts"]

    # ------------------------------------------------
    # 1️⃣ contar en cuántos taxones aparece cada kmer
    # ------------------------------------------------

    kmer_taxon_freq = {}

    for taxon, kmers in kmer_counts.items():
        for kmer in kmers:

            kmer_taxon_freq[kmer] = kmer_taxon_freq.get(kmer, 0) + 1

    # ------------------------------------------------
    # 2️⃣ construir índice filtrando kmers comunes
    # ------------------------------------------------

    kmer_index = {}

    for taxid, (taxon, kmers) in enumerate(kmer_counts.items()):

        for kmer, weight in kmers.items():

            # filtrar kmers no discriminativos
            if kmer_taxon_freq[kmer] > max_taxa_per_kmer:
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
