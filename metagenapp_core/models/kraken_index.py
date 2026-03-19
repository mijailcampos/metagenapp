def build_kraken_index(raw_model, max_taxa_per_kmer=1000):

    k = raw_model["kmer_size"]
    taxa = raw_model["taxonomy_labels"]
    kmer_counts = raw_model["kmer_counts"]  # dict[kmer] -> dict[taxid]

    # ------------------------------------------------
    # 1️⃣ contar en cuántos taxones aparece cada kmer
    # ------------------------------------------------

    kmer_taxon_freq = {}

    for kmer, tax_dict in kmer_counts.items():
        kmer_taxon_freq[kmer] = len(tax_dict)

    # ------------------------------------------------
    # 2️⃣ construir índice correctamente
    # ------------------------------------------------

    kmer_index = {}

    for kmer, tax_dict in kmer_counts.items():

        # filtrar kmers demasiado comunes (ruido)
        if kmer_taxon_freq[kmer] > max_taxa_per_kmer:
            continue

        # opcional: normalizar pesos
        total = sum(tax_dict.values())

        kmer_index[kmer] = {
            taxid: weight / total
            for taxid, weight in tax_dict.items()
        }

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
