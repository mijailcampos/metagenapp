import math


def build_kraken_index(raw_model, max_taxa_per_kmer=100000):
    """
    Construye el índice kraken-lite con pesos TF-IDF.

    Recibe raw_model con kmer_counts YA INVERTIDO:
      kmer_counts: dict[kmer] -> dict[taxid] -> tf_count

    TF  = frecuencia del kmer en el taxon
    IDF = log( N_taxa / df_kmer ) + 1
        donde df_kmer = len(kmer_counts[kmer]) = número de taxa con ese kmer

    Efecto:
      - Kmers diagnósticos (pocos taxa) → IDF alto → peso alto
      - Kmers ubicuos (muchos taxa)     → IDF bajo → peso bajo
      - NO se elimina ningún kmer por IDF — solo se reponderan
    """
    k           = raw_model["kmer_size"]
    taxa        = raw_model["taxonomy_labels"]
    kmer_counts = raw_model["kmer_counts"]  # dict[kmer] -> dict[taxid] -> tf

    N_taxa = len(taxa)

    # ------------------------------------------------
    # 1. construir índice con peso TF-IDF
    #    kmer_counts ya está invertido: kmer -> {taxid: tf}
    # ------------------------------------------------
    kmer_index = {}

    for kmer, taxid_counts in kmer_counts.items():

        df = len(taxid_counts)  # cuántos taxa tienen este kmer

        # filtrar kmers que aparecen en demasiados taxa
        if df > max_taxa_per_kmer:
            continue

        # IDF suavizado: log(N / df) + 1
        # +1 para que kmers únicos no dominen demasiado
        idf = math.log(N_taxa / df) + 1.0

        entry = {}
        for taxid, tf in taxid_counts.items():
            entry[taxid] = tf * idf

        kmer_index[kmer] = entry

    return {
        "k": k,
        "taxonomy": taxa,
        "kmer_index": kmer_index
    }


def invert_kmer_counts_by_taxon(kmer_counts_by_taxon, taxa_labels):
    """
    Input:  kmer_counts_by_taxon: dict[taxon_label] -> dict[kmer] -> weight
            taxa_labels: list de taxon_label en orden taxid (si aplica)
    Output: dict[kmer] -> dict[taxid] -> weight
    """
    kmer_counts = {}
    for taxid, (taxon_label, kmers) in enumerate(kmer_counts_by_taxon.items()):
        for kmer, w in kmers.items():
            d = kmer_counts.get(kmer)
            if d is None:
                d = {}
                kmer_counts[kmer] = d
            d[taxid] = w
    return kmer_counts