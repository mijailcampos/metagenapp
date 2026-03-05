def adapt_naive_model(model):

    # modelo antiguo
    k = model["kmer_size"]
    taxa = model["taxonomy_labels"]
    kmer_counts = model["kmer_counts"]

    kmer_to_id = {}
    postings_taxon = []
    postings_weight = []

    kid = 0

    for kmer, counts in kmer_counts.items():

        kmer_to_id[kmer] = kid

        tids = []
        ws = []

        for taxid, weight in counts.items():
            tids.append(taxid)
            ws.append(weight)

        postings_taxon.append(tids)
        postings_weight.append(ws)

        kid += 1

    return {
        "k": k,
        "taxa": taxa,
        "kmer_to_id": kmer_to_id,
        "postings_taxon": postings_taxon,
        "postings_weight": postings_weight,
        "kmer_counts": kmer_counts
    }
