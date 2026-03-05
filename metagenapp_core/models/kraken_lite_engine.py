from metagenapp_core.utils.kmers import generate_kmers
from metagenapp_core.utils.taxonomy import lca_taxa

MODEL_NAME = "kraken_lite"
MODEL_TYPE = "taxonomy_classifier"


def classify_seq_kraken(seq_id, seq, model):

    k = model["k"]
    kmer_index = model["kmer_index"]
    taxonomy = model["taxonomy"]

    votes = {}

    for kmer in generate_kmers(seq, k):

        taxid = kmer_index.get(kmer)

        if taxid is None:
            continue

        votes[taxid] = votes.get(taxid, 0) + 1

    if not votes:
        return seq_id, "Unclassified"

    sorted_votes = sorted(votes.items(), key=lambda x: x[1], reverse=True)

    best_taxid, best_votes = sorted_votes[0]

    # si el ganador tiene al menos el doble de votos del segundo
    if len(sorted_votes) == 1 or best_votes >= 2 * sorted_votes[1][1]:

        return seq_id, taxonomy[best_taxid]

    # empate fuerte → usar LCA
    top_taxa = [taxonomy[taxid] for taxid,_ in sorted_votes[:3]]

    tax = lca_taxa(top_taxa)

    return seq_id, tax
