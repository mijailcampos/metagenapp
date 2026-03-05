from metagenapp_core.utils.minimizers import generate_minimizers
from metagenapp_core.utils.taxonomy import lca_taxa


def classify_seq_minimizer(seq_id, seq, model):

    k = model["k"]
    m = model["m"]
    index = model["minimizer_index"]
    taxonomy = model["taxonomy"]

    taxa_hits = []

    for minim in generate_minimizers(seq, k, m):

        taxids = index.get(minim)

        if taxids is None:
            continue

        for taxid in taxids:
            taxa_hits.append(taxonomy[taxid])

    if not taxa_hits:
        return seq_id, "Unclassified"

    tax = lca_taxa(taxa_hits)

    return seq_id, tax