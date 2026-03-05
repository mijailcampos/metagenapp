from collections import Counter
from metagenapp_core.utils.kmers import generate_kmers


def classify_seq_v2(seq_id: str, seq: str, model_v2: dict):

    k = model_v2["k"]
    taxa = model_v2["taxa"]
    kmer_to_id = model_v2["kmer_to_id"]
    postings_taxon = model_v2["postings_taxon"]
    postings_weight = model_v2["postings_weight"]

    counts = Counter(generate_kmers(seq, k))

    scores = [0.0] * len(taxa)

    for kmer, c in counts.items():

        kid = kmer_to_id.get(kmer)
        if kid is None:
            continue

        tids = postings_taxon[kid]
        ws = postings_weight[kid]

        for j in range(len(tids)):
            scores[tids[j]] += ws[j] * c

    best_tid = 0
    best_score = scores[0]

    for i in range(1, len(scores)):
        s = scores[i]
        if s > best_score:
            best_score = s
            best_tid = i

    return seq_id, taxa[best_tid]
