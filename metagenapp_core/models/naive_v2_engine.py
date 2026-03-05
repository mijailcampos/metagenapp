from collections import Counter, defaultdict
from metagenapp_core.utils.kmers import generate_kmers


def _prefixes(lineage: str, max_level: int = 6):
    parts = [p for p in lineage.split(";") if p.strip()]
    out = []
    for i in range(1, min(len(parts), max_level) + 1):
        out.append(";".join(parts[:i]) + ";")
    return out


def classify_seq_v2(seq_id: str, seq: str, model_v2: dict):

    k = model_v2["k"]
    kmer_to_id = model_v2["kmer_to_id"]
    kid_to_lca = model_v2.get("kid_to_lca")

    postings_taxon = model_v2.get("postings_taxon")
    postings_weight = model_v2.get("postings_weight")
    taxa = model_v2.get("taxa")

    counts = Counter(generate_kmers(seq, k))

    # ---------------------------
    # Caso v3
    # ---------------------------
    if kid_to_lca is not None:

        scores = defaultdict(float)

        for kmer, c in counts.items():

            kid = kmer_to_id.get(kmer)
            if kid is None:
                continue

            lca = kid_to_lca[kid]
            if not lca:
                continue

            psize = len(postings_taxon[kid]) if postings_taxon is not None else 1

            if psize >= 20:
                continue

            inform = 1.0 / psize

            for pref in _prefixes(lca, max_level=6):
                scores[pref] += inform * c

        if not scores:
            return seq_id, None

        # agrupar por nivel
        by_level = {i: [] for i in range(1,7)}

        for tax,score in scores.items():

            lvl = tax.count(";")

            if lvl > 6:
                lvl = 6

            by_level[lvl].append((tax,score))

        # mejor taxón de cada nivel
        best = {}

        for lvl in range(1,7):
            if by_level[lvl]:
                best[lvl] = max(by_level[lvl], key=lambda x:x[1])

        # fallback dominio
        chosen_tax, chosen_score = best.get(1,(None,0))

        # intentar bajar a filo
        if 2 in best:
            chosen_tax, chosen_score = best[2]

        # intentar bajar niveles
        for lvl in range(3,7):

            if lvl not in best:
                break

            tax,score = best[lvl]

            ratio = score/(chosen_score+1e-9)

            if ratio > 0.85:
                chosen_tax = tax
                chosen_score = score
            else:
                break

        return seq_id, chosen_tax

    # ---------------------------
    # Caso viejo
    # ---------------------------
    if taxa is None or postings_taxon is None or postings_weight is None:
        return seq_id, None

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