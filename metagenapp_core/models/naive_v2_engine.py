import random
from collections import Counter, defaultdict
from metagenapp_core.utils.kmers import generate_kmers


def _prefixes(lineage: str, max_level: int = 6):
    parts = [p for p in lineage.split(";") if p.strip()]
    out = []
    for i in range(1, min(len(parts), max_level) + 1):
        out.append(";".join(parts[:i]) + ";")
    return out


def classify_seq_v2(
    seq_id: str,
    seq: str,
    model_v2: dict,
    psize_max: int = None,
    descent_ratio: float = 0.85,
    confidence: float = None,
    n_bootstrap: int = 100,
):
    # Usar parámetros recomendados del modelo si no se especifican explícitamente.
    # Esto permite que modelos grandes (SILVA) usen umbrales calibrados automáticamente.
    if psize_max is None:
        psize_max = model_v2.get("recommended_psize_max", 20)
    if confidence is None:
        confidence = model_v2.get("recommended_confidence", 0.80)

    k          = model_v2["k"]
    kmer_to_id = model_v2["kmer_to_id"]
    taxa       = model_v2.get("taxa")
    kid_to_lca = model_v2.get("kid_to_lca")

    # ------------------------------------------------------------------
    # Detectar formato del modelo:
    #   csr_v1  — postings en arrays numpy (pt_data / pw_data / pt_indptr)
    #   list_v1 — postings en listas Python (postings_taxon / postings_weight)
    #   legacy  — sin postings, solo kid_to_lca
    # ------------------------------------------------------------------
    is_csr = "pt_data" in model_v2
    if is_csr:
        pt_data   = model_v2["pt_data"]
        pw_data   = model_v2["pw_data"]
        pt_indptr = model_v2["pt_indptr"]
        has_postings = True
    else:
        postings_taxon  = model_v2.get("postings_taxon")
        postings_weight = model_v2.get("postings_weight")
        has_postings = postings_taxon is not None

    counts = Counter(generate_kmers(seq, k))

    # ------------------------------------------------------------------
    # Wang-style bootstrap confidence
    #
    # Activado cuando el modelo tiene postings (cualquier formato) y taxa.
    # NO requiere kid_to_lca.
    #
    # Replica Wang et al. (2007) / Mothur:
    #   1. Para cada k-mer informativo, votar por todos los taxa en postings
    #      con peso uniforme 1/psize.
    #   2. Bootstrap: n_bootstrap iteraciones de muestreo con reemplazo.
    #   3. Confianza = % iteraciones donde el mismo taxón gana en cada nivel.
    #   4. Asignar al nivel más profundo con confianza >= threshold.
    # ------------------------------------------------------------------
    if has_postings and taxa is not None:

        pool     = []   # lista de kids (repetido según count del k-mer en la query)
        kid_data = {}   # kid → lista de (weight, prefixes) — caché

        for kmer, c in counts.items():
            kid = kmer_to_id.get(kmer)
            if kid is None:
                continue

            # Obtener tids según formato
            if is_csr:
                s, e = int(pt_indptr[kid]), int(pt_indptr[kid + 1])
                tids  = pt_data[s:e]
                psize = e - s
            else:
                tids  = postings_taxon[kid]
                psize = len(tids) if tids else 0

            if psize == 0 or psize >= psize_max:
                continue

            pool.extend([kid] * c)

            if kid not in kid_data:
                w = 1.0 / psize
                kid_data[kid] = [
                    (w, _prefixes(taxa[int(tid)], max_level=6))
                    for tid in tids
                ]

        if not pool:
            return seq_id, None

        n = len(pool)

        # --- Paso 2: clasificación completa (sin bootstrap) ---
        taxon_votes = defaultdict(float)
        for kid in pool:
            for w, prefs in kid_data[kid]:
                for pref in prefs:
                    taxon_votes[pref] += w

        by_level = defaultdict(dict)
        for tax, votes in taxon_votes.items():
            lvl = tax.count(";")
            if 1 <= lvl <= 6:
                by_level[lvl][tax] = votes

        top_at_level = {}
        for lvl in range(1, 7):
            if by_level.get(lvl):
                top_at_level[lvl] = max(by_level[lvl], key=lambda x: by_level[lvl][x])

        if not top_at_level:
            return seq_id, None

        # --- Paso 3: bootstrap para estimar confianza ---
        wins = defaultdict(int)
        for _ in range(n_bootstrap):
            bs_votes = defaultdict(float)
            for idx in random.choices(range(n), k=n):
                kid = pool[idx]
                for w, prefs in kid_data[kid]:
                    for pref in prefs:
                        bs_votes[pref] += w

            for lvl, top_tax in top_at_level.items():
                bs_lv = {t: v for t, v in bs_votes.items() if t.count(";") == lvl}
                if bs_lv and max(bs_lv, key=lambda x: bs_lv[x]) == top_tax:
                    wins[lvl] += 1

        # --- Paso 4: nivel más profundo con confianza >= threshold ---
        chosen_tax = None
        for lvl in range(1, 7):
            if lvl not in top_at_level:
                break
            if wins[lvl] / n_bootstrap >= confidence:
                chosen_tax = top_at_level[lvl]
            else:
                break

        return seq_id, chosen_tax

    # ------------------------------------------------------------------
    # Fallback LCA (modelos con kid_to_lca pero sin postings completos)
    # ------------------------------------------------------------------
    if kid_to_lca is not None:

        scores = defaultdict(float)

        for kmer, c in counts.items():
            kid = kmer_to_id.get(kmer)
            if kid is None:
                continue

            lca = kid_to_lca[kid]
            if not lca:
                continue

            if has_postings:
                psize = len(postings_taxon[kid]) if postings_taxon[kid] else 1
            else:
                psize = 1

            if psize >= psize_max:
                continue

            inform = 1.0 / psize

            for pref in _prefixes(lca, max_level=6):
                scores[pref] += inform * c

        if not scores:
            return seq_id, None

        by_level = {i: [] for i in range(1, 7)}
        for tax, score in scores.items():
            lvl = tax.count(";")
            if lvl > 6:
                lvl = 6
            by_level[lvl].append((tax, score))

        best = {}
        for lvl in range(1, 7):
            if by_level[lvl]:
                best[lvl] = max(by_level[lvl], key=lambda x: x[1])

        chosen_tax, chosen_score = best.get(1, (None, 0))
        if 2 in best:
            chosen_tax, chosen_score = best[2]

        for lvl in range(3, 7):
            if lvl not in best:
                break
            tax, score = best[lvl]
            if score / (chosen_score + 1e-9) > descent_ratio:
                chosen_tax = tax
                chosen_score = score
            else:
                break

        return seq_id, chosen_tax

    # ------------------------------------------------------------------
    # Formato legacy (sin postings ni kid_to_lca)
    # ------------------------------------------------------------------
    if taxa is None or not has_postings:
        return seq_id, None

    scores = [0.0] * len(taxa)
    for kmer, c in counts.items():
        kid = kmer_to_id.get(kmer)
        if kid is None:
            continue
        tids = postings_taxon[kid]
        ws   = postings_weight[kid]
        for j in range(len(tids)):
            scores[tids[j]] += ws[j] * c

    best_tid   = 0
    best_score = scores[0]
    for i in range(1, len(scores)):
        s = scores[i]
        if s > best_score:
            best_score = s
            best_tid   = i

    return seq_id, taxa[best_tid]
