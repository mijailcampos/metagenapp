from collections import Counter

from metagenapp_core.utils.kmers import generate_kmers
from metagenapp_core.utils.taxonomy import lca_taxa

MODEL_NAME = "kraken_lite"
MODEL_TYPE = "taxonomy_classifier"

# Threshold para resolución genus dentro del phylum dominante
CONFIDENCE_THRESHOLD = 0.30

# Fracción de votos que el top phylum necesita para considerarse "dominante"
# y restringir la resolución jerárquica a ese phylum
PHYLUM_CONFIDENCE = 0.55

# Fracción mínima del top phylum para retornar a nivel phylum
# (zona gris: no hay phylum dominante pero hay un ganador parcial)
PHYLUM_MIN = 0.25

# Threshold mínimo para LCA fallback
LCA_MIN_FRACTION = 0.05

MAX_LCA_HITS = 10


DNA_MAP = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}


def encode_kmer(seq):
    val = 0
    for b in seq:
        if b not in DNA_MAP:
            return None
        val = (val << 2) | DNA_MAP[b]
    return val


def classify_seq_kraken(seq_id, seq, model):

    seq = seq.upper()

    k = model["k"]
    # k_lca permite usar un k distinto para kmer_lca (ej. k=13 para kmer_index,
    # k=15 para kmer_lca). Si no está presente, usa el mismo k que kmer_index.
    k_lca = model.get("kmer_lca_k", k)
    kmer_index = model["kmer_index"]
    taxonomy = model["taxonomy"]  # lista, no dict
    kmer_lca  = model.get("kmer_lca")  # int_key → lca_str; None en índices viejos

    votes = {}
    # Acumulador LCA: int_key → lca_str, activo solo si el índice lo soporta
    phylum_lca_votes = Counter() if kmer_lca is not None else None

    # ============================================================
    # 1. recolectar votos (taxid) con k=k
    # ============================================================

    for kmer in generate_kmers(seq, k):

        key = encode_kmer(kmer)
        if key is None:
            continue

        # Votos taxid — solo si el k-mer está en el índice filtrado por df
        posting = kmer_index.get(key)
        if posting is not None:
            n_taxa_in_posting = len(posting)
            for taxid, weight in posting.items():
                votes[taxid] = votes.get(taxid, 0) + weight / n_taxa_in_posting

    # ============================================================
    # 1b. recolectar votos phylum LCA con k=k_lca (puede diferir de k)
    # ============================================================

    if phylum_lca_votes is not None:
        for kmer in generate_kmers(seq, k_lca):
            key = encode_kmer(kmer)
            if key is None:
                continue
            lca_str = kmer_lca.get(key)
            if lca_str is not None:
                parts = [p.strip().replace('"', '') for p in lca_str.split(";") if p.strip()]
                if len(parts) >= 2:
                    phylum_lca_votes[";".join(parts[:2]) + ";"] += 1

    # Sin taxid votes ni phylum LCA votes → nada que clasificar
    if not votes and not phylum_lca_votes:
        return seq_id, "Unclassified"

    # Si kmer_lca tiene señal clara de phylum pero kmer_index no tiene hits
    # (kmer_index fue filtrado por df_cap y no hay taxid votes), clasificar
    # directamente por phylum LCA sin pasar por la resolución a genus.
    if not votes and phylum_lca_votes:
        ph_total = sum(phylum_lca_votes.values())
        top_phylum, top_ph_count = phylum_lca_votes.most_common(1)[0]
        top_ph_fraction = top_ph_count / ph_total
        if top_ph_fraction >= PHYLUM_MIN:
            ph_parts = [p for p in top_phylum.split(";") if p]
            if len(ph_parts) > 1:
                return seq_id, top_phylum
        return seq_id, "Unclassified"

    total_vote_sum = sum(votes.values())

    # ============================================================
    # 2. filtrar ruido
    # ============================================================

    votes = {
        taxid: v for taxid, v in votes.items()
        if v / total_vote_sum >= 0.005
    }

    if not votes:
        return seq_id, "Unclassified"

    total_vote_sum = sum(votes.values())
    sorted_votes = sorted(votes.items(), key=lambda x: x[1], reverse=True)

    # ============================================================
    # 3. Phylum vote como árbitro — resuelve la ambigüedad entre
    #    phyla antes de intentar genus-level classification.
    #
    #    Problema anterior: el voto genus cruzaba phyla; Lachnospiraceae
    #    (Firmicutes) ganaba sobre Haemophilus (Proteobacteria) porque
    #    SILVA tiene más secuencias de ese género, aunque el voto de
    #    phylum fuera casi 50/50.
    #
    #    Nuevo flujo:
    #      a) Si un phylum domina (>= PHYLUM_CONFIDENCE):
    #         filtrar votos a ese phylum y resolver genus dentro de él.
    #      b) Si hay un ganador parcial (>= PHYLUM_MIN):
    #         retornar a nivel phylum (conservador pero correcto).
    #      c) Ambiguo (<PHYLUM_MIN): LCA sobre los top hits.
    # ============================================================

    # level_idx: índice en el array de partes de la taxonomía SILVA
    # Formato: Bacteria;Phylum;Class;Order;Family;Genus;Species
    #   1 = Phylum (parts[0:2])
    #   2 = Class  (parts[0:3])
    #   3 = Order  (parts[0:4])
    #   4 = Family (parts[0:5])
    #   5 = Genus  (parts[0:6])

    # ── 3a. agregar votos a nivel phylum ──────────────────────────
    # Estrategia híbrida (v6):
    #   1. Siempre computar taxid-based phylum votes desde kmer_index hits.
    #   2. Si los taxid votes tienen señal clara (top phylum >= PHYLUM_CONFIDENCE),
    #      usarlos como árbitro (comportamiento v3 — preciso para reads con buena
    #      cobertura en kmer_index, ej. faringe).
    #   3. Si la señal taxid es débil o ausente, caer a kmer_lca (cubre reads de
    #      gut cuyas especies no están bien representadas en kmer_index).
    # Esto preserva RMSE~3pp en faringe y mejora gut sin sesgar faringe.
    taxid_phylum_votes = Counter()
    for taxid, v in votes.items():
        if taxid >= len(taxonomy):
            continue
        lineage = taxonomy[taxid]
        if not lineage:
            continue
        parts = [p.strip().replace('"', '') for p in lineage.split(";") if p]
        if len(parts) >= 2:
            ph_key = ";".join(parts[:2]) + ";"
            if ph_key and ph_key != ";":
                taxid_phylum_votes[ph_key] += v

    if taxid_phylum_votes:
        tx_total = sum(taxid_phylum_votes.values())
        _, top_tx_count = taxid_phylum_votes.most_common(1)[0]
        if top_tx_count / tx_total >= PHYLUM_MIN or not phylum_lca_votes:
            # Top phylum tiene >= PHYLUM_MIN (25%) de los votos taxid → señal suficiente
            # para la lógica de PHYLUM_MIN/PHYLUM_CONFIDENCE downstream (v3 behavior).
            # Esto preserva el caso 40-55%: PHYLUM_CONFIDENCE no se cumple pero
            # PHYLUM_MIN sí, y el engine retorna phylum-level correctamente.
            phylum_votes = taxid_phylum_votes
        else:
            # Señal taxid muy débil (< 25% para cualquier phylum) + kmer_lca disponible
            # → usar kmer_lca para reads sin cobertura en kmer_index (ej. gut Bacteroidetes)
            phylum_votes = phylum_lca_votes
    elif phylum_lca_votes:
        phylum_votes = phylum_lca_votes
    else:
        phylum_votes = Counter()

    if not phylum_votes:
        return seq_id, "Unclassified"

    ph_total = sum(phylum_votes.values())
    top_phylum, top_ph_count = phylum_votes.most_common(1)[0]
    top_ph_fraction = top_ph_count / ph_total

    # ── 3b. phylum dominante → resolución jerárquica restringida ──
    if top_ph_fraction >= PHYLUM_CONFIDENCE:

        # Filtrar votos al phylum dominante
        dominant_votes = {}
        for taxid, v in votes.items():
            if taxid >= len(taxonomy):
                continue
            lineage = taxonomy[taxid]
            if not lineage:
                continue
            parts = [p.strip().replace('"', '') for p in lineage.split(";") if p]
            if len(parts) >= 2:
                ph_key = ";".join(parts[:2]) + ";"
                if ph_key == top_phylum:
                    dominant_votes[taxid] = v

        if dominant_votes:
            best_tax = None
            best_conf = 0.0

            for level_idx in [5, 4, 3, 2]:
                level_votes = Counter()
                for taxid, count in dominant_votes.items():
                    lineage = taxonomy[taxid]
                    if not lineage:
                        continue
                    parts = [p.strip().replace('"', '') for p in lineage.split(";") if p]
                    if len(parts) > level_idx:
                        key = ";".join(parts[:level_idx + 1]).strip() + ";"
                        if key and key != ";":
                            level_votes[key] += count

                if not level_votes:
                    continue

                lv_total = sum(level_votes.values())
                top_tax, top_count = level_votes.most_common(1)[0]
                confidence = top_count / lv_total

                if confidence > best_conf:
                    best_conf = confidence
                    best_tax = top_tax

                if confidence >= CONFIDENCE_THRESHOLD:
                    parts = [p for p in top_tax.split(";") if p]
                    if len(parts) > 1:
                        return seq_id, top_tax

            # fallback dentro del phylum dominante
            if best_tax and best_conf >= 0.20:
                parts = [p for p in best_tax.split(";") if p]
                if len(parts) > 1:
                    return seq_id, best_tax

        # no se pudo resolver a genus/family — retornar phylum
        ph_parts = [p for p in top_phylum.split(";") if p]
        if len(ph_parts) > 1:
            return seq_id, top_phylum

    # ── 3c. phylum parcialmente dominante → retornar a nivel phylum
    if top_ph_fraction >= PHYLUM_MIN:
        ph_parts = [p for p in top_phylum.split(";") if p]
        if len(ph_parts) > 1:
            return seq_id, top_phylum

    # ── 3d. phylum ambiguo → LCA ───────────────────────────────────

    # ============================================================
    # 4. fallback LCA para secuencias genuinamente ambiguas
    # ============================================================

    top_taxa = [
        taxonomy[taxid] if taxid < len(taxonomy) else ""
        for taxid, count in sorted_votes[:MAX_LCA_HITS]
        if count / total_vote_sum >= LCA_MIN_FRACTION
    ]

    top_taxa = [t for t in top_taxa if t and t.strip() and t.strip() != ";"]

    if not top_taxa:
        return seq_id, "Unclassified"

    tax = lca_taxa(top_taxa)

    parts = [p for p in tax.split(";") if p]

    if len(parts) <= 1:
        return seq_id, "Unclassified"

    if not tax or tax.strip() == "" or tax.strip() == ";":
        return seq_id, "Unclassified"

    return seq_id, tax
