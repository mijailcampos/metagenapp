import math

DNA_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}


def _encode_kmer(seq):
    """Codifica un k-mer string como entero 2-bit (mismo esquema que el engine)."""
    val = 0
    for b in seq:
        if b not in DNA_MAP:
            return None
        val = (val << 2) | DNA_MAP[b]
    return val


def _extract_phylum_from_lineage(lineage: str) -> str:
    """Extrae el phylum (parts[1]) de un linaje SILVA semicolon-separated."""
    parts = [p.strip().strip('"') for p in lineage.split(";") if p.strip()]
    return parts[1] if len(parts) >= 2 else ""


def build_kraken_index(raw_model, max_taxa_per_kmer=100000, use_phylum_idf=True,
                       raw_model_for_lca=None):
    """
    Construye el índice kraken-lite con pesos TF × taxa_IDF × phylum_IDF.

    Recibe raw_model con kmer_counts YA INVERTIDO:
      kmer_counts: dict[kmer_str] -> dict[taxid] -> tf_count

    Las claves se codifican como enteros 2-bit para que coincidan con
    el esquema de encode_kmer() usado durante la inferencia.

    taxa_IDF  = log( N_taxa / df_kmer ) + 1
        donde df_kmer = número de taxa con ese k-mer.
        Penaliza k-mers ubicuos dentro de un mismo phylum.

    phylum_IDF = log( N_phyla / pf_kmer ) + 1   [solo si use_phylum_idf=True]
        donde pf_kmer = número de phyla distintos que contienen ese k-mer.
        Amplifica k-mers phylum-discriminativos; penaliza k-mers compartidos
        entre phyla que de otro modo sesgan la clasificación.

        Problema que resuelve: un k-mer en 500 Firmicutes + 1 Bacteroidetes
        tiene taxa_IDF bajo (alta df), y al sumar votos acumula ~500x más
        peso en Firmicutes que en Bacteroidetes. El phylum_IDF lo penaliza
        porque aparece en 2 phyla → pf=2 → phylum_IDF bajo. En cambio,
        un k-mer exclusivo de Bacteroidetes tiene pf=1 → phylum_IDF máximo.

    Efecto combinado:
      - K-mers cross-phyla (compartidos entre Firmicutes y Bacteroidetes)
        → taxa_IDF bajo Y phylum_IDF bajo → peso muy bajo
      - K-mers phylum-específicos (ej. solo Bacteroidetes)
        → phylum_IDF alto → peso amplificado aunque aparezcan en muchos taxa
          de ese phylum
      - K-mers ubicuos (en todos los phyla)
        → phylum_IDF ≈ 1.0 → sin amplificación
    """
    k           = raw_model["kmer_size"]
    taxa        = raw_model["taxonomy_labels"]
    kmer_counts = raw_model["kmer_counts"]  # dict[kmer_str] -> dict[taxid] -> tf

    N_taxa = len(taxa)

    # ── Phylum IDF: precomputar phylum de cada taxid ──────────────────────────
    if use_phylum_idf:
        taxid_to_phylum = [_extract_phylum_from_lineage(lin) for lin in taxa]
        N_phyla = max(len({p for p in taxid_to_phylum if p}), 1)
    else:
        taxid_to_phylum = None
        N_phyla = 1

    # ── Construir índice con peso TF × taxa_IDF × phylum_IDF ─────────────────
    kmer_index = {}

    for kmer_str, taxid_counts in kmer_counts.items():

        df = len(taxid_counts)  # cuántos taxa tienen este kmer

        # filtrar kmers que aparecen en demasiados taxa
        if df > max_taxa_per_kmer:
            continue

        key = _encode_kmer(kmer_str)
        if key is None:
            continue  # kmer con bases ambiguas

        # taxa-level IDF: penaliza k-mers comunes entre muchos taxa
        taxa_idf = math.log(N_taxa / df) + 1.0

        # phylum-level IDF: amplifica k-mers phylum-discriminativos
        if use_phylum_idf:
            phyla_with_kmer = {
                taxid_to_phylum[t]
                for t in taxid_counts
                if t < len(taxa) and taxid_to_phylum[t]
            }
            n_phyla_kmer = len(phyla_with_kmer) or 1
            phylum_idf = math.log(N_phyla / n_phyla_kmer) + 1.0
        else:
            phylum_idf = 1.0

        combined_idf = taxa_idf * phylum_idf

        entry = {}
        for taxid, tf in taxid_counts.items():
            entry[taxid] = tf * combined_idf

        kmer_index[key] = entry

    # ── taxid_to_phylum_key: mapa taxid → ph_key para engine (distribución igualitaria) ──
    # El engine usa esta lista para distribuir votos de k-mers igualitariamente
    # entre phyla: cada phylum presente en un posting recibe 1/n_phyla del voto,
    # sin importar cuántos taxa tenga. Esto elimina el sesgo por diversidad de taxa
    # por phylum (Firmicutes ~20K taxa >> Bacteroidetes ~5K taxa en SILVA NR99).
    from collections import Counter as _Counter
    taxid_to_phylum_key = []
    phylum_sizes: dict = _Counter()
    for lineage in taxa:
        parts = [p.strip().strip('"') for p in lineage.split(";") if p.strip()]
        ph_key = ";".join(parts[:2]) + ";" if len(parts) >= 2 else ""
        taxid_to_phylum_key.append(ph_key)
        if ph_key:
            phylum_sizes[ph_key] += 1

    # ── Compute kmer_lca for LCA-based phylum voting in engine ──────────────
    # build_kid_to_lca() calcula el LCA de cada k-mer sobre todos sus taxids.
    # K-mers compartidos entre phyla distintos tienen LCA=Bacteria (depth=1)
    # y NO votan por ningún phylum en el engine → elimina el sesgo por número
    # de taxa por phylum (Firmicutes ~20K taxa >> Bacteroidetes ~5K taxa).
    # kmer_lca: mapa int_key → lca_str para LCA-based phylum voting en el engine.
    #
    # Fuente: raw_model_for_lca (si se proporciona) o raw_model (fallback).
    #
    # Problema de cobertura: con genus_cap=150, el raw_model solo contiene 149K
    # secuencias SILVA. Para cepas de gut (ej. Prevotella copri), SILVA tiene
    # cientos de secuencias pero solo 150 entran en el trainset. Los 13-mers
    # específicos de gut Bacteroidetes pueden NO estar en el trainset → cobertura
    # ~1% para amplicons gut.
    #
    # Solución: construir kmer_lca desde un raw_model SIN genus_cap (todos los
    # 451K secuencias SILVA) → cobertura ~60% para gut → phylum LCA voting funciona.
    # El kmer_index (taxid voting para genus resolution) sigue usando raw_model con
    # genus_cap para evitar sesgos de abundancia en la resolución a nivel genus.
    _lca_source = raw_model_for_lca if raw_model_for_lca is not None else raw_model
    from metagenapp_core.models.train_kmer_postings import build_kid_to_lca as _build_lca
    kmer_lca_str_map = _build_lca(_lca_source)
    kmer_lca: dict = {}
    for kmer_str, lca_str in kmer_lca_str_map.items():
        key = _encode_kmer(kmer_str)
        if key is not None:
            kmer_lca[key] = lca_str
    del kmer_lca_str_map

    return {
        "k": k,
        "taxonomy": taxa,
        "kmer_index": kmer_index,
        "kmer_lca": kmer_lca,          # int_key → lca_lineage_str; para phylum voting
        "use_phylum_idf": use_phylum_idf,
        "taxid_to_phylum_key": taxid_to_phylum_key,
        "phylum_sizes": dict(phylum_sizes),   # diagnóstico
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