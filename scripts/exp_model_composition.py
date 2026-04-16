"""
Experimento: inspeccionar la composición del modelo v4.
Cuenta cuántos taxa de referencia tiene cada filo en el índice.
Un filo sobrerepresentado en el modelo tiende a "ganar" más clasificaciones.

Uso:
    python scripts/exp_model_composition.py
"""

import pickle
from collections import Counter
from metagenapp.metagen_config import NAIVE_MODEL_PATH


def get_phylum(tax):
    parts = [p for p in tax.split(";") if p.strip()]
    return parts[1] if len(parts) >= 2 else parts[0] if parts else "Unknown"


def get_domain(tax):
    parts = [p for p in tax.split(";") if p.strip()]
    return parts[0] if parts else "Unknown"


def main():
    print(f"Cargando modelo {NAIVE_MODEL_PATH} ...", flush=True)
    with open(NAIVE_MODEL_PATH, "rb") as f:
        model = pickle.load(f)

    print(f"\nKeys del modelo: {list(model.keys())}\n")

    taxa = model.get("taxa", [])
    k    = model.get("k")
    print(f"k-mer size: {k}")
    print(f"Total taxa de referencia: {len(taxa):,}")

    # Composición por filo
    phylum_counts = Counter(get_phylum(t) for t in taxa if get_domain(t) == "Bacteria")
    domain_counts = Counter(get_domain(t) for t in taxa)

    print(f"\nDominios en el modelo:")
    for d, c in domain_counts.most_common():
        print(f"  {d:<25} {c:>6,} taxa")

    print(f"\nTop-15 filos en el modelo (taxa de referencia):")
    total_bacteria = sum(phylum_counts.values())
    for phylum, count in phylum_counts.most_common(15):
        pct = 100.0 * count / total_bacteria
        print(f"  {phylum:<30} {count:>6,}  ({pct:.1f}%)")

    # Comparar con lo esperado en faringe
    EXPECTED_FARINGE = {
        "Firmicutes":     46.0,
        "Proteobacteria": 22.0,
        "Bacteroidetes":  15.0,
        "Actinobacteria":  9.0,
        "Fusobacteria":    5.0,
    }

    print("\n── Sobrerepresentación en modelo vs abundancia esperada en faringe ──")
    print(f"  {'Filo':<25} {'% en modelo':>12}  {'% esperado':>12}  {'Ratio':>8}")
    print("  " + "-" * 62)
    for phylum, expected in sorted(EXPECTED_FARINGE.items(), key=lambda x: -x[1]):
        model_pct = 100.0 * phylum_counts.get(phylum, 0) / total_bacteria
        ratio     = model_pct / expected if expected > 0 else float("inf")
        flag      = "⚠️ " if ratio > 2.0 or ratio < 0.5 else "✓ "
        print(f"  {flag}{phylum:<23} {model_pct:>11.1f}%  {expected:>11.1f}%  {ratio:>7.2f}x")

    # Estadísticas del índice invertido
    postings_taxon = model.get("postings_taxon")
    kmer_to_id     = model.get("kmer_to_id")

    if postings_taxon is not None and kmer_to_id is not None:
        psize_dist = Counter()
        for kid in range(len(postings_taxon)):
            psize = len(postings_taxon[kid])
            if psize <= 5:
                psize_dist["1-5"] += 1
            elif psize <= 10:
                psize_dist["6-10"] += 1
            elif psize <= 20:
                psize_dist["11-20"] += 1
            elif psize <= 50:
                psize_dist["21-50"] += 1
            else:
                psize_dist[">50"] += 1

        total_kmers = len(postings_taxon)
        print(f"\n── Distribución de psize en el índice ({total_kmers:,} k-mers únicos) ──")
        for bucket in ["1-5", "6-10", "11-20", "21-50", ">50"]:
            c   = psize_dist.get(bucket, 0)
            pct = 100.0 * c / total_kmers
            bar = "█" * int(pct / 2)
            print(f"  psize {bucket:>6}  {c:>8,}  ({pct:>5.1f}%)  {bar}")


if __name__ == "__main__":
    main()
