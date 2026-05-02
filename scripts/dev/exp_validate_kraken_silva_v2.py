"""
exp_validate_kraken_silva_v2.py
================================
Valida el nuevo índice kraken-lite SILVA v2 (83k taxa, k=13, TF-IDF)
clasificando el dataset faringe_nariz de referencia.

Compara distribución a nivel de phylum contra los valores de referencia
Mothur (gold standard para este dataset).

Uso:
    cd ~/MetagenApp
    python3 scripts/exp_validate_kraken_silva_v2.py
"""

import pickle
from collections import Counter
from multiprocessing import cpu_count

from metagenapp.metagen_config import KRAKEN_INDEX_PATH
from metagenapp_core.models.kraken_lite import classify_kraken_parallel

# ── Dataset de referencia ─────────────────────────────────────────────────────
FASTA = (
    "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta"
)

OUTPUT = "/tmp/kraken_silva_v2_validation.taxonomy"

# Referencia Mothur (gold standard faringe humana)
REFERENCE = {
    "Firmicutes":      46.0,
    "Proteobacteria":  22.0,
    "Bacteroidetes":   15.0,
    "Actinobacteria":   9.0,
    "Fusobacteria":     5.0,
}

# Nombres alternativos SILVA 138 para los mismos phyla
SILVA_ALIASES = {
    "Bacillota":          "Firmicutes",
    "Pseudomonadota":     "Proteobacteria",
    "Bacteroidota":       "Bacteroidetes",
    "Actinomycetota":     "Actinobacteria",
    "Fusobacteriota":     "Fusobacteria",
}

THREADS = max(1, cpu_count() - 4)


def get_phylum(tax):
    if not tax or tax.strip() in ("", "Unclassified"):
        return "Unclassified"
    parts = [p.strip() for p in tax.split(";") if p.strip()]
    phylum = parts[1] if len(parts) >= 2 else parts[0] if parts else "Unclassified"
    return SILVA_ALIASES.get(phylum, phylum)


def print_distribution(counts, reference):
    total = sum(counts.values())
    if total == 0:
        print("  (sin resultados)")
        return

    all_phyla = list(reference.keys()) + [
        p for p in sorted(counts, key=counts.get, reverse=True)
        if p not in reference and p != "Unclassified" and counts[p] / total >= 0.01
    ] + ["Unclassified"]

    for phylum in all_phyla:
        n   = counts.get(phylum, 0)
        val = 100.0 * n / total
        ref = reference.get(phylum)
        if ref is not None:
            delta = val - ref
            flag  = "⚠ " if abs(delta) > 5 else "✓ "
            print(f"  {flag}{phylum:<22} {val:>5.1f}%  ref={ref:.1f}%  Δ={delta:+.1f}pp  (n={n:,})")
        elif val >= 0.5:
            print(f"     {phylum:<22} {val:>5.1f}%  (n={n:,})")


def main():
    print("=" * 60)
    print("  Validación kraken-lite SILVA v2")
    print("=" * 60)
    print(f"  Índice  : {KRAKEN_INDEX_PATH}")
    print(f"  FASTA   : {FASTA}")
    print(f"  Threads : {THREADS}")
    print("=" * 60)

    print("\n[1/3] Cargando índice kraken-lite...")
    with open(KRAKEN_INDEX_PATH, "rb") as f:
        model = pickle.load(f)
    print(f"      K-mers  : {len(model['kmer_index']):,}")
    print(f"      Taxa    : {len(model['taxonomy']):,}")
    print(f"      k       : {model['k']}")

    print(f"\n[2/3] Clasificando secuencias ({THREADS} threads)...")
    classify_kraken_parallel(FASTA, OUTPUT, model, threads=THREADS)
    print(f"      Resultados en: {OUTPUT}")

    print("\n[3/3] Analizando distribución de phyla...")
    counts = Counter()
    with open(OUTPUT) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            phylum = get_phylum(parts[1])
            counts[phylum] += 1

    total = sum(counts.values())
    unclassified = counts.get("Unclassified", 0)

    print(f"\n── Distribución phylum (n={total:,}) ──")
    print_distribution(counts, REFERENCE)

    print("\n── Resumen ──")
    print(f"  Total clasificadas  : {total - unclassified:,}  ({100*(total-unclassified)/total:.1f}%)")
    print(f"  Unclassified        : {unclassified:,}  ({100*unclassified/total:.1f}%)")

    bact = counts.get("Bacteroidetes", 0)
    bact_pct = 100.0 * bact / total
    ref_bact  = REFERENCE["Bacteroidetes"]
    delta     = bact_pct - ref_bact

    print(f"\n── Diagnóstico Bacteroidetes ──")
    print(f"  Resultado  : {bact_pct:.1f}%")
    print(f"  Referencia : {ref_bact:.1f}%")
    print(f"  Δ          : {delta:+.1f}pp")

    if abs(delta) <= 5:
        print("  → BIEN: dentro de ±5pp del gold standard.")
    elif delta > 5:
        print("  → SOBREESTIMADO: el índice aún infla Bacteroidetes.")
    else:
        print("  → SUBESTIMADO: el índice está por debajo del esperado.")

    print()


if __name__ == "__main__":
    main()
