"""
Experimento: efecto del umbral descent_ratio sobre la distribución de filos.
Clasifica los centroides OTU del run faringe con distintos valores de descent_ratio
y reporta la abundancia relativa de cada filo para cada umbral.

Uso:
    python scripts/exp_descent_ratio.py
"""

import pickle
from collections import Counter
from Bio import SeqIO
from metagenapp_core.models.naive_v2_engine import classify_seq_v2
from metagenapp.metagen_config import NAIVE_MODEL_PATH

CENTROIDS  = "/data/results/paper01_faringe_naive_v2/centroids_all.fasta"
RATIOS     = [0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.95, 0.99]

REFERENCE = {
    "Firmicutes":      46.0,
    "Proteobacteria":  22.0,
    "Bacteroidetes":   15.0,
    "Actinobacteria":   9.0,
    "Fusobacteria":     5.0,
}

def get_phylum(tax):
    if tax is None:
        return "Unclassified"
    parts = [p for p in tax.split(";") if p.strip()]
    return parts[1] if len(parts) >= 2 else parts[0] if parts else "Unclassified"


def run_ratio(records, model, ratio):
    phylum_counts = Counter()
    for sid, seq in records:
        _, tax = classify_seq_v2(sid, seq, model, descent_ratio=ratio)
        phylum_counts[get_phylum(tax)] += 1
    return phylum_counts


def phylum_table(counts):
    total = sum(counts.values())
    if total == 0:
        return {}
    return {k: 100.0 * v / total for k, v in counts.most_common()}


def main():
    print(f"Cargando modelo desde {NAIVE_MODEL_PATH} ...", flush=True)
    with open(NAIVE_MODEL_PATH, "rb") as f:
        model = pickle.load(f)

    print(f"Leyendo centroides desde {CENTROIDS} ...", flush=True)
    records = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(CENTROIDS, "fasta")]
    print(f"  → {len(records)} OTU centroides\n", flush=True)

    top_phyla = list(REFERENCE.keys()) + ["Unclassified"]

    header = f"{'ratio':>8}  " + "  ".join(f"{p:>15}" for p in top_phyla)
    print(header)
    print("-" * len(header))

    results = {}
    for ratio in RATIOS:
        counts = run_ratio(records, model, ratio)
        pct    = phylum_table(counts)
        results[ratio] = pct

        row = f"{ratio:>8.2f}  " + "  ".join(
            f"{pct.get(p, 0.0):>14.1f}%" for p in top_phyla
        )
        print(row, flush=True)

    print("\n── Error absoluto en Bacteroidetes vs referencia (15%) ──")
    best_ratio = None
    best_err   = float("inf")
    for ratio, pct in results.items():
        bact = pct.get("Bacteroidetes", 0.0)
        err  = abs(bact - REFERENCE["Bacteroidetes"])
        print(f"  ratio={ratio:.2f}  Bacteroidetes={bact:.1f}%  error={err:.1f}pp")
        if err < best_err:
            best_err   = err
            best_ratio = ratio

    print(f"\n→ Ratio óptimo para Bacteroidetes: {best_ratio:.2f} (error={best_err:.1f}pp)")

    print("\n── Distribución completa para ratio óptimo ──")
    pct_opt = results[best_ratio]
    for phylum, ref in REFERENCE.items():
        val   = pct_opt.get(phylum, 0.0)
        delta = val - ref
        flag  = "⚠️ " if abs(delta) > 5 else "✓ "
        print(f"  {flag}{phylum:<20} naive-v2={val:.1f}%  ref={ref:.1f}%  Δ={delta:+.1f}pp")


if __name__ == "__main__":
    main()
