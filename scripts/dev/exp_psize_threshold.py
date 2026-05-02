"""
Experimento: efecto del umbral psize_max sobre la distribución de filos.
Clasifica los centroides OTU del run faringe con distintos valores de psize_max
y reporta la abundancia relativa de cada filo para cada umbral.

Uso:
    python scripts/exp_psize_threshold.py
"""

import pickle
import sys
from collections import Counter
from Bio import SeqIO
from metagenapp_core.models.naive_v2_engine import classify_seq_v2
from metagenapp.metagen_config import NAIVE_MODEL_PATH

CENTROIDS  = "/data/results/paper01_faringe_naive_v2/centroids_all.fasta"
THRESHOLDS = [5, 8, 10, 15, 20, 30, 50]

# Referencia Mothur/QIIME2 (media de ambos pipelines a nivel de filo)
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


def run_threshold(records, model, psize_max):
    phylum_counts = Counter()
    for sid, seq in records:
        _, tax = classify_seq_v2(sid, seq, model, psize_max=psize_max)
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

    # Filos de interés para el reporte
    top_phyla = list(REFERENCE.keys()) + ["Unclassified"]

    # Cabecera
    header = f"{'psize_max':>10}  " + "  ".join(f"{p:>15}" for p in top_phyla)
    print(header)
    print("-" * len(header))

    results = {}
    for pmax in THRESHOLDS:
        counts = run_threshold(records, model, pmax)
        pct    = phylum_table(counts)
        results[pmax] = pct

        row = f"{pmax:>10}  " + "  ".join(
            f"{pct.get(p, 0.0):>14.1f}%" for p in top_phyla
        )
        print(row, flush=True)

    # Resumen: qué umbral minimiza error vs referencia en Bacteroidetes
    print("\n── Error absoluto en Bacteroidetes vs referencia (15%) ──")
    best_pmax = None
    best_err  = float("inf")
    for pmax, pct in results.items():
        bact = pct.get("Bacteroidetes", 0.0)
        err  = abs(bact - REFERENCE["Bacteroidetes"])
        print(f"  psize_max={pmax:>3}  Bacteroidetes={bact:.1f}%  error={err:.1f}pp")
        if err < best_err:
            best_err  = err
            best_pmax = pmax

    print(f"\n→ Umbral óptimo para Bacteroidetes: psize_max={best_pmax} (error={best_err:.1f}pp)")

    # Verificar que el umbral óptimo no deteriora otros filos
    print("\n── Distribución completa para psize_max óptimo ──")
    pct_opt = results[best_pmax]
    for phylum, ref in REFERENCE.items():
        val = pct_opt.get(phylum, 0.0)
        delta = val - ref
        flag = "⚠️ " if abs(delta) > 5 else "✓ "
        print(f"  {flag}{phylum:<20} naive-v2={val:.1f}%  ref={ref:.1f}%  Δ={delta:+.1f}pp")


if __name__ == "__main__":
    main()
