"""
Experimento v2: clasificar reads QC-filtrados de Mothur (SIN alineamiento) con naive-v2.
Usa stability.trim.contigs.good.fasta — reads limpios, sin gaps de alineamiento.

Si Bacteroidetes baja a ~15% → el sesgo viene del QC diferencial de MetagenApp.
Si Bacteroidetes se mantiene en ~28% → el sesgo viene del modelo v4.

Uso:
    python scripts/exp_classify_mothur_reads.py
"""

import pickle
import random
from collections import Counter
from Bio import SeqIO
from metagenapp_core.models.naive_v2_engine import classify_seq_v2
from metagenapp.metagen_config import NAIVE_MODEL_PATH

# Reads QC-filtrados de Mothur, SIN alineamiento (1,286,486 reads)
MOTHUR_FASTA_RAW = (
    "/data/results/paper01_faringe_mothur/"
    "stability.trim.contigs.good.fasta"
)

# Reads QC-filtrados de MetagenApp (408,889 reads)
METAGEN_FASTA = "/data/results/paper01_faringe_naive_v2/final_clean.fasta"

SAMPLE_SIZE = 50_000  # muestra representativa para velocidad

REFERENCE = {
    "Firmicutes":     46.0,
    "Proteobacteria": 22.0,
    "Bacteroidetes":  15.0,
    "Actinobacteria":  9.0,
    "Fusobacteria":    5.0,
}


def get_phylum(tax):
    if tax is None:
        return "Unclassified"
    parts = [p for p in tax.split(";") if p.strip()]
    return parts[1] if len(parts) >= 2 else parts[0] if parts else "Unclassified"


def reservoir_sample(fasta_path, n):
    """Muestreo reservoir — no carga todo en RAM."""
    reservoir = []
    for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
        seq = str(record.seq).upper()
        if len(reservoir) < n:
            reservoir.append((record.id, seq))
        else:
            j = random.randint(0, i)
            if j < n:
                reservoir[j] = (record.id, seq)
    return reservoir


def classify_records(records, model, label):
    counts = Counter()
    for i, (sid, seq) in enumerate(records):
        _, tax = classify_seq_v2(sid, seq, model)
        counts[get_phylum(tax)] += 1
        if (i + 1) % 10_000 == 0:
            print(f"  [{label}] {i+1:,} / {len(records):,} ...", flush=True)
    print(f"  [{label}] listo: {len(records):,} secuencias", flush=True)
    return counts


def print_distribution(label, counts, reference):
    total = sum(counts.values())
    print(f"\n── {label} (n={total:,}) ──")
    for phylum in list(reference.keys()) + ["Unclassified"]:
        val = 100.0 * counts.get(phylum, 0) / total
        ref = reference.get(phylum)
        if ref is not None:
            delta = val - ref
            flag  = "⚠️ " if abs(delta) > 5 else "✓ "
            print(f"  {flag}{phylum:<20} {val:>5.1f}%   ref={ref:.1f}%  Δ={delta:+.1f}pp")
        else:
            print(f"     {phylum:<20} {val:>5.1f}%")


def main():
    random.seed(42)

    print(f"Cargando modelo {NAIVE_MODEL_PATH} ...", flush=True)
    with open(NAIVE_MODEL_PATH, "rb") as f:
        model = pickle.load(f)

    print(f"\n[1/2] Muestreando {SAMPLE_SIZE:,} reads de Mothur (sin alineamiento) ...", flush=True)
    mothur_records = reservoir_sample(MOTHUR_FASTA_RAW, SAMPLE_SIZE)
    print(f"  Muestra obtenida: {len(mothur_records):,} reads", flush=True)
    counts_mothur = classify_records(mothur_records, model, "Mothur-raw-QC")

    print(f"\n[2/2] Muestreando {SAMPLE_SIZE:,} reads de MetagenApp (final_clean) ...", flush=True)
    metagen_records = reservoir_sample(METAGEN_FASTA, SAMPLE_SIZE)
    print(f"  Muestra obtenida: {len(metagen_records):,} reads", flush=True)
    counts_metagen = classify_records(metagen_records, model, "MetagenApp-QC")

    print_distribution("naive-v2 sobre input Mothur (QC crudo, sin alinear)", counts_mothur, REFERENCE)
    print_distribution("naive-v2 sobre input MetagenApp (final_clean)",       counts_metagen, REFERENCE)

    bact_mothur  = 100.0 * counts_mothur.get("Bacteroidetes", 0)  / sum(counts_mothur.values())
    bact_metagen = 100.0 * counts_metagen.get("Bacteroidetes", 0) / sum(counts_metagen.values())
    diff = bact_metagen - bact_mothur

    print("\n── Diagnóstico ──")
    print(f"  Bacteroidetes con input Mothur:     {bact_mothur:.1f}%")
    print(f"  Bacteroidetes con input MetagenApp: {bact_metagen:.1f}%")
    print(f"  Diferencia atribuible al QC:        {diff:+.1f}pp")

    if bact_mothur < 20.0:
        print("\n→ CONCLUSIÓN: El sesgo viene principalmente del QC diferencial.")
        print("  El motor naive-v2 es correcto. Hay que alinear el preprocesamiento.")
    elif diff < 3.0:
        print("\n→ CONCLUSIÓN: El sesgo viene principalmente del modelo v4.")
        print("  El QC no explica la diferencia. Hay que recalibrar el índice de referencia.")
    else:
        print("\n→ CONCLUSIÓN: Sesgo mixto (QC + modelo). Ambos contribuyen.")
        print(f"  QC aporta ~{diff:.1f}pp, modelo aporta el resto.")


if __name__ == "__main__":
    main()
