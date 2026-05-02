#!/usr/bin/env python3
"""
rebuild_kraken_silva_v4.py
==========================
Construye el índice kraken-lite SILVA v4 con phylum-level IDF.

Diferencia respecto a v3:
  v3: TF × taxa_IDF  (IDF basado en número de taxa con el k-mer)
  v4: TF × taxa_IDF × phylum_IDF  (IDF adicional basado en número de
      phyla distintos que contienen el k-mer)

El problema de v3: k-mers compartidos entre Firmicutes y Bacteroidetes
  reciben IDF bajo por alta df (taxa), pero el mayor número de taxa
  de Firmicutes hace que su voto sea dominante incluso con trainset
  balanceado. Resultado: Bacteroidetes subestimado en gut (1.8% vs 25%).

La corrección: phylum_IDF = log(N_phyla / n_phyla_with_kmer) + 1
  - K-mer en Firmicutes + Bacteroidetes → n_phyla=2 → phylum_IDF bajo
  - K-mer exclusivo de Bacteroidetes   → n_phyla=1 → phylum_IDF alto
  - K-mer en todos los phyla            → n_phyla=N → phylum_IDF ≈ 1.0

Trainset: IDÉNTICO a v3 (max_seqs_per_genus=150, sin phylum cap).
  El phylum cap NO ayuda porque la diversidad de géneros por phylum
  sigue desbalanceada incluso con equal phylum seq counts.

Salida:
  /data/databases/metagenapp_refs/16S/kraken_index_silva_v4.pkl

Uso:
  cd ~/MetagenApp
  python3 scripts/rebuild_kraken_silva_v4.py
"""

import argparse
import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA  = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX    = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
OUTPUT = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v4.pkl"

# ── Parámetros ────────────────────────────────────────────────────────────────
K                   = 13
TRAIN_CAP           = 10_000
MAX_TAXA_PER_KMER   = 10_000    # ajustable con --df-cap; equiv. a psize_max de naive-v2
MAX_SEQS_PER_GENUS  = 150
MAX_SEQS_PER_PHYLUM = None

parser = argparse.ArgumentParser(description="Rebuild kraken-lite SILVA v4")
parser.add_argument("--df-cap", type=int, default=MAX_TAXA_PER_KMER,
    help=f"max_taxa_per_kmer en index build (default: {MAX_TAXA_PER_KMER})")
parser.add_argument("--output", type=str, default=OUTPUT,
    help=f"Ruta de salida (default: {OUTPUT})")
args = parser.parse_args()
MAX_TAXA_PER_KMER = args.df_cap
OUTPUT = args.output

# ── Validación ────────────────────────────────────────────────────────────────
for path in [FASTA, TAX]:
    if not Path(path).exists():
        raise FileNotFoundError(
            f"\nNo se encontró: {path}\n"
            "Corre primero: python3 scripts/build_silva_trainset.py"
        )

Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)

# ── Pipeline ──────────────────────────────────────────────────────────────────
print("=" * 65)
print("  rebuild_kraken_silva_v4.py — MetagenApp")
print("=" * 65)
print(f"  FASTA                : {FASTA}")
print(f"  TAX                  : {TAX}")
print(f"  OUTPUT               : {OUTPUT}")
print(f"  k                    : {K}")
print(f"  train_cap            : {TRAIN_CAP:,}")
print(f"  max_taxa_per_kmer    : {MAX_TAXA_PER_KMER:,}")
print(f"  max_seqs_per_genus   : {MAX_SEQS_PER_GENUS}")
print(f"  max_seqs_per_phylum  : {MAX_SEQS_PER_PHYLUM}  (sin cap — igual que v3)")
print(f"  use_phylum_idf       : True  ← nuevo en v4")
print("=" * 65)

t0 = time.time()

# ── 1. Kmer postings (mismo trainset que v3) ──────────────────────────────────
print("\n[1/3] Training raw kmer postings (genus-balanced, same as v3)...")
raw_model = train_raw_kmer_postings(
    FASTA, TAX,
    k=K,
    max_taxa_per_kmer=TRAIN_CAP,
    max_seqs_per_genus=MAX_SEQS_PER_GENUS,
    max_seqs_per_phylum=MAX_SEQS_PER_PHYLUM,
)
n_kmers = len(raw_model["kmer_counts"])
n_taxa  = len(raw_model["taxonomy_labels"])
print(f"      Taxa indexadas   : {n_taxa:,}")
print(f"      K-mers únicos    : {n_kmers:,}")
print(f"      Tiempo           : {time.time() - t0:.1f}s")

# ── 2. Construir índice con TF × taxa_IDF × phylum_IDF ───────────────────────
t2 = time.time()
print(f"\n[2/3] Building kraken-lite index (phylum_idf=True, max_taxa_per_kmer={MAX_TAXA_PER_KMER:,})...")
index = build_kraken_index(raw_model, max_taxa_per_kmer=MAX_TAXA_PER_KMER, use_phylum_idf=True)
print(f"      K-mers en índice : {len(index['kmer_index']):,}")
print(f"      Taxa en índice   : {len(index['taxonomy']):,}")
print(f"      Tiempo           : {time.time() - t2:.1f}s")

# ── 3. Guardar ────────────────────────────────────────────────────────────────
t3 = time.time()
print(f"\n[3/3] Saving index → {OUTPUT}")
with open(OUTPUT, "wb") as f:
    pickle.dump(index, f)

size_mb = Path(OUTPUT).stat().st_size / 1e6
print(f"      Tamaño archivo : {size_mb:.0f} MB")
print(f"      Tiempo         : {time.time() - t3:.1f}s")

# ── Resumen ───────────────────────────────────────────────────────────────────
total = time.time() - t0
print("\n" + "=" * 65)
print(f"  Índice SILVA v4 listo en {total/60:.1f} min")
print(f"  Trainset: genus_cap={MAX_SEQS_PER_GENUS}  phylum_idf=True")
print(f"  {OUTPUT}")
print("=" * 65)
print()
print("Para usar este índice en MetagenApp, actualizar en metagen_config.py:")
print(f'  KRAKEN_INDEX_PATH = REF_16S_ROOT / "kraken_index_silva_v4.pkl"')
print()
print("Para validar en ambos datasets:")
print("  python3 scripts/exp_validate_kraken_silva_v4.py")
