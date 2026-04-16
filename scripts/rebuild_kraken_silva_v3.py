#!/usr/bin/env python3
"""
rebuild_kraken_silva_v3.py
==========================
Construye el índice kraken-lite SILVA v3 con trainset balanceado por genus.

Diferencia respecto a v2:
  max_seqs_per_genus=50 → limita a 50 secuencias por genus en el trainset,
  evitando que genera sobre-representados en SILVA (Streptococcus, Lachnospiraceae,
  Veillonella) dominen el índice y sesgen la clasificación hacia Firmicutes.

Salida:
  /data/databases/metagenapp_refs/16S/kraken_index_silva_v3.pkl

Uso:
  cd ~/MetagenApp
  python3 scripts/rebuild_kraken_silva_v3.py
"""

import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA  = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX    = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
OUTPUT = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v3.pkl"

K                  = 13
TRAIN_CAP          = 10_000   # cap por k-mer durante entrenamiento
MAX_TAXA_PER_KMER  = 50_000   # cap en índice final
MAX_SEQS_PER_GENUS = 150      # cap de balanceo por genus ← clave de v3 (barrido: 150 minimiza RMSE)

# ── Validación ────────────────────────────────────────────────────────────────
for path in [FASTA, TAX]:
    if not Path(path).exists():
        raise FileNotFoundError(
            f"\nNo se encontró: {path}\n"
            "Corre primero: python3 scripts/build_silva_trainset.py"
        )

Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)

# ── Pipeline ──────────────────────────────────────────────────────────────────
print("=" * 60)
print("  rebuild_kraken_silva_v3.py — MetagenApp")
print("=" * 60)
print(f"  FASTA              : {FASTA}")
print(f"  TAX                : {TAX}")
print(f"  OUTPUT             : {OUTPUT}")
print(f"  k                  : {K}")
print(f"  train_cap          : {TRAIN_CAP:,}")
print(f"  max_taxa_per_kmer  : {MAX_TAXA_PER_KMER:,}")
print(f"  max_seqs_per_genus : {MAX_SEQS_PER_GENUS}")
print("=" * 60)

t0 = time.time()

# ── 1. Kmer postings con cap por genus ────────────────────────────────────────
print("\n[1/3] Training raw kmer postings (genus-balanced)...")
raw_model = train_raw_kmer_postings(
    FASTA, TAX,
    k=K,
    max_taxa_per_kmer=TRAIN_CAP,
    max_seqs_per_genus=MAX_SEQS_PER_GENUS,
)
n_kmers = len(raw_model["kmer_counts"])
n_taxa  = len(raw_model["taxonomy_labels"])
print(f"      Taxa indexadas   : {n_taxa:,}")
print(f"      K-mers únicos    : {n_kmers:,}")
print(f"      Tiempo           : {time.time() - t0:.1f}s")

# ── 2. Construir índice kraken-lite con TF-IDF ────────────────────────────────
t2 = time.time()
print(f"\n[2/3] Building kraken-lite index (max_taxa_per_kmer={MAX_TAXA_PER_KMER:,})...")
index = build_kraken_index(raw_model, max_taxa_per_kmer=MAX_TAXA_PER_KMER)
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
print("\n" + "=" * 60)
print(f"  Índice SILVA v3 listo en {total/60:.1f} min")
print(f"  {OUTPUT}")
print("=" * 60)
print()
print("Para usar este índice en MetagenApp, actualizar en metagen_config.py:")
print(f'  KRAKEN_INDEX_PATH = REF_16S_ROOT / "kraken_index_silva_v3.pkl"')
