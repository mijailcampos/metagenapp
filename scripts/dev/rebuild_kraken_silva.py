#!/usr/bin/env python3
"""
rebuild_kraken_silva.py
=======================
Construye el índice kraken-lite de MetagenApp desde el trainset SILVA completo
(Bacteria + Archaea, ~451,000 secuencias, 83,759 taxa a nivel de especie).

Requiere haber corrido primero:
  python3 scripts/build_silva_trainset.py

Salida:
  /data/databases/metagenapp_refs/16S/kraken_index_silva_v2.pkl

Uso:
  cd ~/MetagenApp
  python3 scripts/rebuild_kraken_silva.py

Nota de memoria:
  Con 256 GB RAM y ~451k secuencias, k=13, max_taxa_per_kmer=10_000,
  el proceso usa ~60-120 GB RAM en pico.
  Monitorear con: watch -n5 free -h
"""

import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index

# ── Paths ──────────────────────────────────────────────────────────────────────
FASTA  = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX    = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
OUTPUT = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v2.pkl"

# K más largo → mejor discriminación a nivel de género
# Con 256 GB RAM, k=13 con max_taxa_per_kmer=10_000 es manejable
K                = 13
# Cap por kmer durante el entrenamiento (limita RAM pico)
TRAIN_CAP        = 10_000
# Cap durante la construcción del índice (solo excluye los más ubicuos)
MAX_TAXA_PER_KMER = 50_000

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
print("  rebuild_kraken_silva.py — MetagenApp")
print("=" * 60)
print(f"  FASTA            : {FASTA}")
print(f"  TAX              : {TAX}")
print(f"  OUTPUT           : {OUTPUT}")
print(f"  k                : {K}")
print(f"  train_cap        : {TRAIN_CAP:,}")
print(f"  max_taxa_per_kmer: {MAX_TAXA_PER_KMER:,}")
print("=" * 60)

t0 = time.time()

# ── 1. Kmer postings ──────────────────────────────────────────────────────────
# train_raw_kmer_postings ya entrega el índice invertido:
#   kmer_counts: dict[kmer_str] -> dict[taxid] -> count
# NO hay que invertir de nuevo.
print("\n[1/3] Training raw kmer postings...")
raw_model = train_raw_kmer_postings(FASTA, TAX, k=K, max_taxa_per_kmer=TRAIN_CAP)
n_kmers = len(raw_model["kmer_counts"])
n_taxa  = len(raw_model["taxonomy_labels"])
print(f"      Taxa indexadas   : {n_taxa:,}")
print(f"      K-mers únicos    : {n_kmers:,}")
print(f"      Tiempo           : {time.time() - t0:.1f}s")

# ── 2. Construir índice kraken-lite con TF-IDF ────────────────────────────────
# Las claves se codifican como enteros 2-bit en build_kraken_index
# para coincidir con encode_kmer() del engine.
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
print(f"  Índice SILVA v2 listo en {total/60:.1f} min")
print(f"  {OUTPUT}")
print("=" * 60)
print()
print("Para usar este índice en MetagenApp, actualizar en metagen_config.py:")
print(f'  KRAKEN_INDEX_PATH = REF_16S_ROOT / "kraken_index_silva_v2.pkl"')
