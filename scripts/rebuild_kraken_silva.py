#!/usr/bin/env python3
"""
rebuild_kraken_silva.py
=======================
Construye el índice kraken-lite de MetagenApp desde el trainset SILVA completo
(Bacteria + Archaea, ~380,000 secuencias).

Requiere haber corrido primero:
  python3 scripts/build_silva_trainset.py

Salida:
  /data/databases/metagenapp_refs/16S/kraken_index_silva_idf_v2.pkl

Uso:
  cd ~/MetagenApp
  python3 scripts/rebuild_kraken_silva.py

Nota de memoria:
  Con ~380k secuencias y k=13 este proceso puede usar 20-60 GB RAM.
  Monitorear con: watch -n5 free -h
"""

import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import (
    build_kraken_index,
    invert_kmer_counts_by_taxon,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
FASTA  = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX    = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
OUTPUT = "/data/databases/metagenapp_refs/16S/kraken_index_silva_idf_v2.pkl"

K                = 13       # igual que kraken_index_v4
MAX_TAXA_PER_KMER = 100_000  # igual que rebuild_kraken_index_v4.py

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
print(f"  FASTA  : {FASTA}")
print(f"  TAX    : {TAX}")
print(f"  OUTPUT : {OUTPUT}")
print(f"  k      : {K}")
print("=" * 60)

t0 = time.time()

# ── 1. Kmer postings ──────────────────────────────────────────────────────────
print("\n[1/4] Training raw kmer postings...")
raw_model = train_raw_kmer_postings(FASTA, TAX, k=K)
print(f"      Taxa indexadas : {len(raw_model['taxonomy_labels']):,}")
print(f"      Tiempo         : {time.time() - t0:.1f}s")

# ── 2. Invertir postings ──────────────────────────────────────────────────────
t1 = time.time()
print("\n[2/4] Inverting kmer postings (kmer → taxon → weight)...")
raw_model["kmer_counts"] = invert_kmer_counts_by_taxon(
    raw_model["kmer_counts"],
    raw_model["taxonomy_labels"],
)
print(f"      Tiempo: {time.time() - t1:.1f}s")

# ── 3. Construir índice kraken-lite ───────────────────────────────────────────
t2 = time.time()
print(f"\n[3/4] Building kraken-lite index (max_taxa_per_kmer={MAX_TAXA_PER_KMER:,})...")
index = build_kraken_index(raw_model, max_taxa_per_kmer=MAX_TAXA_PER_KMER)
print(f"      K-mers en índice : {len(index['kmer_index']):,}")
print(f"      Taxa en índice   : {len(index['taxonomy']):,}")
print(f"      Tiempo           : {time.time() - t2:.1f}s")

# ── 4. Guardar ────────────────────────────────────────────────────────────────
t3 = time.time()
print(f"\n[4/4] Saving index → {OUTPUT}")
with open(OUTPUT, "wb") as f:
    pickle.dump(index, f)

size_mb = Path(OUTPUT).stat().st_size / 1e6
print(f"      Tamaño archivo : {size_mb:.0f} MB")
print(f"      Tiempo         : {time.time() - t3:.1f}s")

# ── Resumen ───────────────────────────────────────────────────────────────────
total = time.time() - t0
print("\n" + "=" * 60)
print(f"  ✅ Índice SILVA listo en {total/60:.1f} min")
print(f"  {OUTPUT}")
print("=" * 60)
print()
print("Para usar este índice en MetagenApp, actualizar en metagen_config.py:")
print(f'  KRAKEN_INDEX = "{OUTPUT}"')
