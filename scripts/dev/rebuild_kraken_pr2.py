#!/usr/bin/env python3
"""
rebuild_kraken_pr2.py
=====================
Construye el índice kraken-lite de MetagenApp desde PR2 v5 (18S, eucariotas).

Requiere:
  - /data/databases/metagenapp_refs/18S/pr2_reference.fasta
  - /data/databases/PR2_v5.1.1_mothur/pr2_v5.1.1.tax

Salida:
  /data/databases/metagenapp_refs/18S/kraken_index_pr2.pkl

Uso:
  cd ~/MetagenApp
  python3 scripts/rebuild_kraken_pr2.py

Nota de memoria:
  PR2 v5 tiene ~240k secuencias. Estima 10-30 GB RAM con k=13.
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
FASTA  = "/data/databases/metagenapp_refs/18S/pr2_reference.fasta"
TAX    = "/data/databases/PR2_v5.1.1_mothur/pr2_v5.1.1.tax"
OUTPUT = "/data/databases/metagenapp_refs/18S/kraken_index_pr2.pkl"

K                 = 13
MAX_TAXA_PER_KMER = 100_000

# ── Validación ────────────────────────────────────────────────────────────────
for path in [FASTA, TAX]:
    if not Path(path).exists():
        raise FileNotFoundError(f"\nNo se encontró: {path}")

Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)

# ── Pipeline ──────────────────────────────────────────────────────────────────
print("=" * 60)
print("  rebuild_kraken_pr2.py — MetagenApp 18S")
print("=" * 60)
print(f"  FASTA  : {FASTA}")
print(f"  TAX    : {TAX}")
print(f"  OUTPUT : {OUTPUT}")
print(f"  k      : {K}")
print("=" * 60)

t0 = time.time()

# ── 1. Kmer postings ──────────────────────────────────────────────────────────
print("\n[1/4] Training raw kmer postings...")
raw_model = train_raw_kmer_postings(FASTA, TAX, k=K, max_taxa_per_kmer=200)
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
print(f"  ✅ Índice PR2 listo en {total/60:.1f} min")
print(f"  {OUTPUT}")
print("=" * 60)
