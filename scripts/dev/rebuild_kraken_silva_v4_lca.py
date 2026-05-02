#!/usr/bin/env python3
"""
rebuild_kraken_silva_v4_lca.py
===============================
Reconstruye el índice kraken-lite SILVA v4 con LCA-based phylum voting.

Diferencia respecto a v4:
  v4     : TF × taxa_IDF × phylum_IDF; phylum arbiter sobre votos taxid
  v4-lca : igual que v4 + kmer_lca en índice; phylum arbiter usa LCA de k-mers
           → k-mers cross-phyla tienen LCA=Bacteria y NO votan por ningún phylum
           → elimina sesgo por número de taxa por phylum (Firmicutes >> Bacteroidetes)

El raw_model se cachea en disco para que re-runs sean rápidos (~2 min vs ~10 min).

Uso:
  cd ~/MetagenApp
  python3 scripts/rebuild_kraken_silva_v4_lca.py           # df-cap default 200
  python3 scripts/rebuild_kraken_silva_v4_lca.py --df-cap 500
  python3 scripts/rebuild_kraken_silva_v4_lca.py --no-cache  # re-entrena
"""

import argparse
import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA      = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX        = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
RAW_CACHE  = "/data/databases/metagenapp_refs/16S/raw_model_silva_genus150.pkl"
OUTPUT     = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v4.pkl"

# ── Parámetros ────────────────────────────────────────────────────────────────
K                   = 13
TRAIN_CAP           = 10_000
MAX_SEQS_PER_GENUS  = 150
DF_CAP_DEFAULT      = 200   # best from df_cap sweep

parser = argparse.ArgumentParser()
parser.add_argument("--df-cap",   type=int, default=DF_CAP_DEFAULT)
parser.add_argument("--output",   type=str, default=OUTPUT)
parser.add_argument("--no-cache", action="store_true", help="Re-entrenar aunque exista caché")
args = parser.parse_args()

print("=" * 65)
print("  rebuild_kraken_silva_v4_lca.py — MetagenApp")
print("=" * 65)
print(f"  FASTA              : {FASTA}")
print(f"  OUTPUT             : {args.output}")
print(f"  k                  : {K}")
print(f"  max_taxa_per_kmer  : {args.df_cap:,}")
print(f"  max_seqs_per_genus : {MAX_SEQS_PER_GENUS}")
print(f"  use_phylum_idf     : True")
print(f"  kmer_lca           : True  ← nuevo; LCA-based phylum voting")
print("=" * 65)

t0 = time.time()

# ── 1. Raw model ───────────────────────────────────────────────────────────────
raw_cache_path = Path(RAW_CACHE)
if not args.no_cache and raw_cache_path.exists():
    print(f"\n[1/3] Cargando raw_model desde caché ({RAW_CACHE})...")
    with open(RAW_CACHE, "rb") as f:
        raw = pickle.load(f)
    print(f"      Taxa : {len(raw['taxonomy_labels']):,}  K-mers: {len(raw['kmer_counts']):,}")
else:
    for path in [FASTA, TAX]:
        if not Path(path).exists():
            raise FileNotFoundError(f"No se encontró: {path}")

    print(f"\n[1/3] Entrenando raw_model (genus_cap={MAX_SEQS_PER_GENUS})...")
    raw = train_raw_kmer_postings(
        FASTA, TAX,
        k=K,
        max_taxa_per_kmer=TRAIN_CAP,
        max_seqs_per_genus=MAX_SEQS_PER_GENUS,
    )
    print(f"      Taxa : {len(raw['taxonomy_labels']):,}  K-mers: {len(raw['kmer_counts']):,}")
    print(f"      Tiempo entrenamiento: {time.time()-t0:.1f}s")

    raw_cache_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"      Guardando caché → {RAW_CACHE}")
    with open(RAW_CACHE, "wb") as f:
        pickle.dump(raw, f)

# ── 2. Construir índice con kmer_lca ──────────────────────────────────────────
t2 = time.time()
print(f"\n[2/3] Construyendo índice (df_cap={args.df_cap:,}, phylum_idf=True, kmer_lca=True)...")
print("      Esto incluye build_kid_to_lca — puede tardar ~2 min...")
index = build_kraken_index(raw, max_taxa_per_kmer=args.df_cap, use_phylum_idf=True)
n_kmers = len(index["kmer_index"])
n_lca   = len(index.get("kmer_lca", {}))
print(f"      K-mers en índice : {n_kmers:,}")
print(f"      K-mers con LCA   : {n_lca:,}  ({100*n_lca/max(n_kmers,1):.1f}%)")
print(f"      Taxa en índice   : {len(index['taxonomy']):,}")
print(f"      Tiempo           : {time.time()-t2:.1f}s")

# ── 3. Guardar ────────────────────────────────────────────────────────────────
t3 = time.time()
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
print(f"\n[3/3] Guardando → {args.output}")
with open(args.output, "wb") as f:
    pickle.dump(index, f)

size_mb = Path(args.output).stat().st_size / 1e6
print(f"      Tamaño : {size_mb:.0f} MB")
print(f"      Tiempo : {time.time()-t3:.1f}s")

total = time.time() - t0
print("\n" + "=" * 65)
print(f"  Índice SILVA v4-LCA listo en {total/60:.1f} min")
print(f"  {args.output}")
print("=" * 65)
print()
print("Para validar:")
print("  python3 scripts/exp_validate_kraken_silva_v4.py")
