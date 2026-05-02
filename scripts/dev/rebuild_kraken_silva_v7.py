#!/usr/bin/env python3
"""
rebuild_kraken_silva_v7.py
==========================
Índice kraken-lite SILVA v7: k=15 (más discriminativo que k=13).

Hipótesis: con k=15 los k-mers son más específicos por phylum.
Los reads de gut Bacteroidetes que con k=13 hacían match a Firmicutes
dejan de hacerlo → taxid votes correctos sin necesidad de kmer_lca override.

Diseño (mismo que v6 mixto pero con k=15):
  kmer_index: raw_model_GENUS (genus_cap=150, k=15) — misma filosofía v3
  kmer_lca  : raw_model_FULL  (sin genus_cap, k=15) — cobertura alta

Uso:
  python3 scripts/rebuild_kraken_silva_v7.py
  python3 scripts/rebuild_kraken_silva_v7.py --no-cache
"""

import argparse
import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index

FASTA           = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX             = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
CACHE_GENUS_K15 = "/data/databases/metagenapp_refs/16S/raw_model_silva_genus150_k15.pkl"
CACHE_FULL_K15  = "/data/databases/metagenapp_refs/16S/raw_model_silva_full_k15.pkl"
OUTPUT          = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v7.pkl"

K                  = 15
TRAIN_CAP          = 10_000
MAX_SEQS_PER_GENUS = 150
DF_CAP             = 50_000   # sin filtrado efectivo (como v3), k=15 ya es discriminativo

parser = argparse.ArgumentParser()
parser.add_argument("--no-cache", action="store_true")
args = parser.parse_args()

print("=" * 65)
print("  rebuild_kraken_silva_v7.py — MetagenApp")
print(f"  k=15 (vs k=13 en versiones anteriores)")
print(f"  kmer_index: genus_cap=150  |  kmer_lca: full SILVA")
print("=" * 65)

t0 = time.time()

# ── 1a. raw_model_GENUS k=15 ──────────────────────────────────────────────────
if not args.no_cache and Path(CACHE_GENUS_K15).exists():
    print(f"\n[1a/4] raw_model_GENUS k=15 — cargando caché...")
    with open(CACHE_GENUS_K15,"rb") as f: raw_genus = pickle.load(f)
else:
    print(f"\n[1a/4] raw_model_GENUS k=15 — entrenando (genus_cap={MAX_SEQS_PER_GENUS})...")
    print(f"       Esto tarda ~5-8 min...")
    t1 = time.time()
    raw_genus = train_raw_kmer_postings(
        FASTA, TAX, k=K, max_taxa_per_kmer=TRAIN_CAP,
        max_seqs_per_genus=MAX_SEQS_PER_GENUS,
    )
    print(f"       Taxa: {len(raw_genus['taxonomy_labels']):,}  K-mers: {len(raw_genus['kmer_counts']):,}")
    print(f"       Tiempo: {time.time()-t1:.1f}s")
    Path(CACHE_GENUS_K15).parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_GENUS_K15,"wb") as f: pickle.dump(raw_genus, f)
    print(f"       Caché → {CACHE_GENUS_K15}")

print(f"  Taxa: {len(raw_genus['taxonomy_labels']):,}  K-mers: {len(raw_genus['kmer_counts']):,}")

# ── 1b. raw_model_FULL k=15 ───────────────────────────────────────────────────
if not args.no_cache and Path(CACHE_FULL_K15).exists():
    print(f"\n[1b/4] raw_model_FULL k=15 — cargando caché...")
    with open(CACHE_FULL_K15,"rb") as f: raw_full = pickle.load(f)
else:
    print(f"\n[1b/4] raw_model_FULL k=15 — entrenando (sin genus_cap, ~451K seqs)...")
    print(f"       Esto tarda ~10-15 min...")
    t1 = time.time()
    raw_full = train_raw_kmer_postings(
        FASTA, TAX, k=K, max_taxa_per_kmer=TRAIN_CAP,
        max_seqs_per_genus=None,
    )
    print(f"       Taxa: {len(raw_full['taxonomy_labels']):,}  K-mers: {len(raw_full['kmer_counts']):,}")
    print(f"       Tiempo: {time.time()-t1:.1f}s")
    Path(CACHE_FULL_K15).parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_FULL_K15,"wb") as f: pickle.dump(raw_full, f)
    print(f"       Caché → {CACHE_FULL_K15}")

print(f"  Taxa: {len(raw_full['taxonomy_labels']):,}  K-mers: {len(raw_full['kmer_counts']):,}")

# ── 2. Verificar cobertura ────────────────────────────────────────────────────
print(f"\n[2/4] Verificando cobertura k={K}...")
try:
    from metagenapp_core.utils.kmers import generate_kmers
    DNA_MAP = {"A":0,"C":1,"G":2,"T":3}
    DATASETS = {
        "gut"    : "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta",
        "faringe": "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta",
    }
    def first_seq(path):
        seq = None
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    if seq is not None: break
                    seq = ""
                elif seq is not None:
                    seq += line.strip().upper()
        return seq

    kc_genus = raw_genus["kmer_counts"]
    kc_full  = raw_full["kmer_counts"]
    for label, path in DATASETS.items():
        seq = first_seq(path)
        if not seq: continue
        total = in_genus = in_full = 0
        for kmer in generate_kmers(seq, K):
            if any(b not in DNA_MAP for b in kmer): continue
            total += 1
            if kmer in kc_genus: in_genus += 1
            if kmer in kc_full:  in_full  += 1
        print(f"       {label:<8}: {total} k-mers  genus={100*in_genus/max(total,1):.1f}%  full={100*in_full/max(total,1):.1f}%")
except Exception as e:
    print(f"       (skip: {e})")

# ── 3. Construir índice mixto (v3-approach kmer_index + kmer_lca de full) ─────
print(f"\n[3/4] Construyendo índice v7...")
print(f"      kmer_index: raw_model_GENUS (df≤{DF_CAP:,}, use_phylum_idf=False)")
print(f"      kmer_lca  : raw_model_FULL")
t3 = time.time()
index = build_kraken_index(
    raw_genus,
    max_taxa_per_kmer=DF_CAP,
    use_phylum_idf=False,
    raw_model_for_lca=raw_full,
)
print(f"      kmer_index: {len(index['kmer_index']):,} k-mers")
print(f"      kmer_lca  : {len(index.get('kmer_lca',{})):,} k-mers")
print(f"      Tiempo: {time.time()-t3:.1f}s")

# ── 4. Guardar ────────────────────────────────────────────────────────────────
print(f"\n[4/4] Guardando → {OUTPUT}")
Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)
with open(OUTPUT,"wb") as f: pickle.dump(index, f)
import os
print(f"      Tamaño: {os.path.getsize(OUTPUT)/1e6:.0f} MB  total: {(time.time()-t0)/60:.1f} min")

print("\n" + "="*65)
print(f"  Índice SILVA v7 (k=15) listo")
print(f"  {OUTPUT}")
print("="*65)
print("\nPara validar:")
print("  python3 scripts/exp_validate_kraken_silva_v7.py")
