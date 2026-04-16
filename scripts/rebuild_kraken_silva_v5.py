#!/usr/bin/env python3
"""
rebuild_kraken_silva_v5.py
==========================
Índice kraken-lite SILVA v5: separa el trainset de phylum LCA (cobertura)
del trainset de resolución taxonómica (precisión).

Problema resuelto:
  Con genus_cap=150, solo ~149K secuencias SILVA entran en el trainset.
  Los 13-mers de amplicons gut (ej. Prevotella copri) pueden no estar en las 150
  secuencias seleccionadas por genus → cobertura ~1.3% → phylum LCA voting inútil.
  Naive-v2 (sin genus_cap, 451K seqs, k=15) logra 60% de cobertura en el mismo query.

Solución v5:
  kmer_lca  → raw_model_FULL (sin genus_cap, 451K seqs) → alta cobertura
  kmer_index → raw_model_GENUS (genus_cap=150, 149K seqs) → sin sesgo por abundancia

Flujo de clasificación en el engine:
  1. Phylum arbiter: LCA de k-mers del query contra raw_model_FULL → cobertura alta,
     sin sesgo por tamaño de phylum (k-mers cross-phyla tienen LCA=Bacteria → no votan).
  2. Genus resolution: taxid votes contra kmer_index genus_capped → sin sesgo Firmicutes.
  3. Si kmer_index no tiene hits (amplicon lejos del trainset genus) pero kmer_lca tiene
     señal → clasificar a nivel phylum (correcto).

Caché:
  raw_model_GENUS → /data/databases/metagenapp_refs/16S/raw_model_silva_genus150.pkl
  raw_model_FULL  → /data/databases/metagenapp_refs/16S/raw_model_silva_full.pkl

Uso:
  python3 scripts/rebuild_kraken_silva_v5.py
  python3 scripts/rebuild_kraken_silva_v5.py --no-cache   # re-entrena todo
  python3 scripts/rebuild_kraken_silva_v5.py --df-cap 500
"""

import argparse
import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA          = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX            = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
CACHE_GENUS    = "/data/databases/metagenapp_refs/16S/raw_model_silva_genus150.pkl"
CACHE_FULL     = "/data/databases/metagenapp_refs/16S/raw_model_silva_full.pkl"
OUTPUT         = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v5.pkl"

# ── Parámetros ────────────────────────────────────────────────────────────────
K                  = 13
TRAIN_CAP          = 10_000
MAX_SEQS_PER_GENUS = 150
DF_CAP_DEFAULT     = 200

parser = argparse.ArgumentParser()
parser.add_argument("--df-cap",   type=int, default=DF_CAP_DEFAULT)
parser.add_argument("--output",   type=str, default=OUTPUT)
parser.add_argument("--no-cache", action="store_true")
args = parser.parse_args()

print("=" * 65)
print("  rebuild_kraken_silva_v5.py — MetagenApp")
print("  kmer_lca: FULL SILVA (sin genus_cap) → cobertura alta")
print("  kmer_index: genus_cap=150 → sin sesgo Firmicutes")
print("=" * 65)

t0 = time.time()

# ── 1a. raw_model_GENUS (genus_cap=150, para kmer_index) ──────────────────────
cache_genus = Path(CACHE_GENUS)
if not args.no_cache and cache_genus.exists():
    print(f"\n[1a/4] raw_model_GENUS — cargando desde caché...")
    with open(CACHE_GENUS, "rb") as f:
        raw_genus = pickle.load(f)
    print(f"       Taxa: {len(raw_genus['taxonomy_labels']):,}  K-mers: {len(raw_genus['kmer_counts']):,}")
else:
    print(f"\n[1a/4] raw_model_GENUS — entrenando (genus_cap={MAX_SEQS_PER_GENUS})...")
    raw_genus = train_raw_kmer_postings(
        FASTA, TAX, k=K, max_taxa_per_kmer=TRAIN_CAP,
        max_seqs_per_genus=MAX_SEQS_PER_GENUS,
    )
    print(f"       Taxa: {len(raw_genus['taxonomy_labels']):,}  K-mers: {len(raw_genus['kmer_counts']):,}")
    cache_genus.parent.mkdir(parents=True, exist_ok=True)
    print(f"       Guardando caché → {CACHE_GENUS}")
    with open(CACHE_GENUS, "wb") as f:
        pickle.dump(raw_genus, f)

# ── 1b. raw_model_FULL (sin genus_cap, para kmer_lca) ─────────────────────────
cache_full = Path(CACHE_FULL)
if not args.no_cache and cache_full.exists():
    print(f"\n[1b/4] raw_model_FULL — cargando desde caché...")
    with open(CACHE_FULL, "rb") as f:
        raw_full = pickle.load(f)
    print(f"       Taxa: {len(raw_full['taxonomy_labels']):,}  K-mers: {len(raw_full['kmer_counts']):,}")
else:
    print(f"\n[1b/4] raw_model_FULL — entrenando (SIN genus_cap, todos los ~451K seqs)...")
    print(f"       Esto puede tardar ~6 min...")
    t1b = time.time()
    raw_full = train_raw_kmer_postings(
        FASTA, TAX, k=K, max_taxa_per_kmer=TRAIN_CAP,
        max_seqs_per_genus=None,   # ← sin cap, usa TODAS las secuencias
    )
    print(f"       Taxa: {len(raw_full['taxonomy_labels']):,}  K-mers: {len(raw_full['kmer_counts']):,}")
    print(f"       Tiempo: {time.time()-t1b:.1f}s")
    cache_full.parent.mkdir(parents=True, exist_ok=True)
    print(f"       Guardando caché → {CACHE_FULL}")
    with open(CACHE_FULL, "wb") as f:
        pickle.dump(raw_full, f)

# ── 2. Verificar cobertura en gut ──────────────────────────────────────────────
print(f"\n[2/4] Verificando cobertura k={K} en gut amplicon...")
try:
    from metagenapp_core.utils.kmers import generate_kmers
    GUT_FASTA = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"
    DNA_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}
    def enc(s):
        v = 0
        for b in s:
            if b not in DNA_MAP: return None
            v = (v << 2) | DNA_MAP[b]
        return v

    seq = None
    with open(GUT_FASTA) as f:
        for line in f:
            if line.startswith(">"):
                if seq: break
                seq = ""
            elif seq is not None:
                seq += line.strip().upper()

    kc_genus = raw_genus["kmer_counts"]
    kc_full  = raw_full["kmer_counts"]
    total, in_genus, in_full = 0, 0, 0
    for kmer in generate_kmers(seq, K):
        if any(b not in DNA_MAP for b in kmer): continue
        total += 1
        if kmer in kc_genus: in_genus += 1
        if kmer in kc_full:  in_full  += 1

    print(f"       Secuencia gut: {len(seq)}bp, {total} k-mers")
    print(f"       Cobertura raw_model_GENUS (kmer_index): {in_genus}/{total} = {100*in_genus/total:.1f}%")
    print(f"       Cobertura raw_model_FULL  (kmer_lca)  : {in_full}/{total} = {100*in_full/total:.1f}%")
except Exception as e:
    print(f"       (skip cobertura: {e})")

# ── 3. Construir índice v5 ────────────────────────────────────────────────────
t3 = time.time()
print(f"\n[3/4] Construyendo índice v5...")
print(f"      kmer_index: raw_model_GENUS (df≤{args.df_cap:,})")
print(f"      kmer_lca  : raw_model_FULL  (todos los k-mers, sin df-filter)")
index = build_kraken_index(
    raw_genus,
    max_taxa_per_kmer=args.df_cap,
    use_phylum_idf=True,
    raw_model_for_lca=raw_full,
)
n_kmers = len(index["kmer_index"])
n_lca   = len(index.get("kmer_lca", {}))
print(f"      kmer_index: {n_kmers:,} k-mers")
print(f"      kmer_lca  : {n_lca:,} k-mers ({100*n_lca/max(len(raw_full['kmer_counts']),1):.1f}% del full trainset)")
print(f"      Tiempo: {time.time()-t3:.1f}s")

# ── 4. Guardar ────────────────────────────────────────────────────────────────
t4 = time.time()
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
print(f"\n[4/4] Guardando → {args.output}")
with open(args.output, "wb") as f:
    pickle.dump(index, f)

size_mb = Path(args.output).stat().st_size / 1e6
print(f"      Tamaño: {size_mb:.0f} MB  ({time.time()-t4:.1f}s)")

print("\n" + "=" * 65)
print(f"  Índice SILVA v5 listo en {(time.time()-t0)/60:.1f} min")
print(f"  {args.output}")
print("=" * 65)
print()
print("Para validar:")
print("  python3 scripts/exp_validate_kraken_silva_v5.py")
