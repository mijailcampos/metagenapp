#!/usr/bin/env python3
"""
rebuild_kraken_silva_v6.py
==========================
Índice kraken-lite SILVA v6: mismo diseño que v5 (kmer_lca desde modelo full
para cobertura alta) pero con max_seqs_per_phylum en raw_model_FULL para
eliminar el sesgo de representación que causó la regresión en faringe.

Problema resuelto en v6 vs v5:
  raw_model_FULL (v5) usa 451K seqs sin ningún cap. SILVA tiene un desbalance
  severo entre phyla: Bacillota 117K seqs, Pseudomonadota 136K seqs vs
  Bacteroidota 43K seqs. Con tantas secuencias de Firmicutes, k-mers de reads
  faringe-Firmicutes también aparecen en otras phyla → su LCA sube a Bacteria
  → dejan de votar por Firmicutes en phylum_lca_votes → Firmicutes colapsa
  de 45.4% a 31.8% y Bacteroidetes se dispara a 31.2% (ref 15%).

Solución v6:
  raw_model_FULL_balanced: max_seqs_per_phylum=20000 (configurable via --ph-cap)
  → Bacillota 20K, Bacteroidota 20K, Pseudomonadota 20K, Actinomycetota 20K
  → k-mers cross-phyla reducidos, phylum_lca_votes más representativos
  → cobertura gut mantenida (20K Bacteroidota >> 149K total genus_cap=150)

Índices:
  kmer_index: raw_model_GENUS (genus_cap=150) — igual que v5, sin cambios
  kmer_lca  : raw_model_FULL_balanced (max_seqs_per_phylum=20000)

Uso:
  python3 scripts/rebuild_kraken_silva_v6.py
  python3 scripts/rebuild_kraken_silva_v6.py --ph-cap 15000
  python3 scripts/rebuild_kraken_silva_v6.py --no-cache
"""

import argparse
import pickle
import time
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index

# ── Paths ──────────────────────────────────────────────────────────────────────
FASTA              = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX                = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
CACHE_GENUS        = "/data/databases/metagenapp_refs/16S/raw_model_silva_genus150.pkl"
OUTPUT_DEFAULT     = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v6.pkl"

# ── Parámetros ─────────────────────────────────────────────────────────────────
K                  = 13
TRAIN_CAP          = 10_000
MAX_SEQS_PER_GENUS = 150
DF_CAP_DEFAULT     = 200
PH_CAP_DEFAULT     = 20_000

parser = argparse.ArgumentParser()
parser.add_argument("--ph-cap",   type=int, default=PH_CAP_DEFAULT,
                    help="max_seqs_per_phylum para raw_model_FULL_balanced (default: 20000)")
parser.add_argument("--df-cap",   type=int, default=DF_CAP_DEFAULT)
parser.add_argument("--output",   type=str, default=OUTPUT_DEFAULT)
parser.add_argument("--no-cache", action="store_true",
                    help="Re-entrenar modelos aunque existan caches")
args = parser.parse_args()

CACHE_BALANCED = (
    f"/data/databases/metagenapp_refs/16S/"
    f"raw_model_silva_balanced_ph{args.ph_cap}.pkl"
)

print("=" * 65)
print("  rebuild_kraken_silva_v6.py — MetagenApp")
print(f"  kmer_lca : FULL SILVA + max_seqs_per_phylum={args.ph_cap}")
print(f"  kmer_index: genus_cap=150 → sin sesgo Firmicutes")
print("=" * 65)

t0 = time.time()

# ── 1a. raw_model_GENUS (genus_cap=150, para kmer_index) — igual que v5 ───────
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

# ── 1b. raw_model_FULL_balanced (ph_cap, para kmer_lca) ──────────────────────
cache_balanced = Path(CACHE_BALANCED)
if not args.no_cache and cache_balanced.exists():
    print(f"\n[1b/4] raw_model_FULL_balanced — cargando desde caché...")
    with open(CACHE_BALANCED, "rb") as f:
        raw_balanced = pickle.load(f)
    print(f"       Taxa: {len(raw_balanced['taxonomy_labels']):,}  K-mers: {len(raw_balanced['kmer_counts']):,}")
else:
    print(f"\n[1b/4] raw_model_FULL_balanced — entrenando (max_seqs_per_phylum={args.ph_cap})...")
    print(f"       Esto puede tardar ~3-5 min...")
    t1b = time.time()
    raw_balanced = train_raw_kmer_postings(
        FASTA, TAX, k=K, max_taxa_per_kmer=TRAIN_CAP,
        max_seqs_per_genus=None,          # sin cap por genus
        max_seqs_per_phylum=args.ph_cap,  # cap por phylum para balancear
    )
    print(f"       Taxa: {len(raw_balanced['taxonomy_labels']):,}  K-mers: {len(raw_balanced['kmer_counts']):,}")
    print(f"       Tiempo: {time.time()-t1b:.1f}s")
    cache_balanced.parent.mkdir(parents=True, exist_ok=True)
    print(f"       Guardando caché → {CACHE_BALANCED}")
    with open(CACHE_BALANCED, "wb") as f:
        pickle.dump(raw_balanced, f)

# ── 2. Verificar cobertura en gut y faringe ────────────────────────────────────
print(f"\n[2/4] Verificando cobertura k={K}...")
try:
    from metagenapp_core.utils.kmers import generate_kmers
    DNA_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}
    DATASETS = {
        "gut"    : "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta",
        "faringe": "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta",
    }

    def first_seq(fasta_path):
        seq = None
        with open(fasta_path) as f:
            for line in f:
                if line.startswith(">"):
                    if seq is not None:
                        break
                    seq = ""
                elif seq is not None:
                    seq += line.strip().upper()
        return seq

    kc_genus    = raw_genus["kmer_counts"]
    kc_balanced = raw_balanced["kmer_counts"]

    for label, path in DATASETS.items():
        seq = first_seq(path)
        if not seq:
            continue
        total, in_genus, in_balanced = 0, 0, 0
        for kmer in generate_kmers(seq, K):
            if any(b not in DNA_MAP for b in kmer):
                continue
            total += 1
            if kmer in kc_genus:
                in_genus += 1
            if kmer in kc_balanced:
                in_balanced += 1
        print(f"       {label:<8}: {total} k-mers  genus={100*in_genus/total:.1f}%  balanced={100*in_balanced/total:.1f}%")
except Exception as e:
    print(f"       (skip cobertura: {e})")

# ── 3. Construir índice v6 ─────────────────────────────────────────────────────
t3 = time.time()
print(f"\n[3/4] Construyendo índice v6...")
print(f"      kmer_index: raw_model_GENUS (df≤{args.df_cap:,})")
print(f"      kmer_lca  : raw_model_FULL_balanced (ph_cap={args.ph_cap:,})")
index = build_kraken_index(
    raw_genus,
    max_taxa_per_kmer=args.df_cap,
    use_phylum_idf=True,
    raw_model_for_lca=raw_balanced,
)
n_kmers = len(index["kmer_index"])
n_lca   = len(index.get("kmer_lca", {}))
print(f"      kmer_index: {n_kmers:,} k-mers")
print(f"      kmer_lca  : {n_lca:,} k-mers")
print(f"      Tiempo: {time.time()-t3:.1f}s")

# ── 4. Guardar ─────────────────────────────────────────────────────────────────
t4 = time.time()
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
print(f"\n[4/4] Guardando → {args.output}")
with open(args.output, "wb") as f:
    pickle.dump(index, f)
size_mb = Path(args.output).stat().st_size / 1e6
print(f"      Tamaño: {size_mb:.0f} MB  ({time.time()-t4:.1f}s)")

print("\n" + "=" * 65)
print(f"  Índice SILVA v6 listo en {(time.time()-t0)/60:.1f} min")
print(f"  {args.output}")
print("=" * 65)
print()
print("Para validar:")
print("  python3 scripts/exp_validate_kraken_silva_v6.py")
