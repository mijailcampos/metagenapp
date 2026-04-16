#!/usr/bin/env python3
"""
exp_kraken_df_cap_sweep.py
===========================
Barrido de max_taxa_per_kmer en build_kraken_index para replicar el efecto
del psize_max de naive-v2 (que descarta k-mers en >N taxa).

Hipótesis: naive-v2 funciona en gut (Bacteroidetes 22%) con los MISMOS 83K taxa
que kraken-lite v2 (Bacteroidetes 3.7%) porque psize_max=20 descarta k-mers
de alto df (cross-phyla, conservados). El equivalente en kraken-lite es bajar
max_taxa_per_kmer en build_kraken_index.

Con 47,904 taxa en v3: psize_max=20 en naive-v2 (1,949 taxa) equivale a
  20/1949 × 47904 ≈ 492 taxa. Barrer: [200, 500, 1000, 2000, 5000, 10000(baseline)].

Para no re-entrenar cada vez, se aprovecha que el raw_model es el mismo
(cambia solo el build_kraken_index).

Uso:
    cd ~/MetagenApp
    python3 scripts/exp_kraken_df_cap_sweep.py
"""

import pickle
import time
from collections import Counter
from multiprocessing import cpu_count
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index
from metagenapp_core.models.kraken_lite import classify_kraken_parallel

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX   = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"
GUT   = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"
FAR   = "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta"

K                  = 13
TRAIN_CAP          = 10_000
MAX_SEQS_PER_GENUS = 150

DF_CAPS = [200, 500, 1000, 2000, 5000, 10_000]   # 10K = baseline (igual a v3)

THREADS = max(1, cpu_count() - 4)

REF_GUT = {"Firmicutes":40,"Bacteroidetes":25,"Actinobacteria":5,"Proteobacteria":5,"Verrucomicrobia":3}
REF_FAR = {"Firmicutes":46,"Proteobacteria":22,"Bacteroidetes":15,"Actinobacteria":9,"Fusobacteria":5}

ALIASES = {
    "Bacillota":"Firmicutes","Pseudomonadota":"Proteobacteria",
    "Bacteroidota":"Bacteroidetes","Actinomycetota":"Actinobacteria",
    "Fusobacteriota":"Fusobacteria","Verrucomicrobiota":"Verrucomicrobia",
}

def get_ph(tax):
    if not tax or tax.strip() in ("","Unclassified"): return "Unclassified"
    p = [x.strip() for x in tax.split(";") if x.strip()]
    ph = p[1] if len(p)>=2 else p[0] if p else "Unclassified"
    return ALIASES.get(ph, ph)

def classify_fasta(fasta, model, tmp):
    classify_kraken_parallel(fasta, tmp, model, threads=THREADS)
    c = Counter()
    with open(tmp) as f:
        for line in f:
            pts = line.strip().split("\t")
            if len(pts)>=2: c[get_ph(pts[1])] += 1
    return c

def rmse(counts, ref):
    total = sum(counts.values()) or 1
    return ((sum((100*counts.get(p,0)/total - r)**2 for p,r in ref.items()))/len(ref))**0.5

print("=" * 65)
print("  Barrido max_taxa_per_kmer — kraken-lite SILVA")
print(f"  Equivalencia: psize_max naive-v2 = 20 → df_cap_kraken ≈ 492")
print(f"  DF caps: {DF_CAPS}")
print("=" * 65)

# ── Entrenar raw_model UNA sola vez ───────────────────────────────────────────
print("\n[1/N] Entrenando raw_model (genus_cap=150, se reutiliza para todos los caps)...")
t0 = time.time()
raw = train_raw_kmer_postings(FASTA, TAX, k=K,
    max_taxa_per_kmer=TRAIN_CAP, max_seqs_per_genus=MAX_SEQS_PER_GENUS)
print(f"  Listo en {time.time()-t0:.0f}s  ({len(raw['taxonomy_labels']):,} taxa, {len(raw['kmer_counts']):,} kmers)")

# ── Resultados ────────────────────────────────────────────────────────────────
results = []

for df_cap in DF_CAPS:
    print(f"\n{'='*50}")
    print(f"  max_taxa_per_kmer = {df_cap:,}")
    print(f"{'='*50}")

    t1 = time.time()
    idx = build_kraken_index(raw, max_taxa_per_kmer=df_cap, use_phylum_idf=True)
    n_kmers = len(idx["kmer_index"])
    print(f"  Índice: {n_kmers:,} kmers  (t={time.time()-t1:.0f}s)")

    gut = classify_fasta(GUT, idx, "/tmp/sweep_df_gut.tax")
    far = classify_fasta(FAR, idx, "/tmp/sweep_df_far.tax")

    gut_total = sum(gut.values()) or 1
    far_total = sum(far.values()) or 1

    bact_gut = 100*gut.get("Bacteroidetes",0)/gut_total
    firm_gut = 100*gut.get("Firmicutes",0)/gut_total
    firm_far = 100*far.get("Firmicutes",0)/far_total
    prot_far = 100*far.get("Proteobacteria",0)/far_total

    gut_r = rmse(gut, REF_GUT)
    far_r = rmse(far, REF_FAR)

    print(f"  GUT: Bact={bact_gut:.1f}% (ref 25%)  Firm={firm_gut:.1f}% (ref 40%)  RMSE={gut_r:.1f}")
    print(f"  FAR: Firm={firm_far:.1f}% (ref 46%)  Prot={prot_far:.1f}% (ref 22%)  RMSE={far_r:.1f}")

    results.append(dict(
        df_cap=df_cap, n_kmers=n_kmers,
        bact_gut=bact_gut, firm_gut=firm_gut,
        firm_far=firm_far, prot_far=prot_far,
        gut_rmse=gut_r, far_rmse=far_r,
        idx=idx,
    ))

# ── Tabla resumen ─────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print(f"  {'df_cap':>8}  {'kmers':>9}  {'Bact_gut%':>10}  {'Firm_gut%':>10}  {'gut_RMSE':>9}  {'Firm_far%':>10}  {'far_RMSE':>9}")
print("─" * 80)
for r in results:
    b_flag = "✓" if abs(r['bact_gut']-25)<=5 else "⚠"
    f_flag = "✓" if abs(r['firm_far']-46)<=5 else "⚠"
    print(f"  {r['df_cap']:>8,}  {r['n_kmers']:>9,}  "
          f"{r['bact_gut']:>8.1f}% {b_flag}  "
          f"{r['firm_gut']:>8.1f}%   "
          f"{r['gut_rmse']:>8.1f}   "
          f"{r['firm_far']:>8.1f}% {f_flag}  "
          f"{r['far_rmse']:>8.1f}")

# ── Mejor cap ─────────────────────────────────────────────────────────────────
# Criterio: minimizar gut_RMSE con faringe_RMSE ≤ 8
candidates = [r for r in results if r['far_rmse'] <= 8.0]
if not candidates:
    candidates = results
best = min(candidates, key=lambda r: r['gut_rmse'])

print(f"\n  Mejor df_cap: {best['df_cap']:,}")
print(f"    gut RMSE      : {best['gut_rmse']:.1f}")
print(f"    Bacteroidetes : {best['bact_gut']:.1f}% (ref 25%)")
print(f"    Firmicutes gut: {best['firm_gut']:.1f}% (ref 40%)")
print(f"    Firmicutes far: {best['firm_far']:.1f}% (ref 46%)")
print(f"    faringe RMSE  : {best['far_rmse']:.1f}")

# ── Guardar mejor índice ──────────────────────────────────────────────────────
OUT = f"/data/databases/metagenapp_refs/16S/kraken_index_silva_v4_dfcap{best['df_cap']}.pkl"
with open(OUT, "wb") as f:
    pickle.dump(best['idx'], f)
size_mb = Path(OUT).stat().st_size / 1e6
print(f"\n  Índice guardado: {OUT}  ({size_mb:.0f} MB)")
print(f"  Rebuild final con:")
print(f"    python3 scripts/rebuild_kraken_silva_v4.py --df-cap {best['df_cap']}")
