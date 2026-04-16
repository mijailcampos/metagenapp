#!/usr/bin/env python3
"""
exp_validate_kraken_silva_v6.py
================================
Valida kraken-lite v3 (baseline) vs v6 (v5 index + hybrid phylum arbiter).

v6 = índice v5 (U→T fix, full-SILVA kmer_lca) + engine híbrido:
  - taxid votes fuertes (>= PHYLUM_CONFIDENCE) → usa taxid (v3 behavior, preciso)
  - taxid votes débiles + kmer_lca disponible → usa kmer_lca (cubre gut reads)

v5 números de referencia (engine anterior, todos los reads usan kmer_lca):
  faringe RMSE=9.9pp  Bacteroidetes=31.2%  (peor que v3)
  gut     RMSE=11.2pp Bacteroidetes=14.5%  (mejor que v3)

Hipótesis: v6 recupera RMSE faringe de v3 (~3pp) y mantiene o mejora gut vs v5.
"""

import pickle
from collections import Counter
from multiprocessing import cpu_count
from pathlib import Path

from metagenapp_core.models.kraken_lite import classify_kraken_parallel

INDEX_V3 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v3.pkl"
INDEX_V6 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v6.pkl"

FARINGE_FASTA = "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta"
GUT_FASTA     = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"

TMP_V3 = "/tmp/kraken_v3_val6.taxonomy"
TMP_V6 = "/tmp/kraken_v6_val6.taxonomy"

THREADS = max(1, cpu_count() - 4)

REF_FARINGE = {"Firmicutes":46.0,"Proteobacteria":22.0,"Bacteroidetes":15.0,
               "Actinobacteria":9.0,"Fusobacteria":5.0}
REF_GUT     = {"Firmicutes":40.0,"Bacteroidetes":25.0,"Actinobacteria":5.0,
               "Proteobacteria":5.0,"Verrucomicrobia":3.0}

ALIASES = {"Bacillota":"Firmicutes","Pseudomonadota":"Proteobacteria",
           "Bacteroidota":"Bacteroidetes","Actinomycetota":"Actinobacteria",
           "Fusobacteriota":"Fusobacteria","Verrucomicrobiota":"Verrucomicrobia"}

def get_phylum(tax):
    if not tax or tax.strip() in ("","Unclassified"): return "Unclassified"
    parts = [p.strip() for p in tax.split(";") if p.strip()]
    ph = parts[1] if len(parts)>=2 else parts[0] if parts else "Unclassified"
    return ALIASES.get(ph, ph)

def classify_and_count(fasta, index, tmp):
    classify_kraken_parallel(fasta, tmp, index, threads=THREADS)
    counts = Counter()
    with open(tmp) as f:
        for line in f:
            pts = line.strip().split("\t")
            if len(pts)>=2: counts[get_phylum(pts[1])] += 1
    return counts

def rmse(counts, ref):
    total = sum(counts.values()) or 1
    return (sum((100*counts.get(p,0)/total - r)**2 for p,r in ref.items())/len(ref))**0.5

# Resultados v5 (engine anterior) para referencia
V5_FARINGE = {"Firmicutes":31.8,"Bacteroidetes":31.2,"Proteobacteria":17.6,
              "Actinobacteria":6.6,"Fusobacteria":6.1}
V5_GUT     = {"Firmicutes":61.1,"Bacteroidetes":14.5,"Actinobacteria":8.4,
              "Proteobacteria":12.3,"Verrucomicrobia":0.2}

def print_comparison(v3, v5_ref, v6, ref, label):
    t3 = sum(v3.values()) or 1
    t6 = sum(v6.values()) or 1
    phyla = list(ref.keys()) + ["Unclassified"]
    print(f"\n  ── {label} ──")
    print(f"  {'phylum':<22}  {'ref':>6}  {'v3':>7}  {'v5*':>7}  {'v6':>7}  {'Δv6':>7}")
    print(f"  {'-'*22}  {'-'*6}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*7}")
    for p in phyla:
        r = ref.get(p)
        p3 = 100*v3.get(p,0)/t3
        p5 = v5_ref.get(p, 0.0)
        p6 = 100*v6.get(p,0)/t6
        if r is not None:
            f6 = "✓" if abs(p6-r)<=5 else "⚠"
            print(f"  {p:<22}  {r:>5.1f}%  {p3:>5.1f}%  {p5:>5.1f}%  {p6:>5.1f}% {f6}  {p6-r:>+6.1f}")
        elif p3>=0.5 or p5>=0.5 or p6>=0.5:
            print(f"  {p:<22}        {p3:>5.1f}%  {p5:>5.1f}%  {p6:>5.1f}%")
    r3 = rmse(v3,ref); r5 = rmse(Counter({k:v for k,v in v5_ref.items()}), ref); r6 = rmse(v6,ref)
    arrow = "↓ MEJOR que v5" if r6<r5 else ("↑ PEOR que v5" if r6>r5 else "= igual v5")
    print(f"\n  RMSE  v3={r3:.1f}pp  v5*={r5:.1f}pp  v6={r6:.1f}pp  {arrow}")
    print(f"  (* v5 = resultados previos con engine original)")

def main():
    print("="*70)
    print("  Validación kraken-lite v3 vs v6 (hybrid phylum arbiter)")
    print("  v6 = nuevo índice (U→T fix, NO phylum_idf, kmer_lca full) + engine híbrido")
    print("  v5* = resultados previos para referencia")
    print("="*70)

    print("\n[1/4] Cargando índices...")
    with open(INDEX_V3,"rb") as f: m3 = pickle.load(f)
    with open(INDEX_V6,"rb") as f: m6 = pickle.load(f)
    print(f"  v3: {len(m3['kmer_index']):,} kmers | kmer_lca: {len(m3.get('kmer_lca',{})):,}")
    print(f"  v6: {len(m6['kmer_index']):,} kmers | kmer_lca: {len(m6.get('kmer_lca',{})):,}")

    print(f"\n[2/4] Clasificando faringe...")
    f3 = classify_and_count(FARINGE_FASTA, m3, TMP_V3)
    f6 = classify_and_count(FARINGE_FASTA, m6, TMP_V6)

    print(f"\n[3/4] Clasificando gut...")
    g3 = classify_and_count(GUT_FASTA, m3, TMP_V3)
    g6 = classify_and_count(GUT_FASTA, m6, TMP_V6)

    print(f"\n[4/4] Resultados...")
    print_comparison(f3, V5_FARINGE, f6, REF_FARINGE, "FARINGE")
    print_comparison(g3, V5_GUT,     g6, REF_GUT,     "GUT/HECES")

    print("\n"+"="*70)
    r3f=rmse(f3,REF_FARINGE); r5f=3.0; r6f=rmse(f6,REF_FARINGE)
    r3g=rmse(g3,REF_GUT);     r5g=11.2; r6g=rmse(g6,REF_GUT)
    print(f"  {'dataset':<12}  {'v3':>6}  {'v5*':>6}  {'v6':>6}  {'v6 vs v3':>10}  {'v6 vs v5':>10}")
    print(f"  {'faringe':<12}  {r3f:>5.1f}   {r5f:>5.1f}   {r6f:>5.1f}   {r6f-r3f:>+8.1f}pp  {r6f-r5f:>+8.1f}pp")
    print(f"  {'gut':<12}  {r3g:>5.1f}   {r5g:>5.1f}   {r6g:>5.1f}   {r6g-r3g:>+8.1f}pp  {r6g-r5g:>+8.1f}pp")

    bact_g3 = 100*g3.get("Bacteroidetes",0)/(sum(g3.values()) or 1)
    bact_g6 = 100*g6.get("Bacteroidetes",0)/(sum(g6.values()) or 1)
    bact_f3 = 100*f3.get("Bacteroidetes",0)/(sum(f3.values()) or 1)
    bact_f6 = 100*f6.get("Bacteroidetes",0)/(sum(f6.values()) or 1)
    print(f"\n  Bacteroidetes gut:     v3={bact_g3:.1f}%  v5*=14.5%  v6={bact_g6:.1f}%  (ref 25%)")
    print(f"  Bacteroidetes faringe: v3={bact_f3:.1f}%  v5*=31.2%  v6={bact_f6:.1f}%  (ref 15%)")

    ok_gut    = abs(bact_g6-25) <= 5
    ok_far    = abs(bact_f6-15) <= 5
    ok_farrmse = r6f <= 5.0
    print()
    if ok_gut and ok_farrmse:
        print("  → AMBOS DATASETS RESUELTOS ✓  — listo para actualizar metagen_config.py")
    elif ok_farrmse and not ok_gut:
        print("  → Faringe recuperada ✓  |  Gut mejorada pero aún lejos de ref")
    elif ok_gut and not ok_farrmse:
        print("  → Gut resuelta ✓  |  Faringe aún sesgada ⚠")
    else:
        print("  → Mejora parcial — continuar diagnóstico")
    print()

if __name__=="__main__": main()
