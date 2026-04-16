#!/usr/bin/env python3
"""
exp_validate_kraken_silva_v5.py
================================
Valida kraken-lite v3 (baseline) vs v5 (U→T fix + full-SILVA kmer_lca + LCA phylum voting)
en faringe y gut.
"""

import pickle
from collections import Counter
from multiprocessing import cpu_count
from pathlib import Path

from metagenapp_core.models.kraken_lite import classify_kraken_parallel

INDEX_V3 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v3.pkl"
INDEX_V5 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v5.pkl"

FARINGE_FASTA = "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta"
GUT_FASTA     = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"

TMP_V3 = "/tmp/kraken_v3_val.taxonomy"
TMP_V5 = "/tmp/kraken_v5_val.taxonomy"

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

def print_comparison(v3, v5, ref, label):
    t3 = sum(v3.values()) or 1
    t5 = sum(v5.values()) or 1
    phyla = list(ref.keys()) + ["Unclassified"]
    print(f"\n  ── {label} ──")
    print(f"  {'phylum':<22}  {'ref':>6}  {'v3':>7}  {'Δv3':>7}  {'v5':>7}  {'Δv5':>7}")
    print(f"  {'-'*22}  {'-'*6}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*7}")
    for p in phyla:
        r = ref.get(p)
        p3 = 100*v3.get(p,0)/t3
        p5 = 100*v5.get(p,0)/t5
        if r is not None:
            f3 = "✓" if abs(p3-r)<=5 else "⚠"
            f5 = "✓" if abs(p5-r)<=5 else "⚠"
            print(f"  {p:<22}  {r:>5.1f}%  {p3:>5.1f}% {f3}  {p3-r:>+6.1f}  {p5:>5.1f}% {f5}  {p5-r:>+6.1f}")
        elif p3>=0.5 or p5>=0.5:
            print(f"  {p:<22}        {p3:>5.1f}%           {p5:>5.1f}%")
    r3 = rmse(v3, ref); r5 = rmse(v5, ref)
    arrow = "↓ MEJOR" if r5<r3 else ("↑ PEOR" if r5>r3 else "= igual")
    print(f"\n  RMSE   v3={r3:.1f}pp  v5={r5:.1f}pp  {arrow}")

def main():
    print("="*65)
    print("  Validación kraken-lite v3 vs v5")
    print("  v5: U→T fix + full-SILVA kmer_lca + LCA phylum voting")
    print("="*65)

    print("\n[1/4] Cargando índices...")
    with open(INDEX_V3,"rb") as f: m3 = pickle.load(f)
    with open(INDEX_V5,"rb") as f: m5 = pickle.load(f)
    print(f"  v3: {len(m3['kmer_index']):,} kmers | kmer_lca: {len(m3.get('kmer_lca',{})):,}")
    print(f"  v5: {len(m5['kmer_index']):,} kmers | kmer_lca: {len(m5.get('kmer_lca',{})):,}")

    print(f"\n[2/4] Clasificando faringe...")
    f3 = classify_and_count(FARINGE_FASTA, m3, TMP_V3)
    f5 = classify_and_count(FARINGE_FASTA, m5, TMP_V5)

    print(f"\n[3/4] Clasificando gut...")
    g3 = classify_and_count(GUT_FASTA, m3, TMP_V3)
    g5 = classify_and_count(GUT_FASTA, m5, TMP_V5)

    print(f"\n[4/4] Resultados...")
    print_comparison(f3, f5, REF_FARINGE, "FARINGE")
    print_comparison(g3, g5, REF_GUT, "GUT/HECES")

    print("\n"+"="*65)
    print(f"  {'dataset':<12}  {'v3':>6}  {'v5':>6}  {'mejora':>8}")
    print(f"  {'faringe':<12}  {rmse(f3,REF_FARINGE):>5.1f}   {rmse(f5,REF_FARINGE):>5.1f}   "
          f"{rmse(f3,REF_FARINGE)-rmse(f5,REF_FARINGE):>+6.1f}pp")
    print(f"  {'gut':<12}  {rmse(g3,REF_GUT):>5.1f}   {rmse(g5,REF_GUT):>5.1f}   "
          f"{rmse(g3,REF_GUT)-rmse(g5,REF_GUT):>+6.1f}pp")

    bact5 = 100*g5.get("Bacteroidetes",0)/(sum(g5.values()) or 1)
    bact3 = 100*g3.get("Bacteroidetes",0)/(sum(g3.values()) or 1)
    print(f"\n  Bacteroidetes gut: v3={bact3:.1f}%  v5={bact5:.1f}%  (ref 25%)")
    if abs(bact5-25)<=5: print(f"  → RESUELTO ✓")
    elif bact5>bact3:    print(f"  → MEJORA ({bact5-bact3:+.1f}pp)")
    else:                print(f"  → SIN MEJORA")
    print()

if __name__=="__main__": main()
