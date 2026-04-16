#!/usr/bin/env python3
"""
exp_validate_kraken_silva_v7.py
================================
Valida kraken-lite v3 (baseline k=13) vs v6 (mixed k=13) vs v7 (k=15)
en faringe y gut.
"""

import pickle
from collections import Counter
from multiprocessing import cpu_count
from pathlib import Path

from metagenapp_core.models.kraken_lite import classify_kraken_parallel

INDEX_V3 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v3.pkl"
INDEX_V6 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v6.pkl"
INDEX_V7 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v7.pkl"

FARINGE_FASTA = "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta"
GUT_FASTA     = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"

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

def print_comparison(v3, v6, v7, ref, label):
    t3 = sum(v3.values()) or 1
    t6 = sum(v6.values()) or 1
    t7 = sum(v7.values()) or 1
    phyla = list(ref.keys()) + ["Unclassified"]
    print(f"\n  ── {label} ──")
    print(f"  {'phylum':<22}  {'ref':>6}  {'v3(k13)':>8}  {'v6(mix)':>8}  {'v7(k15)':>8}  {'Δv7':>7}")
    print(f"  {'-'*22}  {'-'*6}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*7}")
    for p in phyla:
        r = ref.get(p)
        p3 = 100*v3.get(p,0)/t3
        p6 = 100*v6.get(p,0)/t6
        p7 = 100*v7.get(p,0)/t7
        if r is not None:
            f7 = "✓" if abs(p7-r)<=5 else "⚠"
            print(f"  {p:<22}  {r:>5.1f}%  {p3:>6.1f}%   {p6:>6.1f}%   {p7:>6.1f}% {f7}  {p7-r:>+6.1f}")
        elif p3>=0.5 or p6>=0.5 or p7>=0.5:
            print(f"  {p:<22}        {p3:>6.1f}%   {p6:>6.1f}%   {p7:>6.1f}%")
    r3=rmse(v3,ref); r6=rmse(v6,ref); r7=rmse(v7,ref)
    best = min(r3,r6,r7)
    marks = {r3:"v3",r6:"v6",r7:"v7"}
    print(f"\n  RMSE  v3={r3:.1f}pp  v6={r6:.1f}pp  v7={r7:.1f}pp  ← mejor: {marks[best]}")

def main():
    print("="*72)
    print("  Validación kraken-lite v3 vs v6 vs v7")
    print("  v3: k=13 baseline  |  v6: k=13 mixed (v3+kmer_lca)  |  v7: k=15")
    print("="*72)

    print("\n[1/4] Cargando índices...")
    with open(INDEX_V3,"rb") as f: m3 = pickle.load(f)
    with open(INDEX_V6,"rb") as f: m6 = pickle.load(f)
    with open(INDEX_V7,"rb") as f: m7 = pickle.load(f)
    print(f"  v3: {len(m3['kmer_index']):,} kmers | kmer_lca: {len(m3.get('kmer_lca',{})):,}")
    print(f"  v6: {len(m6['kmer_index']):,} kmers | kmer_lca: {len(m6.get('kmer_lca',{})):,}")
    print(f"  v7: {len(m7['kmer_index']):,} kmers | kmer_lca: {len(m7.get('kmer_lca',{})):,}")

    print(f"\n[2/4] Clasificando faringe...")
    f3 = classify_and_count(FARINGE_FASTA, m3, "/tmp/kv3_far.tax")
    f6 = classify_and_count(FARINGE_FASTA, m6, "/tmp/kv6_far.tax")
    f7 = classify_and_count(FARINGE_FASTA, m7, "/tmp/kv7_far.tax")

    print(f"\n[3/4] Clasificando gut...")
    g3 = classify_and_count(GUT_FASTA, m3, "/tmp/kv3_gut.tax")
    g6 = classify_and_count(GUT_FASTA, m6, "/tmp/kv6_gut.tax")
    g7 = classify_and_count(GUT_FASTA, m7, "/tmp/kv7_gut.tax")

    print(f"\n[4/4] Resultados...")
    print_comparison(f3, f6, f7, REF_FARINGE, "FARINGE")
    print_comparison(g3, g6, g7, REF_GUT,     "GUT/HECES")

    print("\n"+"="*72)
    rf3=rmse(f3,REF_FARINGE); rf6=rmse(f6,REF_FARINGE); rf7=rmse(f7,REF_FARINGE)
    rg3=rmse(g3,REF_GUT);     rg6=rmse(g6,REF_GUT);     rg7=rmse(g7,REF_GUT)
    print(f"  {'dataset':<12}  {'v3':>6}  {'v6':>6}  {'v7':>6}")
    print(f"  {'faringe':<12}  {rf3:>5.1f}   {rf6:>5.1f}   {rf7:>5.1f}")
    print(f"  {'gut':<12}  {rg3:>5.1f}   {rg6:>5.1f}   {rg7:>5.1f}")

    bact_g7 = 100*g7.get("Bacteroidetes",0)/(sum(g7.values()) or 1)
    bact_f7 = 100*f7.get("Bacteroidetes",0)/(sum(f7.values()) or 1)
    print(f"\n  Bacteroidetes gut:     v3=4.0%  v6=4.3%  v7={bact_g7:.1f}%  (ref 25%)")
    print(f"  Bacteroidetes faringe: v3=14.0% v6=14.1% v7={bact_f7:.1f}%  (ref 15%)")

    gut_ok    = abs(bact_g7-25) <= 5
    far_ok    = rf7 <= 5.0
    print()
    if gut_ok and far_ok:
        print("  → AMBOS RESUELTOS ✓  — listo para actualizar metagen_config.py con v7")
    elif far_ok and not gut_ok:
        print(f"  → Faringe OK ✓  |  Gut mejorada ({bact_g7:.1f}% vs ref 25%) — continuar diagnóstico")
    elif gut_ok and not far_ok:
        print(f"  → Gut resuelta ✓  |  Faringe aún sesgada ⚠")
    else:
        delta_g = bact_g7 - 4.0
        print(f"  → Bacteroidetes gut: +{delta_g:.1f}pp respecto a v3 — diagnóstico continúa")
    print()

if __name__=="__main__": main()
