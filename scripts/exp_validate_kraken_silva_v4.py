#!/usr/bin/env python3
"""
exp_validate_kraken_silva_v4.py
================================
Valida el índice kraken-lite SILVA v4 (phylum_idf) en dos datasets:
  1. Faringe (ya validado en v3)
  2. Gut/heces (problemático en v3)

Compara contra referencias Mothur/QIIME2.

Uso:
    cd ~/MetagenApp
    python3 scripts/exp_validate_kraken_silva_v4.py
"""

import pickle
from collections import Counter
from multiprocessing import cpu_count
from pathlib import Path

from metagenapp_core.models.kraken_lite import classify_kraken_parallel

# ── Paths ─────────────────────────────────────────────────────────────────────
INDEX_V3 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v3.pkl"
INDEX_V4 = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v4.pkl"

FARINGE_FASTA = "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta"
GUT_FASTA     = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"

TMP_V3 = "/tmp/kraken_v3_val.taxonomy"
TMP_V4 = "/tmp/kraken_v4_val.taxonomy"

THREADS = max(1, cpu_count() - 4)

# ── Referencias gold standard ─────────────────────────────────────────────────
REF_FARINGE = {
    "Firmicutes":      46.0,
    "Proteobacteria":  22.0,
    "Bacteroidetes":   15.0,
    "Actinobacteria":   9.0,
    "Fusobacteria":     5.0,
}
REF_GUT = {
    "Firmicutes":     40.0,
    "Bacteroidetes":  25.0,
    "Actinobacteria":  5.0,
    "Proteobacteria":  5.0,
    "Verrucomicrobia": 3.0,
}

SILVA_ALIASES = {
    "Bacillota":          "Firmicutes",
    "Pseudomonadota":     "Proteobacteria",
    "Bacteroidota":       "Bacteroidetes",
    "Actinomycetota":     "Actinobacteria",
    "Fusobacteriota":     "Fusobacteria",
    "Verrucomicrobiota":  "Verrucomicrobia",
}


def get_phylum(tax: str) -> str:
    if not tax or tax.strip() in ("", "Unclassified"):
        return "Unclassified"
    parts = [p.strip() for p in tax.split(";") if p.strip()]
    phylum = parts[1] if len(parts) >= 2 else parts[0] if parts else "Unclassified"
    return SILVA_ALIASES.get(phylum, phylum)


def classify_and_count(fasta: str, index: dict, tmp: str) -> Counter:
    classify_kraken_parallel(fasta, tmp, index, threads=THREADS)
    counts = Counter()
    with open(tmp) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            counts[get_phylum(parts[1])] += 1
    return counts


def rmse(counts: Counter, reference: dict) -> float:
    total = sum(counts.values()) or 1
    errors = [(100.0 * counts.get(p, 0) / total - r) ** 2 for p, r in reference.items()]
    return (sum(errors) / len(errors)) ** 0.5


def print_comparison(v3_counts: Counter, v4_counts: Counter, reference: dict, label: str):
    total_v3 = sum(v3_counts.values()) or 1
    total_v4 = sum(v4_counts.values()) or 1

    phyla = list(reference.keys()) + ["Unclassified"]
    print(f"\n  ── {label} ──")
    print(f"  {'phylum':<22}  {'ref':>6}  {'v3':>7}  {'Δv3':>7}  {'v4':>7}  {'Δv4':>7}")
    print(f"  {'-'*22}  {'-'*6}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*7}")
    for p in phyla:
        ref = reference.get(p)
        v3  = 100.0 * v3_counts.get(p, 0) / total_v3
        v4  = 100.0 * v4_counts.get(p, 0) / total_v4
        if ref is not None:
            dv3 = v3 - ref
            dv4 = v4 - ref
            f3  = "⚠" if abs(dv3) > 5 else "✓"
            f4  = "⚠" if abs(dv4) > 5 else "✓"
            print(f"  {p:<22}  {ref:>5.1f}%  {v3:>5.1f}% {f3}  {dv3:>+6.1f}  {v4:>5.1f}% {f4}  {dv4:>+6.1f}")
        elif v3 >= 0.5 or v4 >= 0.5:
            print(f"  {p:<22}        {v3:>5.1f}%           {v4:>5.1f}%")

    r3 = rmse(v3_counts, reference)
    r4 = rmse(v4_counts, reference)
    arrow = "↓ MEJOR" if r4 < r3 else ("↑ PEOR" if r4 > r3 else "= igual")
    print(f"\n  RMSE   v3={r3:.1f}pp  v4={r4:.1f}pp  {arrow}")


def main():
    for path in [INDEX_V4, FARINGE_FASTA, GUT_FASTA]:
        if not Path(path).exists():
            print(f"ERROR: no encontrado: {path}")
            if path == INDEX_V4:
                print("  → Corre primero: python3 scripts/rebuild_kraken_silva_v4.py")
            import sys; import sys; sys.exit(1)

    print("=" * 65)
    print("  Validación kraken-lite SILVA v3 vs v4 (phylum_idf)")
    print("=" * 65)
    print(f"  v3 index : {INDEX_V3}")
    print(f"  v4 index : {INDEX_V4}")
    print(f"  Threads  : {THREADS}")
    print("=" * 65)

    print("\n[1/4] Cargando índices...")
    with open(INDEX_V3, "rb") as f:
        model_v3 = pickle.load(f)
    print(f"  v3: {len(model_v3['kmer_index']):,} kmers, {len(model_v3['taxonomy']):,} taxa")

    with open(INDEX_V4, "rb") as f:
        model_v4 = pickle.load(f)
    print(f"  v4: {len(model_v4['kmer_index']):,} kmers, {len(model_v4['taxonomy']):,} taxa")

    print(f"\n[2/4] Clasificando faringe ({FARINGE_FASTA.split('/')[-1]})...")
    faringe_v3 = classify_and_count(FARINGE_FASTA, model_v3, TMP_V3)
    faringe_v4 = classify_and_count(FARINGE_FASTA, model_v4, TMP_V4)

    print(f"\n[3/4] Clasificando gut ({GUT_FASTA.split('/')[-1]})...")
    gut_v3 = classify_and_count(GUT_FASTA, model_v3, TMP_V3)
    gut_v4 = classify_and_count(GUT_FASTA, model_v4, TMP_V4)

    print(f"\n[4/4] Comparación v3 vs v4...")
    print_comparison(faringe_v3, faringe_v4, REF_FARINGE, "FARINGE")
    print_comparison(gut_v3, gut_v4, REF_GUT, "GUT/HECES")

    print("\n" + "=" * 65)
    far_r3 = rmse(faringe_v3, REF_FARINGE)
    far_r4 = rmse(faringe_v4, REF_FARINGE)
    gut_r3 = rmse(gut_v3, REF_GUT)
    gut_r4 = rmse(gut_v4, REF_GUT)

    print(f"  Resumen RMSE (pp)")
    print(f"  {'dataset':<12}  {'v3':>6}  {'v4':>6}  {'mejora':>8}")
    print(f"  {'faringe':<12}  {far_r3:>5.1f}   {far_r4:>5.1f}   {far_r3-far_r4:>+6.1f}pp")
    print(f"  {'gut':<12}  {gut_r3:>5.1f}   {gut_r4:>5.1f}   {gut_r3-gut_r4:>+6.1f}pp")

    gut_bact_v4 = 100.0 * gut_v4.get("Bacteroidetes", 0) / (sum(gut_v4.values()) or 1)
    print(f"\n  Diagnóstico clave — Bacteroidetes gut:")
    print(f"    v3 : {100.0 * gut_v3.get('Bacteroidetes',0)/(sum(gut_v3.values()) or 1):.1f}%  (ref 25%)")
    print(f"    v4 : {gut_bact_v4:.1f}%  (ref 25%)")
    if abs(gut_bact_v4 - 25.0) <= 5:
        print(f"    → RESUELTO: dentro de ±5pp del gold standard")
    elif gut_bact_v4 > 100.0 * gut_v3.get("Bacteroidetes", 0) / (sum(gut_v3.values()) or 1):
        print(f"    → MEJORA pero aún fuera de ±5pp")
    else:
        print(f"    → SIN MEJORA")
    print()


if __name__ == "__main__":
    main()
