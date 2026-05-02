#!/usr/bin/env python3
"""
exp_kraken_gut_phylum_cap.py
============================
Barrido de max_seqs_per_phylum sobre el dataset de heces para encontrar el
cap óptimo que corrija la sub-estimación de Bacteroidetes en gut sin romper
la clasificación de faringe.

El problema:
  v3 (genus cap=150): gut → Bacteroidetes 1.8% (esperado ~25%), Firmicutes +27pp
  Causa: Firmicutes tiene ~5× más géneros que Bacteroidetes en SILVA → domina
         el trainset incluso con cap por genus equitativo.

La solución:
  Cap adicional por phylum (aplicado tras el cap por genus) para que cada
  phylum contribuya a lo sumo N secuencias al trainset.

Valores a barrer: [500, 1000, 2000, 3000, 5000, 8000]
  Con genus_cap=150 el trainset sin phylum_cap tiene Firmicutes ~75K seqs,
  Bacteroidetes ~20K. Un cap de 3000-5000 debería nivelarlos.

Uso:
    cd ~/MetagenApp
    python3 scripts/exp_kraken_gut_phylum_cap.py

Salida:
    Tabla de filos por cap, RMSE global, diagnóstico Bacteroidetes.
    También valida el cap ganador en faringe (no regresión).
"""

import pickle
import sys
import time
from collections import Counter
from pathlib import Path

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import build_kraken_index
from metagenapp_core.models.kraken_lite import classify_kraken_parallel
from multiprocessing import cpu_count

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA  = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta"
TAX    = "/data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax"

GUT_FASTA     = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"
FARINGE_FASTA = "/data/results/runs/ref_runs/faringe_nariz_20260414_1524/final_clean.fasta"

TMP_TAX = "/tmp/sweep_phylum_cap.taxonomy"

# ── Parámetros fijos (idénticos a v3) ─────────────────────────────────────────
K                  = 13
TRAIN_CAP          = 10_000
MAX_TAXA_PER_KMER  = 50_000
MAX_SEQS_PER_GENUS = 150       # se mantiene del v3

# ── Valores a barrer ──────────────────────────────────────────────────────────
PHYLUM_CAPS = [500, 1000, 2000, 3000, 5000, 8000]

THREADS = max(1, cpu_count() - 4)

# ── Referencias gold standard ─────────────────────────────────────────────────
REF_GUT = {
    "Firmicutes":     40.0,   # Bacillota
    "Bacteroidetes":  25.0,   # Bacteroidota
    "Actinobacteria":  5.0,
    "Proteobacteria":  5.0,
    "Verrucomicrobia": 3.0,
}
REF_FARINGE = {
    "Firmicutes":      46.0,
    "Proteobacteria":  22.0,
    "Bacteroidetes":   15.0,
    "Actinobacteria":   9.0,
    "Fusobacteria":     5.0,
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


def classify_fasta(fasta: str, index: dict) -> Counter:
    """Clasifica un FASTA con el índice dado y devuelve Counter de phyla."""
    classify_kraken_parallel(fasta, TMP_TAX, index, threads=THREADS)
    counts = Counter()
    with open(TMP_TAX) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            counts[get_phylum(parts[1])] += 1
    return counts


def rmse(counts: Counter, reference: dict) -> float:
    total = sum(counts.values()) or 1
    errors = []
    for phylum, ref_pct in reference.items():
        val = 100.0 * counts.get(phylum, 0) / total
        errors.append((val - ref_pct) ** 2)
    return (sum(errors) / len(errors)) ** 0.5


def print_distribution(counts: Counter, reference: dict, label: str):
    total = sum(counts.values()) or 1
    phyla = list(reference.keys()) + ["Unclassified"]
    print(f"\n  [{label}]")
    for p in phyla:
        n   = counts.get(p, 0)
        val = 100.0 * n / total
        ref = reference.get(p)
        if ref is not None:
            d   = val - ref
            flag = "⚠ " if abs(d) > 5 else "✓ "
            print(f"    {flag}{p:<22} {val:>5.1f}%  ref={ref:.1f}%  Δ={d:+.1f}pp")
        elif val >= 0.5:
            print(f"       {p:<22} {val:>5.1f}%")


def build_and_eval(phylum_cap):
    print(f"\n{'='*60}")
    print(f"  phylum_cap = {phylum_cap:,}")
    print(f"{'='*60}")

    t0 = time.time()
    raw = train_raw_kmer_postings(
        FASTA, TAX,
        k=K,
        max_taxa_per_kmer=TRAIN_CAP,
        max_seqs_per_genus=MAX_SEQS_PER_GENUS,
        max_seqs_per_phylum=phylum_cap,
    )
    index = build_kraken_index(raw, max_taxa_per_kmer=MAX_TAXA_PER_KMER)
    elapsed = time.time() - t0
    print(f"  Índice construido en {elapsed:.0f}s  "
          f"({len(index['kmer_index']):,} kmers, {len(index['taxonomy']):,} taxa)")

    gut_counts    = classify_fasta(GUT_FASTA, index)
    faringe_counts = classify_fasta(FARINGE_FASTA, index)

    gut_rmse    = rmse(gut_counts, REF_GUT)
    faringe_rmse = rmse(faringe_counts, REF_FARINGE)

    gut_total = sum(gut_counts.values()) or 1
    bact_pct  = 100.0 * gut_counts.get("Bacteroidetes", 0) / gut_total
    firm_pct  = 100.0 * gut_counts.get("Firmicutes", 0) / gut_total

    return {
        "cap":          phylum_cap,
        "gut_rmse":     gut_rmse,
        "faringe_rmse": faringe_rmse,
        "bact_pct":     bact_pct,
        "firm_pct":     firm_pct,
        "gut_counts":   gut_counts,
        "faringe_counts": faringe_counts,
        "index":        index,
    }


def main():
    for path in [FASTA, TAX, GUT_FASTA, FARINGE_FASTA]:
        if not Path(path).exists():
            print(f"ERROR: no se encontró {path}", file=sys.stderr)
            sys.exit(1)

    print("=" * 60)
    print("  Barrido max_seqs_per_phylum — kraken-lite SILVA")
    print(f"  genus_cap={MAX_SEQS_PER_GENUS}  k={K}  threads={THREADS}")
    print("=" * 60)
    print(f"  Gut FASTA    : {GUT_FASTA}")
    print(f"  Faringe FASTA: {FARINGE_FASTA}")
    print(f"  Caps a barrer: {PHYLUM_CAPS}")
    print()

    results = []
    for cap in PHYLUM_CAPS:
        r = build_and_eval(cap)
        results.append(r)
        print(f"  → gut RMSE={r['gut_rmse']:.1f}  "
              f"Bacteroidetes={r['bact_pct']:.1f}%  "
              f"Firmicutes={r['firm_pct']:.1f}%  "
              f"faringe RMSE={r['faringe_rmse']:.1f}")

    # ── Tabla resumen ─────────────────────────────────────────────────────────
    print("\n" + "=" * 75)
    print(f"  {'cap':>8}  {'gut_RMSE':>9}  {'Bacteroidetes%':>15}  "
          f"{'Firmicutes%':>12}  {'faringe_RMSE':>13}")
    print("─" * 75)
    for r in results:
        bact_flag  = "✓" if abs(r['bact_pct'] - 25.0) <= 5 else "⚠"
        far_flag   = "✓" if r['faringe_rmse'] <= 8.0 else "⚠"
        print(f"  {r['cap']:>8,}  {r['gut_rmse']:>9.1f}  "
              f"{r['bact_pct']:>14.1f}% {bact_flag}  "
              f"{r['firm_pct']:>11.1f}%  "
              f"{r['faringe_rmse']:>12.1f} {far_flag}")

    # ── Seleccionar mejor cap ─────────────────────────────────────────────────
    # Criterio: minimizar gut_RMSE con la restricción de faringe_RMSE <= umbral
    FARINGE_RMSE_MAX = 10.0
    candidates = [r for r in results if r['faringe_rmse'] <= FARINGE_RMSE_MAX]
    if not candidates:
        candidates = results  # si todos fallan en faringe, tomar el menos malo

    best = min(candidates, key=lambda r: r['gut_rmse'])

    print("\n" + "=" * 60)
    print(f"  Mejor cap: max_seqs_per_phylum = {best['cap']:,}")
    print(f"    gut RMSE     : {best['gut_rmse']:.1f}")
    print(f"    Bacteroidetes: {best['bact_pct']:.1f}%  (ref 25%,  Δ={best['bact_pct']-25:.1f}pp)")
    print(f"    Firmicutes   : {best['firm_pct']:.1f}%  (ref 40%,  Δ={best['firm_pct']-40:.1f}pp)")
    print(f"    faringe RMSE : {best['faringe_rmse']:.1f}")
    print("=" * 60)

    print_distribution(best['gut_counts'],    REF_GUT,    f"Heces — cap={best['cap']:,}")
    print_distribution(best['faringe_counts'], REF_FARINGE, f"Faringe — cap={best['cap']:,}")

    # ── Guardar el mejor índice ───────────────────────────────────────────────
    OUT = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v4_candidate.pkl"
    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    print(f"\n  Guardando índice candidato → {OUT}")
    with open(OUT, "wb") as f:
        pickle.dump(best['index'], f)
    size_mb = Path(OUT).stat().st_size / 1e6
    print(f"  Tamaño: {size_mb:.0f} MB")
    print(f"\n  Para validar: actualiza KRAKEN_INDEX_PATH en metagen_config.py a:")
    print(f'    "kraken_index_silva_v4_candidate.pkl"')
    print(f"\n  Para rebuild final con este cap:")
    print(f"    python3 scripts/rebuild_kraken_silva_v4.py --phylum-cap {best['cap']}")
    print()


if __name__ == "__main__":
    main()
