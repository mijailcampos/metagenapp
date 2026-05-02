#!/usr/bin/env python3
"""
exp_kraken_lca_diag.py
=======================
Diagnóstico de phylum LCA voting para secuencias gut.

Para una muestra de secuencias del dataset gut, muestra:
  - Cuántos k-mers totales tiene la secuencia
  - Cuántos k-mers están en kmer_lca (con LCA phylum-level)
  - Distribución de phylum_lca_votes
  - A qué phylum las asigna naive-v2 (taxonomy file del run)

Esto permite saber si el problema es:
  a) Cobertura baja en kmer_lca (pocos k-mers votan)
  b) Firmicutes gana incluso con LCA voting
  c) Los votos están dispersos entre muchos phyla

Uso:
  python3 scripts/exp_kraken_lca_diag.py
"""

import pickle
import random
from collections import Counter
from pathlib import Path

from metagenapp_core.utils.kmers import generate_kmers

# ── Paths ─────────────────────────────────────────────────────────────────────
INDEX_V4   = "/data/databases/metagenapp_refs/16S/kraken_index_silva_v4.pkl"
GUT_FASTA  = "/data/results/runs/ref_runs/heces_naive_20260415_0427/final_clean.fasta"
# Taxonomy del run naive-v2 como referencia gold
NAIVE_TAX  = "/data/results/runs/ref_runs/heces_naive_20260415_0427/taxonomy.tsv"

N_SAMPLE = 200
random.seed(42)

DNA_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}

def encode_kmer(seq):
    val = 0
    for b in seq:
        if b not in DNA_MAP:
            return None
        val = (val << 2) | DNA_MAP[b]
    return val

SILVA_ALIASES = {
    "Bacillota": "Firmicutes", "Pseudomonadota": "Proteobacteria",
    "Bacteroidota": "Bacteroidetes", "Actinomycetota": "Actinobacteria",
    "Fusobacteriota": "Fusobacteria", "Verrucomicrobiota": "Verrucomicrobia",
}

def get_phylum(tax: str) -> str:
    if not tax or tax.strip() in ("", "Unclassified"): return "Unclassified"
    parts = [p.strip() for p in tax.split(";") if p.strip()]
    ph = parts[1] if len(parts) >= 2 else parts[0] if parts else "Unclassified"
    return SILVA_ALIASES.get(ph, ph)

# ── Cargar índice ──────────────────────────────────────────────────────────────
print(f"Cargando índice {INDEX_V4}...")
with open(INDEX_V4, "rb") as f:
    model = pickle.load(f)

kmer_index = model["kmer_index"]
kmer_lca   = model.get("kmer_lca", {})
k          = model["k"]

print(f"  kmer_index: {len(kmer_index):,} k-mers")
print(f"  kmer_lca  : {len(kmer_lca):,} k-mers")
print(f"  k={k}")

# ── Cargar taxonomía naive-v2 ─────────────────────────────────────────────────
naive_tax = {}
if Path(NAIVE_TAX).exists():
    with open(NAIVE_TAX) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                naive_tax[parts[0]] = parts[1]
    print(f"  naive-v2 tax: {len(naive_tax):,} secuencias")
else:
    print(f"  AVISO: no se encontró {NAIVE_TAX}  (se omite referencia naive-v2)")

# ── Leer FASTA ────────────────────────────────────────────────────────────────
print(f"\nLeyendo {GUT_FASTA}...")
seqs = {}
sid = None
chunks = []
with open(GUT_FASTA) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if sid: seqs[sid] = "".join(chunks)
            sid = line[1:].split()[0]
            chunks = []
        elif line:
            chunks.append(line.upper())
    if sid: seqs[sid] = "".join(chunks)
print(f"  {len(seqs):,} secuencias")

# ── Muestrear ─────────────────────────────────────────────────────────────────
sample_ids = random.sample(list(seqs.keys()), min(N_SAMPLE, len(seqs)))

# ── Analizar phylum LCA votes por secuencia ───────────────────────────────────
results = []
for sid in sample_ids:
    seq = seqs[sid]
    total_kmers = 0
    in_lca = 0
    phylum_lca_votes = Counter()

    for kmer in generate_kmers(seq, k):
        key = encode_kmer(kmer)
        if key is None: continue
        total_kmers += 1

        lca_str = kmer_lca.get(key)
        if lca_str is not None:
            parts = [p.strip().replace('"', '') for p in lca_str.split(";") if p.strip()]
            if len(parts) >= 2:
                in_lca += 1
                phylum_lca_votes[";".join(parts[:2]) + ";"] += 1

    ph_total = sum(phylum_lca_votes.values())
    top_ph, top_count = phylum_lca_votes.most_common(1)[0] if phylum_lca_votes else ("None", 0)
    top_fraction = top_count / ph_total if ph_total else 0

    naive_ph = get_phylum(naive_tax.get(sid, "")) if naive_tax else "?"

    results.append({
        "sid": sid,
        "len": len(seq),
        "total_kmers": total_kmers,
        "in_lca": in_lca,
        "lca_coverage": in_lca / total_kmers if total_kmers else 0,
        "ph_total": ph_total,
        "top_ph": get_phylum(top_ph) if top_ph != "None" else "None",
        "top_fraction": top_fraction,
        "phylum_lca_votes": dict(phylum_lca_votes.most_common(4)),
        "naive_ph": naive_ph,
    })

# ── Resumen por phylum naive-v2 ────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  Diagnóstico phylum LCA votes — muestras gut")
print("=" * 70)

groups = {}
for r in results:
    groups.setdefault(r["naive_ph"], []).append(r)

for naive_ph in sorted(groups.keys(), key=lambda x: -len(groups[x])):
    grp = groups[naive_ph]
    avg_cov  = sum(r["lca_coverage"] for r in grp) / len(grp)
    avg_frac = sum(r["top_fraction"] for r in grp) / len(grp)
    avg_ph_total = sum(r["ph_total"] for r in grp) / len(grp)

    # ¿cuántas tienen el top phylum correcto (== naive_ph)?
    correct = sum(1 for r in grp if r["top_ph"] == naive_ph)
    top_phyla = Counter(r["top_ph"] for r in grp)

    print(f"\n  naive-v2 phylum: {naive_ph}  (n={len(grp)})")
    print(f"    LCA coverage   : {avg_cov*100:.1f}% de k-mers tienen LCA phylum-level")
    print(f"    ph_total medio : {avg_ph_total:.0f} votos por secuencia")
    print(f"    top fraction   : {avg_frac:.2f} (PHYLUM_MIN=0.25, PHYLUM_CONFIDENCE=0.55)")
    print(f"    top phylum correcto: {correct}/{len(grp)}  ({100*correct/len(grp):.0f}%)")
    print(f"    distribución top phyla votados: {dict(top_phyla.most_common(5))}")

    # Detalles para Bacteroidetes
    if naive_ph == "Bacteroidetes" and grp:
        print(f"\n    Detalle Bacteroidetes:")
        for r in grp[:5]:
            votes_str = "  ".join(f"{get_phylum(k).split(';')[0] if ';' in k else k}:{v}"
                                   for k, v in sorted(r["phylum_lca_votes"].items(),
                                                       key=lambda x: -x[1])[:4])
            print(f"      {r['sid'][:30]}  cov={r['lca_coverage']*100:.0f}%  "
                  f"top={r['top_ph']}({r['top_fraction']:.2f})  votes=[{votes_str}]")

print("\n" + "=" * 70)
print("  Interpretación:")
print("    cov<30%  → k=13 muy específico, pocos k-mers en referencia")
print("    top_frac<0.25 → señal difusa, va a Unclassified")
print("    top_ph correcto pero frac<0.25 → bajar PHYLUM_MIN")
print("    top_ph incorrecto (Firmicutes para Bacteroidetes) → sesgo residual")
print("=" * 70)
