#!/usr/bin/env python3
"""
train_naive_v2.py — SILVA-scale optimized
==========================================
Entrena un modelo naive-v2 para MetagenApp a partir de un trainset
FASTA + taxonomy en formato RDP (compatible con build_homd_trainset.py).

Optimizaciones para trainsets grandes (SILVA 138 NR99, ~510K seqs, ~50K taxa):

  1. Índice invertido directo (sin taxon_kmers intermedio)
     - El script anterior construía un índice forward {taxon → {kmer: count}}
       y lo pivotaba al invertido. Con SILVA esto consume ~60-90 GB en el pico.
     - Ahora se construye el índice invertido directamente durante el scan del FASTA.
     - Ahorro: ~30-50 GB de pico.

  2. Formato CSR (Compressed Sparse Row) por defecto
     - Los postings se guardan como dos arrays numpy int32 planos + offsets int64,
       en lugar de listas de listas Python.
     - Modelo SILVA estimado: ~3 GB vs ~18 GB con listas Python.
     - RAM en inferencia: ~5 GB vs ~20 GB.
     - Usar --no-csr para listas Python (compatible con modelos anteriores).

  3. kmer_counts eliminado del modelo
     - naive_v2_engine.py (Wang bootstrap) no usa kmer_counts.
     - Eliminarlo ahorra ~50% del tamaño del modelo a escala SILVA.

  4. --skip-lca: omite el cálculo de kid_to_lca
     - Wang bootstrap no necesita kid_to_lca.
     - Ahorra ~3 GB de modelo y varios minutos de build.
     - Seguro cuando postings están presentes (todos los modelos v4+).

  5. Monitoreo de RAM (RSS) en cada etapa.

  6. Pickle protocol 5 (más rápido para arrays numpy grandes).

Estructura del modelo generado:

  Formato CSR (default):
    format            — "csr_v1"
    k                 — tamaño de k-mer
    taxa              — list[str]: taxonomías únicas (leaf level)
    kmer_to_id        — dict[str, int]: kmer → ID
    pt_data           — np.int32[N]: taxon IDs de todos los postings concatenados
    pw_data           — np.int32[N]: conteos correspondientes
    pt_indptr         — np.int64[n_kmers+1]: offsets en pt_data por kmer_id
    kid_to_lca        — list[str] | None: LCA por kmer_id (None si --skip-lca)
    max_taxa_per_kmer — int

  Formato lista (--no-csr, backward-compatible con v4):
    format            — "list_v1"
    k, taxa, kmer_to_id, kid_to_lca, max_taxa_per_kmer  (igual que arriba)
    postings_taxon    — list[list[int]]
    postings_weight   — list[list[int]]

Uso típico (modelo SILVA general):

  python3 scripts/train_naive_v2.py \\
      --marker 16S --model-type general \\
      --fasta /data/databases/metagenapp_refs/16S/silva_138_NR99.fasta \\
      --tax   /data/databases/metagenapp_refs/16S/silva_138_NR99.tax \\
      --output /data/databases/metagenapp_refs/16S/naive_model_silva_v1.pkl \\
      --skip-lca

Uso modelo oral (HOMD, pequeño — sin optimizaciones especiales):

  python3 scripts/train_naive_v2.py \\
      --marker 16S --model-type oral \\
      --fasta /data/databases/metagenapp_refs/16S/homd_raw/homd_trainset.fasta \\
      --tax   /data/databases/metagenapp_refs/16S/homd_raw/homd_trainset.tax
"""

import argparse
import gc
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path


# ============================================================
# Utilidades
# ============================================================

def _rss_gb() -> float:
    """RAM del proceso actual (RSS) en GB."""
    import resource
    # ru_maxrss en Linux es kilobytes
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1_048_576


def generate_kmers(seq: str, k: int):
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if "N" not in kmer:
            yield kmer


def compute_lca(tax_strings: list) -> str:
    """
    LCA de un conjunto de strings de taxonomía con formato 'A;B;C;'.
    Ejemplo: ['Bacteria;Firmicutes;', 'Bacteria;Proteobacteria;'] → 'Bacteria;'
    """
    if not tax_strings:
        return ""
    if len(tax_strings) == 1:
        return tax_strings[0]
    split = [t.rstrip(";").split(";") for t in tax_strings]
    lca_parts = []
    for parts in zip(*split):
        if len(set(parts)) == 1:
            lca_parts.append(parts[0])
        else:
            break
    return ";".join(lca_parts) + ";" if lca_parts else ""


def read_fasta(path: Path):
    """
    Genera (seq_id, sequence) para cada entrada del FASTA.
    Convierte U→T automáticamente (trainsets SILVA usan formato RNA).
    """
    seq_id = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip().upper().replace("U", "T"))
    if seq_id is not None:
        yield seq_id, "".join(chunks)


def read_taxonomy(path: Path) -> dict:
    """
    Lee archivo tab-separado: seq_id<TAB>Root;Kingdom;...;
    Devuelve dict {seq_id → taxonomy_string_sin_Root}.
    """
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            seq_id = parts[0].strip()
            tax = parts[1].strip().rstrip(";")
            if tax.startswith("Root;"):
                tax = tax[5:]
            elif tax.startswith("Root"):
                tax = tax[4:].lstrip(";")
            mapping[seq_id] = tax
    return mapping


# ============================================================
# Entrenamiento
# ============================================================

def train(
    fasta_path: Path,
    tax_path: Path,
    k: int = 15,
    min_seq_len: int = 200,
    skip_lca: bool = False,
    use_csr: bool = True,
) -> dict:
    """
    Entrena un modelo naive-v2 desde FASTA + taxonomy.

    Returns:
        dict con todas las keys del modelo, listo para pickle.dump.
    """
    import numpy as np

    # --------------------------------------------------------
    # Leer taxonomía
    # --------------------------------------------------------
    print(f"\n=== Leyendo taxonomía ({tax_path.name}) ===")
    tax_map = read_taxonomy(tax_path)
    print(f"  {len(tax_map):,} entradas  |  RSS: {_rss_gb():.2f} GB")

    # --------------------------------------------------------
    # Fase 1: índice invertido directo
    #
    # kmer_inv[kmer][taxon_id] = count
    #
    # Por qué índice invertido directo en lugar de
    # taxon_kmers + pivot:
    #   - taxon_kmers ocupa ~30-50 GB a escala SILVA
    #   - kmer_inv directo ocupa ~10-15 GB (un solo dict)
    #   - Ahorro de pico: ~30-50 GB
    # --------------------------------------------------------
    print(f"\n=== Fase 1: indexando k-mers (k={k}) — índice invertido directo ===")

    taxon_to_id: dict = {}
    taxa: list = []
    kmer_inv: dict = defaultdict(lambda: defaultdict(int))

    n_seqs = 0
    n_skip = 0
    t0 = time.time()

    for seq_id, seq in read_fasta(fasta_path):
        tax = tax_map.get(seq_id)
        if tax is None or len(seq) < min_seq_len:
            n_skip += 1
            continue

        # Asignar ID de taxón al vuelo (orden de aparición en el FASTA)
        if tax not in taxon_to_id:
            taxon_to_id[tax] = len(taxa)
            taxa.append(tax)
        tid = taxon_to_id[tax]

        for kmer in generate_kmers(seq, k):
            kmer_inv[kmer][tid] += 1

        n_seqs += 1
        if n_seqs % 10_000 == 0:
            elapsed = time.time() - t0
            print(
                f"  {n_seqs:,} seqs | {len(kmer_inv):,} k-mers únicos | "
                f"{len(taxa):,} taxa | {elapsed:.0f}s | RSS {_rss_gb():.2f} GB",
                end="\r",
            )

    elapsed = time.time() - t0
    print(f"\n  {n_seqs:,} seqs procesadas ({elapsed:.1f}s) | {n_skip} omitidas")
    print(f"  {len(taxa):,} taxa  |  {len(kmer_inv):,} k-mers únicos  |  RSS: {_rss_gb():.2f} GB")

    del tax_map
    gc.collect()

    n_kmers = len(kmer_inv)

    # --------------------------------------------------------
    # Fase 2: construir kmer_to_id + postings + LCA
    # --------------------------------------------------------
    print(f"\n=== Fase 2: construyendo postings ({n_kmers:,} k-mers) ===")

    if use_csr:
        model = _build_csr(kmer_inv, taxa, n_kmers, skip_lca)
    else:
        model = _build_list(kmer_inv, taxa, n_kmers, skip_lca)

    model["k"] = k  # set aquí porque _build_* no conocen k

    print(f"\n  max_taxa_per_kmer : {model['max_taxa_per_kmer']}")
    print(f"  RSS final (pre-serialización): {_rss_gb():.2f} GB")
    return model


def _build_csr(kmer_inv, taxa, n_kmers, skip_lca):
    """
    Convierte kmer_inv al formato CSR numpy.
    Mucho más compacto que listas Python para modelos grandes.
    """
    import numpy as np

    kmer_to_id: dict = {}
    pt_data_list: list = []    # taxon IDs  (todos los postings concatenados)
    pw_data_list: list = []    # pesos/conteos correspondientes
    pt_indptr: list = [0]      # pt_indptr[kid] = inicio de postings del k-mer kid
    kid_to_lca: list = [] if not skip_lca else None
    max_taxa_per_kmer = 0

    t0 = time.time()
    for kid, (kmer, tid_counts) in enumerate(kmer_inv.items()):
        kmer_to_id[kmer] = kid

        tids   = list(tid_counts.keys())
        counts = [tid_counts[t] for t in tids]

        pt_data_list.extend(tids)
        pw_data_list.extend(counts)
        pt_indptr.append(len(pt_data_list))

        psize = len(tids)
        if psize > max_taxa_per_kmer:
            max_taxa_per_kmer = psize

        if not skip_lca:
            tax_strs = [taxa[t] + ";" for t in tids]
            kid_to_lca.append(compute_lca(tax_strs))

        if kid % 500_000 == 0 and kid > 0:
            elapsed = time.time() - t0
            print(
                f"  {kid:,}/{n_kmers:,} ({elapsed:.1f}s) | RSS {_rss_gb():.2f} GB",
                end="\r",
            )

    elapsed = time.time() - t0
    print(f"  {n_kmers:,}/{n_kmers:,} k-mers procesados ({elapsed:.1f}s) | RSS {_rss_gb():.2f} GB")

    del kmer_inv
    gc.collect()

    print(f"\n  Convirtiendo a arrays numpy... ", end="", flush=True)
    t_conv = time.time()
    pt_data    = np.array(pt_data_list,  dtype=np.int32)
    pw_data    = np.array(pw_data_list,  dtype=np.int32)
    pt_indptr_arr = np.array(pt_indptr,  dtype=np.int64)
    del pt_data_list, pw_data_list, pt_indptr
    gc.collect()
    print(f"listo ({time.time() - t_conv:.1f}s)")

    total_postings = len(pt_data)
    print(f"  postings totales : {total_postings:,}")
    print(f"  pt_data          : {pt_data.nbytes / 1e9:.3f} GB  ({pt_data.nbytes / 1e6:.1f} MB)")
    print(f"  pw_data          : {pw_data.nbytes / 1e9:.3f} GB  ({pw_data.nbytes / 1e6:.1f} MB)")
    print(f"  pt_indptr        : {pt_indptr_arr.nbytes / 1e6:.1f} MB")
    print(f"  RSS tras conversión: {_rss_gb():.2f} GB")

    return {
        "format":            "csr_v1",
        "taxa":              taxa,
        "kmer_to_id":        kmer_to_id,
        "pt_data":           pt_data,
        "pw_data":           pw_data,
        "pt_indptr":         pt_indptr_arr,
        "kid_to_lca":        kid_to_lca,
        "max_taxa_per_kmer": max_taxa_per_kmer,
    }


def _build_list(kmer_inv, taxa, n_kmers, skip_lca):
    """
    Formato lista Python (backward-compatible con modelos v4 anteriores).
    Usar solo para modelos pequeños (oral, gut con <100K seqs).
    """
    kmer_to_id: dict = {}
    postings_taxon:  list = [None] * n_kmers
    postings_weight: list = [None] * n_kmers
    kid_to_lca: list = ([""] * n_kmers) if not skip_lca else None
    max_taxa_per_kmer = 0

    t0 = time.time()
    for kid, (kmer, tid_counts) in enumerate(kmer_inv.items()):
        kmer_to_id[kmer] = kid

        tids   = list(tid_counts.keys())
        counts = [tid_counts[t] for t in tids]

        postings_taxon[kid]  = tids
        postings_weight[kid] = counts

        psize = len(tids)
        if psize > max_taxa_per_kmer:
            max_taxa_per_kmer = psize

        if not skip_lca:
            tax_strs = [taxa[t] + ";" for t in tids]
            kid_to_lca[kid] = compute_lca(tax_strs)

        if kid % 500_000 == 0 and kid > 0:
            elapsed = time.time() - t0
            print(
                f"  {kid:,}/{n_kmers:,} ({elapsed:.1f}s) | RSS {_rss_gb():.2f} GB",
                end="\r",
            )

    elapsed = time.time() - t0
    print(f"  {n_kmers:,}/{n_kmers:,} k-mers procesados ({elapsed:.1f}s) | RSS {_rss_gb():.2f} GB")

    del kmer_inv
    gc.collect()

    return {
        "format":            "list_v1",
        "taxa":              taxa,
        "kmer_to_id":        kmer_to_id,
        "postings_taxon":    postings_taxon,
        "postings_weight":   postings_weight,
        "kid_to_lca":        kid_to_lca,
        "max_taxa_per_kmer": max_taxa_per_kmer,
    }


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description="Entrena modelo naive-v2 de MetagenApp (SILVA-scale optimized)"
    )
    ap.add_argument("--marker",      default="16S",
                    help="Marcador genético: 16S | 18S")
    ap.add_argument("--model-type",  default="general",
                    help="Tipo de modelo: general | oral | gut | skin | env")
    ap.add_argument("--fasta",       required=True,
                    help="FASTA de entrenamiento")
    ap.add_argument("--tax",         required=True,
                    help="Archivo de taxonomía (RDP format, tab-separado)")
    ap.add_argument("--k",           type=int, default=15,
                    help="Tamaño de k-mer (default: 15)")
    ap.add_argument("--output",      default=None,
                    help="Ruta de salida .pkl (default: resuelve desde config)")
    ap.add_argument("--skip-lca",    action="store_true",
                    help="Omitir kid_to_lca. Wang bootstrap no la necesita. "
                         "Ahorra ~3-5 GB de modelo y varios minutos a escala SILVA.")
    ap.add_argument("--no-csr",      action="store_true",
                    help="Usar formato lista Python en lugar de CSR numpy. "
                         "Compatible con modelos anteriores pero mucho más grande.")
    ap.add_argument("--min-seq-len",   type=int,   default=200,
                    help="Mínimo de longitud de secuencia para incluir (default: 200)")
    ap.add_argument("--confidence",    type=float, default=0.80,
                    help="Umbral de confianza Wang bootstrap recomendado para este modelo "
                         "(default: 0.80 para modelos pequeños; 0.60 recomendado para SILVA)")
    ap.add_argument("--psize-max",     type=int,   default=20,
                    help="psize_max recomendado para este modelo (default: 20)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    tax_path   = Path(args.tax)
    marker     = args.marker.upper()
    model_type = args.model_type.lower()

    # Resolver ruta de salida
    if args.output:
        out_path = Path(args.output)
    else:
        from metagenapp.metagen_config import NAIVE_MODELS
        marker_registry = NAIVE_MODELS.get(marker)
        if marker_registry is None:
            print(f"ERROR: Marker '{marker}' no registrado en NAIVE_MODELS")
            sys.exit(1)
        out_path = marker_registry.get(model_type)
        if out_path is None:
            print(f"ERROR: model-type '{model_type}' no registrado para {marker}")
            sys.exit(1)
        out_path = Path(out_path)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not fasta_path.exists():
        print(f"ERROR: FASTA no encontrado: {fasta_path}"); sys.exit(1)
    if not tax_path.exists():
        print(f"ERROR: taxonomy no encontrado: {tax_path}"); sys.exit(1)

    fmt_label = "lista Python (--no-csr)" if args.no_csr else "CSR numpy (compacto)"
    lca_label = "omitir (--skip-lca)" if args.skip_lca else "calcular"

    print(f"=== MetagenApp — train naive-v2 ===")
    print(f"  Marcador  : {marker}")
    print(f"  Tipo      : {model_type}")
    print(f"  FASTA     : {fasta_path}")
    print(f"  Taxonomy  : {tax_path}")
    print(f"  k         : {args.k}")
    print(f"  Formato   : {fmt_label}")
    print(f"  LCA       : {lca_label}")
    print(f"  Salida    : {out_path}")
    print(f"  RSS inicio: {_rss_gb():.2f} GB")

    t_total = time.time()
    model = train(
        fasta_path,
        tax_path,
        k=args.k,
        min_seq_len=args.min_seq_len,
        skip_lca=args.skip_lca,
        use_csr=not args.no_csr,
    )

    # Guardar parámetros recomendados de inferencia dentro del modelo
    model["recommended_confidence"] = args.confidence
    model["recommended_psize_max"]  = args.psize_max

    print(f"\n=== Guardando modelo ===")
    print(f"  Ruta: {out_path}")
    t_save = time.time()
    with open(out_path, "wb") as f:
        pickle.dump(model, f, protocol=5)
    elapsed_save = time.time() - t_save

    size_gb = out_path.stat().st_size / 1e9
    elapsed_total = time.time() - t_total

    print(f"\n=== Modelo guardado ===")
    print(f"  Formato         : {model.get('format', 'legacy')}")
    print(f"  Taxa            : {len(model['taxa']):,}")
    print(f"  K-mers únicos   : {len(model['kmer_to_id']):,}")
    print(f"  Tamaño en disco : {size_gb:.3f} GB  ({size_gb * 1024:.1f} MB)")
    print(f"  confidence      : {model['recommended_confidence']}")
    print(f"  psize_max       : {model['recommended_psize_max']}")
    print(f"  Tiempo guardado : {elapsed_save:.1f}s")
    print(f"  Tiempo total    : {elapsed_total:.1f}s")
    print(f"\n  Para usar este modelo:")
    print(f"    metagenapp --classifier naive-v2 --model-type {model_type} --marker {marker} ...")


if __name__ == "__main__":
    main()
