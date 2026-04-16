#!/usr/bin/env python3
"""
build_homd_trainset.py
======================
Convierte los archivos crudos de HOMD a un trainset compatible con
el entrenador de naive-v2 de MetagenApp.

Formatos de entrada soportados:
  - FASTA HOMD con taxonomía en el header:
      >HOMD:tax_1234 | Streptococcus mutans UA159 | 16S rRNA | ...
  - Taxonomy file RDP-style:
      seq_id\tRoot;Bacteria;Firmicutes;Bacilli;...
  - Taxonomy file QIIME-style (2 columnas, sep=\t):
      seq_id\tk__Bacteria;p__Firmicutes;c__Bacilli;...

Salida:
  homd_trainset.fasta    — secuencias limpias con IDs simples
  homd_trainset.tax      — id<TAB>Root;Kingdom;Phylum;Class;Order;Family;Genus;Species

Uso:
  python3 scripts/build_homd_trainset.py \
      --fasta /data/databases/metagenapp_refs/16S/homd_raw/homd_16S.fasta \
      --tax   /data/databases/metagenapp_refs/16S/homd_raw/homd_16S.taxonomy \
      --outdir /data/databases/metagenapp_refs/16S/homd_raw/
"""

import argparse
import re
import sys
from pathlib import Path
from collections import Counter


# ============================================================
# Parsers de taxonomía
# ============================================================

def parse_tax_rdp(tax_path: Path) -> dict:
    """
    RDP format:
      seq_id\tRoot;Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;
    """
    mapping = {}
    with open(tax_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            seq_id = parts[0].strip()
            tax    = parts[1].strip().rstrip(";")
            mapping[seq_id] = tax
    return mapping


def parse_tax_qiime(tax_path: Path) -> dict:
    """
    QIIME format:
      seq_id\tk__Bacteria;p__Firmicutes;c__Bacilli;o__...;f__...;g__...;s__...
    Convierte a RDP-style: Root;Bacteria;Firmicutes;Bacilli;...
    """
    prefix_map = {"k__": "", "p__": "", "c__": "", "o__": "", "f__": "", "g__": "", "s__": ""}
    mapping = {}
    with open(tax_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            seq_id = parts[0].strip()
            raw    = parts[1].strip()
            # Strip prefixes k__, p__, etc.
            levels = []
            for level in raw.split(";"):
                level = level.strip()
                for prefix in prefix_map:
                    if level.startswith(prefix):
                        level = level[len(prefix):]
                        break
                if level and level not in ("", "unidentified", "uncultured"):
                    levels.append(level)
            tax = "Root;" + ";".join(levels)
            mapping[seq_id] = tax
    return mapping


def detect_tax_format(tax_path: Path) -> str:
    with open(tax_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "k__" in line or "p__" in line:
                return "qiime"
            return "rdp"
    return "rdp"


# ============================================================
# Parser FASTA — extrae taxonomía del header HOMD si no hay .tax
# ============================================================

def parse_homd_header(header: str) -> tuple[str, str | None]:
    """
    Extrae (seq_id, taxonomy_or_None) del header HOMD.

    Formatos conocidos:
      >HOMD:tax_1234 | Streptococcus mutans UA159 | 16S rRNA gene
      >SEQF1234  Streptococcus mutans  [Firmicutes] 16S
    """
    header = header.lstrip(">").strip()

    # Formato HOMD con pipes
    if "|" in header:
        parts = [p.strip() for p in header.split("|")]
        seq_id = parts[0].split()[0]
        # El nombre de especie suele estar en parts[1]
        species = parts[1] if len(parts) > 1 else ""
        return seq_id, None  # sin lineage completo en el header

    seq_id = header.split()[0]
    return seq_id, None


def read_fasta(fasta_path: Path) -> list[tuple[str, str, str]]:
    """
    Devuelve lista de (seq_id, header, sequence).
    """
    records = []
    current_id = None
    current_header = None
    current_seq = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, current_header, "".join(current_seq)))
                current_header = line[1:].strip()
                current_id     = current_header.split()[0]
                current_seq    = []
            else:
                current_seq.append(line.strip().upper())

    if current_id is not None:
        records.append((current_id, current_header, "".join(current_seq)))

    return records


# ============================================================
# Normalización de taxonomía
# ============================================================

def normalize_taxonomy(tax: str) -> str:
    """
    Asegura que la cadena tenga formato Root;Kingdom;Phylum;...
    Elimina anotaciones de confianza tipo (100) que usa RDP.
    """
    # Elimina confidence scores: "Bacteria(100);Firmicutes(100);..."
    tax = re.sub(r"\(\d+\)", "", tax)
    # Asegura que empiece con Root
    if not tax.startswith("Root"):
        tax = "Root;" + tax.lstrip(";")
    # Elimina nivel vacío al final
    tax = tax.rstrip(";")
    return tax


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(description="Construye trainset HOMD para naive-v2")
    ap.add_argument("--fasta",  required=True, help="FASTA de HOMD")
    ap.add_argument("--tax",    default=None,  help="Archivo de taxonomía (opcional si está en el FASTA)")
    ap.add_argument("--outdir", required=True, help="Directorio de salida")
    ap.add_argument("--min-len", type=int, default=200, help="Longitud mínima de secuencia (default: 200)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    outdir     = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Cargar taxonomía ---
    tax_map = {}
    if args.tax:
        tax_path = Path(args.tax)
        fmt = detect_tax_format(tax_path)
        print(f"Formato taxonomía detectado: {fmt}")
        if fmt == "qiime":
            tax_map = parse_tax_qiime(tax_path)
        else:
            tax_map = parse_tax_rdp(tax_path)
        print(f"  → {len(tax_map):,} entradas en taxonomy file")
    else:
        print("Sin taxonomy file — se intentará extraer del header FASTA")

    # --- Cargar FASTA ---
    print(f"\nLeyendo FASTA: {fasta_path}")
    records = read_fasta(fasta_path)
    print(f"  → {len(records):,} secuencias en FASTA")

    # --- Escribir salida ---
    out_fasta = outdir / "homd_trainset.fasta"
    out_tax   = outdir / "homd_trainset.tax"

    stats = Counter()
    written = 0

    with open(out_fasta, "w") as ff, open(out_tax, "w") as ft:
        for seq_id, header, seq in records:

            # Filtro longitud
            if len(seq) < args.min_len:
                stats["too_short"] += 1
                continue

            # Obtener taxonomía
            tax = tax_map.get(seq_id)
            if tax is None:
                stats["no_taxonomy"] += 1
                continue

            tax = normalize_taxonomy(tax)

            # Verificar al menos 3 niveles (Root;Kingdom;Phylum)
            levels = tax.split(";")
            if len(levels) < 3:
                stats["too_shallow"] += 1
                continue

            ff.write(f">{seq_id}\n{seq}\n")
            ft.write(f"{seq_id}\t{tax}\n")
            stats["written"] += 1
            written += 1

    # --- Resumen ---
    print(f"\n=== Resultado ===")
    print(f"  Escritas:       {stats['written']:>8,} secuencias")
    print(f"  Sin taxonomía:  {stats['no_taxonomy']:>8,}")
    print(f"  Demasiado cortas: {stats['too_short']:>6,} (< {args.min_len} bp)")
    print(f"  Taxonomía superficial: {stats['too_shallow']:>3,}")
    print(f"\n  Salida:")
    print(f"    {out_fasta}")
    print(f"    {out_tax}")

    # Distribución de filos
    print(f"\n=== Distribución de filos ===")
    phyla = Counter()
    with open(out_tax) as ft:
        for line in ft:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            levels = parts[1].split(";")
            phylum = levels[2] if len(levels) > 2 else "?"
            phyla[phylum] += 1
    for ph, cnt in phyla.most_common(15):
        pct = 100 * cnt / written if written else 0
        print(f"  {ph:<35} {cnt:>6,}  ({pct:.1f}%)")

    if written == 0:
        print("\nERROR: No se escribió ninguna secuencia.")
        print("  Verifica que los IDs del FASTA coincidan con los del archivo .taxonomy")
        print("  Ejemplo FASTA ID: >SEQF0001234")
        print("  Ejemplo tax  ID:  SEQF0001234")
        sys.exit(1)

    print("\nSiguiente paso:")
    print(f"  python3 scripts/train_naive_v2.py --marker 16S --model-type oral \\")
    print(f"    --fasta {out_fasta} --tax {out_tax}")


if __name__ == "__main__":
    main()
