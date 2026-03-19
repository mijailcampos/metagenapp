#!/usr/bin/env python3
"""
build_silva_trainset.py
=======================
Filtra SILVA 138.2 SSURef NR99 (Bacteria + Archaea) y genera el trainset
en el formato compatible con MetagenApp / kraken-lite.

Salida:
  /data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.fasta
  /data/databases/metagenapp_refs/trainset_silva/SILVA_NR99_BacArc.tax

Formato .tax (igual que trainset9_032012.pds):
  seq_id<TAB>Dominio;Filo;Clase;Orden;Familia;Género;Especie;

Uso:
  cd ~/MetagenApp
  python3 scripts/build_silva_trainset.py

  # Opcional: otro archivo de entrada
  python3 scripts/build_silva_trainset.py --input /ruta/SILVA.fasta
"""

import sys
import os
import re
import argparse
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
SILVA_FASTA_DEFAULT = Path.home() / "SILVA_138.2_SSURef_NR99_tax_silva.fasta"
OUT_DIR             = Path("/data/databases/metagenapp_refs/trainset_silva")
OUT_FASTA           = OUT_DIR / "SILVA_NR99_BacArc.fasta"
OUT_TAX             = OUT_DIR / "SILVA_NR99_BacArc.tax"

# Dominios a conservar
KEEP_DOMAINS = {"Bacteria", "Archaea"}

# Niveles taxonómicos (igual que trainset9: 7 niveles + ; al final)
TAX_LEVELS = 7

# ── Helpers ────────────────────────────────────────────────────────────────────

def parse_silva_header(header: str):
    """
    Parsea: '>AB000393.1.1510 Bacteria;Pseudomonadota;...'
    Retorna: (seq_id_limpio, tax_raw) o (None, None)
    """
    header = header.lstrip(">").strip()
    parts = header.split(" ", 1)
    if len(parts) < 2:
        return None, None

    seq_id  = parts[0].strip()
    tax_raw = parts[1].strip()

    # AB000393.1.1510 → AB000393_1_1510  (evita puntos en IDs)
    seq_id_clean = seq_id.replace(".", "_")
    return seq_id_clean, tax_raw


def normalize_taxonomy(tax_raw: str) -> str:
    """
    Normaliza taxonomía SILVA al formato trainset9:
      Bacteria;Filo;Clase;Orden;Familia;Género;Especie;
    - Trunca si sobran niveles
    - Rellena con _incertae_sedis si faltan
    - Termina siempre con ';'
    """
    levels = [t.strip() for t in tax_raw.split(";") if t.strip()]
    levels = levels[:TAX_LEVELS]

    while len(levels) < TAX_LEVELS:
        parent = re.sub(r'[^A-Za-z0-9]', '_', levels[-1]) if levels else "unclassified"
        levels.append(f"{parent}_incertae_sedis")

    return ";".join(levels) + ";"


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Build SILVA trainset for MetagenApp")
    parser.add_argument("--input", default=str(SILVA_FASTA_DEFAULT),
                        help="Ruta al SILVA SSURef NR99 .fasta (descomprimido)")
    parser.add_argument("--out-dir", default=str(OUT_DIR),
                        help="Directorio de salida")
    args = parser.parse_args()

    silva_path = Path(args.input)
    out_dir    = Path(args.out_dir)
    out_fasta  = out_dir / "SILVA_NR99_BacArc.fasta"
    out_tax    = out_dir / "SILVA_NR99_BacArc.tax"

    if not silva_path.exists():
        print(f"ERROR: No se encontró {silva_path}", file=sys.stderr)
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  build_silva_trainset.py — MetagenApp")
    print("=" * 60)
    print(f"  Input : {silva_path}  ({silva_path.stat().st_size / 1e6:.0f} MB)")
    print(f"  Output: {out_dir}")
    print(f"  Dominios: {KEEP_DOMAINS}")
    print("=" * 60)

    total_seqs     = 0
    kept_seqs      = 0
    skipped_domain = 0
    skipped_short  = 0
    skipped_header = 0

    current_id  = None
    current_tax = None
    seq_lines   = []

    def flush(fasta_out, tax_out):
        nonlocal kept_seqs, skipped_short
        if not (current_id and current_tax and seq_lines):
            return
        seq = "".join(seq_lines).replace("-", "").replace(".", "").upper()
        if len(seq) < 50:
            skipped_short += 1
            return
        fasta_out.write(f">{current_id}\n{seq}\n")
        tax_out.write(f"{current_id}\t{current_tax}\n")
        kept_seqs += 1

    with open(silva_path, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_fasta, "w") as f_out, \
         open(out_tax,   "w") as t_out:

        for line in fin:
            line = line.rstrip("\n")

            if line.startswith(">"):
                flush(f_out, t_out)
                seq_lines = []
                total_seqs += 1

                if total_seqs % 50_000 == 0:
                    pct = kept_seqs / total_seqs * 100
                    print(f"  [{total_seqs:>7,}] guardadas={kept_seqs:,} ({pct:.1f}%) "
                          f"| omitidas dominio={skipped_domain:,}", flush=True)

                seq_id, tax_raw = parse_silva_header(line)

                if seq_id is None:
                    current_id  = None
                    current_tax = None
                    skipped_header += 1
                    continue

                first_level = tax_raw.split(";")[0].strip()
                if first_level not in KEEP_DOMAINS:
                    current_id  = None
                    current_tax = None
                    skipped_domain += 1
                    continue

                current_id  = seq_id
                current_tax = normalize_taxonomy(tax_raw)

            else:
                if current_id:
                    seq_lines.append(line)

        flush(f_out, t_out)  # última secuencia

    # ── Resumen ────────────────────────────────────────────────────────────────
    print("=" * 60)
    print(f"  Total leídas        : {total_seqs:,}")
    print(f"  Bacteria + Archaea  : {kept_seqs:,}")
    print(f"  Otros dominios      : {skipped_domain:,}")
    print(f"  Secuencias <50 bp   : {skipped_short:,}")
    print(f"  Headers inválidos   : {skipped_header:,}")
    print("=" * 60)
    print(f"  ✓ FASTA → {out_fasta}")
    print(f"  ✓ TAX   → {out_tax}")
    print("=" * 60)
    print()
    print("Siguiente paso:")
    print("  python3 scripts/rebuild_kraken_silva.py")


if __name__ == "__main__":
    main()
