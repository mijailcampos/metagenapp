#!/usr/bin/env python3
"""
Referencia de largo completo derivada de SILVA 138.2 NR99 (Bacteria + Archaea) para el modo ref de MetagenApp.

Independiente de la región amplificada (sirve para cualquier par de primers 16S). La usan:
  - 04_Alignment_vsearch (filtro de identidad mínima + recorte),
  - 12 mapeo de centroides (vsearch) y
  - el alignment fallback de taxonomía,
de modo que las tres cosas salen del mismo archivo y sus IDs siempre coinciden.

Entrada : trainset_silva/SILVA_NR99_BacArc.fasta (ARN) + .tax  (scripts/build_silva_trainset.py)
Salidas (en OUT_DIR):
  silva1382_full_all.fasta / .tax      todas las secuencias (U→T), taxonomía hasta género (para entrenar)
  silva1382_full_derep.fasta / .tax    secuencias únicas; taxonomía mayoritaria (≥90 %) o LCA
  build_report.txt
"""
from collections import defaultdict
from pathlib import Path

from build_v3v4_reference import read_fasta, to_genus, consensus, GENUS_LEVEL, MAJORITY

SRC = Path("/data/databases/metagenapp_refs/trainset_silva")
OUT_DIR = Path("/data/databases/metagenapp_refs/16S_silva1382_full")
MIN_LEN = 900   # SILVA Ref: Bacteria ≥1200, Archaea ≥900


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rep = open(OUT_DIR / "build_report.txt", "w")

    def log(msg):
        print(msg)
        rep.write(msg + "\n")
        rep.flush()

    tax = {}
    with open(SRC / "SILVA_NR99_BacArc.tax") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                tax[p[0]] = to_genus(p[1])

    by_seq = defaultdict(list)
    n_in = n_short = n_notax = n_all = 0
    with open(OUT_DIR / "silva1382_full_all.fasta", "w") as fa, open(OUT_DIR / "silva1382_full_all.tax", "w") as ta:
        for sid, seq in read_fasta(SRC / "SILVA_NR99_BacArc.fasta"):
            n_in += 1
            seq = seq.upper().replace("U", "T")
            if len(seq) < MIN_LEN:
                n_short += 1
                continue
            lin = tax.get(sid)
            if not lin:
                n_notax += 1
                continue
            fa.write(f">{sid}\n{seq}\n")
            ta.write(f"{sid}\t{';'.join(lin)};\n")
            by_seq[seq].append((sid, lin))
            n_all += 1
    log(f"entrada (SILVA 138.2 NR99 BacArc): {n_in:,} | < {MIN_LEN} pb: {n_short:,} | sin taxonomía: {n_notax}")
    log(f"secuencias escritas: {n_all:,}")

    n_conflict = 0
    genera = set()
    with open(OUT_DIR / "silva1382_full_derep.fasta", "w") as fd, open(OUT_DIR / "silva1382_full_derep.tax", "w") as td:
        for seq, members in by_seq.items():
            lins = [lin for _, lin in members]
            distinct = len({tuple(l) for l in lins})
            cons = consensus(lins) if distinct > 1 else lins[0]
            n_conflict += distinct > 1
            sid = members[0][0]
            fd.write(f">{sid}\n{seq}\n")
            td.write(f"{sid}\t{';'.join(cons)};\n")
            if len(cons) >= GENUS_LEVEL:
                genera.add(";".join(cons))
    log(f"secuencias únicas: {len(by_seq):,} | con taxonomías en conflicto (→ consenso ≥{MAJORITY:.0%} o LCA): {n_conflict:,}")
    log(f"géneros distintos en la referencia derep: {len(genera):,}")
    rep.close()


if __name__ == "__main__":
    main()
