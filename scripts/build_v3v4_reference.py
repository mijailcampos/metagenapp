#!/usr/bin/env python3
"""
Referencia única V3–V4 derivada de SILVA 138.2 NR99 (Bacteria + Archaea) para MetagenApp.

Sustituye a 16S/silva_reference.fasta (14,956 secuencias del viejo silva.bacteria de mothur)
en el filtro usearch_global, el mapeo de centroides y el fallback de taxonomía, para que
las tres cosas salgan del mismo archivo.

Entrada : trainset_silva/SILVA_NR99_BacArc.fasta (ARN) + .tax  (scripts/build_silva_trainset.py)
Paso 1  : U→T y extracción in silico con cutadapt, primers enlazados 341F...RC(805R),
          --discard-untrimmed (se conserva solo la región entre primers, como en el pipeline).
Paso 2  : taxonomía truncada a género (6 rangos).
Salidas (en OUT_DIR):
  silva1382_v3v4_all.fasta / .tax      todas las secuencias extraídas (para entrenar el clasificador)
  silva1382_v3v4_derep.fasta / .tax    secuencias únicas; si una secuencia aparece con varias
                                       taxonomías, se guarda la mayoritaria (≥90 %) o, si no hay, el LCA
  build_report.txt
"""
import subprocess, sys
from collections import defaultdict
from pathlib import Path

SRC = Path("/data/databases/metagenapp_refs/trainset_silva")
OUT_DIR = Path("/data/databases/metagenapp_refs/16S_silva1382_v3v4")
CUTADAPT = "/data/envs/edna_env/bin/cutadapt"
PRIMER_F = "CCTACGGGNGGCWGCAG"
PRIMER_R_RC = "GGATTAGATACCCBDGTAGTC"   # reverso complementario de 805R GACTACHVGGGTATCTAATCC
MAX_ERR = 0.2
MIN_LEN = 300
GENUS_LEVEL = 6
MAJORITY = 0.90   # consenso: línea mayoritaria si >= 90 % de las entradas idénticas coinciden (criterio tipo RESCRIPt)


def read_fasta(path):
    sid, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(chunks)
                sid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if sid is not None:
        yield sid, "".join(chunks)


def to_genus(tax):
    levels = [p.strip() for p in tax.rstrip(";").split(";") if p.strip()]
    return levels[:GENUS_LEVEL]


def consensus(lineages):
    """Línea mayoritaria si alcanza MAJORITY; si no, ancestro común (LCA)."""
    from collections import Counter
    top, n = Counter(tuple(l) for l in lineages).most_common(1)[0]
    if n / len(lineages) >= MAJORITY:
        return list(top)
    return lca(lineages)


def lca(lineages):
    out = []
    for ranks in zip(*lineages):
        if len(set(ranks)) != 1:
            break
        out.append(ranks[0])
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rep = open(OUT_DIR / "build_report.txt", "w")

    def log(msg):
        print(msg)
        rep.write(msg + "\n")
        rep.flush()

    dna = OUT_DIR / "tmp_silva_bacarc_dna.fasta"
    n_in = 0
    with open(dna, "w") as fo:
        for sid, seq in read_fasta(SRC / "SILVA_NR99_BacArc.fasta"):
            fo.write(f">{sid}\n{seq.upper().replace('U', 'T')}\n")
            n_in += 1
    log(f"secuencias de entrada (SILVA 138.2 NR99 BacArc): {n_in:,}")

    extracted = OUT_DIR / "tmp_v3v4.fasta"
    cmd = [CUTADAPT, "-g", f"{PRIMER_F}...{PRIMER_R_RC}", "-e", str(MAX_ERR),
           "--discard-untrimmed", "-m", str(MIN_LEN), "-j", "16",
           "-o", str(extracted), str(dna)]
    log("cutadapt: " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log(r.stderr)
        sys.exit(1)
    for line in r.stdout.splitlines():
        if any(w in line for w in ("Total reads processed", "Reads with adapters", "Reads written", "too short")):
            log("  " + line.strip())

    tax = {}
    with open(SRC / "SILVA_NR99_BacArc.tax") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                tax[p[0]] = to_genus(p[1])

    by_seq = defaultdict(list)
    n_all = n_notax = 0
    with open(OUT_DIR / "silva1382_v3v4_all.fasta", "w") as fa, open(OUT_DIR / "silva1382_v3v4_all.tax", "w") as ta:
        for sid, seq in read_fasta(extracted):
            lin = tax.get(sid)
            if not lin:
                n_notax += 1
                continue
            fa.write(f">{sid}\n{seq}\n")
            ta.write(f"{sid}\t{';'.join(lin)};\n")
            by_seq[seq].append((sid, lin))
            n_all += 1
    log(f"extraídas con taxonomía: {n_all:,} (sin taxonomía: {n_notax})")

    n_conflict = 0
    genera = set()
    with open(OUT_DIR / "silva1382_v3v4_derep.fasta", "w") as fd, open(OUT_DIR / "silva1382_v3v4_derep.tax", "w") as td:
        for seq, members in by_seq.items():
            lins = [lin for _, lin in members]
            cons = consensus(lins) if len({tuple(l) for l in lins}) > 1 else lins[0]
            if len({tuple(l) for l in lins}) > 1:
                n_conflict += 1
            sid = members[0][0]
            fd.write(f">{sid}\n{seq}\n")
            td.write(f"{sid}\t{';'.join(cons)};\n")
            if len(cons) >= GENUS_LEVEL:
                genera.add(";".join(cons))
    log(f"secuencias únicas: {len(by_seq):,} | con taxonomías en conflicto (→ consenso ≥{MAJORITY:.0%} o LCA): {n_conflict:,}")
    log(f"géneros distintos en la referencia derep: {len(genera):,}")

    lens = sorted(len(s) for s in by_seq)
    log(f"longitud únicas: min {lens[0]}, mediana {lens[len(lens)//2]}, max {lens[-1]}")
    dna.unlink()
    extracted.unlink()
    rep.close()


if __name__ == "__main__":
    main()
