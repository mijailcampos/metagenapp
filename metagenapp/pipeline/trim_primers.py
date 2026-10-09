"""
trim_primers.py

Recorte de primers via cutadapt sobre contigs ya ensamblados (single-end,
amplicon completo). Quita el primer forward del extremo 5' y el reverso-
complementario del primer reverse del extremo 3', descartando lecturas
donde no se detecten ambos, para consistencia con el resto del benchmark
(QIIME2 usa cutadapt --p-discard-untrimmed sobre las mismas secuencias).
"""

import subprocess
from pathlib import Path
from Bio.Seq import Seq


def _read_count_table(path):
    with open(path) as f:
        header = f.readline()
        rows = {}
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            seq_id = line.split('\t')[0]
            rows[seq_id] = line
    return header, rows


def _read_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def trim_primers_cutadapt(
    fasta_path,
    count_table_path,
    output_fasta,
    output_count,
    primer_f,
    primer_r,
    threads=8,
):
    """
    Recorta primer_f (5') y RC(primer_r) (3') de cada secuencia en
    fasta_path via cutadapt, descarta las que no tengan ambos, y
    sincroniza el count_table con las secuencias sobrevivientes.

    Returns: (output_fasta, error, retained)
    """
    fasta_path = Path(fasta_path)
    count_table_path = Path(count_table_path)
    output_fasta = Path(output_fasta)
    output_count = Path(output_count)

    if not fasta_path.exists():
        return None, f'File not found: {fasta_path}', 0

    primer_r_rc = str(Seq(primer_r).reverse_complement())

    # Adaptador enlazado (linked adapter), anclado en ambos extremos: exige que
    # el amplicon completo empiece con primer_f y termine con RC(primer_r).
    # (probado: -g/-a sueltos sin anclar NO detectan el primer 5' aunque este
    # presente al inicio exacto del contig; y sin adaptador enlazado,
    # --discard-untrimmed solo exige encontrar UNO de los dos, no ambos)
    linked_adapter = f'^{primer_f}...{primer_r_rc}$'

    cmd = [
        'cutadapt',
        '-g', linked_adapter,
        '--discard-untrimmed',
        '-j', str(threads),
        '-o', str(output_fasta),
        str(fasta_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return None, f'cutadapt failed: {result.stderr}', 0

    surviving_ids = set(_read_fasta_ids(output_fasta))

    header, rows = _read_count_table(count_table_path)
    with open(output_count, 'w') as out:
        out.write(header)
        for seq_id in surviving_ids:
            if seq_id in rows:
                out.write(rows[seq_id] + '\n')

    return output_fasta, None, len(surviving_ids)
