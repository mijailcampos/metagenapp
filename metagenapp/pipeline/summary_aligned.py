import os
import re
import pandas as pd
from Bio import SeqIO


def summary_seqs(fasta_path, count_table_path=None):
    """
    Genera un resumen estadístico de las secuencias en un archivo FASTA.
    Calcula longitudes, ambigüedades, homopolímeros y conteos totales.

    Argumentos:
        fasta_path (str): ruta del archivo FASTA alineado.
        count_table_path (str): tabla de conteos (opcional).
    """

    count_data = {}

    # Leer tabla de conteos
    if count_table_path and os.path.exists(count_table_path):
        try:
            count_df = pd.read_csv(count_table_path, sep="\t", index_col=0)
            count_data = count_df.sum(axis=1).to_dict()
        except Exception:
            count_data = {}

    stats = []

    # Procesar archivo FASTA
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        nbases = len(seq)

        ambigs = sum(1 for base in seq if base in "RYWSKMBDHVrywskmbdhv")

        max_poly = 0
        for base_type in "ACGT":
            for match in re.finditer(rf"({base_type}+)", seq, re.IGNORECASE):
                max_poly = max(max_poly, len(match.group(0)))

        count = count_data.get(record.id, 1)
        stats.append([nbases, ambigs, max_poly, count])

    if not stats:
        return {
            "Minimum": {"NBases": 0, "Ambigs": 0, "Polymer": 0, "NumSeqs": 0},
            "Median": {"NBases": 0, "Ambigs": 0, "Polymer": 0, "NumSeqs": 0},
            "Maximum": {"NBases": 0, "Ambigs": 0, "Polymer": 0, "NumSeqs": 0},
            "Mean": {"NBases": 0, "Ambigs": 0, "Polymer": 0, "NumSeqs": 0},
            "# of unique seqs": 0,
            "total # of seqs": 0
        }

    df = pd.DataFrame(
        stats, columns=["NBases", "Ambigs", "Polymer", "NumSeqs"])

    summary = {
        "Minimum": df.min().to_dict(),
        "2.5%-tile": df.quantile(0.025).round().astype(int).to_dict(),
        "25%-tile": df.quantile(0.25).round().astype(int).to_dict(),
        "Median": df.median().round().astype(int).to_dict(),
        "75%-tile": df.quantile(0.75).round().astype(int).to_dict(),
        "97.5%-tile": df.quantile(0.975).round().astype(int).to_dict(),
        "Maximum": df.max().to_dict(),
        "Mean": df.mean().round(1).to_dict(),
        "# of unique seqs": len(df),
        "total # of seqs": df["NumSeqs"].sum()
    }

    return summary
