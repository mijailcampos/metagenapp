import os
import numpy as np
import pandas as pd
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# Rutas estándar del pipeline
INPUT_FASTA = "user_data/outputs/filtered_contigs.fasta"
INPUT_COUNT = "user_data/outputs/filtered_contigs.count_table"

OUTPUT_FASTA = "user_data/outputs/unique_contigs.fasta"
OUTPUT_COUNT = "user_data/outputs/unique_contigs.count_table"


# ============================================================
#   FUNCIÓN PRINCIPAL — GENERAR SECUENCIAS ÚNICAS
# ============================================================
def generate_unique_contigs(
    fasta_path=INPUT_FASTA,
    count_table_path=INPUT_COUNT,
    output_fasta=OUTPUT_FASTA,
    output_count=OUTPUT_COUNT,
):
    """
    Generate unique contigs by collapsing identical sequences and 
    summing their abundances across samples.

    This step eliminates redundant contigs originating from the same biological
    variant and prepares the dataset for downstream alignment, clustering, and 
    taxonomic analysis.

    Parameters
    ----------
    fasta_path : str
        Path to the filtered FASTA file.
    count_table_path : str
        Path to the filtered count table.
    output_fasta : str
        Output FASTA path for unique sequences.
    output_count : str
        Output count table path for unique sequences.

    Returns
    -------
    output_fasta : str
        Path to the final FASTA file with unique contigs.
    output_count : str or None
        Path to the updated count table (None if unavailable).
    num_unique : int
        Number of unique sequences generated.
    """

    # Diccionario: secuencia → lista de IDs originales
    seq_to_ids = defaultdict(list)

    # Diccionario: ID original → SeqRecord (por si se requiere metadatos)
    id_to_record = {}

    # Leer FASTA de entrada
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq)
        seq_to_ids[seq_str].append(record.id)
        id_to_record[record.id] = record

    # Leer tabla de conteos original (si existe)
    original_count_df = None
    sample_names = []

    if count_table_path and os.path.exists(count_table_path):
        try:
            original_count_df = pd.read_csv(
                count_table_path, sep="\t", index_col=0)
            sample_names = original_count_df.columns.tolist()
        except Exception as e:
            print(f"⚠ No se pudo leer count_table: {e}")
            original_count_df = None

    # Preparar estructuras para secuencias únicas
    unique_records = []
    new_count_data = {}

    # Recorrer secuencias únicas
    for index, (seq_str, original_ids) in enumerate(seq_to_ids.items(), 1):

        # Crear nuevo ID formateado
        new_id = f"uniq_{index:05d}"

        # Crear SeqRecord único
        unique_records.append(
            SeqRecord(Seq(seq_str), id=new_id, description="")
        )

        # Sumar abundancias de todos los IDs originales
        if original_count_df is not None:
            summed = np.zeros(len(sample_names), dtype=int)

            for oid in original_ids:
                if oid in original_count_df.index:
                    summed += original_count_df.loc[oid].values

            new_count_data[new_id] = summed.tolist()

    # ---------------------------------------------------------
    # Guardar archivo FASTA de secuencias únicas
    # ---------------------------------------------------------
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    with open(output_fasta, "w") as f:
        SeqIO.write(unique_records, f, "fasta")

    # ---------------------------------------------------------
    # Guardar nueva tabla de conteos (si fue posible)
    # ---------------------------------------------------------
    if original_count_df is not None:
        if new_count_data:
            os.makedirs(os.path.dirname(output_count), exist_ok=True)
            new_df = pd.DataFrame.from_dict(
                new_count_data, orient="index", columns=sample_names
            )
            new_df.index.name = "Sequence_ID"
            new_df.to_csv(output_count, sep="\t")
        else:
            output_count = None
    else:
        output_count = None

    # Regresar rutas y número de secuencias únicas generadas
    return output_fasta, output_count, len(unique_records)
