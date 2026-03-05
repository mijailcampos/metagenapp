# ============================================================
# pipeline/filter_alignment_columns.py  — Biotech version
# ============================================================

import os
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq


def filter_alignment_columns(
    fasta_path,
    output_path="user_data/outputs/aligned_mafft_filtered_columns.fasta",
    vertical=True,
    trump="."
):
    """
    Filter alignment columns from a FASTA multiple sequence alignment.

    This function removes alignment columns in which all sequences contain
    the specified 'trump' character (commonly '.' or '-'), similar to the
    behavior of mothur's filter.seqs when using vertical filtering.

    Parameters
    ----------
    fasta_path : str
        Path to the input aligned FASTA file.
    output_path : str
        Path where the filtered FASTA alignment will be written.
    vertical : bool
        If True, applies vertical filtering: a column is removed if all
        sequences contain the trump character.
    trump : str
        Character that defines a non-informative or removable column.

    Returns
    -------
    tuple
        (output_path, message) if successful,
        (None, error_message) if an error occurs.
    """

    # ------------------------------------------------------------
    # Verificar existencia del archivo de entrada
    # ------------------------------------------------------------
    if not os.path.exists(fasta_path):
        return None, f"Input file not found: {fasta_path}"

    try:
        # ------------------------------------------------------------
        # Leer alineamiento completo
        # ------------------------------------------------------------
        alignment = AlignIO.read(fasta_path, "fasta")
        n_cols = alignment.get_alignment_length()
        n_seqs = len(alignment)

        positions_to_keep = []

        # ------------------------------------------------------------
        # 1. Recorrer columnas y decidir cuáles conservar
        # ------------------------------------------------------------
        for i in range(n_cols):
            column = alignment[:, i]

            # Filtrado vertical: descartar columnas uniformes con 'trump'
            if vertical:
                if all(base == trump for base in column):
                    continue

            positions_to_keep.append(i)

        # ------------------------------------------------------------
        # 2. Reconstruir nuevo alineamiento sin columnas removidas
        # ------------------------------------------------------------
        new_records = []
        for record in alignment:
            new_seq = "".join(record.seq[i] for i in positions_to_keep)
            record.seq = Seq(new_seq)
            new_records.append(record)

        filtered_alignment = MultipleSeqAlignment(new_records)

        # Crear carpeta de salida si es necesario
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Guardar FASTA resultante
        AlignIO.write(filtered_alignment, output_path, "fasta")

        # Mensaje final
        msg = (
            f"{len(positions_to_keep)} columns retained out of {n_cols}. "
            f"{n_cols - len(positions_to_keep)} columns removed based on trump='{trump}'."
        )

        return output_path, msg

    except Exception as e:
        return None, f"Error filtering alignment columns: {e}"
