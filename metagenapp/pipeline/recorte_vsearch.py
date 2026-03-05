import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def recortar_por_alnout(
    aln_path,
    fasta_input,
    fasta_output
):
    """
    Trim FASTA sequences according to alignment coordinates obtained from a
    VSEARCH .alnout file.

    This function parses the alignment output, extracts query-specific
    start/end positions, and trims each sequence to retain only the well-aligned
    region. This is essential for removing overhangs, poorly aligned segments,
    and non-homologous regions prior to clustering or chimera detection.

    Parameters
    ----------
    aln_path : str
        Path to the VSEARCH .alnout alignment file.
    fasta_input : str
        Input FASTA file containing sequences prior to trimming.
    fasta_output : str
        Output FASTA file where trimmed sequences will be saved.

    Returns
    -------
    tuple
        (output_fasta_path, message) on success  
        (None, error_message) on failure
    """

    # Validaciones iniciales
    if not os.path.exists(aln_path):
        return None, "Alignment (.alnout) file not found."

    if not os.path.exists(fasta_input):
        return None, "Input FASTA file not found."

    posiciones_alineamiento = {}
    current_query_id = None

    # ===========================================================
    # 1. Parsear archivo .alnout para obtener coordenadas
    # ===========================================================
    try:
        with open(aln_path, 'r') as f:
            for line in f:

                # Identificar línea de query:  "Query >ID"
                query_match = re.search(r'Query\s+>\s*(\S+)', line)
                if query_match:
                    current_query_id = query_match.group(1)

                    # Inicializar coordenadas mínimas y máximas
                    posiciones_alineamiento[current_query_id] = [
                        float('inf'),
                        float('-inf')
                    ]

                # Identificar coordenadas de alineamiento: "Qry  1 + ... 64"
                coords = re.search(r'Qry\s+(\d+)\s+\S+\s+\S+\s+(\d+)', line)
                if coords and current_query_id:
                    start = int(coords.group(1))
                    end = int(coords.group(2))

                    # Guardar rangos mínimos y máximos
                    posiciones_alineamiento[current_query_id][0] = min(
                        posiciones_alineamiento[current_query_id][0], start
                    )
                    posiciones_alineamiento[current_query_id][1] = max(
                        posiciones_alineamiento[current_query_id][1], end
                    )

    except Exception as e:
        return None, f"Error parsing .alnout file: {e}"

    # ===========================================================
    # 2. Recorte del FASTA de entrada
    # ===========================================================
    recortadas = []
    total = 0
    recortadas_ok = 0

    try:
        for record in SeqIO.parse(fasta_input, "fasta"):
            total += 1

            # Si no existe información de alineamiento para esta secuencia, saltar
            if record.id not in posiciones_alineamiento:
                continue

            ini, fin = posiciones_alineamiento[record.id]

            # Validación de rangos válidos
            if ini <= fin and ini != float('inf') and fin != float('-inf'):
                nueva_seq = record.seq[ini - 1:fin]
                recortadas.append(
                    SeqRecord(
                        nueva_seq,
                        id=record.id,
                        description=""
                    )
                )
                recortadas_ok += 1

        # Guardar FASTA recortado
        os.makedirs(os.path.dirname(fasta_output), exist_ok=True)
        SeqIO.write(recortadas, fasta_output, "fasta")

        # Si ninguna secuencia se recortó correctamente
        if recortadas_ok == 0:
            return None, (
                "No sequences were trimmed. "
                "Check whether FASTA IDs match those in the .alnout file."
            )

        return fasta_output, f"{recortadas_ok} out of {total} sequences successfully trimmed."

    except Exception as e:
        return None, f"Error while trimming and writing FASTA: {e}"
