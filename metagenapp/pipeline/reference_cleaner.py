import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def clean_reference_fasta(
    input_fasta,
    output_fasta,
    remove_gaps=True,
    remove_dots=True,
    uppercase=True
):
    """
    Clean a reference FASTA file by applying several preprocessing steps:

    Operations
    ----------
    - Remove alignment gaps ("-")
    - Remove dot characters (".")
    - Convert sequences to uppercase

    This function standardizes reference databases such as SILVA, PR2,
    Greengenes, RDP, and UNITE before alignment or classification.

    Parameters
    ----------
    input_fasta : str
        Path to the original reference FASTA file.
    output_fasta : str
        Path where the cleaned FASTA will be saved.
    remove_gaps : bool
        Whether to remove gap characters ("-").
    remove_dots : bool
        Whether to remove dot characters (".").
    uppercase : bool
        Whether to convert sequences to uppercase.

    Returns
    -------
    tuple
        (output_fasta, message) if successful  
        (None, error_message) if an error occurred
    """

    # Validación: verificar archivo de entrada
    if not os.path.exists(input_fasta):
        return None, f"File not found: {input_fasta}"

    try:
        # Crear carpeta de salida si no existe
        os.makedirs(os.path.dirname(output_fasta), exist_ok=True)

        cleaned_records = []

        # Procesar cada secuencia del FASTA
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)

            # Eliminar gaps según configuración
            if remove_gaps:
                seq = seq.replace("-", "")

            # Eliminar puntos según configuración
            if remove_dots:
                seq = seq.replace(".", "")

            # Convertir a mayúsculas si está activado
            if uppercase:
                seq = seq.upper()

            cleaned_records.append(
                SeqRecord(
                    Seq(seq),
                    id=record.id,
                    description=record.description
                )
            )

        # Guardar FASTA final limpio
        with open(output_fasta, "w") as out:
            SeqIO.write(cleaned_records, out, "fasta")

        return output_fasta, "Reference FASTA cleaned successfully."

    except Exception as e:
        return None, f"Error cleaning FASTA: {e}"
