import os
from Bio import SeqIO


def extraer_subset_centroides(
    fasta_path,
    output_path,
    num_secuencias=1000
):
    """
    Extrae N secuencias de un FASTA y genera un subset.
    """

    if not os.path.exists(fasta_path):
        return None, "Archivo FASTA de centroides no encontrado."

    try:
        registros = list(SeqIO.parse(fasta_path, "fasta"))

        if len(registros) == 0:
            return None, "El archivo FASTA está vacío."

        # Elegir mínimo entre solicitado y disponible
        subset = registros[:num_secuencias]

        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        with open(output_path, "w") as f_out:
            SeqIO.write(subset, f_out, "fasta")

        return output_path, f"{len(subset)} centroides extraídos correctamente."

    except Exception as e:
        return None, f"Error al extraer centroides: {e}"
