import os
import re
from Bio import SeqIO


# ========================================================
# 1. Convertir FASTA crudo → Formato SINTAX
# ========================================================
def format_reference_sintax(
    input_path,
    output_path
):
    """
    Convert a raw reference FASTA file into SINTAX format.

    Each sequence header is rewritten as:
        >ID tax=Kingdom:1.00,Phylum:1.00,...;

    This reformatted structure is required for SINTAX-compatible
    classifiers (e.g., VSEARCH) and ensures standardized taxonomy
    encoding across the reference database.
    """

    # Validación de existencia del archivo
    if not os.path.exists(input_path):
        return None, f"File not found: {input_path}"

    try:
        # Crear carpeta de salida si no existe
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        with open(output_path, "w") as out_f:
            for record in SeqIO.parse(input_path, "fasta"):

                # Extraer la cadena taxonómica original
                parts = record.description.split()
                tax_string = parts[1] if len(parts) > 1 else "Unclassified"

                # Limpiar caracteres potencialmente conflictivos
                tax_string = re.sub(r'["\'/]', '_', tax_string)

                # Descomponer niveles taxonómicos
                tax_levels = [t for t in tax_string.split(";") if t]

                # Añadir la puntuación SINTAX estándar (1.00)
                scored = ",".join([f"{t}:1.00" for t in tax_levels])

                # Escribir FASTA en formato SINTAX
                out_f.write(f">{record.id} tax={scored};\n{record.seq}\n")

        return output_path, "Reference successfully converted to SINTAX format."

    except Exception as e:
        return None, f"Error during conversion: {e}"


# ========================================================
# 2. Generar tabla taxonómica (.tax)
# ========================================================
def extract_taxonomy_table(
    fasta_sintax_path,
    tax_output_path
):
    """
    Extract a taxonomy table (ID → taxonomy) from a SINTAX-formatted FASTA.

    The output .tax file contains:
        SEQUENCE_ID \\t taxonomy_string
    """

    try:
        with open(tax_output_path, "w") as out:
            for record in SeqIO.parse(fasta_sintax_path, "fasta"):
                header = record.description

                # Extraer la porción de taxonomía desde el header SINTAX
                if "tax=" in header:
                    tax = header.split("tax=")[-1].strip("; ")
                else:
                    tax = "Unclassified"

                out.write(f"{record.id}\t{tax}\n")

        return tax_output_path

    except Exception:
        return None


# ========================================================
# 3. Recortar región (p. ej., 16S V4–V5)
# ========================================================
def trim_reference_region(
    input_fasta,
    output_fasta,
    start,
    end
):
    """
    Trim reference sequences between user-defined nucleotide positions.

    This operation allows extraction of specific rRNA hypervariable
    regions (e.g., V4, V4–V5) to match the target region amplified
    in the sequencing experiment.
    """

    try:
        # Crear carpeta de salida
        os.makedirs(os.path.dirname(output_fasta), exist_ok=True)

        with open(output_fasta, "w") as out:
            for record in SeqIO.parse(input_fasta, "fasta"):
                seq = str(record.seq)

                # Recortar región especificada por el usuario
                trimmed_seq = seq[start - 1:end]

                out.write(f">{record.id}\n{trimmed_seq}\n")

        return output_fasta

    except Exception:
        return None
