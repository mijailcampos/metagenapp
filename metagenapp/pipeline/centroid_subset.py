import os
import pandas as pd
from Bio import SeqIO


# ============================================================
# extraer_subset_centroides
# (Extrae un subconjunto N de centroides en orden)
# ============================================================
def extraer_subset_centroides(
    fasta_path,
    output_path,
    num_secuencias=1000
):
    """
    Extracts a fixed-size subset of centroid sequences.

    Notes for the user:
        - Sequences are taken in their original order.
        - Useful for test or student modes to reduce dataset size.
    """

    # Verificar existencia del archivo
    if not os.path.exists(fasta_path):
        return None, "Centroid FASTA file not found."

    try:
        # Leer todos los centroides
        registros = list(SeqIO.parse(fasta_path, "fasta"))

        # Tomar los primeros N
        subset = registros[:num_secuencias]

        # Crear carpeta si no existe
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Guardar subset
        SeqIO.write(subset, output_path, "fasta")

        if not subset:
            return None, "No sequences were extracted."

        return output_path, f"{len(subset)} centroids extracted successfully."

    except Exception as e:
        return None, f"Error extracting subset: {e}"


# ============================================================
# extraer_mas_abundantes_por_muestra
# (Selecciona los centroides más abundantes por muestra)
# ============================================================
def extraer_mas_abundantes_por_muestra(
    count_table_path,
    fasta_path,
    top_n,
    output_fasta,
    output_count
):
    """
    Extracts the most abundant centroids per sample and
    generates both FASTA and count_table subsets.

    Notes for the user:
        - For each sample, the top N centroid IDs are selected.
        - The union of all selected IDs is used to build FASTA + count_table.
    """

    # Leer count_table
    df = pd.read_csv(count_table_path, sep="\t")
    df.set_index("OTU", inplace=True)

    # Selección de los más abundantes por muestra
    otus_top = set()
    for muestra in df.columns:
        top = df[muestra].sort_values(ascending=False).head(top_n).index
        otus_top.update(top)

    # Convertir FASTA a diccionario rápido
    secuencias = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    # Extraer solo los centroides seleccionados
    subset_seqs = [secuencias[otu] for otu in otus_top if otu in secuencias]
    SeqIO.write(subset_seqs, output_fasta, "fasta")

    # Filtrar count_table para los centroides top
    df_filtrado = df.loc[df.index.isin(otus_top)]
    df_filtrado.reset_index().to_csv(output_count, sep="\t", index=False)

    return output_fasta, output_count, len(subset_seqs)


# ============================================================
# contar_centroides_por_muestra
# (Cuenta centroides no nulos por muestra)
# ============================================================
def contar_centroides_por_muestra(count_table_path):
    """
    Counts how many centroid sequences appear with abundance >0 in each sample.

    Notes for the user:
        - Useful for visualizing sample complexity.
        - Output is sorted from highest to lowest.
    """

    df = pd.read_csv(count_table_path, sep="\t")

    # Verificar columna OTU
    if "OTU" not in df.columns:
        raise ValueError("The count_table does not contain column 'OTU'.")

    df.set_index("OTU", inplace=True)

    # Conteo de centroides con abundancia > 0
    resumen = (df > 0).sum()

    return resumen.sort_values(ascending=False)
