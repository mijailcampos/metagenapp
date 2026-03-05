import pandas as pd


def ordenar_centroides_por_muestra(
    count_table_path,
    salida="user_data/outputs/centroides_ordenados_por_muestra.txt"
):
    """
    Generate a text file containing centroid IDs sorted from highest to lowest
    abundance for each sample.

    This function reads a centroid count_table (OTU vs samples) and produces a 
    plain-text output where centroids are ordered independently for each sample, 
    allowing quick identification of the most abundant OTUs per sample.

    Parameters
    ----------
    count_table_path : str
        Path to the centroid count_table.
    salida : str
        Path where the sorted centroid list will be written.

    Returns
    -------
    str
        Path to the generated output file.
    """

    # Leer tabla de centroides
    df = pd.read_csv(count_table_path, sep="\t", index_col=0)

    # Crear archivo de salida (ordenado por muestra)
    with open(salida, "w") as f:
        for muestra in df.columns:
            # Encabezado por muestra
            f.write(f"# {muestra}\n")

            # Ordenar centroides por abundancia (descendente)
            ordenados = df[muestra].sort_values(ascending=False)

            # Escribir cada ID ordenado
            for centroid_id in ordenados.index:
                f.write(f"{centroid_id}\n")

    return salida
