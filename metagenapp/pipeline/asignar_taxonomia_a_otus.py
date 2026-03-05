def asignar_taxonomia_a_otus(shared_path, taxonomy_path, output_path):
    import pandas as pd

    # Cargar tabla ASV
    df_counts = pd.read_csv(shared_path, sep="\t")

    # Detectar nombre real de la columna de ASV
    posibles_ids = [
        "ASV",
        "Sequence_ID",
        "representative_sequence",
        "Representative_Sequence"
    ]
    
    id_col = None
    for c in df_counts.columns:
        if c in posibles_ids:
            id_col = c
            break

    if id_col is None:
        raise ValueError(
            f"No se encontró columna ASV en {shared_path}. Columnas disponibles: {list(df_counts.columns)}"
        )

    # Cargar taxonomía clasificada
    df_tax = pd.read_csv(
        taxonomy_path,
        sep="\t",
        header=None,
        names=["seq_id", "taxonomy"]
    )

    # Diccionario seq → taxonomy
    tax_dict = dict(zip(df_tax["seq_id"], df_tax["taxonomy"]))

    resultados = []

    for _, row in df_counts.iterrows():

        asv = row[id_col]

        # Obtener taxonomía
        tax = tax_dict.get(asv, "Unclassified")

        # Limpieza de comillas y dobles separadores
        tax = (
            tax.replace('""', '')   # elimina "" repetidas
               .replace('"', '')    # elimina comillas individuales
               .replace(";;", ";")  # arregla separadores
               .strip()
        )

        # Quitar ; final
        if tax.endswith(";"):
            tax = tax[:-1]

        resultados.append([asv, tax])

    # Guardar archivo final
    df_out = pd.DataFrame(resultados, columns=["ASV", "Taxonomía"])
    df_out.to_csv(output_path, sep="\t", index=False)

    return output_path
