def generar_count_table_de_centroides(uc_file, original_count_table, output_count_table):
    import pandas as pd
    from collections import defaultdict

    # Cargar count_table original
    df = pd.read_csv(original_count_table, sep="\t")

    # --- Detectar columna que contiene IDs de secuencia ---
    posibles_columnas = ["OTU", "sequence", "Sequence", "representative",
                         "SeqID", "ID", "id"]

    columna_id = None
    for col in df.columns:
        if col in posibles_columnas:
            columna_id = col
            break

    # Si no encontramos ninguna, asumimos que la PRIMERA columna es la de IDs
    if columna_id is None:
        columna_id = df.columns[0]

    # Renombrar a OTU universalmente
    df = df.rename(columns={columna_id: "OTU"})
    df.set_index("OTU", inplace=True)

    # --- Leer archivo .uc ---
    centroides = defaultdict(list)

    with open(uc_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            tipo = parts[0]
            query = parts[8]
            target = parts[9] if len(parts) > 9 else query

            if tipo == "S":  # Secuencia semilla (centroide)
                centroides[query].append(query)
            elif tipo == "H":  # Secuencia hit
                centroides[target].append(query)

    # --- Construir nueva count_table ---
    nuevas_filas = []

    for centroide, miembros in centroides.items():
        miembros_validos = [m for m in miembros if m in df.index]
        if miembros_validos:
            suma = df.loc[miembros_validos].sum()
            suma.name = centroide
            nuevas_filas.append(suma)

    df_out = pd.DataFrame(nuevas_filas)
    df_out = df_out.reset_index().rename(columns={"index": "OTU"})
    df_out.to_csv(output_count_table, sep="\t", index=False)

    return output_count_table
