import pandas as pd
from collections import defaultdict


def generar_count_table_de_centroides(uc_file, original_count_table, output_count_table):

    # Cargar tabla original de conteos
    df = pd.read_csv(original_count_table, sep="\t")

    # Ajuste de columna
    if "Representative_Sequence" in df.columns:
        df = df.rename(columns={"Representative_Sequence": "OTU"})

    df.set_index("OTU", inplace=True)

    # ----------------------
    # PARSEAR ARCHIVO .UC
    # ----------------------
    centroides = defaultdict(list)

    with open(uc_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            tipo = parts[0]            # S o H
            query = parts[8]           # ID de la secuencia query
            target = parts[9] if len(parts) > 9 else query

            if tipo == "S":
                centroides[query].append(query)

            elif tipo == "H":
                centroides[target].append(query)

    # ----------------------
    # SUMAR CUENTAS POR CENTROIDE
    # ----------------------
    nuevas_filas = []

    for centroide, miembros in centroides.items():

        miembros_validos = [m for m in miembros if m in df.index]

        if miembros_validos:
            suma = df.loc[miembros_validos].sum()
            suma.name = centroide
            nuevas_filas.append(suma)

    df_resultado = pd.DataFrame(nuevas_filas)
    df_resultado = df_resultado.reset_index().rename(columns={"index": "OTU"})

    df_resultado.to_csv(output_count_table, sep="\t", index=False)

    return output_count_table
