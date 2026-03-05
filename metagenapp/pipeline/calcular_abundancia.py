import pandas as pd
import os


def calcular_abundancia_relativa(df_merged, nivel_taxonomico):
    """
    Computes relative abundance per taxonomic group using the ASV summary table.
    The input should be the merged ASV table generated in Step 26.
    """

    # Standard taxonomy columns expected in the merged ASV table
    columnas_tax = [
        "OTU", "Kingdom", "Phylum",
        "Class", "Order", "Family",
        "Genus", "Species"
    ]

    # All columns that are NOT taxonomy-related → assumed to be sample abundance columns
    columnas_muestras = [c for c in df_merged.columns if c not in columnas_tax]

    # Ensure sample columns are numeric
    df_merged[columnas_muestras] = (
        df_merged[columnas_muestras]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
    )

    # Group abundances by the selected taxonomic level
    df_grouped = (
        df_merged.groupby(nivel_taxonomico)[columnas_muestras]
        .sum()
        .reset_index()
    )

    # Compute relative abundance (per sample)
    df_relative = df_grouped.copy()

    for sample in columnas_muestras:
        total = df_relative[sample].sum()
        df_relative[sample] = df_relative[sample] / total if total > 0 else 0

    # Output file path
    output_path = f"outputs/abundancia_relativa_{nivel_taxonomico}.tsv"
    df_relative.to_csv(output_path, sep="\t", index=False)

    return output_path
