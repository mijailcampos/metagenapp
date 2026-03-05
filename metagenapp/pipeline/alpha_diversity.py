# pipeline/alpha_diversity.py
import pandas as pd
import numpy as np


def calcular_alpha_diversidad(otu_table_path: str, sep: str = "\t") -> pd.DataFrame:
    """
    Computes alpha diversity metrics for each sample using the final ASV/OTU table.

    Metrics included:
        • Richness (number of detected ASVs)
        • Shannon index
        • Simpson index
        • Chao1 richness estimator

    Parameters
    ----------
    otu_table_path : str
        Path to the final count table (final_clean.count_table)

    sep : str, optional
        Column separator used in the input file. Default is '\t'.

    Returns
    -------
    pandas.DataFrame
        A table with alpha diversity metrics per sample.
    """

    # Load OTU/ASV table
    otu_df = pd.read_csv(otu_table_path, sep=sep)

    # Detect numeric sample columns
    numeric_cols = otu_df.select_dtypes(include=[np.number]).columns.tolist()

    # Fallback: if numeric detection fails, assume all columns except the first
    if not numeric_cols:
        numeric_cols = otu_df.columns[1:]

    counts = otu_df[numeric_cols]

    resultados = []

    for sample in numeric_cols:

        vec = counts[sample].to_numpy(dtype=float)
        N = vec.sum()                       # total abundance
        riqueza = np.count_nonzero(vec)     # ASV richness

        if N > 0 and riqueza > 0:

            v_pos = vec[vec > 0]
            p = v_pos / N                   # relative abundance

            # Shannon diversity
            shannon = -(p * np.log(p)).sum()

            # Simpson diversity (1 - D)
            simpson = 1.0 - (p ** 2).sum()

            # Singletons and doubletons for Chao1
            f1 = np.sum(vec == 1)
            f2 = np.sum(vec == 2)

            if f2 > 0:
                chao1 = riqueza + (f1 ** 2) / (2.0 * f2)
            else:
                chao1 = np.nan if f1 == 0 else riqueza + (f1 * (f1 - 1)) / 2.0

        else:
            shannon = np.nan
            simpson = np.nan
            chao1 = np.nan

        resultados.append(
            {
                "muestra": sample,
                "riqueza_asvs": int(riqueza),
                "shannon": shannon,
                "simpson": simpson,
                "chao1": chao1,
            }
        )

    alpha_df = pd.DataFrame(resultados)
    alpha_df.attrs["otu_table_path"] = otu_table_path

    return alpha_df
