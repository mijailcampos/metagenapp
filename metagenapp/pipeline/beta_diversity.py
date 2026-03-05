# pipeline/beta_diversity.py

import pandas as pd
import numpy as np
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS


def calcular_beta_diversidad(otu_table_path: str, sep: str = "\t"):
    """
    Computes beta diversity metrics:
      • Bray–Curtis distance matrix
      • PCoA coordinates using metric MDS

    Parameters
    ----------
    otu_table_path : str
        Path to the final ASV count table (final_clean.count_table).
    sep : str
        Field separator used in the input file.

    Returns
    -------
    dist_matrix_df : pd.DataFrame
        Square Bray–Curtis distance matrix (samples × samples).
    pcoa_df : pd.DataFrame
        PCoA coordinates for each sample (PCoA1, PCoA2).
    """

    # ----------------------------------------
    # Load OTU/ASV table
    # ----------------------------------------
    otu_df = pd.read_csv(otu_table_path, sep=sep)

    # Identify numeric columns (samples)
    numeric_cols = otu_df.select_dtypes(include=[np.number]).columns.tolist()
    if not numeric_cols:
        numeric_cols = otu_df.columns[1:]

    # Build sample × feature matrix
    data = otu_df[numeric_cols].T
    data.index.name = "muestra"

    # ----------------------------------------
    # Bray–Curtis distance matrix
    # ----------------------------------------
    dist_condensed = pdist(data.values, metric="braycurtis")
    dist_matrix = squareform(dist_condensed)

    dist_matrix_df = pd.DataFrame(
        dist_matrix,
        index=data.index,
        columns=data.index
    )

    # ----------------------------------------
    # PCoA (classical) using scikit-bio
    # ----------------------------------------
    dm = DistanceMatrix(
        dist_matrix_df.values,
        ids=dist_matrix_df.index.tolist()
    )

    pcoa_res = pcoa(dm)

    pcoa_df = pcoa_res.samples.iloc[:, :2].copy()
    pcoa_df.columns = ["PCoA1", "PCoA2"]
    pcoa_df.index.name = "muestra"

    # ✅ ESTA LÍNEA HABILITA LOS %
    pcoa_df.attrs["variance"] = pcoa_res.proportion_explained.values.tolist()

    return dist_matrix_df, pcoa_df
