import pandas as pd
import numpy as np


def ejecutar_permanova(bray_file, metadata_file, grupo_column):
    """
    Manual PERMANOVA implementation (stable and scikit-bio–free version).

    This function performs a PERMANOVA test by:
      - Computing between-group and within-group sum of squares  
      - Permuting group labels 999 times  
      - Calculating an empirical p-value  
      - Returning F-statistic, p-value, R², and number of permutations

    Notes:
        * Bray–Curtis matrix must be square and aligned with metadata.
        * The metadata file must contain a column named 'sample'.
    """

    # ========================================================
    # 1. Load Bray–Curtis distance matrix
    # ========================================================
    df = pd.read_csv(bray_file, sep="\t", index_col=0)
    dist_matrix = df.values
    samples = df.index.tolist()

    # ========================================================
    # 2. Load metadata and align sample order
    # ========================================================
    meta = pd.read_csv(metadata_file, sep="\t")

    if "sample" not in meta.columns:
        raise ValueError("Metadata must contain a 'sample' column.")

    # Align metadata to match Bray–Curtis order
    meta = meta.set_index("sample").loc[samples]

    grupos = meta[grupo_column].values

    # ========================================================
    # 3. Sum-of-squares calculations
    # ========================================================
    def calcular_ss(dist_matrix, grupos):
        """
        Computes:
            - Total sum of squares
            - Within-group sum of squares
            - Between-group sum of squares
        from the distance matrix and group labels.
        """
        # Total SS (pairwise, distances counted once)
        ss_total = dist_matrix.sum() / 2

        # Within-group SS
        ss_within = 0
        for g in np.unique(grupos):
            idx = np.where(grupos == g)[0]
            if len(idx) > 1:
                sub = dist_matrix[np.ix_(idx, idx)]
                ss_within += sub.sum() / 2

        ss_between = ss_total - ss_within
        return ss_between, ss_within

    # Observed F statistic
    ss_b, ss_w = calcular_ss(dist_matrix, grupos)
    F_obs = ss_b / ss_w if ss_w > 0 else np.inf

    # ========================================================
    # 4. Permutation test (999 permutations)
    # ========================================================
    perm_count = 999
    greater = 0

    for _ in range(perm_count):
        perm = np.random.permutation(grupos)
        ss_b_p, ss_w_p = calcular_ss(dist_matrix, perm)
        F_perm = ss_b_p / ss_w_p if ss_w_p > 0 else np.inf
        if F_perm >= F_obs:
            greater += 1

    # Empirical p-value
    p_value = (greater + 1) / (perm_count + 1)

    # Effect size (R²)
    R2 = ss_b / (ss_b + ss_w)

    # ========================================================
    # 5. Results dictionary
    # ========================================================
    return {
        "F-statistic": float(F_obs),
        "p-value": float(p_value),
        "R2": float(R2),
        "permutations": perm_count
    }
