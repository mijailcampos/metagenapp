# ===============================================
# pipeline/heatmap.py — Biotech PRO Edition
# ===============================================

import pandas as pd
import numpy as np


# ============================================================
# 1. Load OTU/ASV Table
# ============================================================
def load_otu_table(path: str) -> pd.DataFrame:
    """
    Loads an OTU/ASV table (TSV) where:
        - Rows   = OTUs / ASVs
        - Columns = Samples

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by OTU/ASV IDs.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.set_index(df.columns[0])   # First column = OTU/ASV ID
    return df


# ============================================================
# 2. Load Mothur Taxonomy Table
# ============================================================
def load_taxonomy(path: str) -> pd.DataFrame:
    """
    Loads taxonomy in Mothur format:
        OTU \t Kingdom;Phylum;Class;Order;Family;Genus;Species

    Splits taxonomy into structured levels.

    Returns
    -------
    pd.DataFrame
        Taxonomic levels indexed by OTU ID.
    """

    tax = pd.read_csv(path, sep="\t", header=None, names=["OTU", "taxonomy"])
    tax = tax.set_index("OTU")

    # Split taxonomy string into hierarchical levels
    tax_levels = tax["taxonomy"].str.split(";", expand=True)

    level_names = ["Kingdom", "Phylum", "Class", "Order",
                   "Family", "Genus", "Species"]

    # Assign names based on available number of levels
    tax_levels.columns = level_names[:tax_levels.shape[1]]

    return tax_levels


# ============================================================
# 3. Collapse Abundance Table by Taxonomic Level
# ============================================================
def collapse_by_taxon(df_rel: pd.DataFrame,
                      tax_levels: pd.DataFrame,
                      nivel: str,
                      filtrar: bool):
    """
    Collapses relative abundance values by the selected taxonomic level.

    Parameters
    ----------
    df_rel : pd.DataFrame
        Relative abundance table (%), OTUs as rows and samples as columns.
    tax_levels : pd.DataFrame
        Taxonomy table indexed by OTU, with taxonomic ranks as columns.
    nivel : str
        Selected taxonomic rank (e.g., "Phylum", "Genus").
    filtrar : bool
        Whether to remove unassigned taxa.

    Returns
    -------
    pd.DataFrame
        Collapsed abundance table indexed by taxon.
    """

    df_rel = df_rel.copy()

    # Add selected taxon name to each OTU/ASV
    df_rel["taxon"] = tax_levels[nivel].fillna("Unassigned")

    # Optional: remove unassigned OTUs
    if filtrar:
        df_rel = df_rel[df_rel["taxon"] != "Unassigned"]

    # Collapse values by taxon (sum of % abundance)
    df_tax = df_rel.groupby("taxon").sum()

    return df_tax


# ============================================================
# 4. Select Top N Most Abundant Taxa
# ============================================================
def get_top_taxa(df_tax: pd.DataFrame, top_n: int):
    """
    Selects the N most abundant taxa based on mean relative abundance.

    Returns
    -------
    pd.DataFrame
        Table containing only the top-n most abundant taxa.
    """

    taxa = df_tax.mean(axis=1).sort_values(ascending=False).head(top_n).index
    return df_tax.loc[taxa]
