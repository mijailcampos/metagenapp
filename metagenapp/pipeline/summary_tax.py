import pandas as pd


def summary_tax(taxonomy_path, count_table_path, nivel="Phylum"):
    try:
        # -------------------------
        # Load taxonomy
        # -------------------------
        tax_df = pd.read_csv(
            taxonomy_path,
            sep="\t",
            header=None,
            usecols=[0, 1],
            names=["ASV", "Taxonomy"],
            dtype=str,
            engine="python"
        )

        tax_df["ASV"] = tax_df["ASV"].fillna("").astype(str).str.strip()
        tax_df["Taxonomy"] = tax_df["Taxonomy"].fillna("").astype(str).str.strip()
        tax_df = tax_df[tax_df["ASV"] != ""].copy()

        # -------------------------
        # Load count table
        # -------------------------
        count_df = pd.read_csv(
            count_table_path,
            sep="\t",
            dtype={0: str}
        )

        id_col = count_df.columns[0]

        # -------------------------
        # Parse taxonomy safely
        # -------------------------
        def split_levels(tax_string):

            if not tax_string or tax_string == "Unclassified":
                return []

            levels = [x.strip().replace('"','') for x in tax_string.strip(";").split(";") if x.strip()]

            return levels

        tax_df["Levels"] = tax_df["Taxonomy"].apply(split_levels)

        # -------------------------
        # Map taxonomic rank
        # -------------------------
        rank_map = {
            "Kingdom": 0,
            "Phylum": 1,
            "Class": 2,
            "Order": 3,
            "Family": 4,
            "Genus": 5,
            "Species": 6
        }

        idx = rank_map.get(nivel, 1)

        # Extraer nivel taxonómico
        def get_level(levels):
            if len(levels) > idx:
                return levels[idx]
            else:
                return "Unclassified"

        tax_df["Level"] = tax_df["Levels"].apply(get_level)

        # -------------------------
        # Merge taxonomy + counts
        # -------------------------
        merged = pd.merge(
            tax_df[["ASV", "Level"]],
            count_df,
            left_on="ASV",
            right_on=id_col,
            how="inner"
        )

        # -------------------------
        # Summarize counts
        # -------------------------
        sample_cols = [c for c in count_df.columns if c != id_col]

        resumen = merged.groupby("Level")[sample_cols].sum().reset_index()
        resumen = resumen.rename(columns={"Level": nivel})

        return resumen, None

    except Exception as e:
        return None, str(e)