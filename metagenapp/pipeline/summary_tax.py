import pandas as pd


def summary_tax(taxonomy_path, count_table_path, nivel="Phylum"):
    """
    Generates a taxonomic summary table for a chosen taxonomic rank.
    Works with taxonomies from Mothur, VSEARCH SINTAX, PR2, SILVA,
    or Naive Bayes models. Automatically detects delimiters and
    number of levels.
    """

    try:
        # -------------------------
        # 1. Load taxonomy file
        # -------------------------
        tax_df = pd.read_csv(
            taxonomy_path, sep="\t", header=None, names=["OTU", "Taxonomy"]
        )

        # -------------------------
        # 2. Load count table
        # -------------------------
        count_df = pd.read_csv(count_table_path, sep="\t")
        id_column = count_df.columns[0]

        # Standardize IDs
        tax_df["OTU"] = tax_df["OTU"].astype(str).str.strip()
        count_df[id_column] = count_df[id_column].astype(str).str.strip()

        # -------------------------
        # 3. Compute abundance per OTU
        # -------------------------
        count_df["Total"] = count_df.drop(columns=[id_column]).sum(axis=1)
        abundances = dict(zip(count_df[id_column], count_df["Total"]))
        tax_df["Abundance"] = tax_df["OTU"].map(abundances)

        # -------------------------
        # 4. Detect delimiter and split taxonomy
        # -------------------------
        tax_df["delim"] = tax_df["Taxonomy"].apply(
            lambda x: "," if "," in str(x) else ";"
        )

        def split_levels(row):
            raw = str(row["Taxonomy"])
            delim = row["delim"]

            parts = [p.strip() for p in raw.split(delim) if p.strip() != ""]

            # Remove prefixes like "k__", "p__", "g__", etc.
            parts = [p.split("__")[-1] for p in parts]

            return parts

        tax_df["Levels"] = tax_df.apply(split_levels, axis=1)

        # -------------------------
        # 5. Define taxonomic order dynamically
        # -------------------------
        taxonomic_order = ["Kingdom", "Phylum",
                           "Class", "Order", "Family", "Genus"]

        # How many levels does the taxonomy actually have?
        max_levels = tax_df["Levels"].apply(len).max()

        # Extend missing levels with "Unclassified"
        def normalize_levels(levels):
            levels = list(levels)
            while len(levels) < max_levels:
                levels.append("Unclassified")
            return levels

        tax_df["Levels"] = tax_df["Levels"].apply(normalize_levels)

        # -------------------------
        # 6. Select correct index for requested level
        # -------------------------
        if nivel not in taxonomic_order:
            idx = 1  # default to Phylum if unknown
        else:
            idx = taxonomic_order.index(nivel)

        # If the requested index exceeds available levels → use last
        idx = min(idx, max_levels - 1)

        tax_df["Level"] = tax_df["Levels"].apply(lambda x: x[idx])

        # -------------------------
        # 7. Summarize abundances
        # -------------------------
        summary = (
            tax_df.groupby("Level")["Abundance"]
            .sum()
            .sort_values(ascending=False)
            .reset_index()
        )

        return summary, None

    except Exception as e:
        return None, str(e)
