import pandas as pd
from metagenapp.pipeline.analizar_taxonomía import parse_taxonomy


def generar_tabla_resumen_asv(shared_path, taxonomy_path, output_resumen):
    """
    Generates the final ASV summary table containing:
    - Representative ASV IDs
    - Abundances across samples
    - Expanded taxonomy (phylum → genus levels)
    """

    # ======== 1. Load abundance table ========
    df_shared = pd.read_csv(shared_path, sep="\t")

    # Ensure first column is named Representative_Sequence
    if df_shared.columns[0] != "Representative_Sequence":
        df_shared.rename(
            columns={df_shared.columns[0]: "Representative_Sequence"}, inplace=True)

    # ======== 2. Load taxonomy ========
    df_tax = pd.read_csv(
        taxonomy_path,
        sep="\t",
        header=None,
        names=["ASV", "FullTaxonomy"]
    )

    # ======== 3. Merge ASVs with taxonomy ========
    df_merged = df_shared.merge(
        df_tax,
        how="left",
        left_on="Representative_Sequence",
        right_on="ASV"
    )

    # ======== 4. Expand sintax-style taxonomy ========
    expanded = df_merged["FullTaxonomy"].apply(parse_taxonomy)
    tax_df = pd.DataFrame(list(expanded))

    # ======== 5. Integrate taxonomy columns ========
    df_merged = pd.concat([df_merged, tax_df], axis=1)

    # ======== 6. Save final summary table ========
    df_merged.to_csv(output_resumen, sep="\t", index=False)

    return output_resumen
