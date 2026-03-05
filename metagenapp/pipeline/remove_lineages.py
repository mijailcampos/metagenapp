def remove_lineage(
    fasta_path,
    count_table_path,
    taxonomy_path,
    taxon_filter,
    output_fasta,
    output_count,
    output_taxonomy
):
    import pandas as pd
    from Bio import SeqIO

    try:
        # -------------------------
        # Load taxonomy
        # -------------------------
        tax_df = pd.read_csv(
            taxonomy_path,
            sep="\t",
            header=None,
            names=["ASV", "Taxonomy"],
            dtype=str
        )

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
        # Load FASTA
        # -------------------------
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

        # -------------------------
        # Parse taxa to remove
        # -------------------------
        taxa_to_remove = set(taxon_filter.split("-"))

        # -------------------------
        # Filter taxonomy
        # -------------------------
        mask_remove = tax_df["Taxonomy"].apply(
            lambda x: any(taxon in x for taxon in taxa_to_remove)
        )

        tax_df_filt = tax_df.loc[~mask_remove].copy()

        valid_ids = set(tax_df_filt["ASV"])

        # -------------------------
        # Filter count table (NO recalculation)
        # -------------------------
        count_df_filt = count_df[count_df[id_col].isin(valid_ids)].copy()

        # 🔒 CRITICAL: preserve original numeric values
        count_df_filt.reset_index(drop=True, inplace=True)

        # -------------------------
        # Filter FASTA
        # -------------------------
        fasta_filt = [
            seq_dict[asv_id]
            for asv_id in valid_ids
            if asv_id in seq_dict
        ]

        # -------------------------
        # Write outputs
        # -------------------------
        tax_df_filt.to_csv(
            output_taxonomy,
            sep="\t",
            index=False,
            header=False
        )

        count_df_filt.to_csv(
            output_count,
            sep="\t",
            index=False
        )

        with open(output_fasta, "w") as fh:
            SeqIO.write(fasta_filt, fh, "fasta")

        return output_fasta, output_count, output_taxonomy, None

    except Exception as e:
        return None, None, None, str(e)
