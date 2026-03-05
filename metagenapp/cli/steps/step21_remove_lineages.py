import os
from metagenapp.pipeline.remove_lineages import remove_lineage

def run(
    fasta_input,
    count_input,
    taxonomy_input,
    outdir,
    taxa_to_remove
):
    """
    Step 21 — Remove unwanted taxonomic lineages (CLI version)
    """

    output_fasta = os.path.join(outdir, "final_clean.fasta")
    output_count = os.path.join(outdir, "final_clean.count_table")
    output_tax   = os.path.join(outdir, "final_clean.taxonomy")

    fasta_out, count_out, tax_out, msg = remove_lineage(
        fasta_path=fasta_input,
        count_table_path=count_input,
        taxonomy_path=taxonomy_input,
        taxon_filter=taxa_to_remove,
        output_fasta=output_fasta,
        output_count=output_count,
        output_taxonomy=output_tax
    )

    if not fasta_out or not os.path.exists(fasta_out):
        raise RuntimeError(msg or "Lineage removal failed")

    return {
        "fasta": fasta_out,
        "count": count_out,
        "taxonomy": tax_out,
        "message": msg or "Unwanted lineages removed successfully"
    }
