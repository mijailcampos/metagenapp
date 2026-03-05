from metagenapp.pipeline.centroid_subset import (
    extraer_subset_centroides,
    extraer_mas_abundantes_por_muestra
)

def extract_centroids_cli(
    clustered_fasta,
    clustered_count_table,
    mode="test",          # test | student | full
    top_n=None,
    outdir=None
):
    """
    Step 11 — Extract centroid subsets for downstream analysis
    """

    if not outdir:
        raise ValueError("outdir is required")

    if mode not in {"test", "student", "full"}:
        raise ValueError("mode must be: test | student | full")

    # -------------------------
    # 11.0 Subset general
    # -------------------------
    if mode == "test":
        n = 1000
    elif mode == "student":
        n = 10000
    else:
        n = None  # all

    subset_fasta = f"{outdir}/centroids_{mode}.fasta"

    path, msg = extraer_subset_centroides(
        fasta_path=clustered_fasta,
        output_path=subset_fasta,
        num_secuencias=n
    )

    # -------------------------
    # 11.1 Top-N por muestra
    # -------------------------
    top_fasta = None
    top_count = None

    if top_n is not None:
        top_fasta = f"{outdir}/centroids_top{top_n}.fasta"
        top_count = f"{outdir}/centroids_top{top_n}.count_table"

        fasta_file, count_file, total = extraer_mas_abundantes_por_muestra(
            count_table_path=clustered_count_table,
            fasta_path=clustered_fasta,
            top_n=top_n,
            output_fasta=top_fasta,
            output_count=top_count
        )

    return {
        "subset_fasta": path,
        "top_fasta": top_fasta,
        "top_count": top_count,
        "message": msg
    }
