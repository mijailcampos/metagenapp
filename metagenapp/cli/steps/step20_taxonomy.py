import os
from metagenapp.pipeline.clasificacion_tax import classify_naive_por_bloques


def run(
    fasta_input,
    count_input,
    outdir,
    model_path="modelos/naive_model.pkl",
    block_size=10000,
    n_threads=8
):
    try:
        output_taxonomy = os.path.join(outdir, "classification.taxonomy")

        taxonomy_path, msg = classify_naive_por_bloques(
            fasta_path=fasta_input,
            output_path=output_taxonomy,
            modelo_path=model_path,
            block_size=block_size,
            n_threads=n_threads
        )

        if not taxonomy_path or not os.path.exists(taxonomy_path):
            raise RuntimeError(msg or "Naive Bayes classification failed")

        return {
            "taxonomy": taxonomy_path,
            "message": msg or "Naive Bayes classification completed successfully."
        }

    except Exception:
        raise
