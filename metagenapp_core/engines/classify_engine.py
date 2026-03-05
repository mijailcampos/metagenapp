from metagenapp_core.models.registry import MODELS


def classify_sequences(
    fasta_path,
    output_path,
    model="naive_v2",
    threads=4,
    block_size=2000,
):

    if model not in MODELS:
        raise ValueError(f"Unknown model: {model}")

    model_module = MODELS[model]

    return model_module.classify_naive_v2_parallel(
        fasta_path=fasta_path,
        output_path=output_path,
        n_threads=threads,
        block_size=block_size,
    )
