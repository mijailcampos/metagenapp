import os
import platform
import subprocess
import random
from pathlib import Path
from Bio import SeqIO


def count_seqs(fasta_path):
    return sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))


def align_with_mafft(
    input_path,
    output_path,
    threads=8,
    max_full_seqs=50000,
    mode="full",
    logger=None
):
    """
    Align sequences using MAFFT.

    - Blocks full MAFFT for large datasets
    - Forces fast MAFFT mode
    - Enforces STUDENT mode hard limit
    """

    # ============================================================
    # 🧠 Normalizar tipos (Path SIEMPRE)
    # ============================================================
    input_path = Path(input_path)
    output_path = Path(output_path)

    threads_str = str(threads)

    # ============================================================
    # Validación de entrada
    # ============================================================
    if not input_path.exists():
        return None, f"Input FASTA not found: {input_path}"

    # ============================================================
    # Conteo de secuencias
    # ============================================================
    try:
        n_seqs = count_seqs(str(input_path))
    except Exception as e:
        return None, f"Failed to count sequences: {e}"

    if logger:
        logger.info(f"[MODE] {mode}")
        logger.info(f"Centroids detected: {n_seqs}")

    # ============================================================
    # Bloqueo de MAFFT completo en datasets grandes
    # ============================================================
    if n_seqs > max_full_seqs and mode != "student":
        return None, (
            f"MAFFT blocked: {n_seqs} sequences detected. "
            f"Full MAFFT is not viable."
        )

    # ============================================================
    # 🎓 STUDENT MODE — límite duro 10k
    # ============================================================
    if mode == "student" and n_seqs > 10_000:
        MAX_CENTROIDS = 10_000

        records = list(SeqIO.parse(str(input_path), "fasta"))
        records = random.sample(records, MAX_CENTROIDS)

        limited_input = input_path.with_name(
            input_path.stem + "_student10k.fasta"
        )
        SeqIO.write(records, str(limited_input), "fasta")

        input_path = limited_input
        n_seqs = MAX_CENTROIDS

        if logger:
            logger.info(f"Centroids used for MAFFT: {n_seqs}")

    # ============================================================
    # 🛑 Protección REAL basada en número de secuencias
    # ============================================================
    if mode == "student" and n_seqs > 10_000:
        raise RuntimeError(
            f"BUG: MAFFT would run on {n_seqs} sequences in STUDENT mode"
        )

    # ============================================================
    # Construcción del comando MAFFT
    # ============================================================
    system_type = platform.system()
    mafft_core_cmd = f"mafft --retree 2 --maxiterate 0 --thread {threads_str}"

    input_fasta_str = str(input_path)
    output_fasta_str = str(output_path)

    if system_type == "Windows":
        input_mafft = input_fasta_str.replace("D:\\", "/mnt/d/").replace("\\", "/")
        output_mafft = output_fasta_str.replace("D:\\", "/mnt/d/").replace("\\", "/")

        cmd = f"{mafft_core_cmd} '{input_mafft}' > '{output_mafft}'"
        mafft_command = ["wsl", "bash", "-c", cmd]
    else:
        mafft_command = [
            "bash",
            "-c",
            f"{mafft_core_cmd} '{input_fasta_str}' > '{output_fasta_str}'"
        ]

    # ============================================================
    # Ejecución
    # ============================================================
    try:
        subprocess.run(
            mafft_command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if not output_path.exists() or output_path.stat().st_size == 0:
            return None, "MAFFT finished but produced empty output. Aborting."

        return output_path, (
            f"MAFFT alignment completed successfully "
            f"({n_seqs} sequences, fast mode)."
        )

    except subprocess.CalledProcessError as e:
        return None, f"MAFFT failed:\n{e.stderr}"
