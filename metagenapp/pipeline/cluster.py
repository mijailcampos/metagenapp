import os
import subprocess
import multiprocessing


def cluster_seqs(
    fasta_path,
    output_fasta,
    output_uc,
    identity=0.97,
    threads=None,
):
    """
    Cluster trimmed sequences using VSEARCH (cluster_fast).

    Steps
    -----
    1. Sort sequences by length (required by VSEARCH).
    2. Perform greedy clustering at given identity threshold.
    3. Output centroid FASTA and .uc assignment file.

    Parameters
    ----------
    fasta_path : str
        Input FASTA with trimmed, aligned sequences.
    output_fasta : str
        Output FASTA with centroid sequences.
    output_uc : str
        Output UC file with cluster assignments.
    identity : float
        Identity threshold (default: 0.97).
    threads : int or None
        Number of threads (default: all available).

    Returns
    -------
    tuple
        (output_fasta, message) on success
        (None, error_message) on failure
    """

    if not os.path.exists(fasta_path):
        return None, f"FASTA not found: {fasta_path}"

    if threads is None:
        threads = multiprocessing.cpu_count()

    # Temporary sorted FASTA
    sorted_fasta = fasta_path.with_name(
        fasta_path.stem + "_sorted" + fasta_path.suffix
    )


    try:
        # --------------------------------------------------
        # 1. Sort sequences by length
        # --------------------------------------------------
        cmd_sort = [
            "vsearch",
            "--sortbylength", str(fasta_path),
            "--output", str(sorted_fasta),
            "--minseqlength", "1",
            "--threads", str(threads)
        ]

        result_sort = subprocess.run(
            cmd_sort,
            capture_output=True,
            text=True
        )

        if result_sort.returncode != 0:
            return None, f"VSEARCH sort error:\n{result_sort.stderr}"

        # --------------------------------------------------
        # 2. Clustering (cluster_fast)
        # --------------------------------------------------
        cmd_cluster = [
            "vsearch",
            "--cluster_fast", sorted_fasta,
            "--id", str(identity),
            "--centroids", output_fasta,
            "--uc", output_uc,
            "--threads", str(threads)
        ]

        result_cluster = subprocess.run(
            cmd_cluster,
            capture_output=True,
            text=True
        )

        if result_cluster.returncode != 0:
            return None, f"VSEARCH clustering error:\n{result_cluster.stderr}"

        # --------------------------------------------------
        # 3. Basic statistics
        # --------------------------------------------------
        num_in = sum(1 for line in open(fasta_path) if line.startswith(">"))
        num_centroids = sum(
            1 for line in open(output_fasta) if line.startswith(">")
        )

        if num_in == 0:
            return None, "Input FASTA is empty."

        reduction = 100 - (num_centroids / num_in * 100)

        msg = (
            "VSEARCH clustering completed successfully.\n"
            f"Input sequences: {num_in}\n"
            f"Centroids (OTUs): {num_centroids}\n"
            f"Identity threshold: {identity}\n"
            f"Reduction: {reduction:.1f}%"
        )

        return output_fasta, msg

    except Exception as e:
        return None, f"Unexpected clustering error: {e}"

    finally:
        # Clean temporary file
        if os.path.exists(sorted_fasta):
            try:
                os.remove(sorted_fasta)
            except Exception:
                pass
