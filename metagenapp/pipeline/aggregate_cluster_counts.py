"""
Aggregate read counts from unique sequences to their OTU cluster centroids.

After vsearch clustering at 97%, each centroid represents a cluster of similar
sequences. This module uses the vsearch UC file to map every unique sequence to
its centroid, then sums the per-sample read counts from the unique_contigs
count_table. The result is a count table where each row is a centroid with the
correct total read counts — essential for accurate phylum-level abundance
reporting.
"""

import pandas as pd
from pathlib import Path


def aggregate_cluster_counts(
    uc_path,
    count_table_path,
    output_path,
):
    """
    Build a centroid-level count table by aggregating reads via the UC file.

    Args:
        uc_path:           Path to clustered_97.uc (vsearch clustering output)
        count_table_path:  Path to unique_contigs.count_table
        output_path:       Destination path for the aggregated count table

    Returns:
        (output_path_str, n_centroids, total_reads, error_msg)
        On failure: (None, 0, 0, error_msg)
    """
    try:
        # ── 1. Parse UC file: seq_id → centroid_id ──────────────────
        seq_to_centroid = {}
        with open(uc_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                typ = parts[0]
                if typ == "S":
                    # centroid maps to itself
                    seq_to_centroid[parts[8]] = parts[8]
                elif typ == "H":
                    # member maps to its centroid (column 9)
                    seq_to_centroid[parts[8]] = parts[9]

        # ── 2. Load unique-sequence count table ──────────────────────
        ct = pd.read_csv(count_table_path, sep="\t", index_col=0)
        sample_cols = ct.columns.tolist()

        # ── 3. Map each sequence to its centroid ─────────────────────
        ct["_centroid"] = ct.index.map(seq_to_centroid)
        unmapped = ct["_centroid"].isna().sum()
        if unmapped > 0:
            print(f"  ⚠️  {unmapped} sequences not found in UC file — skipped")
        ct = ct.dropna(subset=["_centroid"])

        # ── 4. Aggregate (sum) reads per centroid per sample ─────────
        agg = ct.groupby("_centroid")[sample_cols].sum()
        agg.index.name = "Sequence_ID"

        # ── 5. Write output ──────────────────────────────────────────
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        agg.to_csv(output_path, sep="\t")

        n_centroids = len(agg)
        total_reads = int(agg.values.sum())

        return str(output_path), n_centroids, total_reads, None

    except Exception as exc:
        return None, 0, 0, str(exc)
