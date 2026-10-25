import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from Bio import SeqIO

CHIMERA_MODES = ("pooled", "per-sample")


def remove_chimeras_from_uniques(
    unique_fasta,
    unique_count,
    output_fasta,
    output_count,
    workdir,
    threads=1,
):
    """
    De novo chimera removal on dereplicated (unique) sequences, BEFORE OTU
    clustering, as in the standard VSEARCH/UPARSE workflow.

    Each unique sequence is annotated with its total abundance across samples
    (;size=N, summed from the count table) so that UCHIME can use the abundance
    skew between chimeras and their parents. Sequences classified as chimeric
    ("Y") are removed from both the FASTA and the count table; non-chimeric
    ("N") and borderline ("?") sequences are kept.

    Returns (output_fasta, output_count, stats_dict).
    """
    os.makedirs(workdir, exist_ok=True)
    sized_fasta = os.path.join(workdir, "unique_contigs.sized.fasta")
    uchime_out = os.path.join(workdir, "unique_contigs.uchime")

    df = pd.read_csv(unique_count, sep="\t")
    id_col = df.columns[0]
    sizes = dict(zip(df[id_col], df.iloc[:, 1:].sum(axis=1).astype(int)))

    n_written = 0
    with open(sized_fasta, "w") as out:
        for rec in SeqIO.parse(unique_fasta, "fasta"):
            size = sizes.get(rec.id)
            if size is None:
                raise RuntimeError(f"{rec.id} missing from count table {unique_count}")
            seq = str(rec.seq).replace("-", "").replace(".", "").upper()
            out.write(f">{rec.id};size={size}\n{seq}\n")
            n_written += 1

    cmd = [
        "vsearch",
        "--uchime_denovo", sized_fasta,
        "--uchimeout", uchime_out,
        "--threads", str(threads),
        "--quiet",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"VSEARCH uchime_denovo failed:\n{result.stderr}")

    chimeric = set()
    n_classified = 0
    with open(uchime_out) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            n_classified += 1
            if fields[-1] == "Y":
                chimeric.add(fields[1].split(";size=")[0])

    if n_classified != n_written:
        raise RuntimeError(
            f"UCHIME classified {n_classified} sequences, expected {n_written}"
        )

    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(unique_fasta, "fasta"):
            if rec.id not in chimeric:
                SeqIO.write(rec, out, "fasta")

    keep = ~df[id_col].isin(chimeric)
    df[keep].to_csv(output_count, sep="\t", index=False)

    reads_total = int(sum(sizes.values()))
    reads_chimeric = int(sum(sizes[i] for i in chimeric))
    stats = {
        "uniques_in": n_written,
        "uniques_chimeric": len(chimeric),
        "reads_in": reads_total,
        "reads_chimeric": reads_chimeric,
        "reads_out": reads_total - reads_chimeric,
    }
    return str(output_fasta), str(output_count), stats


def remove_chimeras_per_sample(
    unique_fasta,
    unique_count,
    output_fasta,
    output_count,
    workdir,
    threads=1,
):
    """
    De novo chimera removal run separately in each sample, in parallel: one
    single-threaded `vsearch --uchime_denovo` per sample, up to `threads` at a
    time (uchime_denovo itself is sequential, so the pooled run uses one core).

    UCHIME only looks for parents among the uniques present in the same sample,
    with their abundance in that sample (as mothur chimera.vsearch
    dereplicate=t). A unique flagged "Y" in a sample loses its reads in that
    sample only; it is dropped from the FASTA when no sample keeps any read.

    Same outputs as remove_chimeras_from_uniques(); stats also report how many
    uniques were flagged in some samples but kept in others.
    """
    sample_dir = os.path.join(workdir, "chimera_per_sample")
    os.makedirs(sample_dir, exist_ok=True)

    df = pd.read_csv(unique_count, sep="\t")
    id_col = df.columns[0]
    samples = list(df.columns[1:])
    counts = df[samples].to_numpy(dtype=np.int64)
    row = {sid: i for i, sid in enumerate(df[id_col])}

    # one pass over the FASTA: each unique goes to the samples where it occurs
    sized = [os.path.join(sample_dir, f"sample_{j:04d}.sized.fasta") for j in range(len(samples))]
    handles = [open(p, "w") for p in sized]
    n_written = 0
    try:
        for rec in SeqIO.parse(unique_fasta, "fasta"):
            i = row.get(rec.id)
            if i is None:
                raise RuntimeError(f"{rec.id} missing from count table {unique_count}")
            seq = str(rec.seq).replace("-", "").replace(".", "").upper()
            for j in np.flatnonzero(counts[i]):
                handles[j].write(f">{rec.id};size={counts[i, j]}\n{seq}\n")
            n_written += 1
    finally:
        for h in handles:
            h.close()

    def run_sample(j):
        uchime_out = os.path.join(sample_dir, f"sample_{j:04d}.uchime")
        cmd = [
            "vsearch",
            "--uchime_denovo", sized[j],
            "--uchimeout", uchime_out,
            "--threads", "1",
            "--quiet",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"VSEARCH uchime_denovo failed in sample {samples[j]}:\n{result.stderr}")
        chimeric, n_classified = [], 0
        with open(uchime_out) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                n_classified += 1
                if fields[-1] == "Y":
                    chimeric.append(row[fields[1].split(";size=")[0]])
        expected = int(np.count_nonzero(counts[:, j]))
        if n_classified != expected:
            raise RuntimeError(
                f"UCHIME classified {n_classified} sequences in sample {samples[j]}, expected {expected}"
            )
        return j, chimeric

    kept = counts.copy()
    n_flags = 0
    with ThreadPoolExecutor(max_workers=max(1, min(int(threads), len(samples)))) as pool:
        for j, chimeric in pool.map(run_sample, range(len(samples))):
            kept[chimeric, j] = 0
            n_flags += len(chimeric)

    keep_rows = kept.sum(axis=1) > 0
    flagged_rows = (kept != counts).any(axis=1)
    keep_ids = set(df.loc[keep_rows, id_col])

    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(unique_fasta, "fasta"):
            if rec.id in keep_ids:
                SeqIO.write(rec, out, "fasta")

    df_out = df.copy()
    df_out[samples] = kept
    df_out[keep_rows].to_csv(output_count, sep="\t", index=False)

    reads_total = int(counts.sum())
    reads_out = int(kept.sum())
    stats = {
        "uniques_in": n_written,
        "uniques_chimeric": int((~keep_rows).sum()),
        "uniques_chimeric_some_samples": int((flagged_rows & keep_rows).sum()),
        "sample_flags": n_flags,
        "reads_in": reads_total,
        "reads_chimeric": reads_total - reads_out,
        "reads_out": reads_out,
    }
    return str(output_fasta), str(output_count), stats
