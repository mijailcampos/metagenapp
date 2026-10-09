import os
import subprocess
import pandas as pd
from Bio import SeqIO


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
