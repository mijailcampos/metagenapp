from Bio import SeqIO


def diagnose_alignment(fasta_path, gap_threshold=0.5):
    sequences = list(SeqIO.parse(fasta_path, "fasta"))

    if not sequences:
        return None

    total_seqs = len(sequences)
    alignment_length = len(sequences[0].seq)

    # =============================
    # Mean ungapped length
    # =============================
    ungapped_lengths = [
        str(record.seq).count("-")
        for record in sequences
    ]

    mean_ungapped = sum(
        len(record.seq) - str(record.seq).count("-")
        for record in sequences
    ) / total_seqs

    # =============================
    # Gap ratio per column
    # =============================
    gap_ratios = []

    for i in range(alignment_length):
        gap_count = sum(1 for record in sequences if record.seq[i] == "-")
        gap_ratio = gap_count / total_seqs
        gap_ratios.append(gap_ratio)

    # =============================
    # Detect core region
    # =============================
    core_indices = [
        i for i, ratio in enumerate(gap_ratios)
        if ratio < gap_threshold
    ]

    if not core_indices:
        core_start = None
        core_end = None
        core_width = 0
    else:
        core_start = core_indices[0]
        core_end = core_indices[-1]
        core_width = core_end - core_start + 1

    return {
        "total_seqs": total_seqs,
        "alignment_length": alignment_length,
        "mean_ungapped": round(mean_ungapped, 2),
        "core_start": core_start,
        "core_end": core_end,
        "core_width": core_width
    }


def infer_region(mean_ungapped):
    if mean_ungapped < 300:
        return "Short 16S amplicon (likely V4 region)"
    elif 300 <= mean_ungapped < 450:
        return "Medium 16S amplicon (likely V3–V4 region)"
    elif 450 <= mean_ungapped < 650:
        return "Long 16S amplicon (V1–V4 / V3–V5 range)"
    else:
        return "Very long 16S fragment (check primers)"