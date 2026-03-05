import os
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path

from metagenapp.pipeline.step_tracker import start_step, end_step

# --------------------------------------------------
# Optional Streamlit support (GUI-safe / CLI-safe)
# --------------------------------------------------
try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False


# ==================================================
#   CORE FUNCTION (GUI + CLI)
# ==================================================
def assemble_contigs(
    files_path,
    input_dir,
    output_fasta,
    output_count,
):
    """
    Assemble contigs from paired-end FASTQ files (R1/R2)
    using simple concatenation:
        read1 + reverse_complement(read2)

    Parameters
    ----------
    files_path : str or Path
        Path to input_samples.files
    input_dir : str or Path
        Directory containing FASTQ files
    output_fasta : str or Path
        Output FASTA path
    output_count : str or Path
        Output count table path

    Returns
    -------
    output_fasta, error
    """

    files_path = Path(files_path)
    input_dir = Path(input_dir)
    output_fasta = Path(output_fasta)
    output_count = Path(output_count)

    if not files_path.exists():
        return None, f"File not found: {files_path}"

    count_table_data = defaultdict(lambda: defaultdict(int))

    try:
        with open(files_path) as f:
            for line in f:
                sample, r1_file, r2_file = line.strip().split("\t")

                r1_path = input_dir / r1_file
                r2_path = input_dir / r2_file

                if not r1_path.exists() or not r2_path.exists():
                    print(f"[WARNING] Missing FASTQ for {sample}")
                    continue

                r1_reads = SeqIO.parse(r1_path, "fastq")
                r2_reads = SeqIO.parse(r2_path, "fastq")

                for read1, read2 in zip(r1_reads, r2_reads):
                    contig_seq = str(
                        read1.seq + read2.seq.reverse_complement()
                    )
                    count_table_data[contig_seq][sample] += 1

    except Exception as e:
        return None, f"Error processing FASTQ files: {e}"

    # Ensure output directory exists
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    # Consistent sample order
    all_samples = sorted({
        sample
        for seq_counts in count_table_data.values()
        for sample in seq_counts
    })

    with open(output_fasta, "w") as f_out, open(output_count, "w") as f_count:

        f_count.write("ContigID\t" + "\t".join(all_samples) + "\n")

        for i, (seq, sample_counts) in enumerate(count_table_data.items(), start=1):
            contig_id = f"contig_{i}"
            f_out.write(f">{contig_id}\n{seq}\n")

            row = [str(int(sample_counts.get(s, 0))) for s in all_samples]
            f_count.write(f"{contig_id}\t" + "\t".join(row) + "\n")

    return output_fasta, None


# ==================================================
#   STREAMLIT WRAPPER (GUI ONLY)
# ==================================================
def run_assemble_contigs():
    """
    Streamlit wrapper for assemble_contigs().
    Not used by CLI.
    """

    if not STREAMLIT_AVAILABLE:
        raise RuntimeError("run_assemble_contigs() requires Streamlit")

    st.subheader("Assembling contigs")

    start_step("01_Assemble_contigs")

    # GUI paths (legacy behavior)
    input_dir = Path("user_data/inputs")
    files_path = input_dir / "input_samples.files"
    output_fasta = Path("user_data/outputs/assembled_contigs.fasta")
    output_count = Path("user_data/outputs/assembled_contigs.count_table")

    with st.spinner("Processing R1/R2 pairs and building contigs..."):
        output, error = assemble_contigs(
            files_path=files_path,
            input_dir=input_dir,
            output_fasta=output_fasta,
            output_count=output_count,
        )

    if output:
        end_step(success=True)

        st.success("Contigs assembled successfully.")
        st.code(f"FASTA: {output_fasta}\nCOUNT TABLE: {output_count}")

        try:
            with open(output_fasta) as f:
                preview = "".join(f.readlines()[:10])
                st.text_area("FASTA preview:", preview, height=200)
        except Exception:
            pass
    else:
        end_step(success=False)
        st.error(f"Error assembling contigs:\n{error}")
