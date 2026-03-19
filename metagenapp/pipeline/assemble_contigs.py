import subprocess
import os
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path

from metagenapp.pipeline.step_tracker import start_step, end_step

try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False


def assemble_contigs(
    files_path,
    input_dir,
    output_fasta,
    output_count,
):
    """
    Assemble contigs from paired-end FASTQ files using VSEARCH mergepairs.
    Produces merged FASTQ, converts to FASTA, and generates count_table.
    """

    print("🧬 MetagenApp pipeline started", flush=True)

    files_path = Path(files_path)
    input_dir = Path(input_dir)
    output_fasta = Path(output_fasta)
    output_count = Path(output_count)

    print(f"📂 Input directory: {input_dir}", flush=True)
    print(f"📄 Files list: {files_path}", flush=True)
    print(f"📤 Output FASTA: {output_fasta}", flush=True)
    print(f"📤 Output count table: {output_count}", flush=True)

    if not files_path.exists():
        return None, f"File not found: {files_path}"

    count_table_data = defaultdict(lambda: defaultdict(int))
    merged_fastqs = []

    print("🔬 Starting VSEARCH mergepairs...", flush=True)

    try:
        with open(files_path) as f:
            for line in f:
                sample, r1_file, r2_file = line.strip().split("\t")

                r1_path = input_dir / r1_file
                r2_path = input_dir / r2_file

                if not r1_path.exists() or not r2_path.exists():
                    print(f"[WARNING] Missing FASTQ for {sample}", flush=True)
                    continue

                merged_fastq = output_fasta.parent / f"{sample}_merged.fastq"

                print(f"   ↳ Merging sample: {sample}", flush=True)

                cmd = [
                    "vsearch",
                    "--fastq_mergepairs", str(r1_path),
                    "--reverse", str(r2_path),
                    "--fastqout", str(merged_fastq),
                    "--fastq_minovlen", "20",
                    "--fastq_maxdiffs", "5",
                ]

                subprocess.run(cmd, check=True)

                merged_fastqs.append((sample, merged_fastq))

    except Exception as e:
        return None, f"Error during VSEARCH mergepairs: {e}"

    print(f"✔ VSEARCH merging finished ({len(merged_fastqs)} samples)", flush=True)

    print("📦 Parsing merged FASTQ files...", flush=True)

    try:
        for i, (sample, merged_fastq) in enumerate(merged_fastqs, start=1):

            print(f"   ↳ Reading {merged_fastq.name} ({i}/{len(merged_fastqs)})", flush=True)

            seq_counter = 0

            for record in SeqIO.parse(merged_fastq, "fastq"):
                seq = str(record.seq)
                count_table_data[seq][sample] += 1
                seq_counter += 1

            print(f"      {seq_counter} reads processed", flush=True)

    except Exception as e:
        return None, f"Error parsing merged FASTQ: {e}"

    print("✔ FASTQ parsing finished", flush=True)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    print("🧮 Building count table...", flush=True)

    all_samples = sorted({
        sample
        for seq_counts in count_table_data.values()
        for sample in seq_counts
    })

    print(f"✔ Total unique sequences: {len(count_table_data)}", flush=True)
    print(f"✔ Total samples detected: {len(all_samples)}", flush=True)

    print("💾 Writing assembled FASTA and count table...", flush=True)

    with open(output_fasta, "w") as f_out, open(output_count, "w") as f_count:

        f_count.write("ContigID\t" + "\t".join(all_samples) + "\n")

        for i, (seq, sample_counts) in enumerate(count_table_data.items(), start=1):
            contig_id = f"contig_{i}"

            f_out.write(f">{contig_id}\n{seq}\n")

            row = [str(int(sample_counts.get(s, 0))) for s in all_samples]
            f_count.write(f"{contig_id}\t" + "\t".join(row) + "\n")

            if i % 1000 == 0:
                print(f"   ↳ {i} contigs written", flush=True)

    print("✔ Assemble contigs finished", flush=True)

    return output_fasta, None