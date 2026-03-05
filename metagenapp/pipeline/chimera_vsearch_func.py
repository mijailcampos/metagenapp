import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def detectar_quimeras_vsearch(
    fasta_path,
    count_table_path,
    output_fasta,
    output_count,
    threads=1
):

    # -------------------------------------
    # Input validation
    # -------------------------------------
    if not os.path.exists(fasta_path):
        return None, None, "FASTA file not found."
    if not os.path.exists(count_table_path):
        return None, None, "Count table file not found."

    # =====================================
    # 1. Remove gaps for VSEARCH processing
    # =====================================

    nogap_fasta = os.path.splitext(output_fasta)[0] + "_no_gaps.fasta"

    try:
        with open(nogap_fasta, "w") as out_handle:
            for record in SeqIO.parse(fasta_path, "fasta"):
                seq_no_gaps = str(record.seq).replace("-", "")
                nuevo = SeqRecord(
                    Seq(seq_no_gaps),
                    id=record.id,
                    description=""
                )
                SeqIO.write(nuevo, out_handle, "fasta")

    except Exception as e:
        return None, None, f"Error while removing gaps: {e}"

    # =====================================
    # 2. Run VSEARCH (uchime_denovo)
    # =====================================

    uchime_out = os.path.splitext(output_fasta)[0] + ".uchime"

    cmd = [
        "vsearch",
        "--uchime_denovo", nogap_fasta,
        "--nonchimeras", output_fasta,
        "--uchimeout", uchime_out,
        "--threads", str(threads)
    ]

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )

        if result.stderr:
            print("⚠️ VSEARCH STDERR:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        return None, None, (
            f"VSEARCH failed (exit code {e.returncode}).\n"
            f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
        )
    except FileNotFoundError:
        return None, None, "VSEARCH executable was not found in PATH."

    # =====================================
    # 3. Read non-chimeric sequence IDs
    # =====================================
    if not os.path.exists(output_fasta):
        return None, None, "VSEARCH did not generate an output FASTA file."

    ids_no_chimera = {
        record.id for record in SeqIO.parse(output_fasta, "fasta")
    }

    # =====================================
    # 4. Filter count table
    # =====================================
    df = pd.read_csv(count_table_path, sep="\t")
    df_filtered = df[df.iloc[:, 0].isin(ids_no_chimera)]

    df_filtered.to_csv(output_count, sep="\t", index=False)

    return (
        output_fasta,
        output_count,
        f"Chimera detection completed. {len(ids_no_chimera)} non-chimeric sequences retained."
    )
