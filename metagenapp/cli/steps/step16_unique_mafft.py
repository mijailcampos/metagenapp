from Bio import SeqIO
import os
from metagenapp.pipeline.step_tracker import start_step, end_step


def run(
    fasta_input,
    outdir,
    count_input=None,
    output_fasta_name="unique_mafft.fasta"
):
    start_step("16_Unique_MAFFT")

    try:
        os.makedirs(outdir, exist_ok=True)

        output_fasta = os.path.join(outdir, output_fasta_name)

        records = list(SeqIO.parse(fasta_input, "fasta"))
        SeqIO.write(records, output_fasta, "fasta")

        output_count = None
        if count_input:
            output_count = os.path.join(outdir, "unique_mafft.count_table")
            with open(count_input) as fin, open(output_count, "w") as fout:
                fout.write(fin.read())

        end_step(success=True)

        return {
            "fasta": output_fasta,
            "count": output_count,
            "n_unique": len(records)
        }

    except Exception as e:
        end_step(success=False)
        raise RuntimeError(f"Unique MAFFT step failed: {e}")
