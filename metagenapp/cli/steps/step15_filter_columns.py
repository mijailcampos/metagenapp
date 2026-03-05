from metagenapp.pipeline.filter_alignment_columns import filter_alignment_columns
from metagenapp.pipeline.step_tracker import start_step, end_step
import os


def run(fasta_input, outdir, vertical=True, trump="."):

    start_step("15_Filter_Columns")

    try:
        os.makedirs(outdir, exist_ok=True)
        output_fasta = os.path.join(
            outdir, "aligned_mafft_filtered_columns.fasta"
        )

        out_path, msg = filter_alignment_columns(
            fasta_path=fasta_input,
            output_path=output_fasta,
            vertical=vertical,
            trump=trump
        )

        if not out_path:
            raise RuntimeError(msg)

        end_step(success=True)

        return {
            "output_fasta": out_path,
            "message": msg
        }

    except Exception as e:
        end_step(success=False)
        raise e
