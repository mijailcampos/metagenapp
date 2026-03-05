from metagenapp.pipeline.chimera_vsearch_func import detectar_quimeras_vsearch
from metagenapp.pipeline.step_tracker import start_step, end_step
import os


def run(
    fasta_input,
    count_input,
    outdir,
    output_fasta_name="non_chimeras.fasta",
    output_count_name="non_chimeras.count_table"
):
    start_step("18_ChimeraDetection")

    try:
        os.makedirs(outdir, exist_ok=True)

        output_fasta = os.path.join(outdir, output_fasta_name)
        output_count = os.path.join(outdir, output_count_name)

        out_fasta, out_count, msg = detectar_quimeras_vsearch(
            fasta_path=fasta_input,
            count_table_path=count_input,
            output_fasta=output_fasta,
            output_count=output_count
        )

        if not out_fasta or not os.path.exists(out_fasta):
            raise RuntimeError(msg or "Chimera detection failed")

        end_step(success=True)

        return {
            "fasta": out_fasta,
            "count": out_count,
            "message": msg
        }

    except Exception:
        end_step(success=False)
        raise
