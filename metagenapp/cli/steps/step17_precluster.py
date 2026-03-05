from metagenapp.pipeline.precluster_parallel import pre_cluster_seqs_parallel
from metagenapp.pipeline.step_tracker import start_step, end_step
import os
from multiprocessing import cpu_count


def run(
    fasta_input,
    count_input,
    outdir,
    diffs=2,
    n_cpus=None,
    output_fasta_name="preclustered.fasta",
    output_count_name="preclustered.count_table"
):
    start_step("17_Precluster")

    try:
        os.makedirs(outdir, exist_ok=True)

        if n_cpus is None:
            n_cpus = max(1, cpu_count() - 2)

        output_fasta = os.path.join(outdir, output_fasta_name)
        output_count = os.path.join(outdir, output_count_name)

        fasta_out, count_out, msg = pre_cluster_seqs_parallel(
            fasta_path=fasta_input,
            count_table_path=count_input,
            diffs=diffs,
            output_fasta=output_fasta,
            output_count=output_count,
            n_cpus=n_cpus
        )

        if not fasta_out or not os.path.exists(fasta_out):
            raise RuntimeError(msg or "Precluster failed")

        end_step(success=True)

        return {
            "fasta": fasta_out,
            "count": count_out,
            "message": msg
        }

    except Exception:
        end_step(success=False)
        raise
