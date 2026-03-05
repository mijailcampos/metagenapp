import os
from metagenapp.pipeline.step_tracker import start_step, end_step
from metagenapp.pipeline.generar_otus_desde_asvs import generar_otus_desde_asvs


def run(
    fasta_input,
    count_input,
    outdir,
    identity=0.97,
    threads=8
):
    """
    Step 23 — Generate 97% OTUs from final ASVs
    """

    start_step("23_Generate_OTUs")

    try:
        otu_path, msg = generar_otus_desde_asvs(
            fasta_asv=fasta_input,
            count_table_asv=count_input,
            outdir=outdir,
            identity=identity,
            threads=threads
        )

        if not otu_path or not os.path.exists(otu_path):
            print("⚠ Warning:", msg or "OTU generation issue")
            print("⚠ Continuing pipeline despite missing OTU assignments")
            return {
                "otu_table": None,
                "message": msg or "OTU generation skipped"
            }

        end_step(success=True)

        return {
            "otu_table": otu_path,
            "message": msg or "97% OTUs generated successfully"
        }

    except Exception:
        end_step(success=False)
        raise
