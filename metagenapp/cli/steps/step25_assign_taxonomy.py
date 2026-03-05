# cli/steps/step25_assign_taxonomy.py

import os

from metagenapp.pipeline.step_tracker import start_step, end_step
from metagenapp.pipeline.asignar_taxonomia_a_otus import asignar_taxonomia_a_otus


def run(
    count_input,
    taxonomy_input,
    outdir
):
    """
    Step 25 — Assign taxonomy to final ASVs
    """

    start_step("25_Assign_Taxonomy_to_ASVs")

    try:
        output_tax = os.path.join(outdir, "final_asv.taxonomy")

        resultado = asignar_taxonomia_a_otus(
            count_input,
            taxonomy_input,
            output_tax
        )

        if not resultado or not os.path.exists(resultado):
            raise RuntimeError("ASV taxonomy assignment failed")

        end_step(success=True)

        return {
            "taxonomy": resultado,
            "message": "Taxonomy successfully assigned to ASVs"
        }

    except Exception:
        end_step(success=False)
        raise
