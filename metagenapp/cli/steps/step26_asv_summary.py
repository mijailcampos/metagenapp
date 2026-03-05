# cli/steps/step26_asv_summary.py

import os
from metagenapp.pipeline.step_tracker import start_step, end_step
from metagenapp.pipeline.func_tabla_resumen_asv import generar_tabla_resumen_asv


def run(
    shared_path,
    taxonomy_path,
    outdir,
    output_name="asv_resumen.tsv"
):
    """
    Step 26 — Generate ASV summary table (abundance + taxonomy)
    """

    start_step("26_ASV_Summary_Table")

    try:
        output_path = os.path.join(outdir, output_name)

        result = generar_tabla_resumen_asv(
            shared_path,
            taxonomy_path,
            output_path
        )

        if not result or not os.path.exists(result):
            raise RuntimeError("ASV summary table was not generated")

        end_step(success=True)

        return {
            "summary": result,
            "message": "ASV summary table generated successfully"
        }

    except Exception:
        end_step(success=False)
        raise
