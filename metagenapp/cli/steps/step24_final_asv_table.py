# cli/steps/step24_final_asv_table.py

import os
import shutil
from metagenapp.pipeline.step_tracker import start_step, end_step


def run(
    count_input,
    outdir
):
    """
    Step 24 — Final ASV Table
    Registers and exports the final high-quality ASV count table.
    """

    start_step("24_Final_ASV_Table")

    try:
        if not os.path.exists(count_input):
            raise FileNotFoundError(
                f"Final ASV count table not found: {count_input}"
            )

        # Output path (standardized name)
        output_table = os.path.join(outdir, "final_asv_table.tsv")

        # Copy to final location
        shutil.copy(count_input, output_table)

        end_step(success=True)

        return {
            "asv_table": output_table,
            "message": "Final ASV table exported successfully"
        }

    except Exception:
        end_step(success=False)
        raise
