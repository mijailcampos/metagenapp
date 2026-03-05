import os
from metagenapp.pipeline.step_tracker import start_step, end_step
from metagenapp.pipeline.summary_tax import summary_tax

def run(
    taxonomy_input,
    count_input,
    outdir,
    nivel="Phylum"
):
    start_step("22_TaxonomicSummary")

    try:
        output_summary = os.path.join(
            outdir,
            f"summary_tax_{nivel.lower()}.tsv"
        )

        resumen, error = summary_tax(
            taxonomy_path=taxonomy_input,
            count_table_path=count_input,
            nivel=nivel
        )

        if resumen is None:
            raise RuntimeError(error or "Taxonomic summary failed")

        resumen.to_csv(output_summary, sep="\t", index=False)

        end_step(success=True)

        return {
            "summary": output_summary,
            "message": f"Taxonomic summary generated at level {nivel}"
        }

    except Exception:
        end_step(success=False)
        raise
