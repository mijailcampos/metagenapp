from pathlib import Path

from metagenapp.metagen_config import (
    PR2_REFERENCE,
    SILVA_REFERENCE_RAW,
    SILVA_REFERENCE_ALN
)

# ===============================
# Pipeline core
# ===============================
from metagenapp.pipeline.generate_files import generate_input_files
from metagenapp.pipeline.assemble_contigs import assemble_contigs
from metagenapp.pipeline.filter_contigs import filter_contigs
from metagenapp.pipeline.unique_contigs import generate_unique_contigs
from metagenapp.pipeline.alignment import align_centroids_refmode_professional
from metagenapp.pipeline.recorte_vsearch import recortar_por_alnout
from metagenapp.pipeline.cluster import cluster_seqs
from metagenapp.pipeline.centroid_subset import extraer_subset_centroides
from metagenapp.pipeline.align_with_mafft import align_with_mafft
from metagenapp.pipeline.filtrar_mafft_avanzado import filtrar_mafft_por_posicion
from metagenapp.pipeline.step_tracker import start_step, end_step
from metagenapp.pipeline.alignment_diagnostics import (
    diagnose_alignment,
    infer_region
)
from metagenapp.pipeline.alignment import align_with_vsearch

# ===============================
# CLI steps
# ===============================
from metagenapp.cli.steps.step15_filter_columns import run as run_step15
from metagenapp.cli.steps.step16_unique_mafft import run as run_step16
from metagenapp.cli.steps.step17_precluster import run as run_step17
from metagenapp.cli.steps.step18_chimera import run as run_step18
from metagenapp.cli.steps.step20_taxonomy import run as run_step20
from metagenapp.cli.steps.step21_remove_lineages import run as run_step21
from metagenapp.cli.steps.step22_summary_tax import run as run_step22
from metagenapp.cli.steps.step23_generate_otus import run as run_step23
from metagenapp.cli.steps.step24_final_asv_table import run as run_step24
from metagenapp.cli.steps.step25_assign_taxonomy import run as run_step25
from metagenapp.cli.steps.step26_asv_summary import run as run_step26
from metagenapp.cli.utils.final_count_table import build_final_count_table

def run_pipeline(
    input_dir,
    outdir,
    threads=8,
    mode="student",
    classifier="flat",
    marker="16S",   # 👈 NUEVO
    from_step=None,
    min_length=250,
    max_length=600,
    max_ambigs=0,
    max_poly=8,
    extract_centroids="full",
):


    print("🧬 MetagenApp pipeline started")
    print(f"Input dir: {input_dir}")
    print(f"Output dir: {outdir}")

    pipeline_success = True

    input_dir = Path(input_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    PIPELINE_STEPS = [
        "00_Generate_files",
        "01_Assemble_contigs",
        "02_Filter_contigs",
        "03_Unique_contigs",
        "04_Alignment_vsearch",
        "09_Trim_vsearch",
        "10_Clustering_97",
        "11_Extract_Centroids",
        "12_MAFFT_Alignment",
        "14_MAFFT_Filter",
        "15_Filter_Columns",
        "16_Unique_MAFFT",
        "17_Precluster",
        "18_ChimeraDetection",
        "20_TaxonomicClassification",
        "21_RemoveLineages",
        "22_TaxonomicSummary",
        "23_Generate_OTUs",
        "24_Final_ASV_Table",
        "25_Assign_Taxonomy_to_ASVs",
        "26_ASV_Summary_Table",
    ]

    if from_step:
        if from_step not in PIPELINE_STEPS:
            raise ValueError(
                f"--from-step '{from_step}' inválido.\n"
                f"Opciones válidas:\n" + "\n".join(PIPELINE_STEPS)
            )
        start_index = PIPELINE_STEPS.index(from_step)
        print(f"🔁 Reanudando pipeline desde: {from_step}")
    else:
        start_index = 0

    def should_run(step_name: str) -> bool:
        return PIPELINE_STEPS.index(step_name) >= start_index

    # ============================================================
    # Alignment strategy based on mode
    # ============================================================
    if mode == "ref":
        align_target = "reference"
    else:
        align_target = "centroids"

    # ============================================================
    # Reference selection (RAW vs ALIGNED)
    # ============================================================
    if marker == "18S":
        reference_raw = PR2_REFERENCE
        reference_aln = PR2_REFERENCE  # 👈 si aún no tienes PR2 alineada
    else:
        reference_raw = SILVA_REFERENCE_RAW
        reference_aln = SILVA_REFERENCE_ALN

    # ============================================================
    # Paths (siempre definidos)
    # ============================================================
    files_path = outdir / "input_samples.files"

    assembled_fasta = outdir / "assembled_contigs.fasta"
    assembled_count = outdir / "assembled_contigs.count_table"

    filtered_fasta = outdir / "filtered_contigs.fasta"
    filtered_count = outdir / "filtered_contigs.count_table"

    unique_fasta = outdir / "unique_contigs.fasta"
    unique_count = outdir / "unique_contigs.count_table"

    aligned_aln = outdir / "aligned_contigs.aln"
    trimmed_fasta = outdir / "aligned_trimmed.fasta"

    clustered_fasta = outdir / "clustered_97.fasta"
    clustered_uc = outdir / "clustered_97.uc"

    centroid_fasta = outdir / "centroids_all.fasta"
    aligned_mafft = outdir / "aligned_mafft.fasta"
    filtered_mafft = outdir / "aligned_mafft_filtered.fasta"

    # outputs step15/16/17...
    column_filtered_alignment = outdir / "aligned_mafft_filtered_columns.fasta"
    unique_mafft_fasta = outdir / "unique_mafft.fasta"
    unique_mafft_count = outdir / "unique_mafft.count_table"

    preclustered_fasta = outdir / "preclustered.fasta"
    preclustered_count = outdir / "preclustered.count_table"

    nonchimera_fasta = outdir / "non_chimeras.fasta"
    nonchimera_count = outdir / "non_chimeras.count_table"

    taxonomy_file = outdir / "classification.taxonomy"

    final_fasta = outdir / "final_clean.fasta"
    final_count = outdir / "final_clean.count_table"
    final_tax = outdir / "final_clean.taxonomy"

    summary_tax_file = outdir / "summary_tax_phylum.tsv"
    otu_table = outdir / "otu_table_0_03.tsv"

    final_asv_table = outdir / "final_asv_table.tsv"
    asv_taxonomy = outdir / "final_asv.taxonomy"
    asv_summary = outdir / "asv_resumen.tsv"

    # ============================================================
    # Step 00 — Generate files
    # ============================================================
    if should_run("00_Generate_files"):
        start_step("00_Generate_files")
        generate_input_files(input_dir=input_dir, output_files=files_path)
        end_step(success=True)

    # ============================================================
    # Step 01 — Assemble contigs
    # ============================================================
    if should_run("01_Assemble_contigs"):
        start_step("01_Assemble_contigs")
        _, error = assemble_contigs(
            files_path=files_path,
            input_dir=input_dir,
            output_fasta=assembled_fasta,
            output_count=assembled_count,
        )
        if error:
            end_step(success=False)
            raise RuntimeError(error)
        end_step(success=True)
        print("✔ Assemble contigs finished")

    # ============================================================
    # Step 02 — Filter contigs
    # ============================================================
    if should_run("02_Filter_contigs"):
        start_step("02_Filter_contigs")
        _, _, retained, error = filter_contigs(
            fasta_path=assembled_fasta,
            count_table_path=assembled_count,
            output_fasta=filtered_fasta,
            output_count=filtered_count,
            min_length=min_length,
            max_length=max_length,
            max_ambigs=max_ambigs,
            max_poly=max_poly,
        )

        if error:
            end_step(success=False)
            raise RuntimeError(error)

        if retained == 0:
            end_step(success=False)
            raise RuntimeError(
                "Filter contigs removed ALL sequences. "
                "Check min_length, max_length, max_ambigs, max_poly."
            )

        end_step(success=True)
        print(f"✔ Filter contigs finished — {retained} contigs retained")

    # ============================================================
    # Step 03 — Unique contigs
    # ============================================================
    if should_run("03_Unique_contigs"):
        start_step("03_Unique_contigs")
        try:
            _, _, n_unique = generate_unique_contigs(
                fasta_path=filtered_fasta,
                count_table_path=filtered_count,
                output_fasta=unique_fasta,
                output_count=unique_count,
            )
        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Unique contigs failed: {e}")
        end_step(success=True)
        print(f"✔ Unique contigs finished — {n_unique} sequences")

    # ============================================================
    # Step 04 — Alignment vsearch
    # ============================================================
    if should_run("04_Alignment_vsearch"):
        start_step("04_Alignment_vsearch")

        try:
            out = align_with_vsearch(
                fasta_path=unique_fasta,
                reference_path=reference_raw,
                output_path=aligned_aln,
                threads=threads,
                min_identity=0.70,
            )

            aligned_aln = Path(aligned_aln)

            if not out or not aligned_aln.exists():
                raise RuntimeError("Alignment failed")

        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Alignment failed: {e}")

        end_step(success=True)
        print("✔ Alignment finished")

    # ============================================================
    # Step 05 — Trim
    # ===================================================s=========
    if should_run("09_Trim_vsearch"):
        start_step("09_Trim_vsearch")
        try:
            out, msg = recortar_por_alnout(
                aln_path=aligned_aln,
                fasta_input=unique_fasta,
                fasta_output=trimmed_fasta
            )
            trimmed_fasta = Path(trimmed_fasta)

            if not out or not trimmed_fasta.exists():
                raise RuntimeError(msg)
        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Trim failed: {e}")
        end_step(success=True)
        print(f"✔ Trim finished — {msg}")

    # ============================================================
    # Step 06 — Clustering 97
    # ============================================================
    if should_run("10_Clustering_97"):
        start_step("10_Clustering_97")
        try:
            out, msg = cluster_seqs(
                fasta_path=trimmed_fasta,
                output_fasta=clustered_fasta,
                output_uc=clustered_uc,
                identity=0.97,
                threads=threads
            )
            clustered_fasta = Path(clustered_fasta)

            if not out or not clustered_fasta.exists():
                raise RuntimeError(msg)
        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Clustering failed: {e}")
        end_step(success=True)
        print(msg)

    # ============================================================
    # Step 07 — Extract centroids
    # ============================================================
    if should_run("11_Extract_Centroids"):
        start_step("11_Extract_Centroids")

        if extract_centroids == "test":
            n_centroids = 1000
        elif extract_centroids == "student":
            n_centroids = 10000
        else:
            n_centroids = None  # ALL centroides reales

        centroid_fasta = outdir / (
            f"centroids_{n_centroids}.fasta" if n_centroids else "centroids_all.fasta"
        )

        try:
            centroid_path, msg = extraer_subset_centroides(
                fasta_path=clustered_fasta,   # 👈 FASTA del clustering 97%
                output_path=centroid_fasta,
                num_secuencias=n_centroids
            )
            if not centroid_path:
                end_step(success=False)
                raise RuntimeError(msg)
        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Centroid extraction failed: {e}")

        end_step(success=True)
        label = f"{n_centroids} centroides" if n_centroids else "ALL centroides"
        print(f"✔ Centroid extraction finished ({label})")

    if should_run("12_MAFFT_Alignment"):
        print(f"🧪 MAFFT will run in mode={mode}")


    


    # ============================================================
    # Step 14 — MAFFT positional filter
    # ============================================================
    if mode != "ref" and should_run("14_MAFFT_Filter"):
        start_step("14_MAFFT_Filter")

        filtered_mafft = outdir / "aligned_mafft_filtered.fasta"

        try:

            print("⚙️ Using default positional window (50–650)")

            out, retained, total, msg = filtrar_mafft_por_posicion(
                fasta_path=aligned_mafft,
                output_path=filtered_mafft,
                start=50,
                end=650,
                count_table_path=None,
                output_count_path=None
            )

            filtered_mafft = Path(filtered_mafft)

            if not out or retained == 0 or not filtered_mafft.exists():
                raise RuntimeError("MAFFT positional filter removed ALL centroides")

            print("✔ MAFFT positional filtering finished")
            print(msg)

            end_step(success=True)

        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"MAFFT positional filtering failed: {e}")

    # ============================================================
    # Step 10 — Filter columns (mafft)
    # ============================================================
    # fallback si no se corrió el filtro posicional
    
    filtered_mafft = outdir / "aligned_mafft_filtered.fasta"
    filtered_mafft = Path(filtered_mafft)
    
    if not filtered_mafft.exists():
        filtered_mafft = aligned_mafft

    if mode != "ref" and should_run("15_Filter_Columns"):
        start_step("15_Filter_Columns")
        try:
            result = run_step15(
                fasta_input=filtered_mafft,
                outdir=outdir,
                vertical=True,
                trump="."
            )
            column_filtered_alignment = Path(result["output_fasta"])
        except Exception:
            end_step(success=False)
            raise

        end_step(success=True)
        print("✔ Filter alignment columns finished")

    # ============================================================
    # Step 11 — Unique MAFFT
    # ============================================================
    if mode != "ref" and not column_filtered_alignment.exists():
        raise RuntimeError("Column-filtered alignment missing before Unique MAFFT")

    if mode != "ref" and should_run("16_Unique_MAFFT"):
        start_step("16_Unique_MAFFT")
        try:
            result = run_step16(
                fasta_input=column_filtered_alignment,
                outdir=outdir
            )
            unique_mafft_fasta = Path(result["fasta"])
            unique_mafft_count = filtered_count  # hereda count_table
        except Exception:
            end_step(success=False)
            raise

        end_step(success=True)
        print("✔ Unique MAFFT sequences generated")
        print(f"  Total unique sequences: {result['n_unique']}")

    # ============================================================
    # Step 12 — REF Professional Engine
    # ============================================================
    if should_run("12_MAFFT_Alignment"):
        start_step("12_MAFFT_Alignment")

        try:

            if mode == "ref":

                print("🚀 REF mode (professional): EDLib mapping centroids → reference")

                aligned_mafft = outdir / "centroids_vs_reference.edlib.tsv"

                aligned_path = align_centroids_refmode_professional(
                    centroid_fasta=centroid_fasta,
                    reference_fasta=reference_raw,
                    out_tsv=aligned_mafft,
                    threads=min(threads, 16),
                    auto=True,
                )

                aligned_mafft = Path(aligned_path)

            else:

                # modo clásico MAFFT
                aligned_path, msg = align_with_mafft(
                    input_path=centroid_fasta,
                    output_path=aligned_mafft,
                    threads=threads,
                    mode=mode,
                    logger=None
                )

                if not aligned_path:
                    raise RuntimeError(msg)

                aligned_mafft = Path(aligned_path)

            if not aligned_mafft.exists():
                raise RuntimeError("Alignment output file does not exist")

        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Alignment step failed: {e}")

        end_step(success=True)
        print("✔ Alignment step finished")
        # ------------------------------------------------------------
        # REF mode: skip chimera → use centroids directly
        # ------------------------------------------------------------
        if mode == "ref":
            nonchimera_fasta = centroid_fasta
            nonchimera_count = unique_count

    # ============================================================
    # Step 13 — Chimera detection
    # ============================================================
    if mode == "ref":
        print("⏭ Skipping Chimera Detection in REF mode (professional engine active)")
    elif should_run("18_ChimeraDetection"):
        start_step("18_ChimeraDetection")
        try:
            result = run_step18(
                fasta_input=preclustered_fasta,
                count_input=preclustered_count,
                outdir=outdir
            )
            nonchimera_fasta = Path(result["fasta"])
            nonchimera_count = Path(result["count"])
        except Exception:
            end_step(success=False)
            raise
        end_step(success=True)
        print("✔ Chimera detection completed successfully")

    # REF mode uses centroids directly
    if mode == "ref":
        if not centroid_fasta.exists():
            raise RuntimeError("Centroid FASTA missing for REF mode")

        nonchimera_fasta = centroid_fasta
        nonchimera_count = unique_count   

    # ============================================================
    # Step 20 — Taxonomic classification
    # ============================================================

    if should_run("20_TaxonomicClassification"):
        start_step("20_TaxonomicClassification")

        try:

            taxonomy_file = outdir / "classification.taxonomy"

            if classifier == "flat":

                print("🚀 Kraken-lite classifier activated")

                from metagenapp_core.models.kraken_lite import classify_kraken_parallel
                from metagenapp.metagen_config import KRAKEN_INDEX_PATH
                import pickle

                with open(KRAKEN_INDEX_PATH, "rb") as f:
                    kraken_model = pickle.load(f)

                classify_kraken_parallel(
                    fasta_path=nonchimera_fasta,
                    output_path=taxonomy_file,
                    model=kraken_model,
                    threads=threads
                )

            elif classifier == "pro-engine":

                print("🚀 PRO-ENGINE activated (parallel classifier)")

                from metagenapp.pipeline.clasificacion_tax import classify_naive_por_bloques
                from metagenapp.metagen_config import NAIVE_MODEL_PATH

                classify_naive_por_bloques(
                    fasta_path=nonchimera_fasta,
                    output_path=taxonomy_file,
                    modelo_path=NAIVE_MODEL_PATH,
                    block_size=2000,
                    n_threads=threads
                )

            else:
                raise ValueError(f"Unknown classifier: {classifier}")

        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Taxonomic classification failed: {e}")

        end_step(success=True)
        print("✔ Taxonomic classification completed successfully")
        print(f"  Output: {taxonomy_file}")

    # ============================================================
    # Step 15 — Remove unwanted lineages
    # ============================================================
    if should_run("21_RemoveLineages"):
        start_step("21_RemoveLineages")
        try:

            # 🔬 Diferenciar filtro según marcador
            if marker == "18S":
                # En 18S NO remover Eukaryota ni unknown
                taxa_to_remove = "Chloroplast-Mitochondria-Archaea"
            else:
                # En 16S sí remover eucariotas y unknown
                taxa_to_remove = "Chloroplast-Mitochondria-unknown-Archaea-Eukaryota"

            result = run_step21(
                fasta_input=nonchimera_fasta,
                count_input=nonchimera_count,
                taxonomy_input=taxonomy_file,
                outdir=outdir,
                taxa_to_remove=taxa_to_remove
            )

            final_fasta = Path(result["fasta"])
            final_count = Path(result["count"])
            final_tax = Path(result["taxonomy"])

        except Exception:
            end_step(success=False)
            raise

        end_step(success=True)
        print("✔ Unwanted lineages removed successfully")
        print(f"  Remaining FASTA: {final_fasta}")

    # ============================================================
    # Step 16 — Taxonomic summary
    # ============================================================
    if should_run("22_TaxonomicSummary"):
        start_step("22_TaxonomicSummary")
        try:
            result = run_step22(
                taxonomy_input=final_tax,
                count_input=final_count,
                outdir=outdir,
                nivel="Phylum"
            )
            summary_tax_file = result["summary"]
        except Exception:
            end_step(success=False)
            raise
        end_step(success=True)
        print("✔ Taxonomic summary generated successfully")
        print(f"  Summary file: {summary_tax_file}")

    # ============================================================
    # Step 17 — Generate 97% OTUs
    # ============================================================
    if mode != "ref" and should_run("23_Generate_OTUs"):
        start_step("23_Generate_OTUs")
        try:
            result = run_step23(
                fasta_input=final_fasta,
                count_input=final_count,
                outdir=outdir,
                threads=threads
            )
            otu_table = result["otu_table"]
        except Exception:
            end_step(success=False)
            raise
        end_step(success=True)
        print("✔ 97% OTUs generated successfully")
        print(f"  OTU table: {otu_table}")

    # ============================================================
    # Step 18 — Final ASV table
    # ============================================================
    if should_run("24_Final_ASV_Table"):
        start_step("24_Final_ASV_Table")
        try:
            result = run_step24(count_input=final_count, outdir=outdir)
            final_asv_table = result["asv_table"]
        except Exception:
            end_step(success=False)
            raise
        end_step(success=True)
        print("✔ Final ASV table generated successfully")
        print(f"  ASV table: {final_asv_table}")

    # ============================================================
    # Step 25 — Assign taxonomy to ASVs
    # ============================================================
    if should_run("25_Assign_Taxonomy_to_ASVs"):
        start_step("25_Assign_Taxonomy_to_ASVs")
        try:
            result = run_step25(
                count_input=final_count,
                taxonomy_input=final_tax,
                outdir=outdir
            )
            asv_taxonomy = Path(result["taxonomy"])
        except Exception:
            end_step(success=False)
            raise
        end_step(success=True)
        print("✔ ASV taxonomy assigned successfully")
        print(f"  Output: {asv_taxonomy}")

    # ============================================================
    # Step 26 — ASV Summary Table
    # ============================================================
    if should_run("26_ASV_Summary_Table"):
        start_step("26_ASV_Summary_Table")
        try:
            result = run_step26(
                shared_path=final_count,
                taxonomy_path=final_tax,
                outdir=outdir
            )
            asv_summary = Path(result["summary"])
        except Exception:
            pipeline_success = False
            end_step(success=False)
            raise
        end_step(success=True)
        print("✔ ASV summary table generated successfully")
        print(f"  Summary file: {asv_summary}")

    # ============================================================
    # Pipeline end
    # ============================================================
    if pipeline_success:
        print("🧠💻 Pipeline completed successfully")
        print(f"📁 Results written to: {outdir}")
