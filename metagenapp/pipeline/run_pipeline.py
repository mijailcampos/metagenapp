import os
import time
import resource
from pathlib import Path

from metagenapp.metagen_config import (
    PR2_REFERENCE,
    SILVA_REFERENCE_RAW,
    SILVA_REFERENCE_ALN,
    NAIVE_MODEL_PATH,
    get_naive_model_path,
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
from metagenapp.pipeline.aggregate_cluster_counts import aggregate_cluster_counts

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
    marker="16S",
    model_type="general",
    from_step=None,
    min_length=250,
    max_length=600,
    max_ambigs=0,
    max_poly=8,
    extract_centroids="full",
):


    print()
    print("═" * 54)
    print("  MetagenApp — metabarcoding pipeline")
    print("═" * 54)
    print(f"  Mode:       {mode}")
    print(f"  Marker:     {marker}")
    print(f"  Classifier: {classifier}")
    print(f"  Threads:    {threads}")
    print(f"  Input:      {input_dir}")
    print(f"  Output:     {outdir}")
    print("═" * 54)
    print()

    t_start = time.time()
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

    _merge_results = []   # (sample, pairs, merged, pct, fastq_path)
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
        _, error, _merge_results = assemble_contigs(  # noqa: F841
            files_path=files_path,
            input_dir=input_dir,
            output_fasta=assembled_fasta,
            output_count=assembled_count,
        )
        if error:
            end_step(success=False)
            raise RuntimeError(error)
        end_step(success=True)

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

        print(f"   {retained:,} contigs retained")
        end_step(success=True)

    # ============================================================
    # QA early exit — solo ensamblado + filtrado
    # ============================================================
    if mode == "qa":
        from metagenapp.pipeline.summarize_contigs import summarize_contigs
        from metagenapp.pipeline.qa import _print_summary_table, _print_marker_section, _MARKERS
        from metagenapp.pipeline.step_tracker import emit_qa_result

        print()
        print("═" * 54)
        print("  QA — Resumen de calidad del ensamblado")
        print("═" * 54)

        # Tabla por muestra (de assemble_contigs)
        total_pairs   = sum(r[1] for r in _merge_results)
        total_merged  = sum(r[2] for r in _merge_results)
        overall_pct   = (total_merged / total_pairs * 100) if total_pairs else 0

        qa_data = {
            "retained": retained,
            "filters": {
                "min_length": min_length,
                "max_length": max_length,
                "max_ambigs": max_ambigs,
                "max_poly":   max_poly,
            },
            "assembly": {
                "samples": [
                    {"sample": r[0], "pairs": r[1], "merged": r[2], "pct": round(r[3], 1)}
                    for r in _merge_results
                ],
                "total_pairs":  total_pairs,
                "total_merged": total_merged,
                "overall_pct":  round(overall_pct, 1),
            },
            "stats": {},
            "stats_table": {},
            "markers": [],
        }

        try:
            summary = summarize_contigs(
                fasta_path=filtered_fasta,
                count_table_path=filtered_count,
            )
            _print_summary_table(summary)
            median_bp = float(summary["Median"]["NBases"])
            min_bp    = float(summary["Minimum"]["NBases"])
            mean_bp   = float(summary["Mean"]["NBases"])
            _print_marker_section(median_bp)

            qa_data["stats"] = {
                "total_assembled":  int(summary.get("total sequences", 0)),
                "unique_contigs":   int(summary.get("# of contigs", 0)),
                "median_bp":  median_bp,
                "min_bp":     min_bp,
                "mean_bp":    mean_bp,
            }
            # Tabla completa Min/Median/Mean × columnas
            cols = ["Start", "End", "NBases", "Ambigs", "Polymer", "NumSeqs"]
            qa_data["stats_table"] = {
                stat: {col: summary[stat][col] for col in cols}
                for stat in ("Minimum", "Median", "Mean")
                if stat in summary
            }
            qa_data["markers"] = [
                {"region": name, "primers": primers, "min_bp": lo, "max_bp": hi,
                 "match": lo <= median_bp <= hi}
                for name, primers, lo, hi in _MARKERS
            ]
        except Exception as e:
            print(f"   (resumen no disponible: {e})")

        print()
        print(f"  Filtros aplicados: min={min_length}bp  max={max_length}bp  ambigs≤{max_ambigs}  poly≤{max_poly}")
        print("═" * 54)

        emit_qa_result(qa_data)
        return

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
        print(f"   {n_unique:,} unique sequences")
        end_step(success=True)

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

    # ============================================================
    # Step 05 — Trim
    # ============================================================
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

        label = f"{n_centroids:,}" if n_centroids else "all"
        print(f"   {label} centroids extracted")
        end_step(success=True)

        # default path for downstream steps
        nonchimera_fasta = centroid_fasta
        nonchimera_count = unique_count

        # ── Aggregate read counts from unique sequences to centroids ──
        # unique_count has one row per unique sequence — not per centroid.
        # We use the UC file to sum reads into their cluster centroids,
        # producing the correct read-weighted count table for all
        # downstream steps (taxonomy summary, ASV table, etc.).
        centroids_count = outdir / "centroids_all.count_table"
        out, n_c, total_r, err = aggregate_cluster_counts(
            uc_path=clustered_uc,
            count_table_path=unique_count,
            output_path=centroids_count,
        )
        if err:
            raise RuntimeError(f"Cluster count aggregation failed: {err}")
        nonchimera_count = Path(out)
        print(f"   {n_c:,} centroids, {total_r:,} reads aggregated")

        if classifier == "kraken-lite":
            preclustered_fasta = centroid_fasta
            preclustered_count = nonchimera_count


    



    # ============================================================
    # Step 12 — Alignment (MAFFT or REF engine)
    # ============================================================

    if classifier != "kraken-lite" and should_run("12_MAFFT_Alignment"):

        start_step("12_MAFFT_Alignment")

        try:

            if mode == "ref":

                print("   EDLib mapping centroids → reference")

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

    else:
        print("\n   [kraken-lite] skipping MAFFT alignment")

    # ============================================================
    # Step 14 — MAFFT positional filter
    # ============================================================
    if classifier != "kraken-lite" and mode != "ref" and should_run("14_MAFFT_Filter"):
        start_step("14_MAFFT_Filter")

        filtered_mafft = outdir / "aligned_mafft_filtered.fasta"

        try:

            # ----------------------------------------------------
            # Detectar longitud del alineamiento
            # ----------------------------------------------------
            from Bio import SeqIO

            with open(aligned_mafft) as f:
                first_seq = next(SeqIO.parse(f, "fasta"))
                aln_len = len(first_seq.seq)

            # ----------------------------------------------------
            # Ventana adaptativa
            # ----------------------------------------------------
            start = max(10, int(aln_len * 0.05))
            end = int(aln_len * 0.95)

            print(f"   window {start}–{end} / {aln_len} bp")

            out, retained, total, msg = filtrar_mafft_por_posicion(
                fasta_path=aligned_mafft,
                output_path=filtered_mafft,
                start=start,
                end=end,
                count_table_path=None,
                output_count_path=None
            )

            filtered_mafft = Path(filtered_mafft)

            # ----------------------------------------------------
            # Protección: evitar eliminar todo
            # ----------------------------------------------------
            if retained == 0 or not filtered_mafft.exists():
                print("⚠️ Positional filter removed all sequences — reverting to original alignment")
                filtered_mafft = aligned_mafft
                retained = total

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

    if classifier != "kraken-lite" and mode != "ref" and should_run("15_Filter_Columns"):
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

    # ============================================================
    # Step 16 — Unique MAFFT
    # ============================================================
    if classifier != "kraken-lite" and mode != "ref" and not column_filtered_alignment.exists():
        raise RuntimeError("Column-filtered alignment missing before Unique MAFFT")

    if classifier != "kraken-lite" and mode != "ref" and should_run("16_Unique_MAFFT"):
        start_step("16_Unique_MAFFT")
        try:
            result = run_step16(
                fasta_input=column_filtered_alignment,
                outdir=outdir
            )
            unique_mafft_fasta = Path(result["fasta"])
            unique_mafft_count = nonchimera_count  # centroid-level count table
        except Exception:
            end_step(success=False)
            raise

        print(f"   {result['n_unique']:,} unique sequences")
        end_step(success=True)


    # ============================================================
    # Step 13 — Chimera detection
    # ============================================================

    if classifier == "kraken-lite":
        print("\n   [kraken-lite] skipping chimera detection")

    elif mode == "ref" and classifier == "naive-v2":
        # En ref mode + naive-v2 los centroides no pasan por MAFFT,
        # pero sí necesitan chimera removal antes de clasificar.
        if should_run("18_ChimeraDetection"):
            start_step("18_ChimeraDetection")
            try:
                result = run_step18(
                    fasta_input=centroid_fasta,
                    count_input=nonchimera_count,
                    outdir=outdir
                )
                nonchimera_fasta = Path(result["fasta"])
                nonchimera_count = Path(result["count"]) if result["count"] else None
                end_step(success=True)
            except Exception as e:
                end_step(success=False)
                raise RuntimeError(str(e))

    elif mode == "ref":
        pass  # ref mode skips chimera (professional engine handles it)

    else:

        if should_run("18_ChimeraDetection"):

            start_step("18_ChimeraDetection")

            try:

                result = run_step18(
                    fasta_input=unique_mafft_fasta,
                    count_input=unique_mafft_count,
                    outdir=outdir
                )

                nonchimera_fasta = Path(result["fasta"])
                nonchimera_count = Path(result["count"]) if result["count"] else None

                end_step(success=True)

            except Exception as e:

                end_step(success=False)
                raise RuntimeError(str(e))
            
    # ============================================================
    # Step 20 — Taxonomic classification
    # ============================================================

    if should_run("20_TaxonomicClassification"):
        start_step("20_TaxonomicClassification")

        try:

            taxonomy_file = outdir / "classification.taxonomy"

            # ----------------------------------------------------
            # FLAT (simple taxonomy assigner)
            # ----------------------------------------------------
            if classifier == "flat":

                print("🧠 Naive Bayes classifier activated")

                from metagenapp.pipeline.clasificacion_tax import classify_naive_por_bloques


                classify_naive_por_bloques(
                    fasta_path=centroid_path,
                    output_path=taxonomy_file,
                    modelo_path=NAIVE_MODEL_PATH,
                    block_size=10000,
                    n_threads=threads
                )

            # ----------------------------------------------------
            # NAIVE V2
            # ----------------------------------------------------
            elif classifier == "naive-v2":
                naive_model_path = get_naive_model_path(marker=marker, model_type=model_type)
                print(f"🚀 Naive V2 classifier activated [{model_type}]")
                print(f"   Model: {naive_model_path}")
                import pickle, gc
                from Bio import SeqIO
                from metagenapp_core.models.naive_v2_engine import classify_seq_v2
                with open(naive_model_path, "rb") as f:
                    model_v2 = pickle.load(f)
                with open(taxonomy_file, "w") as fout:
                    for i, record in enumerate(SeqIO.parse(nonchimera_fasta, "fasta")):
                        seq_id, taxon = classify_seq_v2(record.id, str(record.seq).upper(), model_v2)
                        fout.write(f"{seq_id}\t{taxon or 'Unclassified'}\n")
                        if (i+1) % 1000 == 0:
                            print(f"  ↳ {i+1} sequences classified")
                            gc.collect()

            # ----------------------------------------------------
            # KRAKEN-LITE
            # ----------------------------------------------------
            elif classifier == "kraken-lite":

                from metagenapp_core.models.kraken_lite import classify_kraken_parallel
                from metagenapp.metagen_config import KRAKEN_INDEX_PATH, KRAKEN_INDEX_PATH_18S
                import pickle

                kraken_index = KRAKEN_INDEX_PATH_18S if marker == "18S" else KRAKEN_INDEX_PATH

                if not os.path.exists(kraken_index):
                    if marker == "18S":
                        raise RuntimeError(
                            f"Índice PR2 no encontrado: {kraken_index}\n"
                            "  Constrúyelo con: python3 scripts/rebuild_kraken_pr2.py"
                        )
                    raise RuntimeError(f"Kraken index not found: {kraken_index}")

                with open(kraken_index, "rb") as f:
                    kraken_model = pickle.load(f)

                classify_kraken_parallel(
                    fasta_path=nonchimera_fasta,
                    output_path=taxonomy_file,
                    model=kraken_model,
                    threads=threads
                )
                
            # ----------------------------------------------------
            # METASPECIES (SSI+ANI species-level)
            # ----------------------------------------------------
            elif classifier == "metaspecies":

                print("🔬 MetaSpecies classifier activated (SSI+ANI)")

                from metagenapp.pipeline.metaspecies_classify import classify_fasta_to_taxonomy

                DEFAULT_REF = "/data/projects/metaspecies/data/ref_index_silva"
                ref_dir = os.environ.get("METASPECIES_REF_DIR", DEFAULT_REF)

                def _ms_progress(current, total, msg=""):
                    if msg:
                        print(f"  ↳ {msg}")

                result_path, error = classify_fasta_to_taxonomy(
                    fasta_path=str(nonchimera_fasta),
                    ref_dir=ref_dir,
                    output_path=str(taxonomy_file),
                    progress_callback=_ms_progress,
                )

                if not result_path:
                    raise RuntimeError(f"MetaSpecies classification failed: {error}")

            # ----------------------------------------------------
            # PRO ENGINE
            # ----------------------------------------------------
            elif classifier == "pro-engine":

                print("🚀 PRO-ENGINE activated")

                from metagenapp.pipeline.clasificacion_tax import classify_pro_engine

                classify_pro_engine(
                    fasta_path=nonchimera_fasta,
                    output_path=taxonomy_file,
                    threads=threads
                )

            else:
                raise ValueError(f"Unknown classifier: {classifier}")

        except Exception as e:
            end_step(success=False)
            raise RuntimeError(f"Taxonomic classification failed: {e}")

        end_step(success=True)

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

    # ============================================================
    # Pipeline end
    # ============================================================
    if pipeline_success:
        elapsed   = time.time() - t_start
        ram_mb    = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

        elapsed_m = int(elapsed) // 60
        elapsed_s = int(elapsed) % 60

        def _fsize(p):
            try:
                b = Path(p).stat().st_size
                if b >= 1_000_000:
                    return f"{b/1_000_000:.1f} MB"
                if b >= 1_000:
                    return f"{b/1_000:.1f} KB"
                return f"{b} B"
            except Exception:
                return "—"

        def _table_dims(p):
            try:
                with open(p) as fh:
                    header = fh.readline()
                    n_cols = len(header.split("\t"))
                    n_otus = sum(1 for _ in fh)
                # count_table: cols = ID + samples (no "total" column)
                n_samples = n_cols - 1
                return f"{n_otus:,} OTUs  ×  {n_samples} muestras"
            except Exception:
                return "—"

        _CY = "\033[96m"   # cian neón
        _MG = "\033[95m"   # magenta neón
        _GR = "\033[92m"   # verde
        _BL = "\033[94m"   # azul
        _RS = "\033[0m"    # reset

        def _dirsize(p):
            try:
                total = sum(f.stat().st_size for f in Path(p).rglob("*") if f.is_file())
                if total >= 1_000_000_000:
                    return f"{total/1_000_000_000:.2f} GB"
                return f"{total/1_000_000:.1f} MB"
            except Exception:
                return "—"

        print()
        print(f"║{_MG}▓▓░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▓▓{_RS}║")
        print(f"║  {_CY}██████╗  ██████╗ ███╗  ██╗███████╗ ██╗{_RS}            ║")
        print(f"║  {_CY}██╔══██╗██╔═══██╗████╗ ██║██╔════╝ ██║{_RS}            ║")
        print(f"║  {_CY}██║  ██║██║   ██║██╔██╗██║█████╗   ██║{_RS}            ║")
        print(f"║  {_CY}██║  ██║██║   ██║██║╚████║██╔══╝   ╚═╝{_RS}            ║")
        print(f"║  {_CY}██████╔╝╚██████╔╝██║ ╚███║███████╗ ██╗{_RS}            ║")
        print(f"║  {_CY}╚═════╝  ╚═════╝ ╚═╝  ╚══╝╚══════╝ ╚═╝{_RS}            ║")
        print(f"║{_MG}▓▓░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▓▓{_RS}║")
        print()
        print("  ── RESULTADOS ──────────────────────────────────────")
        print(f"  {_GR}📂 Carpeta      {outdir}{_RS}")
        print(f"  📊 Tabla ASV    {_fsize(final_asv_table)}  ·  {_table_dims(final_asv_table)}")
        print(f"  🏷  Taxonomía    {_fsize(asv_taxonomy)}")
        print(f"  📋 Resumen      {_fsize(asv_summary)}")
        print()
        print("  ── MÉTRICAS ─────────────────────────────────────────")
        print(f"  ⏱  Tiempo total  {_BL}{elapsed_m}m {elapsed_s:02d}s{_RS}")
        print(f"  🧠 RAM máxima    {ram_mb:.0f} MB  ({ram_mb/1024:.2f} GB)")
        print(f"  💾 Carpeta       {_dirsize(outdir)}")
        print("  ─────────────────────────────────────────────────────")
        print()
