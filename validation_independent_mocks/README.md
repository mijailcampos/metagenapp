# Independent mock-community validation: evidence archive

This directory documents the evaluation of MetagenApp, mothur and QIIME 2 on ten ZymoBIOMICS
mock-community libraries that were not used during the development of MetagenApp
(sets A, D and E), reported in the manuscript as Results §9, Table 8 and Methods §1.2 and §6.3.

## Timeline (UTC)

| Date and time | Event | Record |
|---|---|---|
| 2026-11-11 19:51 | Analysis plan frozen and published | `PREREGISTRATION.md`, commit e353b80 |
| 2026-11-11 20:00 | First FASTQ file downloaded from ENA (MD5 checked) | `results/manifest_download.tsv` |
| 2026-11-11 20:26 – 2026-11-12 00:56 | Pre-specified runs: MetagenApp v1.1, QIIME 2, mothur | `logs/metagenapp_v11.log`, `logs/qiime2.log`, `logs/mothur.log` |
| 2026-11-12 | Post hoc reanalysis A: tag v1.1 with only `trim_primers.py` replaced | `logs/posthoc_*.log` |
| 2026-11-12 19:20 | Amendment 1 published | `AMENDMENT_1_v1.2.3_reevaluation.md`, commit 14c6a1c |
| 2026-11-12 19:22 – 19:38 | Post hoc reanalysis B: released MetagenApp v1.2.3 | `logs/reeval_v123_*.log` |

Both post hoc reanalyses were run after the pre-specified results were known. They are not an independent
validation of the corrected implementation.

## Contents

- `PREREGISTRATION.md`: the pre-specified plan (unchanged since e353b80).
- `AMENDMENT_1_v1.2.3_reevaluation.md`: the amendment declaring the v1.2.3 re-evaluation before it was run.
- `scripts/`: download, run and evaluation scripts, as executed. Paths refer to the analysis server.
  - `download.sh`, `input_map.tsv`, `qiime2_manifest.tsv`: data download and library naming.
  - `run_metagenapp_v11.sh`, `run_qiime2.sh`, `run_mothur.sh`, `mothur_validation.batch`: pre-specified runs.
  - `run_posthoc_primerfix.sh`, `eval_posthoc.py`: reanalysis A, including the regression on the mock communities of Table 7.
  - `trim_primers_ESTRICTO_min_overlap_len.py`, `LEEME_estricto.txt`: an earlier, stricter version of the corrected step (minimum overlap = full primer length), superseded before the final reanalysis A. Its notes are in Spanish.
  - `run_reeval_v123.sh`, `eval_reeval_v123.py`: reanalysis B (v1.2.3).
  - `eval_quick.py`, `eval_q2.py`, `eval_mothur.py`: per-pipeline metrics.
  - `tabla_resultados_manuscrito.py`, `tabla_v123_amendment1.py`: the values of Table 8.
  - `regression_table7.py`: compares reanalysis A with Table 7.
  - `otu_abundance_check.py`: abundance distribution of the OTUs.
  - `qiime2_secondary_per_library.py`: per-library secondary Bray–Curtis dissimilarity for QIIME 2.
  - `sensibilidad_salmonella_3herramientas.py`: Enterobacteriaceae-to-Salmonella reassignment applied to all three pipelines.
- `logs/`: complete run logs. `metagenapp_v11_FALLIDA_enlaces_rotos.log` is a first attempt that failed because of a scripting error in the input links ("fallida, enlaces rotos" = "failed, broken links"). No reads were processed in that attempt.
- `results/`: final tables and outputs.
  - Per-pipeline feature tables with taxonomy: `metagenapp_v11/`, `metagenapp_posthoc/` (reanalysis A), `metagenapp_v123/` (reanalysis B), `qiime2/`, `mothur/`, plus `eval_qiime2_resumen.tsv` and `eval_mothur_resumen.tsv`.
  - `tabla_resultados_manuscrito.tsv`, `tabla_v123_amendment1.txt`: metrics per set.
  - `pares_reales_fastq.tsv`: read pairs counted in the downloaded FASTQ files. These counts are the retention denominator.
  - `regression_table7/`: reanalysis A outputs for the mock communities of Table 7.
  - `salmonella_check/`: VSEARCH comparison of the QIIME 2 Enterobacteriaceae ASVs with the official ZymoBIOMICS reference sequences.
  - `Zymo_D6306_16S_theoretical.tsv`: expected 16S rRNA gene copy composition.

## Deviations from the plan

1. **Post hoc reanalyses.** After the pre-specified evaluation, the primer-trimming step was corrected. Reanalysis A used v1.1 with only that step replaced; reanalysis B used v1.2.3, as declared in Amendment 1. Both are reported separately from the pre-specified results.
2. **Metric scripts.** The plan names `mock_metrics.py`. The metrics were calculated instead with the scripts in `scripts/`, using the same definitions and formulae. They reproduce Table 7 exactly (`regression_table7.py`).
3. **Per-library dissimilarities.** The manuscript table gives mean ± SD per set. Per-library values, from which the minimum and maximum can be read, are given in the supplementary material.
4. **Retention denominator.** Retention is relative to the read pairs counted in the FASTQ files (`pares_reales_fastq.tsv`). The ENA `read_count` field in `manifest_download.tsv` counts R1 and R2 separately for the ERR runs of set D.
5. **Additional sensitivity analysis.** This analysis was not pre-specified. The Enterobacteriaceae-to-Salmonella reassignment, pre-specified for QIIME 2 only, was also applied to MetagenApp and mothur. It did not change their values.

## Not included

- Raw reads: available from ENA/SRA (BioProjects PRJNA1069904, PRJEB88026 and PRJNA781121; run accessions in `scripts/input_map.tsv`).
- Intermediate files of the runs (about 157 GB, mostly mothur alignments).
- The ZymoBIOMICS reference sequences used in `salmonella_check/`, which are distributed by Zymo Research.

Some file names and log messages are in Spanish.
