# Pre-specified plan: independent mock-community validation of MetagenApp v1.1

Status: FROZEN on 2026-11-11. This document will be committed to the public GitHub
repository (github.com/mijailcampos/metagenapp) **before** the full read sets are
downloaded and before any pipeline is run. The commit date is the public record that the
datasets, parameters, metrics and reporting rules below were fixed before results were seen.

## 1. Purpose

MetagenApp v1.1 was finalized after preliminary analyses of the benchmark dataset
(PRJNA1051577) and of the mock datasets used in the manuscript (PRJNA667146, PRJNA659641).
Its accuracy on those mocks may therefore be optimistic. This validation evaluates the
three pipelines, with versions and parameters frozen beforehand, on mock-community
libraries that played no role in developing or configuring any of them.

## 2. What was seen before freezing (disclosure)

To check technical compatibility only, the first 1,000 read pairs of one library per
candidate project were inspected for read length and for the presence of the 341F/805R
primer sequences. No pipeline was run, no taxonomy was assigned and no composition was
examined on any of the libraries below.

## 3. Datasets

All libraries are Illumina MiSeq paired-end V3–V4 amplicons generated with 341F
(CCTACGGGNGGCWGCAG) and 805R (GACTACHVGGGTATCTAATCC; called 785R in some sources, same
sequence). All are ZymoBIOMICS standards (8 bacterial strains).

| Set | BioProject | Material | Runs | Libraries (run accessions) |
|---|---|---|---|---|
| A | PRJNA1069904 (Finland; PLoS ONE 2025, doi:10.1371/journal.pone.0317513) | ZymoBIOMICS Microbial Community Standard (D6300, cells; DNA extraction included) | 1 run, shared with natural samples and negative controls | SRR27790137 (ZY1), SRR27790136 (ZY2) |
| D | PRJEB88026 (Linnaeus University) | ZymoBIOMICS Microbial Community DNA Standard | 5 MiSeq runs (flow cells KKBHR, KLC7M, KJJFC, KV2C6, KV3TP) | ERR14820732, ERR14820828, ERR14820924, ERR14821020, ERR14821116 |
| E | PRJNA781121 (Linnaeus University) | ZymoBIOMICS Microbial Community DNA Standard | released 2021; run structure to be recorded from metadata | SRR16997123, SRR16997307, SRR16997311 |

Sets D and E come from the same laboratory and are not fully independent of each other;
results are reported per set as well as pooled.

Every library listed is included. A library is excluded only for a technical reason fixed
here: fewer than 1,000 raw read pairs, or a corrupted download that cannot be repaired
after re-download (checked against ENA MD5). Exclusions apply to all three pipelines and
are reported.

### Optional set F (decided by rule, not by results)

PRJEB74474 (Lausanne) contains ATCC mock libraries (MOCK_ATCC_1_2_3 ×2, MOCK_7_8 ×2) whose
product and composition are not documented publicly. Set F is included only if, by
**25 November 2026**, the authors confirm (i) the exact product and lot,
(ii) its expected composition, and (iii) the library protocol, and the reads are
compatible with the unchanged parameters. Otherwise F is not analysed. If included, F is
analysed with exactly the same parameters and metrics, and reported as a separate set.

## 4. Expected composition

Genus-level expected profile = theoretical **16S rRNA gene copy** composition from the
manufacturer's official sheet ("16S Only" column), the same reference type as Table 7
(file `Zymo_D6306_16S_theoretical.tsv`, genera named as in SILVA 138.2):

Pseudomonas 4.2, Escherichia-Shigella 10.1, Salmonella 10.4, Limosilactobacillus 18.4,
Enterococcus 9.9, Staphylococcus 15.5, Listeria 14.1, Bacillus 17.4 (%).

Checked before freezing: the "16S Only" column of the official D6300 instruction manual
(Ver. 1.2.0, revised 24 January 2023, Table 1) gives the same values as the D6305/D6306
sheet used for Table 7. Cell (A) and DNA (D, E) standards therefore share the expected 16S
profile, but A also includes DNA-extraction bias; this is stated when interpreting A.
Some strains were replaced by similar strains from lot ZRC190633 onwards; the species and
genus composition did not change, so the genus-level profile applies to all lots.

## 5. Software and parameters (frozen)

Same server (biodyssey), 16 threads, same reference files as the manuscript.

- **MetagenApp v1.1** as released (tag v1.1; Zenodo doi:10.5281/zenodo.23325861), with
  reference files from doi:10.5281/zenodo.23326733.
  `--mode ref --marker 16S --classifier wang --primer-f CCTACGGGNGGCWGCAG
  --primer-r GACTACHVGGGTATCTAATCC --min-length 400 --max-length 440 --max-ambigs 0
  --max-poly 8 -t 16`
- **mothur v.1.48.5**, batch identical to `mocks/aachen/mothur/mothur_aachen_seed.batch`
  (SILVA SEED v138.2 alignment, positions 6428–23440; classification against SILVA 138.2
  NR99), with only input/output directories changed.
- **QIIME 2 2024.10.1**: cutadapt trim-paired with the same primers and
  `--p-discard-untrimmed`; DADA2 denoise-paired `--p-trunc-len-f 280 --p-trunc-len-r 220`
  (benchmark values; in the inspected reads, median lengths after primer removal are
  ≥ 284 bp for R1 and ≥ 257 bp for R2); classify-sklearn with
  `silva-138.2-99-nb-classifier.qza` (Zenodo record of the benchmark outputs).

No parameter is changed after seeing results. If a pipeline fails technically on a
library or set, the failure is reported as a result; it is not fixed by tuning.

## 6. Metrics (same code as Table 7)

Computed with `mock_metrics.py` (formula unchanged), per pipeline and per set (A, D, E),
and pooled over all libraries:

1. Expected genera detected (pooled, and in every library).
2. Percentage of reads assigned outside the expected genera, including reads of expected
   members classified only above genus.
3. Bray–Curtis dissimilarity between observed and expected genus profiles (unexpected
   reads count as disagreement), pooled and per library (mean, SD, min, max).
4. Number of features (OTUs/ASVs) and final read retention.

Secondary, reported alongside and never in place of metric 3: Bray–Curtis for QIIME 2
counting Enterobacteriaceae-only reads as Salmonella, as in Table 7.

No rarefaction. No performance thresholds are set: results are reported descriptively and
compared with the variation among libraries within each set.

## 7. Reporting commitments

- All results are reported in the manuscript (new table or rows in Table 7, Methods,
  Results, Discussion), including results unfavourable to MetagenApp.
- The text states that these libraries were not used during development and that this
  plan was committed before analysis (commit hash cited).
- The original limitation (benchmark and original mocks used for both development and
  evaluation) remains declared.
- Any deviation from this plan is listed explicitly with its reason.
- Outputs, scripts and logs are archived with the release that accompanies the manuscript.
