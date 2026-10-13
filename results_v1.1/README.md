# results_v1.1 — statistics of the MetagenApp v1.1 benchmark (manuscript)

Scripts and outputs behind the statistics and figure of the 36-library benchmark
reported in the manuscript (MetagenApp v1.1 vs mothur v.1.48.5 vs QIIME 2 2024.10.1).

The input is the final output of the three pipelines, deposited on Zenodo:
**[10.5281/zenodo.23349702](https://doi.org/10.5281/zenodo.23349702)**
(`benchmark_outputs_v1.1.tar.gz`, 2.7 MB). Its README describes the content, units
(MetagenApp and mothur: 97% OTUs; QIIME 2: DADA2 ASVs) and the feature/read counts of
Tables 2 and 7.

## Reproduce

```bash
cd results_v1.1
bash run_all.sh          # downloads and checks the tarball if data/ is missing (~5 min)
```

Requirements: R ≥ 4.3 with `vegan` (tested with R 4.3.3, vegan 2.6-8) and Python 3 with
`pandas`. Set `RSCRIPT=` / `PYTHON=` to use specific interpreters.

Run from the Zenodo tarball on 13 Oct 2026, every file in `tables/`, `rarefaction/` and
`figures/` was byte-identical to the one used for the manuscript. The logs in `logs/` are
from that run.

## Scripts (in run order)

| # | Script | Output |
|---|---|---|
| 1 | `rarefaction_analysis.R` | `rarefaction/rarefaction_curves.png` |
| 2 | `rarefied_diversity_all36.R` | `rarefaction/alpha_*_rarefied_all36.csv` (alpha diversity, per-sample mean of 1,000 rarefactions to the minimum depth) |
| 3 | `alpha_diversity_stats_all36.R` | `tables/alpha_diversity_paired_stats_all36.csv` (Friedman + Holm-corrected Wilcoxon) |
| 4 | `beta_diversity_pcoa_permanova.R` | `tables/beta_diversity_summary.txt`, `tables/pcoa_coords_*.csv`, `figures/Figure_benchmark_beta_diversity.png` (Bray–Curtis PCoA, PERMANOVA, betadisper, Mantel, Procrustes) |
| 5 | `pcoa_variance.R` | `tables/pcoa_variance_explained.csv` |
| 6 | `taxonomic_composition.R` | `tables/taxonomic_composition_phylum.csv` |
| 7 | `permanova_strata.R` | `logs/07_permanova_strata.log` (combined PERMANOVA with free and within-sample restricted permutations) |
| 8 | `site_permanova.R` | `logs/08_site_permanova.log` (site and carrier PERMANOVA per pipeline) |
| 9 | `genus_resolution.py` | `tables/genus_resolution_table6.tsv` (Table 6) |

Random seeds are fixed in the scripts (19760620; one seed per pipeline for the rarefactions),
so rarefied values and permutation p-values are reproducible.

The feature and read counts of Table 2 (36 libraries) and Table 7 (mock communities) are
read directly from the deposited outputs; see the correspondence table in the Zenodo README.
