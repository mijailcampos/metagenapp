#!/usr/bin/env Rscript
# Recompute just the % variance explained per PCoA axis (PC1/PC2) per pipeline.
# Same distance matrices as beta_diversity_pcoa_permanova.R, just re-deriving
# two numbers per pipeline that weren't saved as text the first time.

suppressPackageStartupMessages({
  library(vegan)
  library(ape)
})

RESULTS_MA <- "data/samples36/metagenapp"
RESULTS_MO <- "data/samples36/mothur"
RESULTS_QI <- "data/samples36/qiime2"
TABLES     <- "tables"

ma_raw <- read.delim(file.path(RESULTS_MA, "final_asv_table.tsv"), row.names = 1, check.names = FALSE)
colnames(ma_raw) <- sub("_.+", "", colnames(ma_raw))

mo_raw <- read.delim(
  file.path(RESULTS_MO, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared"),
  check.names = FALSE
)
rownames(mo_raw) <- mo_raw$Group
mo_raw <- mo_raw[, -(1:3)]

qi_raw <- read.delim(file.path(RESULTS_QI, "feature-table.tsv"),
                     skip = 1, row.names = 1, check.names = FALSE)
colnames(qi_raw) <- sub("_S[0-9]+$", "", colnames(qi_raw))
qi_raw <- t(as.matrix(qi_raw))

common_samples <- Reduce(intersect, list(colnames(ma_raw), rownames(mo_raw), rownames(qi_raw)))
common_samples <- sort(common_samples)

rel_ab <- function(m) sweep(m, 1, rowSums(m), "/")
ma_rel <- rel_ab(t(ma_raw[, common_samples]))
mo_rel <- rel_ab(mo_raw[common_samples, ])
qi_rel <- rel_ab(qi_raw[common_samples, ])

variance_explained <- function(rel_mat, label) {
  d <- vegdist(rel_mat, method = "bray")
  pc <- pcoa(d)
  var_pct <- pc$values$Relative_eig * 100
  cat(sprintf("%-12s PC1=%.1f%% PC2=%.1f%%\n", label, var_pct[1], var_pct[2]))
  data.frame(pipeline = label, PC1_pct = round(var_pct[1], 1), PC2_pct = round(var_pct[2], 1))
}

out <- rbind(
  variance_explained(ma_rel, "MetagenApp"),
  variance_explained(mo_rel, "mothur"),
  variance_explained(qi_rel, "QIIME2")
)
write.csv(out, file.path(TABLES, "pcoa_variance_explained.csv"), row.names = FALSE)
cat("Guardado:", file.path(TABLES, "pcoa_variance_explained.csv"), "\n")

cat("\n===== sessionInfo =====\n"); print(sessionInfo())
