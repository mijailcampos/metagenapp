#!/usr/bin/env Rscript
# Taxonomic composition (phylum level) per pipeline, for the ECCB2026 poster.
# Reuses the exact taxonomy-parsing logic already validated in
# beta_diversity_pcoa_permanova.R (same parse_phylum_* functions, same file paths).
# Read-only aggregation, no new pipeline computation.

suppressPackageStartupMessages(library(dplyr))

RESULTS_MA <- "data/samples36/metagenapp"
RESULTS_MO <- "data/samples36/mothur"
RESULTS_QI <- "data/samples36/qiime2"
TABLES     <- "tables"

# ── Load feature tables (same as beta diversity script) ────────────────────
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
cat("Common samples:", length(common_samples), "\n")

rel_ab <- function(m) sweep(m, 1, rowSums(m), "/")

# ── Taxonomy parsing (identical to beta_diversity_pcoa_permanova.R) ────────
parse_phylum_ma <- function(tax_str) {
  parts <- strsplit(tax_str, ";")[[1]]
  if (length(parts) >= 2) trimws(parts[2]) else "Unclassified"
}
parse_phylum_qi <- function(tax_str) {
  parts <- strsplit(tax_str, ";")[[1]]
  p <- grep("^\\s*p__", parts, value = TRUE)
  if (length(p) > 0) sub("^\\s*p__", "", trimws(p[1])) else "Unclassified"
}
parse_phylum_mo <- function(tax_str) {
  parts <- strsplit(tax_str, ";")[[1]]
  if (length(parts) >= 2) sub("\\([0-9]+\\)", "", parts[2]) else "Unclassified"
}

agg_to_phylum <- function(counts, phylum_vec) {
  ph <- phylum_vec[rownames(counts)]
  ph[is.na(ph)] <- "Unclassified"
  agg <- rowsum(counts, ph)
  rel_ab(t(agg))
}

ma_tax <- read.delim(file.path(RESULTS_MA, "final_asv.taxonomy"), header = TRUE, stringsAsFactors = FALSE)
colnames(ma_tax) <- c("feature", "taxonomy")
ma_phy_vec <- setNames(sapply(ma_tax$taxonomy, parse_phylum_ma), ma_tax$feature)
ma_mat <- as.matrix(ma_raw[, common_samples])
ma_ph_rel <- agg_to_phylum(ma_mat, ma_phy_vec)

qi_tax <- read.delim(file.path(RESULTS_QI, "taxonomy.tsv"), stringsAsFactors = FALSE)
colnames(qi_tax)[1:2] <- c("feature","taxonomy")
qi_phy_vec <- setNames(sapply(qi_tax$taxonomy, parse_phylum_qi), qi_tax$feature)
qi_mat <- t(qi_raw[common_samples, ])
qi_ph_rel <- agg_to_phylum(qi_mat, qi_phy_vec)

mo_tax <- read.delim(
  file.path(RESULTS_MO, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"),
  stringsAsFactors = FALSE
)
colnames(mo_tax) <- c("feature","size","taxonomy")
mo_phy_vec <- setNames(sapply(as.character(mo_tax$taxonomy), parse_phylum_mo), mo_tax$feature)
mo_mat <- t(as.matrix(mo_raw[common_samples, ]))
mo_ph_rel <- agg_to_phylum(mo_mat, mo_phy_vec)

# ── Mean relative abundance per phylum per pipeline (averaged across 36 samples) ─
mean_ab <- function(ph_rel) sort(colMeans(ph_rel), decreasing = TRUE)

ma_mean <- mean_ab(ma_ph_rel)
mo_mean <- mean_ab(mo_ph_rel)
qi_mean <- mean_ab(qi_ph_rel)

cat("\n== Top phyla per pipeline (mean relative abundance) ==\n")
cat("-- MetagenApp --\n"); print(round(head(ma_mean, 10) * 100, 2))
cat("-- mothur --\n"); print(round(head(mo_mean, 10) * 100, 2))
cat("-- QIIME2 --\n"); print(round(head(qi_mean, 10) * 100, 2))

# ── Global top 8 phyla (by max mean abundance across the 3 pipelines), rest -> Other ─
all_phyla <- union(union(names(ma_mean), names(mo_mean)), names(qi_mean))
all_phyla <- setdiff(all_phyla, "Unclassified")
max_ab <- sapply(all_phyla, function(p) {
  max(c(if (p %in% names(ma_mean)) ma_mean[p] else 0,
        if (p %in% names(mo_mean)) mo_mean[p] else 0,
        if (p %in% names(qi_mean)) qi_mean[p] else 0))
})
top8 <- names(sort(max_ab, decreasing = TRUE))[1:8]
cat("\nTop 8 filos globales (elegidos por abundancia maxima en cualquier pipeline):\n")
print(top8)

build_row <- function(mean_vec, pipeline) {
  vals <- sapply(top8, function(p) if (p %in% names(mean_vec)) mean_vec[p] else 0)
  unclass <- if ("Unclassified" %in% names(mean_vec)) mean_vec["Unclassified"] else 0
  other <- 1 - sum(vals) - unclass
  data.frame(pipeline = pipeline, phylum = c(top8, "Unclassified", "Other"),
             mean_rel_abundance = c(vals, unclass, max(other, 0)))
}

out <- rbind(
  build_row(ma_mean, "MetagenApp"),
  build_row(mo_mean, "mothur"),
  build_row(qi_mean, "QIIME2")
)
write.csv(out, file.path(TABLES, "taxonomic_composition_phylum.csv"), row.names = FALSE)
cat("\nGuardado:", file.path(TABLES, "taxonomic_composition_phylum.csv"), "\n")

cat("\n===== sessionInfo =====\n"); print(sessionInfo())
