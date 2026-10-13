#!/usr/bin/env Rscript
# Paired statistical tests for alpha diversity across pipelines (36 matched samples).
# Friedman test (omnibus, paired, 3 groups) + post-hoc pairwise Wilcoxon signed-rank
# (Holm-corrected within each metric's 3 comparisons), for each of the 4 metrics.
# Base R only (stats package) — no new dependencies. Read-only on already-computed
# rarefied per-sample CSVs, no new pipeline computation.

RARE_DIR <- "rarefaction"
TABLES <- "tables"

ma <- read.csv(file.path(RARE_DIR, "alpha_metagenapp_rarefied_all36.csv"), stringsAsFactors = FALSE)
mo <- read.csv(file.path(RARE_DIR, "alpha_mothur_rarefied_all36.csv"), stringsAsFactors = FALSE)
qi <- read.csv(file.path(RARE_DIR, "alpha_qiime2_rarefied_all36.csv"), stringsAsFactors = FALSE)

ma$sample_id <- sub("_.+", "", ma$sample)
mo$sample_id <- mo$sample
qi$sample_id <- sub("_S[0-9]+$", "", qi$sample)

metrics <- c("richness_obs", "chao1", "shannon", "invsimpson")
common <- Reduce(intersect, list(ma$sample_id, mo$sample_id, qi$sample_id))
common <- sort(common)
cat("Matched samples:", length(common), "\n\n")

results <- data.frame()

for (m in metrics) {
  mat <- cbind(
    MetagenApp = ma[[m]][match(common, ma$sample_id)],
    mothur     = mo[[m]][match(common, mo$sample_id)],
    QIIME2     = qi[[m]][match(common, qi$sample_id)]
  )
  rownames(mat) <- common

  fr <- friedman.test(mat)
  cat(sprintf("== %s ==\n", m))
  cat(sprintf("Friedman: chi2=%.3f, df=%d, p=%.5f\n", fr$statistic, fr$parameter, fr$p.value))

  pairs <- list(c("MetagenApp","mothur"), c("MetagenApp","QIIME2"), c("mothur","QIIME2"))
  raw_p <- sapply(pairs, function(pr) {
    wilcox.test(mat[, pr[1]], mat[, pr[2]], paired = TRUE, exact = FALSE)$p.value
  })
  adj_p <- p.adjust(raw_p, method = "holm")

  for (i in seq_along(pairs)) {
    cat(sprintf("  Wilcoxon %s vs %s: raw p=%.5f, holm-adjusted p=%.5f\n",
                pairs[[i]][1], pairs[[i]][2], raw_p[i], adj_p[i]))
    results <- rbind(results, data.frame(
      metric = m, comparison = paste(pairs[[i]], collapse = " vs "),
      friedman_chi2 = round(as.numeric(fr$statistic), 3),
      friedman_p = signif(fr$p.value, 4),
      wilcoxon_p_raw = signif(raw_p[i], 4),
      wilcoxon_p_holm = signif(adj_p[i], 4)
    ))
  }
  cat("\n")
}

write.csv(results, file.path(TABLES, "alpha_diversity_paired_stats_all36.csv"), row.names = FALSE)
cat("Guardado:", file.path(TABLES, "alpha_diversity_paired_stats_all36.csv"), "\n")

cat("\n===== sessionInfo =====\n"); print(sessionInfo())
