#!/usr/bin/env Rscript
# Genus-level taxonomic resolution statistics for the three pipelines
# Metrics: % classified, % ambiguous, % unclassified — by features AND by reads

suppressPackageStartupMessages({
  library(biomformat)
  library(dplyr)
})

RESULTS <- "/data/results/runs/01_benchmarking/benchmark_controlled/results_20260418_0205"

# ── helpers ───────────────────────────────────────────────────────────────────

label_genus <- function(g) {
  # Returns "classified" | "ambiguous" | "unclassified"
  if (is.na(g) || g == "" || nchar(trimws(g)) == 0) return("unclassified")
  g <- trimws(g)
  ambig_patterns <- c("unclassified", "Unclassified", "uncultured", "Uncultured",
                       "Incertae Sedis", "incertae sedis", "metagenome",
                       "ambiguous", "unknown", "Unknown")
  if (any(sapply(ambig_patterns, function(p) grepl(p, g, fixed = TRUE))))
    return("ambiguous")
  return("classified")
}

summarise_resolution <- function(labels, reads) {
  # labels: character vector ("classified"|"ambiguous"|"unclassified")
  # reads:  numeric vector of total reads per feature
  n_feat   <- length(labels)
  n_reads  <- sum(reads)
  tbl_feat <- table(labels) / n_feat * 100
  tbl_read <- tapply(reads, labels, sum) / n_reads * 100
  # fill missing categories
  for (cat in c("classified","ambiguous","unclassified")) {
    if (is.null(tbl_feat[cat]) || is.na(tbl_feat[cat])) tbl_feat[cat] <- 0
    if (is.null(tbl_read[cat]) || is.na(tbl_read[cat])) tbl_read[cat] <- 0
  }
  list(feat = tbl_feat, read = tbl_read,
       n_feat = n_feat, n_reads = n_reads)
}

print_summary <- function(name, s) {
  cat(sprintf("\n── %s (%d features | %s reads) ────────────────\n",
              name, s$n_feat, format(s$n_reads, big.mark = ",")))
  for (cat in c("classified","ambiguous","unclassified")) {
    cat(sprintf("  %-14s  features: %5.1f%%   reads: %5.1f%%\n",
                cat, s$feat[cat], s$read[cat]))
  }
}

# ── 1. MetagenApp ─────────────────────────────────────────────────────────────
ma_tax <- read.delim(file.path(RESULTS, "metagenapp/final_asv.taxonomy"),
                     stringsAsFactors = FALSE)
colnames(ma_tax) <- c("feature","taxonomy")

ma_tbl <- read.delim(file.path(RESULTS, "metagenapp/final_asv_table.tsv"),
                     row.names = 1, check.names = FALSE)

# genus = position 6 (1-indexed), semicolon-separated
ma_tax$genus <- sapply(ma_tax$taxonomy, function(x) {
  p <- strsplit(x, ";")[[1]]
  if (length(p) >= 6) trimws(p[6]) else ""
})

# reads per feature = row sum
ma_reads <- rowSums(ma_tbl[ma_tax$feature, ], na.rm = TRUE)

ma_labels <- sapply(ma_tax$genus, label_genus)
s_ma <- summarise_resolution(ma_labels, ma_reads)
print_summary("MetagenApp", s_ma)

# genus counts (classified)
cat("  Top 10 genera (by features):\n")
top_ma <- sort(table(ma_tax$genus[ma_labels == "classified"]), decreasing = TRUE)[1:10]
print(top_ma)

# ── 2. mothur ─────────────────────────────────────────────────────────────────
mo_tax <- read.delim(
  file.path(RESULTS, "mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"),
  stringsAsFactors = FALSE
)
colnames(mo_tax) <- c("feature","size","taxonomy")

# genus = position 6, strip confidence "(nn)"
mo_tax$genus <- sapply(mo_tax$taxonomy, function(x) {
  p <- strsplit(x, ";")[[1]]
  raw <- if (length(p) >= 6) p[6] else ""
  trimws(sub("\\([0-9]+\\)", "", raw))
})

# reads per OTU = Size column
mo_reads <- mo_tax$size

# _unclassified suffix → ambiguous (propagated from family level)
mo_labels <- sapply(mo_tax$genus, function(g) {
  if (grepl("_unclassified$", g, ignore.case = TRUE)) return("ambiguous")
  label_genus(g)
})
s_mo <- summarise_resolution(mo_labels, mo_reads)
print_summary("mothur", s_mo)

# ── 3. QIIME 2 ────────────────────────────────────────────────────────────────
qi_tax <- read.delim(file.path(RESULTS, "qiime2/taxonomy_export/taxonomy.tsv"),
                     stringsAsFactors = FALSE)
colnames(qi_tax)[1:3] <- c("feature","taxonomy","confidence")

biom_obj <- read_biom(file.path(RESULTS, "qiime2/table_export/feature-table.biom"))
qi_mat   <- as(biom_data(biom_obj), "matrix")
qi_reads <- rowSums(qi_mat)  # reads per ASV

# genus: field with g__ prefix
qi_tax$genus <- sapply(qi_tax$taxonomy, function(x) {
  p <- strsplit(x, ";")[[1]]
  gf <- grep("^\\s*g__", p, value = TRUE)
  if (length(gf) == 0) return("")
  g <- trimws(sub("^\\s*g__", "", gf[1]))
  if (g == "" || g == "__") return("")
  g
})

qi_reads_named <- qi_reads[qi_tax$feature]
qi_reads_named[is.na(qi_reads_named)] <- 0

qi_labels <- sapply(qi_tax$genus, label_genus)
s_qi <- summarise_resolution(qi_labels, qi_reads_named)
print_summary("QIIME 2", s_qi)

# ── 4. Summary table ──────────────────────────────────────────────────────────
cat("\n\n══ TABLA RESUMEN (% de reads) ══════════════════════════════════════\n")
cat(sprintf("%-14s  %12s  %12s  %12s\n",
            "Pipeline", "Classified", "Ambiguous", "Unclassified"))
cat(strrep("-", 56), "\n")
for (nm in c("MetagenApp","mothur","QIIME 2")) {
  s <- list(MetagenApp=s_ma, mothur=s_mo, `QIIME 2`=s_qi)[[nm]]
  cat(sprintf("%-14s  %11.1f%%  %11.1f%%  %11.1f%%\n",
              nm, s$read["classified"], s$read["ambiguous"], s$read["unclassified"]))
}

cat("\n══ TABLA RESUMEN (% de features) ══════════════════════════════════\n")
cat(sprintf("%-14s  %12s  %12s  %12s\n",
            "Pipeline", "Classified", "Ambiguous", "Unclassified"))
cat(strrep("-", 56), "\n")
for (nm in c("MetagenApp","mothur","QIIME 2")) {
  s <- list(MetagenApp=s_ma, mothur=s_mo, `QIIME 2`=s_qi)[[nm]]
  cat(sprintf("%-14s  %11.1f%%  %11.1f%%  %11.1f%%\n",
              nm, s$feat["classified"], s$feat["ambiguous"], s$feat["unclassified"]))
}

# ── 5. Genus richness (number of unique genera classified) ────────────────────
cat("\n══ RIQUEZA DE GÉNEROS (n géneros únicos clasificados) ══════════════\n")
cat("MetagenApp:", length(unique(ma_tax$genus[ma_labels == "classified"])), "\n")
cat("mothur:    ", length(unique(mo_tax$genus[mo_labels == "classified"])), "\n")
cat("QIIME 2:   ", length(unique(qi_tax$genus[qi_labels == "classified"])), "\n")
