#!/usr/bin/env Rscript
# Beta diversity benchmark: PCoA Bray-Curtis + Mantel + Procrustes + PERMANOVA
# Run: results_20260418_0205

suppressPackageStartupMessages({
  library(vegan)
  library(ape)
  library(biomformat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

RESULTS <- "/data/results/runs/01_benchmarking/benchmark_controlled/results_20260418_0205"
FIGURES <- "/home/mijailcc/MetagenApp/user_data/figures"
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)

COLS <- c(MetagenApp = "#4E79A7", mothur = "#F28E2B", QIIME2 = "#59A14F")

# ── 1. Load feature tables ─────────────────────────────────────────────────────

# MetagenApp: rows=features, cols=samples
ma_raw <- read.delim(file.path(RESULTS, "metagenapp/final_asv_table.tsv"),
                     row.names = 1, check.names = FALSE)
colnames(ma_raw) <- sub("_.+", "", colnames(ma_raw))   # Bustos10_S310_L001 -> Bustos10

# mothur: rows=samples (label, Group, numOtus, OTUs...)
mo_raw <- read.delim(
  file.path(RESULTS, "mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared"),
  check.names = FALSE
)
rownames(mo_raw) <- mo_raw$Group                         # Bustos1..36
mo_raw <- mo_raw[, -(1:3)]                               # drop label, Group, numOtus

# QIIME2: BIOM
biom_obj <- read_biom(file.path(RESULTS, "qiime2/table_export/feature-table.biom"))
qi_raw   <- as(biom_data(biom_obj), "matrix")
colnames(qi_raw) <- sub("_S[0-9]+$", "", colnames(qi_raw))  # Bustos1_S297 -> Bustos1
qi_raw <- t(qi_raw)                                       # rows=samples

# ── 2. Harmonise samples & relative abundance ──────────────────────────────────

# Keep only samples present in all three pipelines
common_samples <- Reduce(intersect,
                          list(colnames(ma_raw), rownames(mo_raw), rownames(qi_raw)))
common_samples <- sort(common_samples)
cat("Common samples:", length(common_samples), "\n")

ma <- t(ma_raw[, common_samples])  # rows=samples
mo <- mo_raw[common_samples, ]
qi <- qi_raw[common_samples, ]

# Relative abundance (proportional)
rel_ab <- function(m) sweep(m, 1, rowSums(m), "/")
ma_rel <- rel_ab(ma)
mo_rel <- rel_ab(mo)
qi_rel <- rel_ab(qi)

# ── 3. Bray-Curtis distance matrices ──────────────────────────────────────────

bc_ma <- vegdist(ma_rel, method = "bray")
bc_mo <- vegdist(mo_rel, method = "bray")
bc_qi <- vegdist(qi_rel, method = "bray")

# ── 4. PCoA ───────────────────────────────────────────────────────────────────

run_pcoa <- function(d, label) {
  pc  <- pcoa(d)
  var <- pc$values$Relative_eig * 100
  df  <- data.frame(
    sample   = rownames(as.matrix(d)),
    PC1      = pc$vectors[, 1],
    PC2      = pc$vectors[, 2],
    pipeline = label,
    stringsAsFactors = FALSE
  )
  list(df = df, var1 = round(var[1], 1), var2 = round(var[2], 1))
}

pcoa_ma <- run_pcoa(bc_ma, "MetagenApp")
pcoa_mo <- run_pcoa(bc_mo, "mothur")
pcoa_qi <- run_pcoa(bc_qi, "QIIME 2")

# ── 5. Figure: 3-panel PCoA ───────────────────────────────────────────────────

# Compute common axis limits for visual comparability
all_pc1 <- c(pcoa_ma$df$PC1, pcoa_mo$df$PC1, pcoa_qi$df$PC1)
all_pc2 <- c(pcoa_ma$df$PC2, pcoa_mo$df$PC2, pcoa_qi$df$PC2)
xlim <- range(all_pc1) * 1.12
ylim <- range(all_pc2) * 1.12

make_panel <- function(pcr, title, color) {
  ggplot(pcr$df, aes(PC1, PC2)) +
    geom_point(color = color, size = 2.5, alpha = 0.85) +
    stat_ellipse(color = color, linetype = "dashed", linewidth = 0.5, level = 0.95) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      title = title,
      x = paste0("PC1 (", pcr$var1, "%)"),
      y = paste0("PC2 (", pcr$var2, "%)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

p_ma <- make_panel(pcoa_ma, "MetagenApp", COLS["MetagenApp"])
p_mo <- make_panel(pcoa_mo, "mothur",     COLS["mothur"])
p_qi <- make_panel(pcoa_qi, "QIIME 2",   COLS["QIIME2"])

fig <- p_ma + p_mo + p_qi +
  plot_annotation(
    caption = "Bray-Curtis PCoA. Elipse de confianza 95%. Abundancias relativas; n = 36 muestras por pipeline.",
    theme   = theme(plot.caption = element_text(size = 9, hjust = 0))
  )

ggsave(file.path(FIGURES, "Figure_benchmark_beta_diversity.png"),
       fig, width = 12, height = 4.2, dpi = 300)
cat("Figura guardada: Figure_benchmark_beta_diversity.png\n")

# ── 6. Mantel tests (concordance entre pipelines) ─────────────────────────────

cat("\n── Mantel tests (Bray-Curtis) ──────────────────────────────────────\n")

mant_ma_mo <- mantel(bc_ma, bc_mo, permutations = 9999)
mant_ma_qi <- mantel(bc_ma, bc_qi, permutations = 9999)
mant_mo_qi <- mantel(bc_mo, bc_qi, permutations = 9999)

print_mantel <- function(label, m) {
  cat(sprintf("%-28s  r = %.3f  p = %.4f\n", label, m$statistic, m$signif))
}
print_mantel("MetagenApp vs mothur", mant_ma_mo)
print_mantel("MetagenApp vs QIIME2", mant_ma_qi)
print_mantel("mothur vs QIIME2",     mant_mo_qi)

# ── 7. Procrustes tests ────────────────────────────────────────────────────────

cat("\n── Procrustes (protest) ────────────────────────────────────────────\n")

prot_ma_mo <- protest(pcoa_ma$df[, c("PC1","PC2")], pcoa_mo$df[, c("PC1","PC2")],
                      permutations = 9999)
prot_ma_qi <- protest(pcoa_ma$df[, c("PC1","PC2")], pcoa_qi$df[, c("PC1","PC2")],
                      permutations = 9999)
prot_mo_qi <- protest(pcoa_mo$df[, c("PC1","PC2")], pcoa_qi$df[, c("PC1","PC2")],
                      permutations = 9999)

print_prot <- function(label, p) {
  cat(sprintf("%-28s  M² = %.4f  p = %.4f\n", label, p$ss, p$signif))
}
print_prot("MetagenApp vs mothur", prot_ma_mo)
print_prot("MetagenApp vs QIIME2", prot_ma_qi)
print_prot("mothur vs QIIME2",     prot_mo_qi)

# ── 8. PERMANOVA combinado (nivel filo, pipeline como factor) ──────────────────
# Agrega cada pipeline a filo → une tablas → testa si pipeline explica varianza

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

# Helper: rows=features, cols=samples → aggregate by phylum → rows=samples, cols=phyla
agg_to_phylum <- function(counts, phylum_vec) {
  # counts: matrix rows=features, cols=samples
  # phylum_vec: named vector feature→phylum (same names as rownames(counts))
  ph <- phylum_vec[rownames(counts)]
  ph[is.na(ph)] <- "Unclassified"
  agg <- rowsum(counts, ph)          # rows=phyla, cols=samples
  rel_ab(t(agg))                     # rows=samples, cols=phyla
}

# MetagenApp
ma_tax <- read.delim(file.path(RESULTS, "metagenapp/final_asv.taxonomy"),
                     header = TRUE, stringsAsFactors = FALSE)
colnames(ma_tax) <- c("feature", "taxonomy")
ma_phy_vec <- setNames(sapply(ma_tax$taxonomy, parse_phylum_ma), ma_tax$feature)
ma_mat <- as.matrix(ma_raw[, common_samples])   # rows=features, cols=samples
ma_ph_rel <- agg_to_phylum(ma_mat, ma_phy_vec)

# QIIME2
qi_tax <- read.delim(file.path(RESULTS, "qiime2/taxonomy_export/taxonomy.tsv"),
                     stringsAsFactors = FALSE)
colnames(qi_tax)[1:2] <- c("feature","taxonomy")
qi_phy_vec <- setNames(sapply(qi_tax$taxonomy, parse_phylum_qi), qi_tax$feature)
qi_mat <- t(qi_raw[common_samples, ])            # rows=features, cols=samples
qi_ph_rel <- agg_to_phylum(qi_mat, qi_phy_vec)

# mothur: cols = OTU, Size, Taxonomy
mo_tax <- read.delim(
  file.path(RESULTS, "mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"),
  stringsAsFactors = FALSE
)
colnames(mo_tax) <- c("feature","size","taxonomy")
mo_phy_vec <- setNames(sapply(as.character(mo_tax$taxonomy), parse_phylum_mo), mo_tax$feature)
mo_mat <- t(as.matrix(mo[common_samples, ]))     # rows=features (OTUs), cols=samples
mo_ph_rel <- agg_to_phylum(mo_mat, mo_phy_vec)

# Align columns (phyla present in all three)
common_phyla <- Reduce(intersect,
                        list(colnames(ma_ph_rel), colnames(mo_ph_rel), colnames(qi_ph_rel)))
common_phyla <- common_phyla[common_phyla != "Unclassified"]
cat("\nFilos en común para PERMANOVA:", length(common_phyla), "\n")

# Build combined matrix — drop samples with zero total at phylum level
build_block <- function(ph_rel, phyla, tag) {
  m <- ph_rel[, phyla, drop = FALSE]
  keep <- rowSums(m) > 0
  if (!all(keep)) {
    cat(sprintf("  [%s] dropping %d zero-sum samples\n", tag, sum(!keep)))
  }
  m[keep, ]
}
ma_block <- build_block(ma_ph_rel, common_phyla, "MetagenApp")
mo_block <- build_block(mo_ph_rel, common_phyla, "mothur")
qi_block <- build_block(qi_ph_rel, common_phyla, "QIIME2")

# Keep only samples present in all three blocks
shared <- Reduce(intersect, list(rownames(ma_block), rownames(mo_block), rownames(qi_block)))
cat("  Samples in combined PERMANOVA:", length(shared), "\n")

combined <- rbind(ma_block[shared,], mo_block[shared,], qi_block[shared,])
pipeline_grp <- factor(rep(c("MetagenApp","mothur","QIIME2"), each = length(shared)))
sample_grp   <- factor(rep(shared, times = 3))

bc_combined <- vegdist(combined, method = "bray")

cat("\n── PERMANOVA combinado (phylum-level, Bray-Curtis) ─────────────────\n")
perm_res <- adonis2(bc_combined ~ pipeline_grp + sample_grp,
                    permutations = 9999, by = "margin")
print(perm_res)

# betadisper (homogeneidad de dispersión)
cat("\n── Betadisper (dispersión por pipeline) ────────────────────────────\n")
bd <- betadisper(bc_combined, pipeline_grp)
print(permutest(bd, permutations = 9999))
cat("\nDistancias medias al centroide:\n")
print(round(tapply(bd$distances, pipeline_grp, mean), 4))
