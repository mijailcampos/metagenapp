#!/usr/bin/env Rscript
# Beta diversity benchmark: PCoA Bray-Curtis + Mantel + Procrustes + PERMANOVA + betadisper
# Adapted (7 ago 2026) from Figure_benchmark_beta_diversity.R to point at the corrected,
# consistent runs in results_v2 (SILVA 138.2 uniform, primer trimming fixed).
# Differences from the original script:
#   - paths point to data/samples36/{metagenapp,mothur,qiime2} (Zenodo 10.5281/zenodo.23349702)
#   - QIIME2 table read from feature-table.tsv instead of feature-table.biom (avoids
#     the biomformat dependency, which is not installed in this env)
#   - 3-panel PCoA figure built with base grid instead of patchwork (not installed)
#   - PERMANOVA/Mantel/Procrustes/betadisper results also written to a text summary
#     in results_v2/tables/ for pulling exact numbers into the poster

suppressPackageStartupMessages({
  library(vegan)
  library(ape)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})
set.seed(19760620)  # agregado 30 sep: permutaciones reproducibles (betadisper, etc.)

RESULTS_MA <- "data/samples36/metagenapp"
RESULTS_MO <- "data/samples36/mothur"
RESULTS_QI <- "data/samples36/qiime2"
FIGURES    <- "figures"
TABLES     <- "tables"
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)

COLS <- c(MetagenApp = "#4E79A7", mothur = "#F28E2B", QIIME2 = "#59A14F")

summary_lines <- c()
log_line <- function(...) {
  txt <- sprintf(...)
  cat(txt, "\n")
  summary_lines <<- c(summary_lines, txt)
}

# ── 1. Load feature tables ─────────────────────────────────────────────────

# MetagenApp: rows=features, cols=samples
ma_raw <- read.delim(file.path(RESULTS_MA, "final_asv_table.tsv"),
                     row.names = 1, check.names = FALSE)
colnames(ma_raw) <- sub("_.+", "", colnames(ma_raw))   # Bustos10_S310_L001 -> Bustos10

# mothur: rows=samples (label, Group, numOtus, OTUs...)
mo_raw <- read.delim(
  file.path(RESULTS_MO, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared"),
  check.names = FALSE
)
rownames(mo_raw) <- mo_raw$Group                         # Bustos1..36
mo_raw <- mo_raw[, -(1:3)]                                # drop label, Group, numOtus

# QIIME2: TSV export (skip the "# Constructed from biom file" comment line)
qi_raw <- read.delim(file.path(RESULTS_QI, "feature-table.tsv"),
                     skip = 1, row.names = 1, check.names = FALSE)
colnames(qi_raw) <- sub("_S[0-9]+$", "", colnames(qi_raw))  # Bustos1_S297 -> Bustos1
qi_raw <- t(as.matrix(qi_raw))                              # rows=samples

# ── 2. Harmonise samples & relative abundance ──────────────────────────────

common_samples <- Reduce(intersect,
                          list(colnames(ma_raw), rownames(mo_raw), rownames(qi_raw)))
common_samples <- sort(common_samples)
log_line("Common samples: %d", length(common_samples))

ma <- t(ma_raw[, common_samples])  # rows=samples
mo <- mo_raw[common_samples, ]
qi <- qi_raw[common_samples, ]

rel_ab <- function(m) sweep(m, 1, rowSums(m), "/")
ma_rel <- rel_ab(ma)
mo_rel <- rel_ab(mo)
qi_rel <- rel_ab(qi)

# ── 3. Bray-Curtis distance matrices ────────────────────────────────────────

bc_ma <- vegdist(ma_rel, method = "bray")
bc_mo <- vegdist(mo_rel, method = "bray")
bc_qi <- vegdist(qi_rel, method = "bray")

# ── 4. PCoA ──────────────────────────────────────────────────────────────────

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

cat(sprintf("PCoA var MetagenApp: PC1=%.1f%% PC2=%.1f%%\n", pcoa_ma$var1, pcoa_ma$var2))
cat(sprintf("PCoA var mothur:     PC1=%.1f%% PC2=%.1f%%\n", pcoa_mo$var1, pcoa_mo$var2))
cat(sprintf("PCoA var QIIME2:     PC1=%.1f%% PC2=%.1f%%\n", pcoa_qi$var1, pcoa_qi$var2))
write.csv(pcoa_ma$df, file.path(TABLES, "pcoa_coords_metagenapp.csv"), row.names = FALSE)
write.csv(pcoa_mo$df, file.path(TABLES, "pcoa_coords_mothur.csv"), row.names = FALSE)
write.csv(pcoa_qi$df, file.path(TABLES, "pcoa_coords_qiime2.csv"), row.names = FALSE)


write.csv(pcoa_ma$df, file.path(TABLES, "pcoa_coords_metagenapp.csv"), row.names = FALSE)
write.csv(pcoa_mo$df, file.path(TABLES, "pcoa_coords_mothur.csv"), row.names = FALSE)
write.csv(pcoa_qi$df, file.path(TABLES, "pcoa_coords_qiime2.csv"), row.names = FALSE)

# ── 5. Figure: 3-panel PCoA (base grid, no patchwork) ───────────────────────

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
p_qi <- make_panel(pcoa_qi, "QIIME 2",    COLS["QIIME2"])

png(file.path(FIGURES, "Figure_benchmark_beta_diversity.png"),
    width = 12, height = 4.4, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3, heights = unit(c(1, 0.12), "null"))))
print(p_ma, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_mo, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_qi, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:3))
grid.text("Bray-Curtis PCoA. 95% confidence ellipse. Relative abundance; n = 36 samples per pipeline.",
          gp = gpar(fontsize = 9), just = "left", x = unit(0.02, "npc"))
popViewport()
dev.off()
log_line("Figura guardada: Figure_benchmark_beta_diversity.png")

# ── 6. Mantel tests (concordancia entre pipelines) ──────────────────────────

log_line("\n== Mantel tests (Bray-Curtis) ==")
mant_ma_mo <- mantel(bc_ma, bc_mo, permutations = 9999)
mant_ma_qi <- mantel(bc_ma, bc_qi, permutations = 9999)
mant_mo_qi <- mantel(bc_mo, bc_qi, permutations = 9999)

print_mantel <- function(label, m) log_line("%-28s  r = %.3f  p = %.4f", label, m$statistic, m$signif)
print_mantel("MetagenApp vs mothur", mant_ma_mo)
print_mantel("MetagenApp vs QIIME2", mant_ma_qi)
print_mantel("mothur vs QIIME2",     mant_mo_qi)

# ── 7. Procrustes tests ──────────────────────────────────────────────────────

log_line("\n== Procrustes (protest) ==")
prot_ma_mo <- protest(pcoa_ma$df[, c("PC1","PC2")], pcoa_mo$df[, c("PC1","PC2")], permutations = 9999)
prot_ma_qi <- protest(pcoa_ma$df[, c("PC1","PC2")], pcoa_qi$df[, c("PC1","PC2")], permutations = 9999)
prot_mo_qi <- protest(pcoa_mo$df[, c("PC1","PC2")], pcoa_qi$df[, c("PC1","PC2")], permutations = 9999)

print_prot <- function(label, p) log_line("%-28s  M2 = %.4f  p = %.4f", label, p$ss, p$signif)
print_prot("MetagenApp vs mothur", prot_ma_mo)
print_prot("MetagenApp vs QIIME2", prot_ma_qi)
print_prot("mothur vs QIIME2",     prot_mo_qi)

# ── 8. PERMANOVA combinado (nivel filo, pipeline como factor) ────────────────

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

# MetagenApp
ma_tax <- read.delim(file.path(RESULTS_MA, "final_asv.taxonomy"),
                     header = TRUE, stringsAsFactors = FALSE)
colnames(ma_tax) <- c("feature", "taxonomy")
ma_phy_vec <- setNames(sapply(ma_tax$taxonomy, parse_phylum_ma), ma_tax$feature)
ma_mat <- as.matrix(ma_raw[, common_samples])
ma_ph_rel <- agg_to_phylum(ma_mat, ma_phy_vec)

# QIIME2
qi_tax <- read.delim(file.path(RESULTS_QI, "taxonomy.tsv"),
                     stringsAsFactors = FALSE)
colnames(qi_tax)[1:2] <- c("feature","taxonomy")
qi_phy_vec <- setNames(sapply(qi_tax$taxonomy, parse_phylum_qi), qi_tax$feature)
qi_mat <- t(qi_raw[common_samples, ])
qi_ph_rel <- agg_to_phylum(qi_mat, qi_phy_vec)

# mothur
mo_tax <- read.delim(
  file.path(RESULTS_MO, "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"),
  stringsAsFactors = FALSE
)
colnames(mo_tax) <- c("feature","size","taxonomy")
mo_phy_vec <- setNames(sapply(as.character(mo_tax$taxonomy), parse_phylum_mo), mo_tax$feature)
mo_mat <- t(as.matrix(mo[common_samples, ]))
mo_ph_rel <- agg_to_phylum(mo_mat, mo_phy_vec)

common_phyla <- Reduce(intersect,
                        list(colnames(ma_ph_rel), colnames(mo_ph_rel), colnames(qi_ph_rel)))
common_phyla <- common_phyla[common_phyla != "Unclassified"]
log_line("\nFilos en comun para PERMANOVA: %d", length(common_phyla))

build_block <- function(ph_rel, phyla, tag) {
  m <- ph_rel[, phyla, drop = FALSE]
  keep <- rowSums(m) > 0
  if (!all(keep)) log_line("  [%s] dropping %d zero-sum samples", tag, sum(!keep))
  m[keep, ]
}
ma_block <- build_block(ma_ph_rel, common_phyla, "MetagenApp")
mo_block <- build_block(mo_ph_rel, common_phyla, "mothur")
qi_block <- build_block(qi_ph_rel, common_phyla, "QIIME2")

shared <- Reduce(intersect, list(rownames(ma_block), rownames(mo_block), rownames(qi_block)))
log_line("  Samples in combined PERMANOVA: %d", length(shared))

combined <- rbind(ma_block[shared,], mo_block[shared,], qi_block[shared,])
pipeline_grp <- factor(rep(c("MetagenApp","mothur","QIIME2"), each = length(shared)))
sample_grp   <- factor(rep(shared, times = 3))

bc_combined <- vegdist(combined, method = "bray")

log_line("\n== PERMANOVA combinado (phylum-level, Bray-Curtis) ==")
perm_res <- adonis2(bc_combined ~ pipeline_grp + sample_grp, permutations = 9999, by = "margin")
summary_lines <- c(summary_lines, capture.output(print(perm_res)))
print(perm_res)

log_line("\n== Betadisper (dispersion por pipeline) ==")
bd <- betadisper(bc_combined, pipeline_grp)
pmt <- permutest(bd, permutations = 9999)
summary_lines <- c(summary_lines, capture.output(print(pmt)))
print(pmt)
log_line("\nDistancias medias al centroide:")
centroid_dist <- round(tapply(bd$distances, pipeline_grp, mean), 4)
summary_lines <- c(summary_lines, capture.output(print(centroid_dist)))
print(centroid_dist)

writeLines(summary_lines, file.path(TABLES, "beta_diversity_summary.txt"))
cat("\nResumen guardado en:", file.path(TABLES, "beta_diversity_summary.txt"), "\n")

cat("\n===== sessionInfo =====\n"); print(sessionInfo())
