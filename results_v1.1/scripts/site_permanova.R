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
FIGURES    <- "site_permanova"
TABLES     <- "site_permanova"
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

# ── #2/#9 checklist (30 sep): ¿cada pipeline recupera faringe vs nariz? ──
# Mismas matrices Bray-Curtis a nivel feature que la PCoA/Mantel (n = 36).
# Metadatos: SRA runinfo PRJNA1051577 (LibraryName P/N + PC/IC/NC).
meta <- read.csv("data/samples36/sample_site_carrier.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
run_site <- function(bc, tag) {
  lab <- attr(bc, "Labels"); md <- meta[lab, ]
  stopifnot(!any(is.na(md$site)))
  cat(sprintf("\n==================== %s (n = %d; faringe %d, nariz %d) ====================\n",
      tag, length(lab), sum(md$site=="pharynx"), sum(md$site=="nose")))
  set.seed(19760620)
  cat("\n-- PERMANOVA sitio (faringe vs nariz) --\n")
  print(adonis2(bc ~ site, data = md, permutations = 9999))
  set.seed(19760620)
  cat("\n-- betadisper por sitio --\n")
  bd <- betadisper(bc, factor(md$site)); print(permutest(bd, permutations = 9999))
  print(round(tapply(bd$distances, md$site, mean), 4))
  for (s in c("pharynx","nose")) {
    keep <- md$site == s
    sub <- as.dist(as.matrix(bc)[keep, keep])
    set.seed(19760620)
    cat(sprintf("\n-- PERMANOVA tipo de portador dentro de %s (n = %d) --\n", s, sum(keep)))
    print(adonis2(sub ~ carrier, data = md[keep, ], permutations = 9999))
  }
}
run_site(bc_ma, "MetagenApp"); run_site(bc_mo, "mothur"); run_site(bc_qi, "QIIME2")
cat("\n===== sessionInfo =====\n"); print(sessionInfo())
