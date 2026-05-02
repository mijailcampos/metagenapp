########################################################
# Figure — Alpha diversity boxplots (3 panels)
# Benchmark controlado: MetagenApp vs mothur vs QIIME 2
#
# Rscript scripts/Figure_benchmark_alpha_diversity.R
########################################################

library(tidyverse)

BENCH   <- "/data/results/runs/01_benchmarking/benchmark_controlled/results_20260418_0205"
OUT_DIR <- file.path(normalizePath("~"), "MetagenApp", "user_data", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Funciones de diversidad
# ============================================================

chao1 <- function(counts) {
  counts <- counts[counts > 0]
  obs <- length(counts)
  n1  <- sum(counts == 1)
  n2  <- sum(counts == 2)
  if (n2 == 0) obs + n1*(n1-1)/2 else obs + n1^2 / (2*n2)
}

shannon <- function(counts) {
  counts <- counts[counts > 0]
  p <- counts / sum(counts)
  -sum(p * log(p))
}

inv_simpson <- function(counts) {
  counts <- counts[counts > 0]
  p <- counts / sum(counts)
  1 / sum(p^2)
}


# ============================================================
# 1. MetagenApp
# ============================================================

meta_tbl <- read.table(
  file.path(BENCH, "metagenapp/final_asv_table.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE
)
sample_cols_meta <- setdiff(colnames(meta_tbl), "Sequence_ID")

meta_alpha <- tibble(
  Sample   = sample_cols_meta,
  Pipeline = "MetagenApp",
  Chao1    = sapply(sample_cols_meta, function(s) chao1(meta_tbl[[s]])),
  Shannon  = sapply(sample_cols_meta, function(s) shannon(meta_tbl[[s]])),
  InvSimp  = sapply(sample_cols_meta, function(s) inv_simpson(meta_tbl[[s]]))
)


# ============================================================
# 2. mothur
# ============================================================

mothur_shared <- read.table(
  file.path(BENCH, "mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared"),
  header = TRUE, sep = "\t", check.names = FALSE
)
otu_cols <- colnames(mothur_shared)[-(1:3)]

mothur_alpha <- mothur_shared %>%
  select(Group, all_of(otu_cols)) %>%
  rowwise() %>%
  mutate(
    Pipeline = "mothur",
    Sample   = Group,
    Chao1    = chao1(c_across(all_of(otu_cols))),
    Shannon  = shannon(c_across(all_of(otu_cols))),
    InvSimp  = inv_simpson(c_across(all_of(otu_cols)))
  ) %>%
  ungroup() %>%
  select(Sample, Pipeline, Chao1, Shannon, InvSimp)


# ============================================================
# 3. QIIME 2
# ============================================================

q2_tbl <- read.table(
  "/tmp/q2table.tsv",
  header = TRUE, sep = "\t", comment.char = "", skip = 1,
  check.names = FALSE
)
colnames(q2_tbl)[1] <- "FeatureID"
sample_cols_q2 <- setdiff(colnames(q2_tbl), "FeatureID")

q2_alpha <- tibble(
  Sample   = sample_cols_q2,
  Pipeline = "QIIME 2",
  Chao1    = sapply(sample_cols_q2, function(s) chao1(q2_tbl[[s]])),
  Shannon  = sapply(sample_cols_q2, function(s) shannon(q2_tbl[[s]])),
  InvSimp  = sapply(sample_cols_q2, function(s) inv_simpson(q2_tbl[[s]]))
)


# ============================================================
# 4. Combinar y formato largo para facets
# ============================================================

pipeline_order <- c("MetagenApp", "mothur", "QIIME 2")

alpha_long <- bind_rows(meta_alpha, mothur_alpha, q2_alpha) %>%
  mutate(Pipeline = factor(Pipeline, levels = pipeline_order)) %>%
  pivot_longer(cols = c(Chao1, Shannon, InvSimp),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
    "Chao1"   = "Chao1",
    "Shannon" = "Shannon (H)",
    "InvSimp" = "Inv. Simpson"
  )) %>%
  mutate(Metric = factor(Metric, levels = c("Chao1", "Shannon (H)", "Inv. Simpson")))


# ============================================================
# 5. Paleta por pipeline
# ============================================================

pipeline_colors <- c(
  "MetagenApp" = "#4E79A7",
  "mothur"     = "#F28E2B",
  "QIIME 2"    = "#59A14F"
)


# ============================================================
# 6. Figura
# ============================================================

p <- ggplot(alpha_long, aes(x = Pipeline, y = Value, fill = Pipeline)) +
  geom_boxplot(width = 0.55, outlier.shape = 21, outlier.size = 1.8,
               outlier.fill = "white", linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.0, alpha = 0.45, color = "grey30") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pipeline_colors) +
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.1)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", size = 11),
    strip.background   = element_rect(fill = "grey95", color = NA),
    axis.text.x        = element_text(face = "bold", size = 10, color = "black"),
    axis.text.y        = element_text(size = 9, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor   = element_blank(),
    legend.position    = "none"
  )


# ============================================================
# 7. Guardar
# ============================================================

out_png <- file.path(OUT_DIR, "Figure_benchmark_alpha_diversity.png")
out_pdf <- file.path(OUT_DIR, "Figure_benchmark_alpha_diversity.pdf")

ggsave(out_png, plot = p, width = 9, height = 4, dpi = 300)
ggsave(out_pdf, plot = p, width = 9, height = 4)

# Resumen estadístico
cat("\n=== Resumen alpha diversity (mediana) ===\n")
alpha_long %>%
  group_by(Pipeline, Metric) %>%
  summarise(median = round(median(Value), 2),
            mean   = round(mean(Value), 2),
            min    = round(min(Value), 2),
            max    = round(max(Value), 2),
            .groups = "drop") %>%
  print(n = 30)

cat("\nFigura guardada en:", out_png, "\n")
