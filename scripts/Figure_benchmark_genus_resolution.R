#!/usr/bin/env Rscript
# Figure: genus-level resolution comparison (stacked bar)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

FIGURES <- "/home/mijailcc/MetagenApp/user_data/figures"
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)

# ── Data (from compute_genus_resolution.R) ────────────────────────────────────

# By reads (post alignment-fallback for MetagenApp)
reads_df <- data.frame(
  pipeline     = c("MetagenApp", "mothur", "QIIME 2"),
  Classified   = c(91.2, 93.6, 97.7),
  Ambiguous    = c( 3.8,  6.4,  0.3),
  Unclassified = c( 5.0,  0.0,  2.0)
)

# By features (post alignment-fallback for MetagenApp)
feat_df <- data.frame(
  pipeline     = c("MetagenApp", "mothur", "QIIME 2"),
  Classified   = c(80.0, 83.5, 94.6),
  Ambiguous    = c(11.2, 16.5,  1.0),
  Unclassified = c( 8.8,  0.0,  4.4)
)

# ── reshape to long ───────────────────────────────────────────────────────────

to_long <- function(df, metric_label) {
  df %>%
    pivot_longer(-pipeline, names_to = "category", values_to = "pct") %>%
    mutate(metric = metric_label)
}

plot_df <- bind_rows(
  to_long(reads_df, "Por lecturas"),
  to_long(feat_df,  "Por features")
) %>%
  mutate(
    pipeline = factor(pipeline, levels = c("MetagenApp", "mothur", "QIIME 2")),
    category = factor(category, levels = c("Classified", "Ambiguous", "Unclassified")),
    metric   = factor(metric, levels = c("Por lecturas", "Por features"))
  )

# ── colours ───────────────────────────────────────────────────────────────────

CAT_COLS <- c(
  Classified   = "#4E79A7",
  Ambiguous    = "#F28E2B",
  Unclassified = "#BAB0AC"
)

PIPE_COLS <- c(MetagenApp = "#4E79A7", mothur = "#F28E2B", `QIIME 2` = "#59A14F")

# ── figure ────────────────────────────────────────────────────────────────────

p <- ggplot(plot_df, aes(x = pipeline, y = pct, fill = category)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(pct >= 2, paste0(pct, "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3.2, color = "white", fontface = "bold"
  ) +
  facet_wrap(~ metric) +
  scale_fill_manual(
    values = CAT_COLS,
    labels = c("Classified" = "Clasificado a género",
               "Ambiguous"  = "Ambiguo / propagado",
               "Unclassified" = "No clasificado a género")
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0), limits = c(0, 102)) +
  labs(
    x    = NULL,
    y    = "Porcentaje (%)",
    fill = NULL,
    caption = "Clasificado: género definido sin ambigüedad. Ambiguo: etiquetas propagadas o 'Incertae Sedis'.\nNo clasificado: sin asignación a nivel género."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    strip.background   = element_rect(fill = "grey92"),
    strip.text         = element_text(face = "bold"),
    plot.caption       = element_text(size = 8, hjust = 0, color = "grey40")
  )

ggsave(file.path(FIGURES, "Figure_benchmark_genus_resolution.png"),
       p, width = 8, height = 5, dpi = 300)
cat("Figura guardada: Figure_benchmark_genus_resolution.png\n")

# ── genus richness annotation ─────────────────────────────────────────────────
cat("\nRiqueza de géneros únicos clasificados:\n")
cat("  MetagenApp: 609  (post alignment-fallback)\n")
cat("  mothur:     572\n")
cat("  QIIME 2:    154\n")
