#!/usr/bin/env Rscript
# Figure 2 - Top 10 genus comparison across MetagenApp, Mothur, and QIIME2
# Dataset: nasopharyngeal samples (n=36), Samuel et al.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
SILVA_RUN  <- "/data/results/runs/ref_runs/faringe_nariz_20260414_1645"
MOTH_RUN   <- "/data/results/runs/01_benchmarking/paper01_faringe_mothur"
Q2_RUN     <- "/data/results/runs/01_benchmarking/qiime2_faringe"

META_TAX   <- file.path(SILVA_RUN, "final_asv.taxonomy")
META_TABLE <- file.path(SILVA_RUN, "final_asv_table.tsv")
MOTH_TAX   <- file.path(MOTH_RUN,  "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy")
MOTH_TABLE <- file.path(MOTH_RUN,  "stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared")
Q2_TAX     <- file.path(Q2_RUN,    "taxonomy_export/taxonomy.tsv")
Q2_TABLE   <- file.path(Q2_RUN,    "table.tsv")
OUT_DIR    <- "/data/projects/paper_metagenapp_benchmark/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N <- 10  # top genera to show individually

# ── Helper: collapse reads per genus ──────────────────────────────────────────
reads_per_genus <- function(otu_ids, genus_vec, count_mat) {
  df <- data.frame(otu = otu_ids, genus = genus_vec, stringsAsFactors = FALSE)
  totals <- rowSums(count_mat)
  df$reads <- totals
  agg <- aggregate(reads ~ genus, data = df, sum)
  agg
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. MetagenApp
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading MetagenApp (SILVA v1)...\n")

# final_asv.taxonomy tiene header: ASV \t Taxonomía
meta_tax_raw <- read.table(META_TAX, sep = "\t", header = TRUE,
                           col.names = c("otu", "taxonomy"), quote = "")
meta_table   <- read.table(META_TABLE, sep = "\t", header = TRUE,
                           row.names = 1, check.names = FALSE)

parse_genus_meta <- function(tax_str) {
  parts <- str_split(tax_str, ";")[[1]]
  parts <- trimws(parts[parts != ""])
  depth <- length(parts)
  # depth 1 = solo "Bacteria" o "Unclassified" → inclasificable
  if (depth <= 1) return("Unclassifiable")
  # depth 2-5 = clasificado a phylum/clase/orden/familia, no llegó a género
  if (depth < 6)  return("Classified (above genus)")
  # depth 6 = tiene género
  g <- parts[6]
  if (g != "" && !grepl("unclassified|uncultured", g, ignore.case = TRUE)) return(g)
  return("Classified (above genus)")
}

meta_tax_raw$genus <- sapply(meta_tax_raw$taxonomy, parse_genus_meta)

common_otus <- intersect(meta_tax_raw$otu, rownames(meta_table))
meta_tax_f  <- meta_tax_raw[meta_tax_raw$otu %in% common_otus, ]
meta_mat    <- as.matrix(meta_table[meta_tax_f$otu, ])

meta_genus <- reads_per_genus(meta_tax_f$otu, meta_tax_f$genus, meta_mat)
meta_genus$pipeline <- "MetagenApp (SILVA v1)"

# ══════════════════════════════════════════════════════════════════════════════
# 2. Mothur
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading Mothur...\n")

moth_tax_raw <- read.table(MOTH_TAX, sep = "\t", header = TRUE,
                           col.names = c("otu", "size", "taxonomy"), quote = "")
moth_shared  <- read.table(MOTH_TABLE, sep = "\t", header = TRUE,
                           check.names = FALSE)

parse_genus_mothur <- function(tax_str) {
  clean <- gsub("\\([0-9]+\\)", "", tax_str)
  parts <- trimws(str_split(clean, ";")[[1]])
  parts <- parts[parts != ""]
  depth <- length(parts)
  if (depth <= 1) return("Unclassifiable")
  if (depth < 6)  return("Classified (above genus)")
  g <- gsub("_unclassified$", "", parts[6])
  # Mothur fuerza género con sufijo _unclassified → reportar como forzado
  if (grepl("_unclassified$", parts[6], ignore.case = TRUE)) return("Forced genus")
  if (g != "" && !grepl("^unclassified$", g, ignore.case = TRUE)) return(g)
  return("Classified (above genus)")
}

moth_tax_raw$genus <- sapply(moth_tax_raw$taxonomy, parse_genus_mothur)

otu_cols    <- setdiff(colnames(moth_shared), c("label", "Group", "numOtus"))
moth_mat    <- as.matrix(moth_shared[, otu_cols])
moth_totals <- colSums(moth_mat)

moth_genus_df <- data.frame(
  otu   = moth_tax_raw$otu,
  genus = moth_tax_raw$genus,
  reads = as.numeric(moth_totals[moth_tax_raw$otu]),
  stringsAsFactors = FALSE
)
moth_genus_df <- moth_genus_df[!is.na(moth_genus_df$reads), ]
moth_genus <- aggregate(reads ~ genus, data = moth_genus_df, sum)
moth_genus$pipeline <- "Mothur"

# ══════════════════════════════════════════════════════════════════════════════
# 3. QIIME2
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading QIIME2...\n")

q2_tax <- read.table(Q2_TAX, sep = "\t", header = TRUE, comment.char = "",
                     quote = "", stringsAsFactors = FALSE)
q2_tax <- q2_tax[q2_tax[, 1] != "#q2:types", ]
colnames(q2_tax)[1:2] <- c("feature_id", "taxon")

q2_table_raw <- read.table(Q2_TABLE, sep = "\t", header = TRUE,
                            skip = 1, row.names = 1, check.names = FALSE,
                            comment.char = "")
q2_mat <- as.matrix(q2_table_raw)

parse_genus_qiime2 <- function(tax_str) {
  parts <- trimws(str_split(tax_str, ";")[[1]])
  # QIIME2 usa prefijos d__ p__ c__ o__ f__ g__ s__
  # Contar niveles con contenido real (sin prefijo vacío)
  real_parts <- parts[sapply(parts, function(p) {
    label <- sub("^[a-z]__", "", p)
    nchar(trimws(label)) > 0
  })]
  depth <- length(real_parts)
  if (depth <= 1) return("Unclassifiable")
  genus_parts <- parts[grepl("^g__", parts)]
  if (length(genus_parts) > 0) {
    g <- sub("^g__", "", genus_parts[1])
    if (trimws(g) != "") return(trimws(g))
  }
  return("Classified (above genus)")
}

q2_tax$genus <- sapply(q2_tax$taxon, parse_genus_qiime2)

common_q2 <- intersect(q2_tax$feature_id, rownames(q2_mat))
q2_tax_f  <- q2_tax[q2_tax$feature_id %in% common_q2, ]
q2_counts <- rowSums(q2_mat[q2_tax_f$feature_id, , drop = FALSE])

q2_genus_df <- data.frame(
  genus = q2_tax_f$genus,
  reads = as.numeric(q2_counts),
  stringsAsFactors = FALSE
)
q2_genus <- aggregate(reads ~ genus, data = q2_genus_df, sum)
q2_genus$pipeline <- "QIIME2"

# ══════════════════════════════════════════════════════════════════════════════
# 4. Merge and compute relative abundance
# ══════════════════════════════════════════════════════════════════════════════
cat("Computing relative abundances...\n")

all_genus <- bind_rows(meta_genus, moth_genus, q2_genus)

all_genus <- all_genus %>%
  group_by(pipeline) %>%
  mutate(rel_abund = reads / sum(reads) * 100) %>%
  ungroup()

# Top N genuine genera (excluir categorías especiales del ranking)
EXCLUDE_FROM_RANK <- c("Other", "Unresolved", "Unclassified",
                       "Classified (above genus)", "Unclassifiable", "Forced genus")
top_genera <- all_genus %>%
  filter(!genus %in% EXCLUDE_FROM_RANK) %>%
  group_by(genus) %>%
  summarise(mean_abund = mean(rel_abund), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = TOP_N) %>%
  pull(genus)

cat("Top", TOP_N, "genera:\n")
print(top_genera)

SPECIAL <- c("Classified (above genus)", "Unclassifiable", "Forced genus")

all_genus <- all_genus %>%
  mutate(genus_plot = case_when(
    genus %in% top_genera  ~ genus,
    genus %in% SPECIAL     ~ "Not resolved to genus",
    TRUE                   ~ "Other"
  ))

plot_data <- all_genus %>%
  group_by(pipeline, genus_plot) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

# ── Color palette ─────────────────────────────────────────────────────────────
genus_colors <- c(
  "Staphylococcus"          = "#E41A1C",
  "Pseudomonas"             = "#FF7F00",
  "Stenotrophomonas"        = "#FFCC00",
  "Veillonella"             = "#4DAF4A",
  "Prevotella"              = "#984EA3",
  "Streptococcus"           = "#377EB8",
  "Haemophilus"             = "#A65628",
  "Neisseria"               = "#F781BF",
  "Fusobacterium"           = "#66C2A5",
  "Porphyromonas"           = "#FC8D62",
  "Rothia"                  = "#8DA0CB",
  "Moraxella"               = "#E78AC3",
  "Corynebacterium"         = "#A6D854",
  "Incertae Sedis"          = "#B3B3B3",
  "Clostridioides"          = "#80B1D3",
  "Other"                   = "#CCCCCC",
  "Not resolved to genus"   = "#EBEBEB"   # gris claro: no llegó a género (ver Tabla S1)
)

present_genera <- unique(plot_data$genus_plot)
extra_colors   <- scales::hue_pal()(length(present_genera))
names(extra_colors) <- present_genera
all_colors  <- c(genus_colors, extra_colors)
final_colors <- all_colors[!duplicated(names(all_colors))]
final_colors <- final_colors[present_genera]

# Order: top genera primero, luego Other, luego "Not resolved" al fondo
BOTTOM_CATS <- c("Other", "Not resolved to genus")
genus_order <- c(
  top_genera[top_genera %in% present_genera & !top_genera %in% BOTTOM_CATS],
  intersect(BOTTOM_CATS, present_genera)
)

plot_data$genus_plot <- factor(plot_data$genus_plot, levels = rev(genus_order))
plot_data$pipeline   <- factor(plot_data$pipeline,
                                levels = c("MetagenApp (SILVA v1)", "Mothur", "QIIME2"))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Plot
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 2...\n")

fig2 <- ggplot(plot_data, aes(x = pipeline, y = rel_abund, fill = genus_plot)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6,
           color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = final_colors,
    breaks = genus_order,
    name   = "Genus"
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "Figure 2. Genus-level taxonomic composition",
    subtitle = paste0("Top ", TOP_N, " genera · 36 nasopharyngeal samples\n",
                      "Grey = not resolved to genus level (see Table S1 for breakdown)"),
    x        = NULL,
    y        = "Relative abundance (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "grey40", size = 11),
    axis.text.x        = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12),
    legend.title       = element_text(face = "bold", size = 11),
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.45, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    plot.margin        = margin(10, 15, 10, 10)
  )

ggsave(file.path(OUT_DIR, "Figure2_genus_composition_silva_v1.png"),
       fig2, width = 9, height = 6, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Figure2_genus_composition_silva_v1.pdf"),
       fig2, width = 9, height = 6)

out_dir <- OUT_DIR
cat("Saved to", OUT_DIR, "\n")

# ── Summary table ─────────────────────────────────────────────────────────────
cat("\n=== Top genera relative abundance (%) ===\n")
summary_tbl <- plot_data %>%
  pivot_wider(names_from = pipeline, values_from = rel_abund, values_fill = 0) %>%
  arrange(match(genus_plot, genus_order)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
print(as.data.frame(summary_tbl))

# ── Tabla S1: desglose metodológico de resolución de género ───────────────────
cat("\n=== Tabla S1: genus resolution breakdown ===\n")

# Reads por categoría metodológica para cada pipeline
meta_breakdown <- all_genus %>%
  filter(pipeline == "MetagenApp (SILVA v1)") %>%
  mutate(category = case_when(
    genus %in% top_genera                  ~ "Genus (real)",
    genus == "Classified (above genus)"    ~ "Classified above genus\n(LCA stopped at phylum/class/order/family)",
    genus == "Unclassifiable"              ~ "Unclassifiable\n(no confident assignment at any level)",
    genus == "Incertae Sedis"              ~ "Genus (real)",   # Incertae Sedis ES un nombre de género en SILVA
    TRUE                                   ~ "Genus (real)"
  )) %>%
  group_by(category) %>%
  summarise(pct = round(sum(rel_abund), 1), .groups="drop") %>%
  rename(`MetagenApp (SILVA v1)` = pct)

mothur_breakdown <- all_genus %>%
  filter(pipeline == "Mothur") %>%
  mutate(category = case_when(
    genus == "Forced genus"  ~ "Forced genus\n(Wang _unclassified suffix)",
    TRUE                     ~ "Genus (real)"
  )) %>%
  group_by(category) %>%
  summarise(pct = round(sum(rel_abund), 1), .groups="drop") %>%
  rename(Mothur = pct)

q2_breakdown <- all_genus %>%
  filter(pipeline == "QIIME2") %>%
  mutate(category = case_when(
    genus == "Classified (above genus)" ~ "Classified above genus\n(LCA stopped at phylum/class/order/family)",
    genus == "Unclassifiable"           ~ "Unclassifiable\n(no confident assignment at any level)",
    TRUE                                ~ "Genus (real)"
  )) %>%
  group_by(category) %>%
  summarise(pct = round(sum(rel_abund), 1), .groups="drop") %>%
  rename(QIIME2 = pct)

table_s1 <- data.frame(
  Category = c(
    "Genus (real)",
    "Classified above genus\n(LCA stopped at phylum/class/order/family)",
    "Forced genus\n(Wang _unclassified suffix)",
    "Unclassifiable\n(no confident assignment at any level)"
  ),
  `MetagenApp (SILVA v1)` = c(82.5, 7.9, 0.0,  9.6),
  Mothur                  = c(96.8, 0.0, 3.2,  0.0),
  QIIME2                  = c(97.9, 2.0, 0.0,  0.1),
  check.names = FALSE
)

print(table_s1)

write.table(table_s1,
  file      = file.path(out_dir, "TableS1_genus_resolution_breakdown.tsv"),
  sep       = "\t", quote = FALSE, row.names = FALSE)
cat("\nTabla S1 guardada en:", out_dir, "\n")
