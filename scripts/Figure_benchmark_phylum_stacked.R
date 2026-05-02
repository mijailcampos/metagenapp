########################################################
# Figure — Phylum-level 100% stacked bar chart
# Benchmark controlado: MetagenApp vs mothur vs QIIME 2
#
# Dataset: faringe/nariz Samuel (n=36)
# Ejecutar desde el raíz del proyecto MetagenApp:
#   Rscript scripts/Figure_benchmark_phylum_stacked.R
########################################################

library(tidyverse)

BENCH <- "/data/results/runs/01_benchmarking/benchmark_controlled/results_20260418_0205"
OUT_DIR <- file.path(normalizePath("~"), "MetagenApp", "user_data", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. MetagenApp
# ============================================================

meta_df <- read.table(file.path(BENCH, "metagenapp/summary_tax_phylum.tsv"),
                      header = TRUE, sep = "\t", check.names = FALSE)
meta_counts <- meta_df %>%
  pivot_longer(-Phylum, names_to = "Sample", values_to = "Count") %>%
  group_by(Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Pipeline = "MetagenApp")


# ============================================================
# 2. mothur
# ============================================================

mothur_tax_file <- file.path(BENCH, "mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy")
mothur_shared_file <- file.path(BENCH, "mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared")

mothur_tax_raw <- read.table(mothur_tax_file, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
mothur_tax_raw$Phylum <- sapply(mothur_tax_raw$Taxonomy, function(t) {
  t <- gsub("\\(\\d+\\)", "", t)
  parts <- strsplit(t, ";")[[1]]
  if (length(parts) >= 3) trimws(parts[3]) else trimws(parts[length(parts)])
})

class_to_phylum <- c(
  "Gammaproteobacteria"   = "Proteobacteria",
  "Betaproteobacteria"    = "Proteobacteria",
  "Alphaproteobacteria"   = "Proteobacteria",
  "Epsilonproteobacteria" = "Proteobacteria",
  "Deltaproteobacteria"   = "Proteobacteria",
  "Bacilli"               = "Firmicutes",
  "Clostridia"            = "Firmicutes",
  "Erysipelotrichia"      = "Firmicutes",
  "Negativicutes"         = "Firmicutes",
  "Bacteroidia"           = "Bacteroidetes",
  "Flavobacteria"         = "Bacteroidetes",
  "Sphingobacteria"       = "Bacteroidetes",
  "Fusobacteria"          = "Fusobacteria",
  "Fusobacteriia"         = "Fusobacteria",
  "Actinobacteria"        = "Actinobacteria",
  "Spirochaetia"          = "Spirochaetes",
  "Spirochaetes"          = "Spirochaetes"
)
mothur_tax_raw$Phylum <- ifelse(
  mothur_tax_raw$Phylum %in% names(class_to_phylum),
  class_to_phylum[mothur_tax_raw$Phylum],
  mothur_tax_raw$Phylum
)
otu_phylum <- setNames(mothur_tax_raw$Phylum, mothur_tax_raw$OTU)

mothur_shared <- read.table(mothur_shared_file, header = TRUE, sep = "\t",
                              check.names = FALSE)
otu_cols <- colnames(mothur_shared)[-(1:3)]

mothur_counts <- mothur_shared %>%
  select(all_of(otu_cols)) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(everything(), names_to = "OTU", values_to = "Count") %>%
  mutate(Phylum = otu_phylum[OTU]) %>%
  group_by(Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Pipeline = "mothur")


# ============================================================
# 3. QIIME 2
# ============================================================

qiime_table <- read.table(
  "/tmp/q2table.tsv",
  header = TRUE, sep = "\t", comment.char = "", skip = 1,
  check.names = FALSE
)
colnames(qiime_table)[1] <- "FeatureID"

qiime_tax <- read.table(
  file.path(BENCH, "qiime2/taxonomy_export/taxonomy.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  comment.char = ""
)
colnames(qiime_tax)[1:2] <- c("FeatureID", "Taxon")

qiime <- left_join(qiime_table, qiime_tax, by = "FeatureID")
qiime$Phylum <- trimws(gsub("p__", "",
                       stringr::str_extract(qiime$Taxon, "p__[^;]+")))
qiime$Phylum[is.na(qiime$Phylum) | qiime$Phylum == ""] <- "Unclassified"

sample_cols_q <- colnames(qiime_table)[-1]
qiime_counts <- qiime %>%
  pivot_longer(cols = all_of(sample_cols_q),
               names_to = "Sample", values_to = "Count") %>%
  group_by(Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Pipeline = "QIIME 2")


# ============================================================
# 4. Normalizar nombres de filos (GTDB → nomenclatura clásica)
# ============================================================

phylum_recode <- c(
  "Bacillota"        = "Firmicutes",
  "Pseudomonadota"   = "Proteobacteria",
  "Actinomycetota"   = "Actinobacteria",
  "Actinobacteriota" = "Actinobacteria",
  "Bacteroidota"     = "Bacteroidetes",
  "Fusobacteriota"   = "Fusobacteria",
  "Spirochaetota"    = "Spirochaetes",
  "Patescibacteria"  = "Patescibacteria",
  "Campilobacterota" = "Proteobacteria",
  "Mycoplasmatota"   = "Tenericutes",
  "Firmicutes_A"     = "Firmicutes"
)

recode_phyla <- function(df) {
  df %>%
    mutate(Phylum = recode(Phylum, !!!phylum_recode)) %>%
    mutate(Phylum = if_else(is.na(Phylum) | Phylum == "", "Unclassified", Phylum))
}

combined <- bind_rows(
  recode_phyla(meta_counts),
  recode_phyla(mothur_counts),
  recode_phyla(qiime_counts)
) %>%
  group_by(Pipeline, Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop")


# ============================================================
# 5. Top filos + "Other"
# ============================================================

# Calcular abundancia relativa por pipeline primero
combined_rel <- combined %>%
  group_by(Pipeline) %>%
  mutate(RelAbund = Count / sum(Count)) %>%
  ungroup()

# Filos con >= 1% de abundancia media entre pipelines
top_phyla <- combined_rel %>%
  filter(!Phylum %in% c("Unclassified", "Other")) %>%
  group_by(Phylum) %>%
  summarise(mean_ab = mean(RelAbund)) %>%
  filter(mean_ab >= 0.01) %>%
  pull(Phylum)

plot_data <- combined_rel %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla, Phylum, "Other")) %>%
  group_by(Pipeline, Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Pipeline) %>%
  mutate(RelAbund = Count / sum(Count)) %>%
  ungroup()


# ============================================================
# 6. Orden
# ============================================================

phylum_order <- plot_data %>%
  filter(!Phylum %in% c("Other", "Unclassified")) %>%
  group_by(Phylum) %>%
  summarise(mean_ab = mean(RelAbund)) %>%
  arrange(desc(mean_ab)) %>%
  pull(Phylum)

phylum_levels <- c(phylum_order, "Other", "Unclassified")
pipeline_order <- c("MetagenApp", "mothur", "QIIME 2")

plot_data <- plot_data %>%
  mutate(
    Phylum   = factor(Phylum,   levels = phylum_levels),
    Pipeline = factor(Pipeline, levels = pipeline_order)
  )


# ============================================================
# 7. Paleta
# ============================================================

palette <- c(
  "Firmicutes"      = "#F28E2B",
  "Proteobacteria"  = "#4E79A7",
  "Bacteroidetes"   = "#76B7B2",
  "Actinobacteria"  = "#E15759",
  "Fusobacteria"    = "#59A14F",
  "Spirochaetes"    = "#B07AA1",
  "Patescibacteria" = "#FF9DA7",
  "Tenericutes"     = "#9C755F",
  "Other"           = "#BAB0AC"
)


# ============================================================
# 8. Figura
# ============================================================

p <- ggplot(plot_data, aes(x = Pipeline, y = RelAbund, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.65, color = "white", linewidth = 0.3) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(
    x    = NULL,
    y    = "Relative abundance (%)",
    fill = "Phylum"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x        = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y        = element_text(size = 10, color = "black"),
    axis.title.y       = element_text(size = 11),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    legend.title       = element_text(face = "bold", size = 10),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.45, "cm"),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor   = element_blank()
  )


# ============================================================
# 9. Guardar
# ============================================================

out_png <- file.path(OUT_DIR, "Figure_benchmark_phylum_stacked.png")
out_pdf <- file.path(OUT_DIR, "Figure_benchmark_phylum_stacked.pdf")

ggsave(out_png, plot = p, width = 7, height = 5.5, dpi = 300)
ggsave(out_pdf, plot = p, width = 7, height = 5.5)

cat("\n=== Abundancias relativas por pipeline ===\n")
print(
  plot_data %>%
    arrange(Pipeline, desc(RelAbund)) %>%
    mutate(RelAbund = scales::percent(RelAbund, accuracy = 0.1)) %>%
    select(Pipeline, Phylum, RelAbund),
  n = 60
)
cat("\nFigura guardada en:", out_png, "\n")
