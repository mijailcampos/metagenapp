########################################################
# Figure 1 — Phylum-level 100% stacked bar chart
# Comparación: MetagenApp (general) vs MetagenApp (oral)
#              vs Mothur vs QIIME2
#
# Dataset: faringe/nariz Samuel
# Ejecutar desde el raíz del proyecto MetagenApp:
#   Rscript scripts/Figure1_phylum_stacked.R
########################################################

library(tidyverse)

OUT_DIR <- file.path(normalizePath("~"), "MetagenApp", "user_data", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. MetagenApp — modelo general (v4)
# ============================================================

load_metagenapp <- function(path, pipeline_name) {
  df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
  df %>%
    pivot_longer(-Phylum, names_to = "Sample", values_to = "Count") %>%
    group_by(Phylum) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    mutate(Pipeline = pipeline_name)
}

metagen_general <- load_metagenapp(
  "/data/results/paper01_faringe_naivev2_v4/summary_tax_phylum.tsv",
  "MetagenApp (general)"
)

metagen_oral <- load_metagenapp(
  "/data/results/paper01_faringe_oral_v1/summary_tax_phylum.tsv",
  "MetagenApp (oral)"
)


# ============================================================
# 2. Mothur — shared table + cons.taxonomy
# ============================================================

mothur_tax_file <- "/data/results/paper01_faringe_mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"
mothur_shared_file <- "/data/results/paper01_faringe_mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared"

mothur_tax_raw <- read.table(mothur_tax_file, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
# Extraer filo (nivel 3 en cadena Root;Kingdom;Phylum;...)
mothur_tax_raw$Phylum <- sapply(mothur_tax_raw$Taxonomy, function(t) {
  t <- gsub("\\(\\d+\\)", "", t)          # eliminar confidence scores
  parts <- strsplit(t, ";")[[1]]
  if (length(parts) >= 3) parts[3] else parts[length(parts)]
})
otu_phylum <- setNames(mothur_tax_raw$Phylum, mothur_tax_raw$OTU)

# Mothur usa clases de Proteobacteria y de Firmicutes — unificar a filo
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
  "Actinobacteria"        = "Actinobacteria",
  "Spirochaetia"          = "Spirochaetes",
  "Spirochaetes"          = "Spirochaetes"
)
mothur_tax_raw$Phylum <- ifelse(
  mothur_tax_raw$Phylum %in% names(class_to_phylum),
  class_to_phylum[mothur_tax_raw$Phylum],
  mothur_tax_raw$Phylum
)
otu_phylum_merged <- setNames(mothur_tax_raw$Phylum, mothur_tax_raw$OTU)

mothur_shared <- read.table(mothur_shared_file, header = TRUE, sep = "\t",
                              check.names = FALSE)
otu_cols <- colnames(mothur_shared)[-(1:3)]

mothur_counts <- mothur_shared %>%
  select(all_of(otu_cols)) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(everything(), names_to = "OTU", values_to = "Count") %>%
  mutate(Phylum = otu_phylum_merged[OTU]) %>%
  group_by(Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Pipeline = "Mothur")


# ============================================================
# 3. QIIME2 — feature table + taxonomy
# ============================================================

qiime_table <- read.table(
  "/data/projects/qiime2_faringe/table_export/feature-table.tsv",
  header = TRUE, sep = "\t", comment.char = "", skip = 1,
  check.names = FALSE
)
colnames(qiime_table)[1] <- "FeatureID"

qiime_tax <- read.table(
  "/data/projects/qiime2_faringe/taxonomy_export/taxonomy.tsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
colnames(qiime_tax)[1] <- "FeatureID"

qiime <- left_join(qiime_table, qiime_tax, by = "FeatureID")
qiime$Phylum <- gsub("p__", "",
                      stringr::str_extract(qiime$Taxon, "p__[^;]+"))
qiime$Phylum[is.na(qiime$Phylum) | qiime$Phylum == ""] <- "Unclassified"

sample_cols_q <- colnames(qiime_table)[-1]

qiime_counts <- qiime %>%
  pivot_longer(cols = all_of(sample_cols_q),
               names_to = "Sample", values_to = "Count") %>%
  group_by(Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Pipeline = "QIIME2")


# ============================================================
# 4. Combinar y normalizar nombres de filos
# ============================================================

# Tabla de estandarización de nombres (GTDB → clásico)
phylum_recode <- c(
  "Bacillota"        = "Firmicutes",
  "Pseudomonadota"   = "Proteobacteria",
  "Actinomycetota"   = "Actinobacteria",
  "Bacteroidota"     = "Bacteroidetes",
  "Fusobacteriota"   = "Fusobacteria",
  "Cyanobacteriota"  = "Cyanobacteria",
  "Spirochaetota"    = "Spirochaetes",
  "Patescibacteria"  = "Patescibacteria",
  "Mycoplasmatota"   = "Tenericutes",
  "Chloroflexota"    = "Chloroflexi",
  "Verrucomicrobiota"= "Verrucomicrobia",
  "Synergistota"     = "Synergistetes",
  "Campilobacterota" = "Proteobacteria",  # Epsilonproteobacteria en QIIME2
  "Actinobacteriota" = "Actinobacteria",
  "Firmicutes_A"     = "Firmicutes",
  "Ascomycota"       = "Unclassified"     # hongo = contaminante
)

recode_phyla <- function(df) {
  df %>%
    mutate(Phylum = recode(Phylum, !!!phylum_recode)) %>%
    mutate(Phylum = if_else(Phylum == "" | is.na(Phylum), "Unclassified", Phylum))
}

combined <- bind_rows(
  recode_phyla(metagen_general),
  recode_phyla(metagen_oral),
  recode_phyla(mothur_counts),
  recode_phyla(qiime_counts)
) %>%
  group_by(Pipeline, Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop")


# ============================================================
# 5. Seleccionar top filos y agrupar resto como "Other"
# ============================================================

TOP_N <- 7

top_phyla <- combined %>%
  filter(Phylum != "Unclassified") %>%
  group_by(Phylum) %>%
  summarise(total = sum(Count)) %>%
  arrange(desc(total)) %>%
  slice_head(n = TOP_N) %>%
  pull(Phylum)

plot_data <- combined %>%
  mutate(Phylum = if_else(Phylum %in% c(top_phyla, "Unclassified"),
                           Phylum, "Other")) %>%
  group_by(Pipeline, Phylum) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Pipeline) %>%
  mutate(RelAbund = Count / sum(Count)) %>%
  ungroup()


# ============================================================
# 6. Orden de filos (mayor → menor abundancia media)
# ============================================================

phylum_order <- plot_data %>%
  filter(!Phylum %in% c("Other", "Unclassified")) %>%
  group_by(Phylum) %>%
  summarise(mean_ab = mean(RelAbund)) %>%
  arrange(desc(mean_ab)) %>%
  pull(Phylum)

phylum_levels <- c(phylum_order, "Other", "Unclassified")

pipeline_order <- c("MetagenApp (general)", "MetagenApp (oral)",
                    "Mothur", "QIIME2")

plot_data <- plot_data %>%
  mutate(
    Phylum   = factor(Phylum,   levels = phylum_levels),
    Pipeline = factor(Pipeline, levels = pipeline_order)
  )


# ============================================================
# 7. Paleta de colores
# ============================================================

palette <- c(
  "Proteobacteria"  = "#4E79A7",
  "Firmicutes"      = "#F28E2B",
  "Bacteroidetes"   = "#76B7B2",
  "Actinobacteria"  = "#E15759",
  "Fusobacteria"    = "#59A14F",
  "Cyanobacteria"   = "#EDC948",
  "Spirochaetes"    = "#B07AA1",
  "Patescibacteria" = "#FF9DA7",
  "Tenericutes"     = "#9C755F",
  "Other"           = "#BAB0AC",
  "Unclassified"    = "#D3D3D3"
)


# ============================================================
# 8. Figura
# ============================================================

p <- ggplot(plot_data, aes(x = Pipeline, y = RelAbund, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.72, color = "white", linewidth = 0.25) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(
    title    = "Phylum-level relative abundance across pipelines",
    subtitle = "Nasopharyngeal microbiome - Faringe dataset (n=36 samples)",
    x        = NULL,
    y        = "Relative abundance (%)",
    fill     = "Phylum"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 10),
    axis.text.x     = element_text(face = "bold", size = 11),
    axis.text.y     = element_text(size = 10),
    legend.title    = element_text(face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )


# ============================================================
# 9. Guardar
# ============================================================

ggsave(file.path(OUT_DIR, "Figure1_phylum_stacked.png"),
       plot = p, width = 8, height = 5.5, dpi = 300)
ggsave(file.path(OUT_DIR, "Figure1_phylum_stacked.pdf"),
       plot = p, width = 8, height = 5.5)

cat("\n=== Abundancias relativas por pipeline ===\n")
print(
  plot_data %>%
    arrange(Pipeline, desc(RelAbund)) %>%
    mutate(RelAbund = scales::percent(RelAbund, accuracy = 0.1)) %>%
    select(Pipeline, Phylum, RelAbund),
  n = 50
)

cat("\nFigura guardada en:", OUT_DIR, "\n")
