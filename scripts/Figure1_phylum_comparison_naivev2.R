########################################################
# Figure 1: Phylum comparison across pipelines
# MetagenApp: naive-v2 classifier, SILVA 138 NR99 (v1 model)
########################################################

library(tidyverse)

# Genus resolution labels for annotation (% of OTUs classified to genus level)
genus_resolution <- c(
  "MetagenApp (SILVA v1)" = 91.1,
  "Mothur"               = 100.0,
  "QIIME2"               = 98.0
)

metagen <- read.table(
  "/data/results/runs/ref_runs/faringe_nariz_20260502_v1/summary_tax_phylum.tsv",
  header=TRUE, sep="\t", check.names=FALSE
)

metagen_long <- metagen %>%
  pivot_longer(-Phylum, names_to="Sample", values_to="Count") %>%
  group_by(Phylum) %>%
  summarise(Count=sum(Count), .groups="drop") %>%
  mutate(Pipeline="MetagenApp (SILVA v1)")

mothur_raw <- read.table(
  "/data/results/runs/01_benchmarking/paper01_faringe_mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.rdp.wang.tax.summary",
  header=TRUE, sep="\t", check.names=FALSE
)

mothur <- mothur_raw %>%
  filter(taxlevel==2) %>%
  select(taxon, total) %>%
  rename(Phylum=taxon, Count=total) %>%
  mutate(Pipeline="Mothur")

qiime_table <- read.table(
  "/data/results/runs/01_benchmarking/qiime2_faringe/table.tsv",
  header=TRUE, sep="\t", comment.char="", check.names=FALSE, skip=1
)
taxonomy <- read.table(
  "/data/results/runs/01_benchmarking/qiime2_faringe/taxonomy_export/taxonomy.tsv",
  header=TRUE, sep="\t"
)
colnames(qiime_table)[1] <- "FeatureID"
colnames(taxonomy)[1]    <- "FeatureID"
qiime <- left_join(qiime_table, taxonomy, by="FeatureID")
qiime$Phylum <- gsub("p__","", stringr::str_extract(qiime$Taxon,"p__[^;]+"))
qiime$Phylum[is.na(qiime$Phylum)] <- "Unclassified"
qiime_long <- qiime %>%
  pivot_longer(cols=-c(FeatureID,Taxon,Phylum), names_to="Sample", values_to="Count") %>%
  group_by(Phylum) %>%
  summarise(Count=sum(Count), .groups="drop") %>%
  mutate(Pipeline="QIIME2")

combined <- bind_rows(metagen_long, mothur, qiime_long)

# Mapa SILVA 138 → nomenclatura clásica
# SILVA 138 introdujo nombres modernos basados en la clasificación NCBI/GTDB;
# se mapean a los nombres clásicos para comparabilidad con Mothur y QIIME2.
silva_to_classic <- c(
  "Pseudomonadota"      = "Proteobacteria",
  "Bacillota"           = "Firmicutes",
  "Actinomycetota"      = "Actinobacteria",
  "Bacteroidota"        = "Bacteroidetes",
  "Fusobacteriota"      = "Fusobacteria",
  "Cyanobacteriota"     = "Cyanobacteria",
  "Chloroflexota"       = "Chloroflexi",
  "Planctomycetota"     = "Planctomycetes",
  "Verrucomicrobiota"   = "Verrucomicrobia",
  "Chlamydiota"         = "Chlamydiae",
  "Spirochaetota"       = "Spirochaetes",
  "Nitrospirota"        = "Nitrospirae",
  "Myxococcota"         = "Proteobacteria",   # delta-Proteobacteria reclasificadas
  "Campylobacterota"    = "Proteobacteria",   # epsilon-Proteobacteria reclasificadas
  "Bdellovibrionota"    = "Proteobacteria",   # delta-Proteobacteria reclasificadas
  "Patescibacteria"     = "Other",            # filo nuevo SILVA 138, sin equivalente clásico
  "Deinococcota"        = "Deinococcus-Thermus",
  "Acidobacteriota"     = "Acidobacteria",
  "Fibrobacterota"      = "Fibrobacteres",
  "Synergistota"        = "Synergistetes",
  "Elusimicrobiota"     = "Elusimicrobia",
  "Gemmatimonadota"     = "Gemmatimonadetes",
  "Firmicutes_A"        = "Firmicutes",
  "Proteobacteria_A"    = "Proteobacteria"
)
combined$Phylum <- recode(combined$Phylum, !!!silva_to_classic)
combined$Phylum[combined$Phylum=="" | is.na(combined$Phylum)] <- "Unclassified"

data_all <- combined %>%
  group_by(Pipeline) %>%
  mutate(RelAbund=Count/sum(Count))

top_phyla <- data_all %>%
  filter(Phylum!="Unclassified") %>%
  group_by(Phylum) %>%
  summarise(total=sum(RelAbund)) %>%
  arrange(desc(total)) %>%
  slice(1:5) %>%
  pull(Phylum)

data_all$Phylum <- ifelse(data_all$Phylum %in% top_phyla, data_all$Phylum, "Other")

plot_data <- data_all %>%
  group_by(Pipeline, Phylum) %>%
  summarise(RelAbund=sum(RelAbund), .groups="drop")

phylum_order <- plot_data %>%
  group_by(Phylum) %>%
  summarise(total=sum(RelAbund)) %>%
  arrange(desc(total)) %>%
  pull(Phylum)

plot_data$Phylum <- factor(plot_data$Phylum, levels=phylum_order)
plot_data$Pipeline <- factor(plot_data$Pipeline,
  levels=c("MetagenApp (SILVA v1)","Mothur","QIIME2"))

cat("\n---- DEBUG: abundancias relativas por pipeline ----\n")
print(plot_data %>% arrange(Pipeline, desc(RelAbund)) %>%
  mutate(RelAbund=scales::percent(RelAbund, accuracy=0.1)))

phylum_colors <- c(
  "Firmicutes"       = "#F28E2B",
  "Bacteroidetes"    = "#76B7B2",
  "Proteobacteria"   = "#4E79A7",
  "Actinobacteria"   = "#E15759",
  "Fusobacteria"     = "#59A14F",
  "Other"            = "#B07AA1",
  "Unclassified"     = "#D3D3D3"
)

# Etiquetas de genus resolution debajo del nombre del pipeline
genus_labels <- data.frame(
  Pipeline = names(genus_resolution),
  label    = paste0("Genus res.: ", genus_resolution, "%")
)
genus_labels$Pipeline <- factor(genus_labels$Pipeline,
  levels=c("MetagenApp (SILVA v1)","Mothur","QIIME2"))

p <- ggplot(plot_data, aes(x=Pipeline, y=RelAbund, fill=Phylum)) +
  geom_bar(stat="identity", width=0.7, color="black", linewidth=0.3) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1),
                     expand=expansion(mult=c(0, 0.08))) +
  scale_fill_manual(values=phylum_colors) +
  geom_text(data=genus_labels,
            aes(x=Pipeline, y=-0.04, label=label),
            inherit.aes=FALSE, size=3.2, color="grey30") +
  labs(title="Phylum-level taxonomic composition across pipelines",
       subtitle="36 pharyngeal microbiome samples · SILVA 138 names mapped to classical nomenclature",
       x="Pipeline", y="Relative abundance", fill="Phylum") +
  theme_classic(base_size=14) +
  theme(legend.position   = "right",
        plot.title        = element_text(hjust=0.5, face="bold"),
        plot.subtitle     = element_text(hjust=0.5, color="grey40", size=10),
        axis.text.x       = element_text(face="bold"),
        plot.margin       = margin(10, 10, 30, 10))

out_dir <- "/data/projects/paper_metagenapp_benchmark/figures"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

ggsave(file.path(out_dir,"Figure1_phylum_comparison_naivev2.png"), plot=p, width=8, height=5, dpi=300)
ggsave(file.path(out_dir,"Figure1_phylum_comparison_naivev2.pdf"), plot=p, width=8, height=5)

cat("\nFigura generada en:", out_dir, "\n")
