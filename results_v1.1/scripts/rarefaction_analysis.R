library(vegan)

# ---- Cargar tablas: filas=OTUs/features, columnas=muestras -> transponer a muestras x OTUs ----

# MetagenApp
metagenapp_raw <- read.delim("data/samples36/metagenapp/final_asv_table.tsv", row.names=1, check.names=FALSE)
metagenapp_mat <- t(as.matrix(metagenapp_raw))
metagenapp_mat[is.na(metagenapp_mat)] <- 0
storage.mode(metagenapp_mat) <- "double"

# QIIME2 (biom-tsv tiene una linea de comentario extra)
qiime_raw <- read.delim("data/samples36/qiime2/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
qiime_mat <- t(as.matrix(qiime_raw))
qiime_mat[is.na(qiime_mat)] <- 0
storage.mode(qiime_mat) <- "double"

# Mothur (.shared: filas=muestras, cols 1-3 son label/Group/numOtus)
mothur_raw <- read.delim("data/samples36/mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared", check.names=FALSE)
rownames(mothur_raw) <- mothur_raw$Group
mothur_mat <- as.matrix(mothur_raw[, -(1:3)])
storage.mode(mothur_mat) <- "double"

cat("Dimensiones (muestras x features):\n")
cat("  MetagenApp:", dim(metagenapp_mat), "\n")
cat("  Mothur:    ", dim(mothur_mat), "\n")
cat("  QIIME2:    ", dim(qiime_mat), "\n")

# ---- Curvas de rarefaccion — una figura con 3 paneles ----
png("rarefaction/rarefaction_curves.png", width=1800, height=650, res=130)
par(mfrow=c(1,3))

rarecurve(metagenapp_mat, step=500, label=FALSE, col="steelblue", main="MetagenApp", xlab="Lecturas", ylab="Features observados")
rarecurve(mothur_mat, step=500, label=FALSE, col="darkorange", main="Mothur", xlab="Lecturas", ylab="OTUs observados")
rarecurve(qiime_mat, step=500, label=FALSE, col="forestgreen", main="QIIME2", xlab="Lecturas", ylab="ASVs observados")

dev.off()
cat("Guardado: rarefaction/rarefaction_curves.png\n")

# ---- Profundidad minima compartida para rarefy ----
min_depth <- min(rowSums(metagenapp_mat), rowSums(mothur_mat), rowSums(qiime_mat))
cat("Profundidad minima compartida (para rarefy):", min_depth, "\n")

saveRDS(list(metagenapp=metagenapp_mat, mothur=mothur_mat, qiime2=qiime_mat, min_depth=min_depth),
        "rarefaction/abundance_matrices.rds")
cat("Guardado: rarefaction/abundance_matrices.rds\n")
cat("\n"); print(sessionInfo())
