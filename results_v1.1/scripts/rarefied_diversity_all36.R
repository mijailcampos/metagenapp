# 7 oct 2026: alfa rarefactada = media POR MUESTRA de 1,000 rarefacciones a min_depth, con semilla propia por pipeline
# (antes: una sola rarefacción y un solo set.seed compartido en secuencia MA→mothur→QIIME 2, de modo que cambiar MA
# alteraba los valores de mothur/QIIME 2). Convergencia 100/500/1,000: rng_check/convergence.R. La muestra sigue siendo
# la unidad experimental: las pruebas pareadas usan los 36 valores resumidos por pipeline (sin pseudorreplicar).
library(vegan)
d <- readRDS("rarefaction/abundance_matrices.rds")
min_depth <- min(sapply(d[c("metagenapp","mothur","qiime2")], function(m) min(rowSums(m))))
cat("n por pipeline:", sapply(d[c("metagenapp","mothur","qiime2")], nrow), " min_depth =", min_depth, "\n")
N_RAREFACTIONS <- 1000
SEEDS <- c(metagenapp = 19760620, mothur = 19760621, qiime2 = 19760622)
one <- function(rare) cbind(richness_obs=specnumber(rare), chao1=estimateR(rare)["S.chao1", ],
                            shannon=diversity(rare, index="shannon"), invsimpson=diversity(rare, index="invsimpson"))
compute_metrics <- function(mat, min_depth, seed) {
  mat <- mat[rowSums(mat) >= min_depth, , drop=FALSE]
  set.seed(seed)
  reps <- replicate(N_RAREFACTIONS, one(rrarefy(mat, sample=min_depth)), simplify="array")
  avg <- apply(reps, c(1, 2), mean)
  data.frame(sample=rownames(mat), avg, row.names=NULL)
}
for (p in names(SEEDS)) {
  a <- compute_metrics(d[[p]], min_depth, SEEDS[[p]])
  write.csv(a, sprintf("rarefaction/alpha_%s_rarefied_all36.csv", p), row.names=FALSE)
  cat(p, "medianas:", sprintf("%s=%.2f", names(a)[-1], sapply(a[-1], median)), "\n")
}
cat("\n===== sessionInfo =====\n"); print(sessionInfo())
