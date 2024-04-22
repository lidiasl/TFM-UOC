#####################################################################
# GRAPHS UMAP
# Lidia Sainz
#####################################################################

library(data.table)
library(tidyr)
library(caret)
library(ggplot2)
library(forcats)
library(DESeq2)
library(edgeR)
library(stats)
library(ggfortify)
library(stringr)
library(DESeq2)
library(dplyr)
library(umap)
library(car)
library(gridExtra)
library(grid)


tablesDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Tables"
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Read training data
######################################
setwd(tablesDirectory)
rnaseq_counts_training <- read.csv(file = "rnaseq_counts_training.csv")
log_rnaseq_training <- read.csv(file = "log_rnaseq_training.csv")
counts.vsd <- read.csv(file = "counts_vsd.csv")
counts.norm <- read.csv(file = "counts_norm.csv")
counts.norm.log <- read.csv(file = "counts_norm_log.csv")
# Filtered
rnaseq_counts_training.f <- read.csv(file = "rnaseq_counts_training_f.csv")
log_rnaseq_training.f <- read.csv(file = "log_rnaseq_training_f.csv")
counts.vsd.f <- read.csv(file = "counts_vsd_f_after.csv")
counts.norm.f <- read.csv(file = "counts_norm_f_after.csv")
counts.norm.f.log <- read.csv(file = "counts_norm_f_after_log.csv")
# phenoData
phenoData_training <- read.csv(file = "phenoData_training.csv")
######################################


rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(rnaseq_counts_training) <- rnaseq_counts_training$X
rnaseq_counts_training <- rnaseq_counts_training[, colnames(rnaseq_counts_training)!="X"]

rownames(rnaseq_counts_training.f) <- rnaseq_counts_training.f$X
rnaseq_counts_training.f <- rnaseq_counts_training.f[, colnames(rnaseq_counts_training.f)!="X"]

rownames(log_rnaseq_training) <- log_rnaseq_training$X
log_rnaseq_training <- log_rnaseq_training[, colnames(log_rnaseq_training)!="X"]

rownames(log_rnaseq_training.f) <- log_rnaseq_training.f$X
log_rnaseq_training.f <- log_rnaseq_training.f[, colnames(log_rnaseq_training.f)!="X"]

rownames(counts.vsd) <- counts.vsd$X
counts.vsd <- counts.vsd[, colnames(counts.vsd)!="X"]

rownames(counts.vsd.f) <- counts.vsd.f$X
counts.vsd.f <- counts.vsd.f[, colnames(counts.vsd.f)!="X"]

rownames(counts.norm) <- counts.norm$X
counts.norm <- counts.norm[, colnames(counts.norm)!="X"]

rownames(counts.norm.f) <- counts.norm.f$X
counts.norm.f <- counts.norm.f[, colnames(counts.norm.f)!="X"]

rownames(counts.norm.log) <- counts.norm.log$X
counts.norm.log <- counts.norm.log[, colnames(counts.norm.log)!="X"]

rownames(counts.norm.f.log) <- counts.norm.f.log$X
counts.norm.f.log <- counts.norm.f.log[, colnames(counts.norm.f.log)!="X"]

######################################

phenoData_training$race[phenoData_training$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"

# Cargar la paleta de colores turbo
library(viridisLite)
colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)
names(colors_blind) <- unique(phenoData_training$primary_disease_set)
colors_blind <- list(primary_disease_set = colors_blind)

library(pheatmap)

table(colnames(rnaseq_counts_training)==phenoData_training$sample)

rnaseq_counts_training_top <- rnaseq_counts_training[order(rowVars(as.matrix(rnaseq_counts_training)), decreasing = TRUE)[1:5000],]

setwd(resultsDirectory)
pdf(file = "Pheatmap2.pdf", width = 80, height = 80)

pheatmap(as.matrix(rnaseq_counts_training_top), 
         annotation_col = data.frame(anno = phenoData_training$primary_disease_set,
                                     row.names = phenoData_training$sample),
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_colors = colors_blind)
dev.off()

colnames(as.matrix(rnaseq_counts_training))
annotation_col = as.factor(phenoData_training$primary_disease_set)



### CLUSTER
library(dendextend)
library(circlize)

# Distance matrix
d <- dist(as.matrix(rnaseq_counts_training_top))
hc <- as.dendrogram(hclust(d))
dend <- as.dendrogram(hc)
plot(dend)

circlize_dendrogram(hc,
                    labels_track_height = NA,
                    dend_track_height = 0.5,
                    facing = "inside")
