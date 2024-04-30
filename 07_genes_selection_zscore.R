#####################################################################
# zscore from log2(fpkm+1) genes list
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
# Filtered
rnaseq_counts_training.f <- read.csv(file = "rnaseq_counts_training_f.csv")
# phenoData
phenoData_training <- read.csv(file = "phenoData_training.csv")

genes_list <- readRDS(file="genes_list.rds")

fpkm_matrix <- read.csv(file="fpkm_matrix.csv")
log_fpkm_matrix <- read.csv(file="log_fpkm_matrix.csv")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(rnaseq_counts_training) <- rnaseq_counts_training$X
rnaseq_counts_training <- rnaseq_counts_training[, colnames(rnaseq_counts_training)!="X"]

rownames(rnaseq_counts_training.f) <- rnaseq_counts_training.f$X
rnaseq_counts_training.f <- rnaseq_counts_training.f[, colnames(rnaseq_counts_training.f)!="X"]

rownames(fpkm_matrix) <- fpkm_matrix$X
fpkm_matrix <- fpkm_matrix[, colnames(fpkm_matrix)!="X"]

rownames(log_fpkm_matrix) <- log_fpkm_matrix$X
log_fpkm_matrix <- log_fpkm_matrix[, colnames(log_fpkm_matrix)!="X"]

######################################

zscore_log_fpkm_matrix <- scale(log_fpkm_matrix, center = TRUE, scale = TRUE)


zscore_log_fpkm_matrix <- zscore_log_fpkm_matrix[, colnames(zscore_log_fpkm_matrix) %in% phenoData_training$sample]
zscore_log_fpkm_matrix <- as.data.frame(zscore_log_fpkm_matrix)

######################################
for(i in 1:length(genes_list)){
  genes_list[[i]] <- genes_list[[i]][-i]
}

genes_list_0.05 <- lapply(genes_list, function(x){
  lapply(x, function(y){
    rownames(y[y$padj <0.05 & y$log2FoldChange > 0,])
  })
})

genes_union <- lapply(genes_list_0.05, function(x){
  Reduce(unique, x)
})

genes_intersect <- lapply(genes_list_0.05, function(x){
  Reduce(intersect, x)
})

######################################
# Union ZSCORE LOG FPKM
PCA_plot_union_zscore_log_fpkm <- list()

for(i in 1:length(genes_union)){
  genes_union_i <- genes_union[[i]]
  zscore_log_fpkm_matrix_i <- zscore_log_fpkm_matrix[rownames(zscore_log_fpkm_matrix) %in% genes_union_i,]
  PCA_zscore_log_fpkm_i <- prcomp(t(as.matrix(zscore_log_fpkm_matrix_i)), 
                       scale. = FALSE, 
                       center = FALSE)
  PCA_zscore_log_fpkm_i_plot <- autoplot(PCA_zscore_log_fpkm_i, data = phenoData_training, colour = "primary_disease_set",
                              main = as.character(names(genes_union)[i]), size=1)
  PCA_plot_union_zscore_log_fpkm[[i]] <- PCA_zscore_log_fpkm_i_plot
  
}

setwd(resultsDirectory)
pdf(file = "PCAs_zscore_log_fpkm_union.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_union_zscore_log_fpkm)){
  print(PCA_plot_union_zscore_log_fpkm[[i]])
}
dev.off()


######################################
# INtersect ZSCORE LOG FPKM

PCA_plot_intersect_zscore_log_fpkm <- list()

for(i in 1:length(genes_intersect)){
  genes_intersect_i <- genes_intersect[[i]]
  zscore_log_fpkm_matrix_i <- zscore_log_fpkm_matrix[rownames(zscore_log_fpkm_matrix) %in% genes_intersect_i,]
  PCA_zscore_log_fpkm_i <- prcomp(t(as.matrix(zscore_log_fpkm_matrix_i)), 
                       scale. = FALSE, 
                       center = FALSE)
  PCA_zscore_log_fpkm_i_plot <- autoplot(PCA_zscore_log_fpkm_i, data = phenoData_training, colour = "primary_disease_set",
                              main = as.character(names(genes_intersect)[i]), size=1)
  PCA_plot_intersect_zscore_log_fpkm[[i]] <- PCA_zscore_log_fpkm_i_plot
  
}


setwd(resultsDirectory)
pdf(file = "PCAs_zscore_log_fpkm_intersect.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_intersect_zscore_log_fpkm)){
  print(PCA_plot_intersect_zscore_log_fpkm[[i]])
}
dev.off()