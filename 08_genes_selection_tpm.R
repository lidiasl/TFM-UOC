#####################################################################
# TPM genes list
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
log_tpm_matrix <- read.csv(file = "log_tpm_matrix.csv")
tpm_matrix <- read.csv(file = "tpm_matrix.csv")
# phenoData
phenoData_training <- read.csv(file = "phenoData_training.csv")

genes_list <- readRDS(file="genes_list.rds")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(tpm_matrix) <- tpm_matrix$X
tpm_matrix <- tpm_matrix[, colnames(tpm_matrix)!="X"]

rownames(log_tpm_matrix) <- log_tpm_matrix$X
log_tpm_matrix <- log_tpm_matrix[, colnames(log_tpm_matrix)!="X"]

######################################
#CHECK:
table(colnames(tpm_matrix)==phenoData_training$sample)
phenoData_training$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
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

tpm_matrix <- tpm_matrix[, colnames(tpm_matrix) %in% phenoData_training$sample]
tpm_matrix <- as.data.frame(tpm_matrix)

log_tpm_matrix <- log_tpm_matrix[, colnames(log_tpm_matrix) %in% phenoData_training$sample]
log_tpm_matrix <- as.data.frame(log_tpm_matrix)


# Union TPM
PCA_plot_union_tpm <- list()

for(i in 1:length(genes_union)){
  genes_union_i <- genes_union[[i]]
  tpm_matrix_i <- tpm_matrix[rownames(tpm_matrix) %in% genes_union_i,]
  PCA_tpm_i <- prcomp(t(as.matrix(tpm_matrix_i)), 
                       scale. = FALSE, 
                       center = FALSE)
  PCA_tpm_i_plot <- autoplot(PCA_tpm_i, data = phenoData_training, colour = "primary_disease_set",
                              main = as.character(names(genes_union)[i]), size=1)
  PCA_plot_union_tpm[[i]] <- PCA_tpm_i_plot
  
}

setwd(resultsDirectory)
pdf(file = "PCAs_tpm_union.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_union_tpm)){
  print(PCA_plot_union_tpm[[i]])
}
dev.off()
