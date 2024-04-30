#####################################################################
# FPKM genes list
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
log_fpkm_matrix <- read.csv(file = "log_fpkm_matrix.csv")
fpkm_matrix <- read.csv(file = "fpkm_matrix.csv")
# phenoData
phenoData_training <- read.csv(file = "phenoData_training.csv")

genes_list <- readRDS(file="genes_list.rds")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(fpkm_matrix) <- fpkm_matrix$X
fpkm_matrix <- fpkm_matrix[, colnames(fpkm_matrix)!="X"]

rownames(log_fpkm_matrix) <- log_fpkm_matrix$X
log_fpkm_matrix <- log_fpkm_matrix[, colnames(log_fpkm_matrix)!="X"]

######################################
#CHECK:
table(colnames(log_fpkm_matrix)==phenoData_training$sample)
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

############################################################################
################################### FPKM ###################################
############################################################################
fpkm_matrix <- fpkm_matrix[, colnames(fpkm_matrix) %in% phenoData_training$sample]
fpkm_matrix <- as.data.frame(fpkm_matrix)

log_fpkm_matrix <- log_fpkm_matrix[, colnames(log_fpkm_matrix) %in% phenoData_training$sample]
log_fpkm_matrix <- as.data.frame(log_fpkm_matrix)

genes_union_all <- Reduce(unique, genes_union)
PCA_genes_union_all <- prcomp(t(as.matrix(log_fpkm_matrix[rownames(log_fpkm_matrix) %in% genes_union_all,])), 
                     scale. = FALSE, 
                     center = FALSE)
autoplot(PCA_genes_union_all, data = phenoData_training, colour = "primary_disease_set", size=1)

genes_intersect_all <- Reduce(unique, genes_intersect)
PCA_genes_intersect_all <- prcomp(t(as.matrix(log_fpkm_matrix[rownames(log_fpkm_matrix) %in% genes_intersect_all,])), 
                              scale. = FALSE, 
                              center = FALSE)
autoplot(PCA_genes_intersect_all, data = phenoData_training, colour = "primary_disease_set", size=1)

# Union FPKM
PCA_plot_union_fpkm <- list()

for(i in 1:length(genes_union)){
  genes_union_i <- genes_union[[i]]
  fpkm_matrix_i <- fpkm_matrix[rownames(fpkm_matrix) %in% genes_union_i,]
  PCA_fpkm_i <- prcomp(t(as.matrix(fpkm_matrix_i)), 
                       scale. = FALSE, 
                       center = FALSE)
  PCA_fpkm_i_plot <- autoplot(PCA_fpkm_i, data = phenoData_training, colour = "primary_disease_set",
                              main = as.character(names(genes_union)[i]), size=1)
  PCA_plot_union_fpkm[[i]] <- PCA_fpkm_i_plot
  
}

setwd(resultsDirectory)
pdf(file = "PCAs_fpkm_union.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_union_fpkm)){
  print(PCA_plot_union_fpkm[[i]])
}
dev.off()


#Union LOG FPKM 
PCA_plot_union_log_fpkm <- list()

for(i in 1:length(genes_union)){
  genes_union_i <- genes_union[[i]]
  log_fpkm_matrix_i <- log_fpkm_matrix[rownames(log_fpkm_matrix) %in% genes_union_i,]
  PCA_log_fpkm_i <- prcomp(t(as.matrix(log_fpkm_matrix_i)), 
                         scale. = FALSE, 
                         center = FALSE)
  PCA_log_fpkm_i_plot <- autoplot(PCA_log_fpkm_i, data = phenoData_training, colour = "primary_disease_set",
                                main = as.character(names(genes_union)[i]), size=1)
  PCA_plot_union_log_fpkm[[i]] <- PCA_log_fpkm_i_plot
  
}


setwd(resultsDirectory)
pdf(file = "PCAs_log_fpkm_union.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_union_log_fpkm)){
  print(PCA_plot_union_log_fpkm[[i]])
}
dev.off()


# Intersect FPKM
PCA_plot_intersect_fpkm <- list()

for(i in 1:length(genes_intersect)){
  genes_intersect_i <- genes_intersect[[i]]
  fpkm_matrix_i <- fpkm_matrix[rownames(fpkm_matrix) %in% genes_intersect_i,]
  PCA_fpkm_i <- prcomp(t(as.matrix(fpkm_matrix_i)), 
                         scale. = FALSE, 
                         center = FALSE)
  PCA_fpkm_i_plot <- autoplot(PCA_fpkm_i, data = phenoData_training, colour = "primary_disease_set",
                              main = as.character(names(genes_intersect)[i]), size=1)
  PCA_plot_intersect_fpkm[[i]] <- PCA_fpkm_i_plot
  
}


setwd(resultsDirectory)
pdf(file = "PCAs_fpkm_intersect.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_intersect_fpkm)){
  print(PCA_plot_intersect_fpkm[[i]])
}
dev.off()

# Intersect LOG FPKM
PCA_plot_intersect_log_fpkm <- list()

for(i in 1:length(genes_intersect)){
  genes_intersect_i <- genes_intersect[[i]]
  log_fpkm_matrix_i <- log_fpkm_matrix[rownames(log_fpkm_matrix) %in% genes_intersect_i,]
  PCA_log_fpkm_i <- prcomp(t(as.matrix(log_fpkm_matrix_i)), 
                           scale. = FALSE, 
                           center = FALSE)
  PCA_log_fpkm_i_plot <- autoplot(PCA_log_fpkm_i, data = phenoData_training, colour = "primary_disease_set",
                                  main = as.character(names(genes_intersect)[i]), size=1)
  PCA_plot_intersect_log_fpkm[[i]] <- PCA_log_fpkm_i_plot
  
}


setwd(resultsDirectory)
pdf(file = "PCAs_log_fpkm_intersect.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_intersect_log_fpkm)){
  print(PCA_plot_intersect_log_fpkm[[i]])
}
dev.off()

############################################################################
################################### FPKM ###################################
############################################################################