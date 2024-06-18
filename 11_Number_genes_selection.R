#####################################################################
# Number genes selection
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/DESeq2"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Read training data
######################################
setwd(tablesDirectory)
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

names_primary_disease_set <- unique(phenoData_training_final$primary_disease_set)

######################################
tables_list <- list()
# Read tables
setwd(resultsDirectory)
for (i in 1:length(names_primary_disease_set)){
  tables_list[[i]] <- read.csv(file=paste0(names_primary_disease_set[i], "_padj_log2FC.csv"))
}

names(tables_list) <- names_primary_disease_set

df_union <- as.data.frame(cbind(tables_list[[1]]$padj, tables_list[[1]]$log2FC))
for(i in 1:length(names_primary_disease_set)){
  df_union <- cbind(df_union, tables_list[[i]]$Union)
}

colnames(df_union) <- c("padj","log2FC",names_primary_disease_set)

df_union[,3:ncol(df_union)] <- round(df_union[,3:ncol(df_union)])

setwd(resultsDirectory)
pdf(file = "Genes_selection.pdf", width = 8, height = 4)
for (i in 1:length(names_primary_disease_set)){
  plot(tables_list[[i]]$Union, type="b", main=names_primary_disease_set[i], ylab="Nº genes", xlab="Indice combinación")
}
dev.off()


setwd(resultsDirectory)
pdf(file = "Genes_selection_Intersect.pdf", width = 8, height = 4)
for (i in 1:length(names_primary_disease_set)){
  plot(tables_list[[i]]$Intersect, type="b", main=names_primary_disease_set[i], ylab="Nº genes", xlab="Indice combinación")
}
dev.off()


setwd(resultsDirectory)
pdf(file = "Genes_selection_Intersect_24.pdf", width = 8, height = 4)
for (i in 1:length(names_primary_disease_set)){
  plot(tables_list[[i]]$Intersect_24, type="b", main=names_primary_disease_set[i], ylab="Nº genes", xlab="Indice combinación")
}
dev.off()