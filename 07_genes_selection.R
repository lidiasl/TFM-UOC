#####################################################################
# DESeq2 select genes list
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

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(rnaseq_counts_training) <- rnaseq_counts_training$X
rnaseq_counts_training <- rnaseq_counts_training[, colnames(rnaseq_counts_training)!="X"]

rownames(rnaseq_counts_training.f) <- rnaseq_counts_training.f$X
rnaseq_counts_training.f <- rnaseq_counts_training.f[, colnames(rnaseq_counts_training.f)!="X"]

######################################
#CHECK:
table(colnames(rnaseq_counts_training.f)==phenoData_training$sample)
table(colnames(rnaseq_counts_training)==phenoData_training$sample)
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

PCA_plot_union <- list()

for(i in 1:length(genes_union)){
  genes_union_i <- genes_union[[i]]
  rnaseq_counts_training_i <- rnaseq_counts_training[rownames(rnaseq_counts_training) %in% genes_union_i,]
  PCA_rnaseq_i <- prcomp(t(as.matrix(rnaseq_counts_training_i)), 
                       scale. = TRUE, 
                       center = TRUE)
  PCA_rnaseq_i_plot <- autoplot(PCA_rnaseq_i, data = phenoData_training, colour = "primary_disease_set",
                              main = as.character(names(genes_union)[i]), size=1)
  PCA_plot_union[[i]] <- PCA_rnaseq_i_plot
  
}

setwd(resultsDirectory)
pdf(file = "PCAs_union.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_union)){
  print(PCA_plot_union[[i]])
}
dev.off()


PCA_plot_intersect <- list()

for(i in 1:length(genes_intersect)){
  genes_intersect_i <- genes_intersect[[i]]
  rnaseq_counts_training_i <- rnaseq_counts_training[rownames(rnaseq_counts_training) %in% genes_intersect_i,]
  PCA_rnaseq_i <- prcomp(t(as.matrix(rnaseq_counts_training_i)), 
                         scale. = TRUE, 
                         center = TRUE)
  PCA_rnaseq_i_plot <- autoplot(PCA_rnaseq_i, data = phenoData_training, colour = "primary_disease_set",
                                main = as.character(names(genes_intersect)[i]), size=1)
  PCA_plot_intersect[[i]] <- PCA_rnaseq_i_plot
  
}

setwd(resultsDirectory)
pdf(file = "PCAs_intersect.pdf", width = 10, height = 8)
for (i in 1:length(PCA_plot_intersect)){
  print(PCA_plot_intersect[[i]])
}
dev.off()



genes_list_0.01 <- lapply(genes_list, function(x){
  lapply(x, function(y){
    rownames(y[y$padj <0.01 & y$log2FoldChange > 0,])
  })
})


genes_union_0.01 <- lapply(genes_list_0.01, function(x){
  Reduce(unique, x)
})

genes_union_in_class_plot <- sapply(genes_union_0.01, function(x){
  list(length(x))
})
genes_union_in_class_plot <- as.data.frame(t(as.data.frame(genes_union_in_class_plot)))
colnames(genes_union_in_class_plot) <- "number"
genes_union_in_class_plot$class <- rownames(genes_union_in_class_plot)
genes_union_in_class_plot$clase <- set_subset$spanish[match(genes_union_in_class_plot$class,set_subset$set)]

ggplot(data=genes_union_in_class_plot, aes(x=clase, y=number)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=number), vjust=-0.3, size=3.5)+
  theme_minimal()+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 12, color="black"), 
        axis.text.x = element_text(size = 10, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 10, color="black"))+
  ggtitle("Número de genes por clase (padj <0.01 & log2FoldChange > 0)")+
  xlab("Tumor primario") + ylab("Número de genes")
