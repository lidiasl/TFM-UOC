#####################################################################
# Samples study Primary Tumor and selection of samples
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
library(MGFR)


tablesDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Tables"
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Read data
######################################
setwd(tablesDirectory)

phenoData_train_primary <- read.csv(file = "phenoData_train_primary.csv")
phenoData_train_primary <- phenoData_train_primary[, colnames(phenoData_train_primary) != "X"]
rownames(phenoData_train_primary) <- phenoData_train_primary$sample

log_rnaseq_primary_train <- read.csv(file = "log_rnaseq_primary_train.csv")
rownames(log_rnaseq_primary_train) <- log_rnaseq_primary_train$X
log_rnaseq_primary_train <- log_rnaseq_primary_train[,colnames(log_rnaseq_primary_train)!="X"]

rnaseq_counts_primary_train <- read.csv(file = "rnaseq_counts_primary_train.csv")
rownames(rnaseq_counts_primary_train) <- rnaseq_counts_primary_train$X
rnaseq_counts_primary_train <- rnaseq_counts_primary_train[,colnames(rnaseq_counts_primary_train)!="X"]


#CHECK:
table(colnames(rnaseq_counts_primary_train)==phenoData_train_primary$sample)
table(colnames(log_rnaseq_primary_train)==phenoData_train_primary$sample)

data(ref.mat)

mat_MGFR <- rnaseq_counts_primary_train
colnames(mat_MGFR) <- phenoData_train_primary$primary_disease_set

print(Sys.time())
res.list <- getMarkerGenes.rnaseq(mat_MGFR, class.vec=colnames(mat_MGFR), 
                                  samples2compare="all", annotate=TRUE, gene.ids.type="ensembl", score.cutoff=1)
print(Sys.time())

names(res.list)


# adrenocortical vs B_cell_lymphoma
httr::GET("http://cran.r-project.org/Rlogo.jpg", config = httr::config(connecttimeout = 120))

print(Sys.time())
res.list1 <- getMarkerGenes.rnaseq(mat_MGFR, class.vec=colnames(mat_MGFR), 
                                  samples2compare=c("adrenocortical", "B_cell_lymphoma"), annotate=TRUE, gene.ids.type="ensembl", score.cutoff=1)
print(Sys.time())

names(res.list1)


df_foo <- mat_MGFR[, colnames(mat_MGFR) %in% c("adrenocortical", "B_cell_lymphoma")]

colnames(df_foo) <- sapply(colnames(df_foo), function(x){
  str_split(x[1], "[.]")[[1]][1]
})

print(Sys.time())
res.list2 <- getMarkerGenes.rnaseq(df_foo, 
                                   class.vec=colnames(df_foo), 
                                   samples2compare="all", annotate=TRUE, gene.ids.type="ensembl", score.cutoff=1)
print(Sys.time())
