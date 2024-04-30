#####################################################################
# Transformations ML
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
# Read data
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

phenoData_training$enfermedad_primaria <- set_subset$spanish[match(phenoData_training$primary_disease_set, set_subset$set)]

#Sexo
#ifelse(test, yes, no)
phenoData_training$sexo <- ifelse(phenoData_training$sex=="MALE", "HOMBRE", "MUJER")


#Raza
phenoData_training$raza[phenoData_training$race=="WHITE"] <- "BLANCO"
phenoData_training$raza[phenoData_training$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "NEGRO"
phenoData_training$raza[phenoData_training$race=="ASIAN"] <- "ASIÁTICO"


#Edad
phenoData_training$edad <- phenoData_training$age


######################################

setwd(tablesDirectory)
gtf <- read.table('gencode.v23.chr_patch_hapl_scaff.basic.annotation.gtf', header = FALSE, sep = '\t')

gtf$gene_id_mod <- sapply(gtf$V9, function(x){
  sub(".*gene_id ([^.]+).*", "\\1", x)
})

gtf.f <- gtf[gtf$V3 == "gene",]

gtf.f <- gtf.f[gtf.f$gene_id_mod %in% rownames(rnaseq_counts_training),]

featureLength <- (gtf.f$V5 - gtf.f$V4)+1

gtf.f$featureLength <- featureLength
gtf.f <- gtf.f[match(rownames(rnaseq_counts_training), gtf$gene_id_mod),]

#install.packages("DGEobj.utils")
library(DGEobj.utils)


fpkm_matrix <- convertCounts(as.matrix(rnaseq_counts_training),
                             unit       = "fpkm",
                             geneLength = featureLength,
                             log        = FALSE,
                             normalize  = "none")

tpm_matrix <- convertCounts(as.matrix(rnaseq_counts_training),
                             unit       = "tpm",
                             geneLength = featureLength,
                             log        = FALSE,
                             normalize  = "none")

log_fpkm_matrix <- log2(fpkm_matrix+0.001)
log_tpm_matrix <- log2(tpm_matrix+0.001)

zscore_log_fpkm <- scale(log_fpkm_matrix)
zscore_log_tpm <- scale(log_tpm_matrix)

######################################
setwd(tablesDirectory)
write.csv(fpkm_matrix, file="fpkm_matrix.csv")
write.csv(log_fpkm_matrix, file="log_fpkm_matrix.csv")
write.csv(zscore_log_fpkm, file = "zscore_log_fpkm.csv")

write.csv(tpm_matrix, file="tpm_matrix.csv")
write.csv(log_tpm_matrix, file="log_tpm_matrix.csv")
write.csv(zscore_log_tpm, file = "zscore_log_tpm.csv")
######################################
