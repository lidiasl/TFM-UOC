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


tablesDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Tables"
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"


######################################
# Read data
######################################
setwd(tablesDirectory)
rnaseq_counts_primary_select <- read.csv(file = "rnaseq_counts_primary_select.csv")
log_rnaseq_primary_select <- read.csv(file = "log_rnaseq_primary_select.csv")
phenoData_select_primary <- read.csv(file = "phenoData_select_primary.csv")

rownames(rnaseq_counts_primary_select) <- rnaseq_counts_primary_select$X
rnaseq_counts_primary_select <- rnaseq_counts_primary_select[, colnames(rnaseq_counts_primary_select)!="X"]

rownames(log_rnaseq_primary_select) <- log_rnaseq_primary_select$X
log_rnaseq_primary_select <- log_rnaseq_primary_select[, colnames(log_rnaseq_primary_select)!="X"]

rownames(phenoData_select_primary) <- phenoData_select_primary$X
phenoData_select_primary <- phenoData_select_primary[, colnames(phenoData_select_primary)!="X"]

set.seed(42)
phenoData_training_index <-createDataPartition(
  phenoData_select_primary$primary_disease_set,
  times = 1,
  p = 0.8,
  list = TRUE
)

phenoData_training <- phenoData_select_primary[phenoData_training_index[["Resample1"]],]

phenoData_training_split <- split(x = phenoData_training, f = phenoData_training$primary_disease_set)


races_list <- c("ASIAN", 
                "BLACK_OR_AFRICAN_AMERICAN", 
                "WHITE")

setwd(resultsDirectory)
pdf(file = "samples_trainig_race.pdf", width = 15, height = 15)
par(mfrow=c(5,6))
for (i in 1:length(names(phenoData_training_split))){
  counts <- table(phenoData_training_split[[i]]$race)
  counts.df <- as.data.frame(cbind(races_list, rep(0, length(races_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$races_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$races_list
  barplot(vector, 
          main = names(phenoData_training_split)[i], 
          ylab = "Nº muestras", 
          col = c("#287D8EFF", "#FDE725FF", "#482677FF"),
          names.arg = c("ASIAN" = "Asia", 
                        "BLACK OR AFRICAN AMERICAN" = "Black", 
                        "WHITE" = "White"),
          cex.axis = 1,
          cex.lab = 1.4,
          cex.main = 1.5,
          cex.names =1.2,
          las = 2)
}

dev.off()


sex_list <- c("MALE", 
              "FEMALE")


setwd(resultsDirectory)
pdf(file = "samples_training_sex.pdf", width = 15, height = 15)
par(mfrow=c(5,6))
for (i in 1:length(names(phenoData_training_split))){
  counts <- table(phenoData_training_split[[i]]$sex)
  counts.df <- as.data.frame(cbind(sex_list, rep(0, length(sex_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$sex_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$sex_list
  barplot(vector, 
          main = names(phenoData_training_split)[i], 
          ylab = "Nº muestras", 
          col = c("deepskyblue3", "coral1"),
          cex.axis = 1,
          cex.lab = 1.4,
          cex.main = 1.5,
          cex.names =1.2,
          las = 2)
}

dev.off()


###### TEST
phenoData_test <- phenoData_select_primary[!rownames(phenoData_select_primary)%in% rownames(phenoData_training),]

#Check
table(rownames(phenoData_training)%in% rownames(phenoData_test))

######################################
# Write select primary
######################################
setwd(tablesDirectory)
write.csv(phenoData_training, file = "phenoData_training.csv")
write.csv(phenoData_test, file = "phenoData_test.csv")
######################################

