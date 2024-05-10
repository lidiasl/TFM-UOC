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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Script 03"
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

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

set.seed(42)
phenoData_training_index <-createDataPartition(
  phenoData_select_primary$primary_disease_set,
  times = 1,
  p = 0.8,
  list = TRUE
)

phenoData_training <- phenoData_select_primary[phenoData_training_index[["Resample1"]],]

phenoData_training_split <- split(x = phenoData_training, f = phenoData_training$primary_disease_set)

phenoData_training_split_español <- phenoData_training_split
names(phenoData_training_split_español) <- set_subset$spanish[match(names(phenoData_training_split_español), set_subset$set)]

races_list <- c("BLANCO", "NEGRO_O_AFROAMERICANO", "ASIÁTICO")


setwd(resultsDirectory)
pdf(file = "muestras_training_raza.pdf", width = 8, height = 12)
par(mfrow=c(6,5))
for (i in 1:length(names(phenoData_training_split_español))){
  counts <- table(phenoData_training_split_español[[i]]$raza)
  counts.df <- as.data.frame(cbind(races_list, rep(0, length(races_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$races_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$races_list
  barplot(vector, 
          main = names(phenoData_training_split_español)[i], 
          ylab = "Nº muestras", 
          col = c("#482677FF", "#FDE725FF", "#287D8EFF"),
          names.arg = c("BLANCO"="BLANCO", 
                        "NEGRO_O_AFROAMERICANO"="NEGRO", 
                        "ASIÁTICO"="ASIÁTICO"),
          cex.axis = 1,
          cex.lab = 1,
          cex.main = 1,
          cex.names =0.8,
          las = 2)
}

dev.off()



sex_list <- c("HOMBRE", 
              "MUJER")

setwd(resultsDirectory)
pdf(file = "muestras_training_sexo.pdf", width = 8, height = 12)
par(mfrow=c(5,6))
for (i in 1:length(names(phenoData_training_split_español))){
  counts <- table(phenoData_training_split_español[[i]]$sexo)
  counts.df <- as.data.frame(cbind(sex_list, rep(0, length(sex_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$sex_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$sex_list
  barplot(vector, 
          main = names(phenoData_training_split_español)[i], 
          ylab = "Nº muestras", 
          col = c("deepskyblue3", "coral1"),
          cex.axis = 1,
          cex.lab = 1,
          cex.main = 1,
          cex.names =0.8,
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


rnaseq_counts_training <- rnaseq_counts_primary_select[,colnames(rnaseq_counts_primary_select) %in% phenoData_training$sample]
rnaseq_counts_training <- rnaseq_counts_training[,match(phenoData_training$sample, colnames(rnaseq_counts_training))]

rnaseq_counts_test <- rnaseq_counts_primary_select[,colnames(rnaseq_counts_primary_select) %in% phenoData_test$sample]
rnaseq_counts_test <- rnaseq_counts_test[,match(phenoData_test$sample, colnames(rnaseq_counts_test))]

log_rnaseq_training <- log_rnaseq_primary_select[, colnames(log_rnaseq_primary_select) %in% phenoData_training$sample]
log_rnaseq_training <- log_rnaseq_training[, match(phenoData_training$sample, colnames(log_rnaseq_training))]

log_rnaseq_test <- log_rnaseq_primary_select[, colnames(log_rnaseq_primary_select) %in% phenoData_test$sample]
log_rnaseq_test <- log_rnaseq_test[, match(phenoData_test$sample, colnames(log_rnaseq_test))]

######################################
# Write training and test rnaseq
######################################
setwd(tablesDirectory)
write.csv(rnaseq_counts_training, file = "rnaseq_counts_training.csv")
write.csv(rnaseq_counts_test, file = "rnaseq_counts_test.csv")
write.csv(log_rnaseq_training, file = "log_rnaseq_training.csv")
write.csv(log_rnaseq_test, file = "log_rnaseq_test.csv")
######################################