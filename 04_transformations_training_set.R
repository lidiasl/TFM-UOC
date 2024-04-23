#####################################################################
# Transformations data training set
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
phenoData_training <- read.csv(file = "phenoData_training.csv")
phenoData_test <- read.csv(file = "phenoData_test.csv")


rownames(rnaseq_counts_primary_select) <- rnaseq_counts_primary_select$X
rnaseq_counts_primary_select <- rnaseq_counts_primary_select[, colnames(rnaseq_counts_primary_select)!="X"]

rownames(log_rnaseq_primary_select) <- log_rnaseq_primary_select$X
log_rnaseq_primary_select <- log_rnaseq_primary_select[, colnames(log_rnaseq_primary_select)!="X"]

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(phenoData_test) <- phenoData_test$X
phenoData_test <- phenoData_test[, colnames(phenoData_test)!="X"]

#Check
table(rownames(rnaseq_counts_primary_select)==rownames(log_rnaseq_primary_select))
# TRUE 
#19063 


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

table(phenoData_training$sample==colnames(rnaseq_counts_training))
#TRUE 
#783 

# DESeq2
dds.train <- DESeqDataSetFromMatrix(countData = rnaseq_counts_training, 
                                    colData = phenoData_training, 
                                    design = ~ 0 + primary_disease_set)
#vsd sin filtrar genes
library("vsn")
vsd <- vst(dds.train, blind=FALSE)

counts.vsd <- assay(vsd)
counts.vsd <- as.data.frame(counts.vsd)

#Check
table(colnames(counts.vsd)==phenoData_training$sample)
#TRUE 
#783 

#Norm sin filtrar genes
dds.train <- estimateSizeFactors(dds.train)
counts.norm <- counts(dds.train, normalized=T)
counts.norm <- as.data.frame(counts.norm)

# Norm log transform
counts.norm.log <- log2(counts.norm+1)

######################################
# Write counts vsd and norm
######################################
setwd(tablesDirectory)
write.csv(counts.vsd, file = "counts_vsd.csv")
write.csv(counts.norm, file = "counts_norm.csv")
write.csv(counts.norm.log, file = "counts_norm_log.csv")
######################################


#####################################
# Findings genes to filter
######################################

#While it is not necessary to pre-filter low count genes before running the DESeq2 functions, 
#there are two reasons which make pre-filtering useful: by removing rows in which there are
#very few reads, we reduce the memory size of the dds data object, and we increase the speed 
#of count modeling within DESeq2. It can also improve visualizations, as features with no information 
#for differential expression are not plotted in dispersion plots or MA-plots.

# Filtering by expression
# Filtering out lowly expressed genes
smallestGroupSize <- min(table(phenoData_training$primary_disease_set)) #29

nrow(rnaseq_counts_training) #19063

keep.train <- rowSums(rnaseq_counts_training >= 10) >= smallestGroupSize

table(keep.train)
# keep 17750 genes
# Remove 1313 genes

# Filter by expression after transformations

rnaseq_counts_training.f <- rnaseq_counts_training[keep.train,]
log_rnaseq_training.f <- log_rnaseq_training[keep.train,]
counts.vsd.f.after <- counts.vsd[keep.train,]
counts.norm.f.after <- counts.norm[keep.train,]

# Norm log transform
counts.norm.f.after.log <- log2(counts.norm.f.after+1)

######################################
# Write training filtering after
######################################
setwd(tablesDirectory)
write.csv(rnaseq_counts_training.f, file = "rnaseq_counts_training_f.csv")
write.csv(log_rnaseq_training.f, file = "log_rnaseq_training_f.csv")
write.csv(counts.vsd.f.after, file = "counts_vsd_f_after.csv")
write.csv(counts.norm.f.after, file = "counts_norm_f_after.csv")
write.csv(counts.norm.f.after.log, file = "counts_norm_f_after_log.csv")
######################################

######################################
# FILTERING before VSD y norm
######################################


dds.train.f <- DESeqDataSetFromMatrix(countData = rnaseq_counts_training.f, 
                                    colData = phenoData_training, 
                                    design = ~ 0 + primary_disease_set)

vsd.f.pre <- vst(dds.train.f, blind=FALSE)
counts.vsd.f.pre <- assay(vsd.f.pre)
counts.vsd.f.pre <- as.data.frame(counts.vsd.f.pre)

pdf(file = "meanSdPlot.pdf", width = 6, height = 4)
meanSdPlot(assay(vsd))
dev.off()

pdf(file = "meanSdPlot_f.pdf", width = 6, height = 4)
meanSdPlot(assay(vsd.f.pre))
dev.off()


dds.train.f <- estimateSizeFactors(dds.train.f)
counts.norm.f.pre <- counts(dds.train.f, normalized=T)
counts.norm.f.pre <- as.data.frame(counts.norm.f.pre)


######################################
# Write training filtering pre
######################################
setwd(tablesDirectory)
write.csv(counts.vsd.f.pre, file = "counts_vsd_f_pre.csv")
write.csv(counts.norm.f.pre, file = "counts_norm_f_pre.csv")
######################################


# Study counts norm f after and f pre
table(rownames(counts.norm.f.pre)==rownames(counts.norm.f.after))

df_foo <- cbind(counts.norm.f.pre$TCGA_OR_A5J6_01, counts.norm.f.after$TCGA_OR_A5J6_01)
df_foo <- as.data.frame(df_foo)
rownames(df_foo) <- rownames(counts.norm.f.pre)
colnames(df_foo) <- c("pre", "after")

ggplot(df_foo, aes(x = pre, y = after)) + geom_point()
table(df_foo$pre==df_foo$after)
#TRUE 
#17750 


######################################
setwd(tablesDirectory)
phenoData_training <- read.csv(file = "phenoData_training.csv")
phenoData_test <- read.csv(file = "phenoData_test.csv")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(phenoData_test) <- phenoData_test$X
phenoData_test <- phenoData_test[, colnames(phenoData_test)!="X"]

# Poner las enfermedades en español

phenoData_training$enfermedad_primaria <- set_subset$spanish[match(phenoData_training$primary_disease_set, set_subset$set)]

phenoData_test$enfermedad_primaria <- set_subset$spanish[match(phenoData_test$primary_disease_set, set_subset$set)]

#Sexo
#ifelse(test, yes, no)
phenoData_training$sexo <- ifelse(phenoData_training$sex=="MALE", "HOMBRE", "MUJER")

phenoData_test$sexo <- ifelse(phenoData_test$sex=="MALE", "HOMBRE", "MUJER")

#Raza
phenoData_training$raza[phenoData_training$race=="WHITE"] <- "BLANCO"
phenoData_training$raza[phenoData_training$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "NEGRO"
phenoData_training$raza[phenoData_training$race=="ASIAN"] <- "ASIÁTICO"

phenoData_test$raza[phenoData_test$race=="WHITE"] <- "BLANCO"
phenoData_test$raza[phenoData_test$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "NEGRO"
phenoData_test$raza[phenoData_test$race=="ASIAN"] <- "ASIÁTICO"

#Edad
phenoData_training$edad <- phenoData_training$age

phenoData_test$edad <- phenoData_test$age

######################################
setwd(tablesDirectory)
write.csv(phenoData_training, file = "phenoData_training.csv")
write.csv(phenoData_test, file = "phenoData_test.csv")
######################################