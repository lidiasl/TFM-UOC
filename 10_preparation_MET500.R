#####################################################################
# Transform data MET500 ML
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/MET500 Preparation"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

#####################################################################
setwd(tablesDirectory)
zscore_log_fpkm.genes_final <- read.csv(file = "zscore_log_fpkm.genes_final.csv")
fpkm_matrix_training <- read.csv(file="fpkm_matrix.csv")

rownames(zscore_log_fpkm.genes_final) <- zscore_log_fpkm.genes_final$X
zscore_log_fpkm.genes_final <- zscore_log_fpkm.genes_final[, colnames(zscore_log_fpkm.genes_final)!="X"]

rownames(fpkm_matrix_training) <- fpkm_matrix_training$X
fpkm_matrix_training <- fpkm_matrix_training[, colnames(fpkm_matrix_training)!="X"]


setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

# Load FPKM for MET500
FPKM_MET500 <- fread("~/Desktop/LIDIA/TCGA_classification/Data/M.mx.txt", sep = "\t", stringsAsFactors = F)
genenames_MET500 <- separate(data = FPKM_MET500, col = sample, into = c("GeneName", "version"), remove = T)
FPKM_MET500 <- FPKM_MET500[,-1]
FPKM_MET500 <- as.data.frame(FPKM_MET500)
rownames(FPKM_MET500) <- genenames_MET500$GeneName

#####################################################################

table(rownames(FPKM_MET500) %in% rownames(fpkm_matrix_training))
#FALSE  TRUE 
#1916 19063 

FPKM_MET500.f <- FPKM_MET500[rownames(FPKM_MET500) %in% rownames(fpkm_matrix_training), ]

#####################################################################
# Samples
setwd(dataDirectory)
phenodata_MET500_1 <- fread(file="MET500_geneExpression_M.meta.plus.txt")
phenodata_MET500_2 <- fread(file="phenodata_complete_MET500.csv")


table(phenodata_MET500_1$Sample_id == colnames(FPKM_MET500.f))
# TRUE: 868

# Muestras:
muestras_polyA <- phenodata_MET500_1$Sample_id[grep("poly", phenodata_MET500_1$run.id)]
muestras_capt <- phenodata_MET500_1$Sample_id[grep("capt", phenodata_MET500_1$run.id)]

muestras_source_polyA <- phenodata_MET500_1$sample_source[grep("poly", phenodata_MET500_1$run.id)]
muestras_source_capt <- phenodata_MET500_1$sample_source[grep("capt", phenodata_MET500_1$run.id)]

table(muestras_source_polyA %in% muestras_source_capt)
#FALSE  TRUE 
#12   372 

muestras_capt_outofpolyA <- muestras_source_capt[!muestras_source_capt %in% muestras_source_polyA]
muestras_capt_outofpolyA <- phenodata_MET500_1$Sample_id[phenodata_MET500_1$sample_source %in% muestras_capt_outofpolyA]

muestras_totales <- c(muestras_polyA, muestras_capt_outofpolyA)


# Filtering
FPKM_MET500.f <- FPKM_MET500.f[, colnames(FPKM_MET500.f) %in% muestras_totales]
phenodata_MET500.f <- phenodata_MET500_1[match(colnames(FPKM_MET500.f), phenodata_MET500_1$run.id),]

## Prepare tissue
table(is.na(phenodata_MET500.f$tissue)) #FALSE
table(is.na(phenodata_MET500.f$biopsy_tissue)) #FALSE

table(phenodata_MET500.f$tissue=="")
#FALSE  TRUE 
#463    33 
table(phenodata_MET500.f$biopsy_tissue=="")
#FALSE 
#496

table(phenodata_MET500.f$sample_type) #tumor

phenodata_MET500.f <- phenodata_MET500.f[phenodata_MET500.f$tissue!="",]
phenodata_MET500.f <- as.data.frame(phenodata_MET500.f)

table(phenodata_MET500.f$tissue==phenodata_MET500.f$biopsy_tissue)

#####################################################################

setwd(tablesDirectory)
write.csv(phenodata_MET500.f, file="phenodata_MET500.f.csv")
write.csv(FPKM_MET500.f, file="FPKM_MET500.f.csv")


#setwd(dataDirectory)
#write.csv(df_tissue_biopsy_tissue, file="df_tissue_biopsy_tissue.csv")
#####################################################################

phenodata_MET500_clean <- phenodata_MET500.f

phenodata_MET500_clean <- phenodata_MET500_clean[phenodata_MET500_clean$tissue!=phenodata_MET500_clean$biopsy_tissue,] #434 samples

table(phenodata_MET500_clean$tissue %in% c("gall_bladder", "parotid", "soft_tissue", "lymph_node", "bone_marrow", "other"))
#FALSE  TRUE 
#387    47 

phenodata_MET500_clean <- phenodata_MET500_clean[!phenodata_MET500_clean$tissue %in% c("gall_bladder", "parotid", "soft_tissue", "lymph_node", "bone_marrow", "other"),] #387

table(phenodata_MET500_clean$tissue, phenodata_MET500_clean$biopsy_tissue)

FPKM_MET500_clean <- FPKM_MET500.f[,colnames(FPKM_MET500.f) %in% phenodata_MET500_clean$Sample_id]
phenodata_MET500_clean <- phenodata_MET500_clean[match(colnames(FPKM_MET500_clean), phenodata_MET500_clean$run.id),]

#####################################################################
setwd(tablesDirectory)
write.csv(phenodata_MET500_clean, file="phenodata_MET500_clean.csv")
write.csv(FPKM_MET500_clean, file="FPKM_MET500_clean.csv")
#####################################################################

log_FPKM_MET500 <- log2(FPKM_MET500.f+0.001)
zscore_log_FPKM_MET500 <- scale(log_FPKM_MET500)

log_FPKM_MET500_clean <- log_FPKM_MET500[, colnames(log_FPKM_MET500) %in% phenodata_MET500_clean$Sample_id]
zscore_log_FPKM_MET500_clean <- zscore_log_FPKM_MET500[, colnames(zscore_log_FPKM_MET500) %in% phenodata_MET500_clean$Sample_id]

table(colnames(log_FPKM_MET500_clean)==phenodata_MET500_clean$Sample_id)
table(colnames(zscore_log_FPKM_MET500_clean)==phenodata_MET500_clean$Sample_id)
#####################################################################
setwd(tablesDirectory)
write.csv(log_FPKM_MET500, file="log_FPKM_MET500.csv")
write.csv(zscore_log_FPKM_MET500, file="zscore_log_FPKM_MET500.csv")
write.csv(log_FPKM_MET500_clean, file="log_FPKM_MET500_clean.csv")
write.csv(zscore_log_FPKM_MET500_clean, file="zscore_log_FPKM_MET500_clean.csv")
#####################################################################


setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file="df_tissue_biopsy_tissue.csv", sep=",")
df_tissue_biopsy_tissue <- df_tissue_biopsy_tissue[,colnames(df_tissue_biopsy_tissue)!="X"]

phenodata_MET500_clean$tissue_set <- df_tissue_biopsy_tissue$set[match(phenodata_MET500_clean$tissue, df_tissue_biopsy_tissue$tissue)]
phenodata_MET500_clean$tejido <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$tissue, df_tissue_biopsy_tissue$tissue)]

phenodata_MET500_clean$biopsy_tissue_set <- df_tissue_biopsy_tissue$set[match(phenodata_MET500_clean$biopsy_tissue, df_tissue_biopsy_tissue$tissue)]
phenodata_MET500_clean$biopsia_tejido <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$biopsy_tissue, df_tissue_biopsy_tissue$tissue)]

#####################################################################
setwd(tablesDirectory)
write.csv(phenodata_MET500_clean, file="phenodata_MET500_clean.csv")
#####################################################################


#####################################################################
############# GRAPHS
#####################################################################

library(viridisLite)
colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

ggplot(data=phenodata_MET500_clean, aes(x=tejido, fill=biopsia_tejido))+
  geom_bar()+ theme(axis.text.x = element_text(size=10, angle=90)) + ggtitle("Muestras MET500 (limpias)") +
  xlab("Tejido") +
  scale_fill_manual(name = "Tejido biopsia",values=c(colors_blind[c(1,2,3,6)],"black", colors_blind[c(9,12)], "grey",colors_blind[c(15,17,19)], "red", colors_blind[27]))

