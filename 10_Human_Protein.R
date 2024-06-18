#####################################################################
# Genes of Human Protein
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
library(fgsea)


# Directories
tablesDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Tables"
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Human Protein"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
setwd(tablesDirectory)
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")


rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

######################################
setwd(dataDirectory)
#### Adrenocortical
# https://www.proteinatlas.org/humanproteome/tissue/adrenal+gland#protein_expression_of_genes_elevated_in_adrenal_gland
genes_adrenocortical_enhanced <- fread(file = "tissue_category_rna_adrenal_enchanced.tsv")
genes_adrenocortical_enhanced_all <- genes_adrenocortical_enhanced$Ensembl[genes_adrenocortical_enhanced$`RNA tissue distribution` %in% 
                                                                             c("Detected in all","Detected in many")]
genes_adrenocortical_enhanced <- genes_adrenocortical_enhanced$Ensembl

#### B_cell_lymphoma
#

#### bladder
#https://www.proteinatlas.org/humanproteome/tissue/urinary+bladder
genes_bladder_enhanced <- fread(file = "tissue_category_rna_urinary_enhanced.tsv")
genes_bladder_enhanced_all <- genes_bladder_enhanced$Ensembl[genes_bladder_enhanced$`RNA tissue distribution`%in% 
                                                               c("Detected in all","Detected in many")]
genes_bladder_enhanced <- genes_bladder_enhanced$Ensembl


#### brain
#https://www.proteinatlas.org/humanproteome/brain/human+brain
genes_brain_enhanced <- fread(file = "tissue_category_rna_brain_Tissue_enhanced.tsv")
genes_brain_enhanced_all <- genes_brain_enhanced$Ensembl[genes_brain_enhanced$`RNA tissue distribution`%in% 
                                                           c("Detected in all","Detected in many")]
genes_brain_enhanced <- genes_brain_enhanced$Ensembl


#### breast
#https://www.proteinatlas.org/humanproteome/tissue/breast
genes_breast_enhanced <- fread(file = "tissue_category_rna_breast_Tissue_enhanced.tsv")
genes_breast_enhanced_all <- genes_breast_enhanced$Ensembl[genes_breast_enhanced$`RNA tissue distribution`%in% 
                                                             c("Detected in all","Detected in many")]
genes_breast_enhanced <- genes_breast_enhanced$Ensembl


#### cervical
#https://www.proteinatlas.org/humanproteome/tissue/cervix#protein_expression_of_genes_elevated_in_cervix
genes_cervical_enhanced <- fread(file = "tissue_category_rna_cervix_Tissue_enhanced.tsv")
genes_cervical_enhanced_all <- genes_cervical_enhanced$Ensembl[genes_cervical_enhanced$`RNA tissue distribution`%in% 
                                                                 c("Detected in all","Detected in many")]
genes_cervical_enhanced <- genes_cervical_enhanced$Ensembl


#### cholangiocarcinoma
#

#### colon
#https://www.proteinatlas.org/humanproteome/tissue/intestine
genes_colon_enhanced <- fread(file = "tissue_category_rna_intestine_Tissue_enhanced.tsv")
genes_colon_enhanced_all <- genes_colon_enhanced$Ensembl[genes_colon_enhanced$`RNA tissue distribution`%in% 
                                                           c("Detected in all","Detected in many")]
genes_colon_enhanced <- genes_colon_enhanced$Ensembl


#### esophageal
#https://www.proteinatlas.org/humanproteome/tissue/esophagus
genes_esophageal_enhanced <- fread(file = "tissue_category_rna_esophagus_Tissue_enhanced.tsv")
genes_esophageal_enhanced_all <- genes_esophageal_enhanced$Ensembl[genes_esophageal_enhanced$`RNA tissue distribution`%in% 
                                                                     c("Detected in all","Detected in many")]
genes_esophageal_enhanced <- genes_esophageal_enhanced$Ensembl


#### head_neck
#

#### kidney
#https://www.proteinatlas.org/humanproteome/tissue/kidney
genes_kidney_enhanced <- fread(file = "tissue_category_rna_kidney_Tissue_enhanced.tsv")
genes_kidney_enhanced_all <- genes_kidney_enhanced$Ensembl[genes_kidney_enhanced$`RNA tissue distribution`%in% 
                                                             c("Detected in all","Detected in many")]
genes_kidney_enhanced <- genes_kidney_enhanced$Ensembl


#### leukemia
#

#### liver
#https://www.proteinatlas.org/humanproteome/tissue/liver
genes_liver_enhanced <- fread(file = "tissue_category_rna_liver_Tissue_enhanced.tsv")
genes_liver_enhanced_all <- genes_liver_enhanced$Ensembl[genes_liver_enhanced$`RNA tissue distribution`%in% 
                                                           c("Detected in all","Detected in many")]
genes_liver_enhanced <- genes_liver_enhanced$Ensembl


#### lung
#https://www.proteinatlas.org/humanproteome/tissue/lung
genes_lung_enhanced <- fread(file = "tissue_category_rna_lung_Tissue_enhanced.tsv")
genes_lung_enhanced_all <- genes_lung_enhanced$Ensembl[genes_lung_enhanced$`RNA tissue distribution`%in% 
                                                         c("Detected in all","Detected in many")]
genes_lung_enhanced <- genes_lung_enhanced$Ensembl


#### mesothelioma
#

#### ovarian
#https://www.proteinatlas.org/humanproteome/tissue/ovary
genes_ovarian_enhanced <- fread(file = "tissue_category_rna_ovary_Tissue_enhanced.tsv")
genes_ovarian_enhanced_all <- genes_ovarian_enhanced$Ensembl[genes_ovarian_enhanced$`RNA tissue distribution`%in% 
                                                               c("Detected in all","Detected in many")]
genes_ovarian_enhanced <- genes_ovarian_enhanced$Ensembl


#### pancreas
#https://www.proteinatlas.org/humanproteome/tissue/pancreas
genes_pancreas_enhanced <- fread(file = "tissue_category_rna_pancreas_Tissue_enhanced.tsv")
genes_pancreas_enhanced_all <- genes_pancreas_enhanced$Ensembl[genes_pancreas_enhanced$`RNA tissue distribution`%in% 
                                                                 c("Detected in all","Detected in many")]
genes_pancreas_enhanced <- genes_pancreas_enhanced$Ensembl


#### paraganglioma
#

#### prostate
#https://www.proteinatlas.org/humanproteome/tissue/prostate
genes_prostate_enhanced <- fread(file = "tissue_category_rna_prostate_Tissue_enhanced.tsv")
genes_prostate_enhanced_all <- genes_prostate_enhanced$Ensembl[genes_prostate_enhanced$`RNA tissue distribution`%in% 
                                                                 c("Detected in all","Detected in many")]
genes_prostate_enhanced <- genes_prostate_enhanced$Ensembl


#### sarcoma
#

#### skin
#https://www.proteinatlas.org/humanproteome/tissue/skin
genes_skin_enhanced <- fread(file = "tissue_category_rna_skin_Tissue_enhanced.tsv")
genes_skin_enhanced_all <- genes_skin_enhanced$Ensembl[genes_skin_enhanced$`RNA tissue distribution`%in% 
                                                         c("Detected in all","Detected in many")]
genes_skin_enhanced <- genes_skin_enhanced$Ensembl


#### stomach
#https://www.proteinatlas.org/humanproteome/tissue/stomach
genes_stomach_enhanced <- fread(file = "tissue_category_rna_stomach_Tissue_enhanced.tsv")
genes_stomach_enhanced_all <- genes_stomach_enhanced$Ensembl[genes_stomach_enhanced$`RNA tissue distribution`%in% 
                                                               c("Detected in all","Detected in many")]
genes_stomach_enhanced <- genes_stomach_enhanced$Ensembl


#### testicles
#https://www.proteinatlas.org/humanproteome/tissue/testis
genes_testicles_enhanced <- fread(file = "tissue_category_rna_testis_Tissue_enhanced.tsv")
genes_testicles_enhanced_all <- genes_testicles_enhanced$Ensembl[genes_testicles_enhanced$`RNA tissue distribution`%in% 
                                                                   c("Detected in all","Detected in many")]
genes_testicles_enhanced <- genes_testicles_enhanced$Ensembl


#### thymoma
#

#### thyroid
#https://www.proteinatlas.org/humanproteome/tissue/thyroid+gland
genes_thyroid_enhanced <- fread(file = "tissue_category_rna_thyroid_enhanced.tsv")
genes_thyroid_enhanced_all <- genes_thyroid_enhanced$Ensembl[genes_thyroid_enhanced$`RNA tissue distribution`%in% 
                                                               c("Detected in all","Detected in many")]
genes_thyroid_enhanced <- genes_thyroid_enhanced$Ensembl


#### uterine

#### uveal

###################################################

genes_union_enhanced <- c(genes_adrenocortical_enhanced,
                          genes_bladder_enhanced,
                          genes_brain_enhanced,
                          genes_breast_enhanced,
                          genes_cervical_enhanced,
                          genes_colon_enhanced,
                          genes_esophageal_enhanced,
                          genes_kidney_enhanced,
                          genes_liver_enhanced,
                          genes_lung_enhanced,
                          genes_ovarian_enhanced,
                          genes_pancreas_enhanced,
                          genes_prostate_enhanced,
                          genes_skin_enhanced,
                          genes_stomach_enhanced,
                          genes_testicles_enhanced,
                          genes_thyroid_enhanced)

genes_union_enhanced <- unique(genes_union_enhanced)

table(is.na(genes_union_enhanced))
#FALSE 
#4137 

########################################
genes_union_enhanced_all <- c(genes_adrenocortical_enhanced_all,
                          genes_bladder_enhanced_all,
                          genes_brain_enhanced_all,
                          genes_breast_enhanced_all,
                          genes_cervical_enhanced_all,
                          genes_colon_enhanced_all,
                          genes_esophageal_enhanced_all,
                          genes_kidney_enhanced_all,
                          genes_liver_enhanced_all,
                          genes_lung_enhanced_all,
                          genes_ovarian_enhanced_all,
                          genes_pancreas_enhanced_all,
                          genes_prostate_enhanced_all,
                          genes_skin_enhanced_all,
                          genes_stomach_enhanced_all,
                          genes_testicles_enhanced_all,
                          genes_thyroid_enhanced_all)

genes_union_enhanced_all <- unique(genes_union_enhanced_all)

table(is.na(genes_union_enhanced_all))

genes_union_enhanced_all_list <- list(genes_adrenocortical_enhanced_all,
                              genes_bladder_enhanced_all,
                              genes_brain_enhanced_all,
                              genes_breast_enhanced_all,
                              genes_cervical_enhanced_all,
                              genes_colon_enhanced_all,
                              genes_esophageal_enhanced_all,
                              genes_kidney_enhanced_all,
                              genes_liver_enhanced_all,
                              genes_lung_enhanced_all,
                              genes_ovarian_enhanced_all,
                              genes_pancreas_enhanced_all,
                              genes_prostate_enhanced_all,
                              genes_skin_enhanced_all,
                              genes_stomach_enhanced_all,
                              genes_testicles_enhanced_all,
                              genes_thyroid_enhanced_all)
names(genes_union_enhanced_all_list) <- c("genes_adrenocortical_enhanced_all",
                                          "genes_bladder_enhanced_all",
                                          "genes_brain_enhanced_all",
                                          "genes_breast_enhanced_all",
                                          "genes_cervical_enhanced_all",
                                          "genes_colon_enhanced_all",
                                          "genes_esophageal_enhanced_all",
                                          "genes_kidney_enhanced_all",
                                          "genes_liver_enhanced_all",
                                          "genes_lung_enhanced_all",
                                          "genes_ovarian_enhanced_all",
                                          "genes_pancreas_enhanced_all",
                                          "genes_prostate_enhanced_all",
                                          "genes_skin_enhanced_all",
                                          "genes_stomach_enhanced_all",
                                          "genes_testicles_enhanced_all",
                                          "genes_thyroid_enhanced_all")
###################################################
setwd(resultsDirectory)
saveRDS(genes_union_enhanced, file="genes_union_enhanced.rds")
saveRDS(genes_union_enhanced_all, file="genes_union_enhanced_all.rds")
saveRDS(genes_union_enhanced_all_list, file="genes_union_enhanced_all_list.rds")
###################################################

genes_union_enhanced_list <- list(genes_adrenocortical_enhanced,
                               genes_bladder_enhanced,
                               genes_brain_enhanced,
                               genes_breast_enhanced,
                               genes_cervical_enhanced,
                               genes_colon_enhanced,
                               genes_esophageal_enhanced,
                               genes_kidney_enhanced,
                               genes_liver_enhanced,
                               genes_lung_enhanced,
                               genes_ovarian_enhanced,
                               genes_pancreas_enhanced,
                               genes_prostate_enhanced,
                               genes_skin_enhanced,
                               genes_stomach_enhanced,
                               genes_testicles_enhanced,
                               genes_thyroid_enhanced)

names(genes_union_enhanced_list) <- c("adrenocortical","bladder","brain","breast","cervical",
                                      "colon","esophageal","kidney","liver","lung","ovarian","pancreas","prostate",
                                      "skin", "stomach","testicles","thyroid")

setwd(tablesDirectory)
padj_log2FC <- read.csv(file="padj_log2FC.csv")
padj_log2FC <- padj_log2FC[,colnames(padj_log2FC)!="X"]


setwd(tablesDirectory)
genes_list_combinations <- readRDS(file="genes_list_all_combinations_padj_log2FC.rds")


#### UNION
genes_included_union <- data.frame()

for (i in 1:nrow(padj_log2FC)){
  genes_union <- lapply(genes_list_combinations[[i]], function(x){
    Reduce(union,x)
  })
  
  vector <- c()
  
  for (j in 1:length(names(genes_union))){
    if (names(genes_union)[j] %in% names(genes_union_enhanced_list)){
      value <- sum(genes_union_enhanced_list[[names(genes_union)[j]]] %in% genes_union[[j]]==TRUE)/length(genes_union_enhanced_list[[names(genes_union)[j]]])*100
      vector <- c(vector, value)
    }
  }
  
  genes_included_union <- rbind(genes_included_union, vector)
}

colnames(genes_included_union) <- names(genes_union)[names(genes_union) %in% names(genes_union_enhanced_list)]
genes_included_union <- cbind(padj_log2FC,genes_included_union)

setwd(resultsDirectory)
write.csv(genes_included_union, file="genes_included_union.csv")

###################################################
# Study top 1000 genes in each tissue
setwd("~/Desktop/LIDIA/TCGA_classification/TFM/Results/FINAL MODELS")
genes_top1000 <- readRDS(file="genes_top1000.rds")

inclusion_list <- list()
for (i in 1:length(genes_union_enhanced_list)){
  inclusion_list[[i]] <- sum(genes_top1000 %in% genes_union_enhanced_list[[i]])
}

names(inclusion_list) <- names(genes_union_enhanced_list)
inclusion_df <- as.data.frame(inclusion_list)
inclusion_df <- as.data.frame(t(inclusion_df))
colnames(inclusion_df) <- "numero"
inclusion_df$enfermedad_primaria <- rownames(inclusion_df)
num_genes_inclusion <- ggplot(inclusion_df, aes(x = enfermedad_primaria, y = numero)) +
  geom_bar(stat = "identity", fill = "blue", alpha=0.5) +
  geom_text(aes(label = numero),angle = 90, vjust = 0.5, hjust=0.5) +
  labs(title = "Número de genes del conjunto de predictores incluidos en el HP",
       x = "Tumor primario",
       y = "Número de genes") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12), 
    plot.title = element_text(size = 15) 
  )

###################################################
# Do the same with top variables in RF
importance <- read.csv(file="importance_RF.csv")
importance <- importance[,c("X","MeanDecreaseGini")]
importance_top100 <- importance[order(importance$MeanDecreaseGini, decreasing = TRUE)[1:100],]


inclusion_RF_list <- list()
for (i in 1:length(genes_union_enhanced_list)){
  inclusion_RF_list[[i]] <- sum(importance_top100$X %in% genes_union_enhanced_list[[i]])
}

names(inclusion_RF_list) <- names(genes_union_enhanced_list)
inclusion_RF_df <- as.data.frame(inclusion_RF_list)
inclusion_RF_df <- as.data.frame(t(inclusion_RF_df))
colnames(inclusion_RF_df) <- "numero"
inclusion_RF_df$enfermedad_primaria <- rownames(inclusion_RF_df)
num_genes_inclusion_RF <- ggplot(inclusion_RF_df, aes(x = enfermedad_primaria, y = numero)) +
  geom_bar(stat = "identity", fill = "blue", alpha=0.5) +
  geom_text(aes(label = numero),angle = 90, vjust = 0.5, hjust=0.5) +
  labs(title = "Número de genes con mayor Mean Decrease Gini (top 100)",
       x = "Tumor primario",
       y = "Número de genes") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12), 
    plot.title = element_text(size = 15) 
  )
