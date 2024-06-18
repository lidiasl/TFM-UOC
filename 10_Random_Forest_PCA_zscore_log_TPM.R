#####################################################################
# Genes PCA TPM in Random Forest
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Random Forest Prueba"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
setwd(tablesDirectory)
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")
zscore_log_tpm <- read.csv(file="zscore_log_tpm.csv")

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_log_tpm) <- zscore_log_tpm$X
zscore_log_tpm <- zscore_log_tpm[, colnames(zscore_log_tpm)!="X"]

genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")
setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file = "df_tissue_biopsy_tissue.csv", header = TRUE)

setwd("~/Desktop/LIDIA/TCGA_classification/TFM/Results/Human Protein")
genes_union_enhanced <- readRDS(file="genes_union_enhanced.rds")
######################################


genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

genes_list_0.0001_5 <- Reduce(c,genes_list_0.0001_5) 
genes_list_0.0001_5 <- unique(genes_list_0.0001_5) #7339
length(genes_list_0.0001_5) #7339

table(genes_union_enhanced %in% genes_list_0.0001_5)
#FALSE  TRUE 
#1663  2474 

genes_intersect <- genes_union_enhanced[genes_union_enhanced %in% genes_list_0.0001_5]


zscore_log_tpm.f <- zscore_log_tpm[rownames(zscore_log_tpm) %in% genes_intersect,]
zscore_log_tpm.f <- zscore_log_tpm.f[,colnames(zscore_log_tpm.f) %in% phenoData_training_final$sample,]
zscore_log_tpm.f <- t(zscore_log_tpm.f)
zscore_log_tpm.f <- as.data.frame(zscore_log_tpm.f)

table(rownames(zscore_log_tpm.f)==phenoData_training_final$sample)
zscore_log_tpm.f$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_log_tpm.f), phenoData_training_final$sample)]


######################################
####### PCA
######################################

set.seed(42)
PCA <- prcomp(as.matrix(zscore_log_tpm.f[,colnames(zscore_log_tpm.f)!="enfermedad_primaria"]), 
              scale. = FALSE, 
              center = FALSE)

library(factoextra)
fviz_eig(PCA, addlabels = FALSE, labelsize=0.1,ncp=30)+ 
  geom_hline(yintercept = 0.1, color = "red") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
var <- get_pca_var(PCA)


transformacion_PCA <- PCA$rotation

# First 100 components

PCA_training <- PCA$x
table(rownames(PCA_training) == phenoData_training_final$sample)

PCA_training <- as.data.frame(PCA_training)
PCA_training <- PCA_training[, 1:100]

PCA_training$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(PCA_training), phenoData_training_final$sample)]
PCA_training$enfermedad_primaria <- as.factor(PCA_training$enfermedad_primaria)


####
tunegrid <- expand.grid(mtry=seq(3,16,2))
ntree_list <- seq(500,1500,100)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)

bestTune_list<- list()

library(doParallel)
registerDoParallel(cores = 8)
set.seed(42) 
toc <- Sys.time()
for (i in 1:length(ntree_list)){
  model_rf <- caret::train(enfermedad_primaria ~ .,
                           data = PCA_training,
                           method = "rf", # this will use the randomForest::randomForest function
                           metric = "Accuracy", # which metric should be optimized for classification
                           # options to be passed to randomForest
                           tuneGrid = tunegrid,
                           trControl = control,
                           ntree = ntree_list[i],
                           keep.forest = TRUE,
                           importance = TRUE)
  bestTune_list[[i]] <- model_rf[["results"]]
  print(i)
}


tic <- Sys.time()

tic-toc #Time difference of 30.07904 mins


for (i in 1:length(ntree_list)){
  bestTune_list[[i]]$ntree <- ntree_list[i]
}

bestTune_df <- Reduce(rbind, bestTune_list)
bestTune_df$ntree <- as.factor(bestTune_df$ntree)


#### 
library(randomForest)
set.seed(42)
model_rf <- randomForest(x = as.matrix(PCA_training[, colnames(PCA_training)!="enfermedad_primaria"]),
                         y=PCA_training$enfermedad_primaria,
                         mtry=3,
                         ntree=1400,
                         keep.forest = TRUE,
                         importance = TRUE)
predRF <- predict(model_rf)
tab_cm <- confusionMatrix(predRF,PCA_training$enfermedad_primaria)
#Accuracy : 0.9645
#Kappa : 0.9631

######################################
########## RF VALIDATION COHORT
######################################

setwd(tablesDirectory)
zscore_log_tpm_test <- read.csv(file="zscore_log_tpm_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_log_tpm_test) <- zscore_log_tpm_test$X
zscore_log_tpm_test <- zscore_log_tpm_test[, colnames(zscore_log_tpm_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]
phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)

zscore_log_tpm_test <- zscore_log_tpm_test[rownames(zscore_log_tpm_test) %in% genes_intersect,]
zscore_log_tpm_test <- as.data.frame(t(zscore_log_tpm_test))

table(colnames(zscore_log_tpm_test)==rownames(transformacion_PCA))

zscore_log_TPM_test_PCA <- as.matrix(zscore_log_tpm_test) %*% transformacion_PCA
zscore_log_TPM_test_PCA <- as.data.frame(zscore_log_TPM_test_PCA)
zscore_log_TPM_test_PCA <- zscore_log_TPM_test_PCA[, 1:100]

phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)


pred_val <- predict(model_rf, zscore_log_TPM_test_PCA)
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)
#Accuracy : 0.9101
#Kappa : 0.9066 



######################################
setwd(tablesDirectory)
zscore_log_tpm_MET500 <- read.csv(file="zscore_log_tpm_MET500_clean.csv")
phenodata_MET500_clean <- read.csv(file="phenodata_MET500_clean.csv")

rownames(zscore_log_tpm_MET500) <- zscore_log_tpm_MET500$X
zscore_log_tpm_MET500 <- zscore_log_tpm_MET500[, colnames(zscore_log_tpm_MET500)!="X"]
colnames(zscore_log_tpm_MET500) <- str_replace_all(colnames(zscore_log_tpm_MET500), "[.]", "-")
table(colnames(zscore_log_tpm_MET500) %in% phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

zscore_log_tpm_MET500 <- zscore_log_tpm_MET500[, match(phenodata_MET500_clean$Sample_id, colnames(zscore_log_tpm_MET500))]
table(colnames(zscore_log_tpm_MET500)==phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

phenodata_MET500_clean$enfermedad_primaria <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$tissue_set, df_tissue_biopsy_tissue$set)]
phenodata_MET500_clean$enfermedad_primaria <- as.factor(phenodata_MET500_clean$enfermedad_primaria)

zscore_log_tpm_MET500 <- zscore_log_tpm_MET500[rownames(zscore_log_tpm_MET500) %in% genes_intersect,]
zscore_log_tpm_MET500 <- as.data.frame(t(zscore_log_tpm_MET500))

zscore_log_TPM_MET500_PCA <- as.matrix(zscore_log_tpm_MET500) %*% transformacion_PCA
zscore_log_TPM_MET500_PCA <- as.data.frame(zscore_log_TPM_MET500_PCA)
zscore_log_TPM_MET500_PCA <- zscore_log_TPM_MET500_PCA[, 1:100]

predRF_MET500 <- predict(model_rf, zscore_log_TPM_MET500_PCA)
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)
#Accuracy : 0.0258
#Kappa : 0.0065
