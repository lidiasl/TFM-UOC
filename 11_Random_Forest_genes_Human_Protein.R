#####################################################################
# Genes of Human Protein in Random Forest
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
zscore_variable_log_fpkm <- read.csv(file="zscore_variable_log_fpkm.csv")

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_variable_log_fpkm) <- zscore_variable_log_fpkm$X
zscore_variable_log_fpkm <- zscore_variable_log_fpkm[, colnames(zscore_variable_log_fpkm)!="X"]

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


zscore_variable_log_fpkm.f <- zscore_variable_log_fpkm[rownames(zscore_variable_log_fpkm) %in% genes_intersect,]
zscore_variable_log_fpkm.f <- zscore_variable_log_fpkm.f[,colnames(zscore_variable_log_fpkm.f) %in% phenoData_training_final$sample,]
zscore_variable_log_fpkm.f <- t(zscore_variable_log_fpkm.f)
zscore_variable_log_fpkm.f <- as.data.frame(zscore_variable_log_fpkm.f)

table(rownames(zscore_variable_log_fpkm.f)==phenoData_training_final$sample)
zscore_variable_log_fpkm.f$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_variable_log_fpkm.f), phenoData_training_final$sample)]

######################################
#### RANDOM FOREST
######################################
zscore_variable_log_fpkm.f$enfermedad_primaria <- as.factor(zscore_variable_log_fpkm.f$enfermedad_primaria)

####
tunegrid <- expand.grid(mtry=seq(35,410,25))
ntree_list <- seq(500,1500,100)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)

bestTune_list<- list()

library(doParallel)
registerDoParallel(cores = 16)
set.seed(42) 
toc <- Sys.time()
for (i in 1:length(ntree_list)){
  model_rf <- caret::train(enfermedad_primaria ~ .,
                           data = zscore_variable_log_fpkm.f,
                           method = "rf", 
                           metric = "Accuracy", 
                           tuneGrid = tunegrid,
                           trControl = control,
                           ntree = ntree_list[i],
                           keep.forest = TRUE,
                           importance = TRUE)
  bestTune_list[[i]] <- model_rf[["results"]]
  print(i)
}


tic <- Sys.time()

tic-toc #Time difference of 16.1135 hours


for (i in 1:length(ntree_list)){
  bestTune_list[[i]]$ntree <- ntree_list[i]
}

bestTune_df <- Reduce(rbind, bestTune_list)
bestTune_df$ntree <- as.factor(bestTune_df$ntree)
ggplot(bestTune_df, aes(x=mtry, y=Kappa, col=ntree)) + geom_line()
ggplot(bestTune_df, aes(x=mtry, y=Accuracy, col=ntree)) + geom_line()

######################################
setwd(resultsDirectory)
write.csv(bestTune_df, file="bestTune_df_Human_Protein.csv")


#### 
library(randomForest)
set.seed(42)
model_rf <- randomForest(x = as.matrix(zscore_variable_log_fpkm.f[, colnames(zscore_variable_log_fpkm.f)!="enfermedad_primaria"]),
                         y=zscore_variable_log_fpkm.f$enfermedad_primaria,
                         mtry=160,
                         ntree=1400,
                         keep.forest = TRUE,
                         importance = TRUE)
predRF <- predict(model_rf)
tab_cm <- confusionMatrix(predRF,zscore_variable_log_fpkm.f$enfermedad_primaria)
#Accuracy : 0.9526 
#Kappa : 0.9508

accuracy_train <- sum(diag(tab_cm$table))/sum(tab_cm$table)
sensibilidad_train <- mean(tab_cm$byClass[,1])
especificidad_train <- mean(tab_cm$byClass[,2])
F1_score_train <- mean(tab_cm$byClass[,7])
kappa_train <- 0.9508
precision_train <- mean(tab_cm$byClass[,5])

######################################
setwd(resultsDirectory)
saveRDS(tab_cm, file="tab_cm_train_RF_Human_Protein.rds")

######################################
########## LASSO VALIDATION
######################################

setwd(tablesDirectory)
zscore_variable_log_FPKM_test <- read.csv(file="zscore_variable_log_FPKM_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_variable_log_FPKM_test) <- zscore_variable_log_FPKM_test$X
zscore_variable_log_FPKM_test <- zscore_variable_log_FPKM_test[, colnames(zscore_variable_log_FPKM_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]
phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)

zscore_variable_log_FPKM_test <- zscore_variable_log_FPKM_test[rownames(zscore_variable_log_FPKM_test) %in% genes_intersect,]
zscore_variable_log_FPKM_test <- as.data.frame(t(zscore_variable_log_FPKM_test))

table(colnames(zscore_variable_log_FPKM_test)==colnames(zscore_variable_log_fpkm.f)[colnames(zscore_variable_log_fpkm.f)!="enfermedad_primaria"])

pred_val <- predict(model_rf, zscore_variable_log_FPKM_test)
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)

#Accuracy : 0.9101
#Kappa : 0.9066 

accuracy_val <- sum(diag(tab_cm_VAL$table))/sum(tab_cm_VAL$table)
sensibilidad_val <- mean(tab_cm_VAL$byClass[,1])
especificidad_val <- mean(tab_cm_VAL$byClass[,2])
F1_score_val <- mean(tab_cm_VAL$byClass[,7])
kappa_val <- 0.9066   
precision_val <- mean(tab_cm_VAL$byClass[,5])
######################################
setwd(resultsDirectory)
saveRDS(tab_cm_VAL, file="tab_cm_VAL_RF_Human_Protein.rds")


######################################
setwd(tablesDirectory)
zscore_variable_log_FPKM_MET500 <- read.csv(file="zscore_variable_log_FPKM_MET500_clean.csv")
phenodata_MET500_clean <- read.csv(file="phenodata_MET500_clean.csv")

rownames(zscore_variable_log_FPKM_MET500) <- zscore_variable_log_FPKM_MET500$X
zscore_variable_log_FPKM_MET500 <- zscore_variable_log_FPKM_MET500[, colnames(zscore_variable_log_FPKM_MET500)!="X"]
colnames(zscore_variable_log_FPKM_MET500) <- str_replace_all(colnames(zscore_variable_log_FPKM_MET500), "[.]", "-")
table(colnames(zscore_variable_log_FPKM_MET500) %in% phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

zscore_variable_log_FPKM_MET500 <- zscore_variable_log_FPKM_MET500[, match(phenodata_MET500_clean$Sample_id, colnames(zscore_variable_log_FPKM_MET500))]
table(colnames(zscore_variable_log_FPKM_MET500)==phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

phenodata_MET500_clean$enfermedad_primaria <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$tissue_set, df_tissue_biopsy_tissue$set)]
phenodata_MET500_clean$enfermedad_primaria <- as.factor(phenodata_MET500_clean$enfermedad_primaria)

zscore_variable_log_FPKM_MET500 <- zscore_variable_log_FPKM_MET500[rownames(zscore_variable_log_FPKM_MET500) %in% genes_intersect,]
zscore_variable_log_FPKM_MET500 <- as.data.frame(t(zscore_variable_log_FPKM_MET500))

table(colnames(zscore_variable_log_FPKM_MET500)==colnames(zscore_variable_log_fpkm.f)[(colnames(zscore_variable_log_fpkm.f)!="enfermedad_primaria")])
zscore_variable_log_FPKM_MET500 <- zscore_variable_log_FPKM_MET500[, match(colnames(zscore_variable_log_fpkm.f)[(colnames(zscore_variable_log_fpkm.f)!="enfermedad_primaria")], colnames(zscore_variable_log_FPKM_MET500))]

predRF_MET500 <- predict(model_rf, zscore_variable_log_FPKM_MET500)
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)

#Accuracy : 0.2248
#Kappa : 0.0381

tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)
tejidos_MET500 <- paste("Class:",tejidos_MET500)

accuracy_MET500 <- sum(diag(tab_cm_MET500$table))/sum(tab_cm_MET500$table)
sensibilidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,1])
especificidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,2])
F1_score_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,7])
kappa_MET500 <- 0.0381    
precision_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,5])

######################################
setwd(resultsDirectory)
saveRDS(tab_cm_MET500, file="tab_cm_MET500_RF_Human_Protein.rds")



