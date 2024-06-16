#####################################################################
# Genes of pvalue adj < 0.0001 & log2FC>3 & top 1000 in SVM
# zscore_log_FPKM
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/FINAL MODELS"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
setwd(tablesDirectory)
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")
zscore_log_fpkm <- read.csv(file="zscore_log_fpkm.csv")

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_log_fpkm) <- zscore_log_fpkm$X
zscore_log_fpkm <- zscore_log_fpkm[, colnames(zscore_log_fpkm)!="X"]

setwd(resultsDirectory)
genes_top1000 <- readRDS(file="genes_top1000.rds")

setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file = "df_tissue_biopsy_tissue.csv", header = TRUE)
######################################

zscore_log_fpkm.f <- zscore_log_fpkm[rownames(zscore_log_fpkm) %in% genes_top1000,]
zscore_log_fpkm.f <- zscore_log_fpkm.f[,colnames(zscore_log_fpkm.f) %in% phenoData_training_final$sample,]

# Top variables
zscore_log_fpkm.var <- zscore_log_fpkm.f

zscore_log_fpkm.var <- t(zscore_log_fpkm.var)
zscore_log_fpkm.var <- as.data.frame(zscore_log_fpkm.var)

table(rownames(zscore_log_fpkm.var)==phenoData_training_final$sample)
zscore_log_fpkm.var$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_log_fpkm.var), phenoData_training_final$sample)]

######################################
#### SVM
######################################

zscore_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_log_fpkm.var$enfermedad_primaria)

####
tunegrid <- expand.grid(C=c(0.001, 0.01, 0.1, 1, 10),
                        sigma=c(0.001, 0.01, 0.1, 1, 10))

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)



library(doParallel)
registerDoParallel(cores = 16)
set.seed(42) 
toc <- Sys.time()
model_SVM <- caret::train(enfermedad_primaria ~ .,
                                 data = zscore_log_fpkm.var,
                                 method = "svmRadial", # this will use the svmRadial function
                                 metric = "Accuracy", # which metric should be optimized for classification
                                 tuneGrid = tunegrid,
                                 trControl = control)

tic <- Sys.time()

tic-toc #Time difference of 13.78551 mins

bestTune_list <- model_SVM$results
bestTune_list$C <- as.numeric(bestTune_list$C)
bestTune_list$sigma <- as.factor(bestTune_list$sigma)
str(bestTune_list)


ggplot(bestTune_list, aes(x=C, y=Kappa, col=sigma)) + geom_line()
ggplot(bestTune_list, aes(x=C, y=Accuracy, col=sigma)) + geom_line()

#Best tune
model_SVM$bestTune
#sigma  C
#0.001 10


######################################
predSVM <- predict(model_SVM)
tab_cm <- confusionMatrix(predSVM,zscore_log_fpkm.var$enfermedad_primaria)
######################################
setwd(resultsDirectory)
write.csv(bestTune_list, file="bestTune_SVM_0.0001_5_top1000_zscore_log_FPKM.csv")


#### 

setwd(resultsDirectory)
save(model_SVM, file = "model_SVM.RData")
#Accuracy : 1
#Kappa : 1   

accuracy_train <- sum(diag(tab_cm$table))/sum(tab_cm$table)
sensibilidad_train <- mean(tab_cm$byClass[,1])
especificidad_train <- mean(tab_cm$byClass[,2])
F1_score_train <- mean(tab_cm$byClass[,7])
kappa_train <- 1
precision_train <- mean(tab_cm$byClass[,5])
NPV_train <- mean(tab_cm$byClass[,4])


######################################
########## PROBAR SVM SOBRE VALIDATION
######################################

setwd(tablesDirectory)
zscore_log_FPKM_test <- read.csv(file="zscore_log_fpkm_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_log_FPKM_test) <- zscore_log_FPKM_test$X
zscore_log_FPKM_test <- zscore_log_FPKM_test[, colnames(zscore_log_FPKM_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]
phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)

zscore_log_FPKM_test <- zscore_log_FPKM_test[rownames(zscore_log_FPKM_test) %in% colnames(zscore_log_fpkm.var),]
zscore_log_FPKM_test <- as.data.frame(t(zscore_log_FPKM_test))

table(colnames(zscore_log_FPKM_test)==colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"])
zscore_log_FPKM_test <- zscore_log_FPKM_test[, match(colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"], colnames(zscore_log_FPKM_test))]

pred_val <- predict(model_SVM, zscore_log_FPKM_test)
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)

#Accuracy : 0.9259   
#Kappa : 0.9231  

accuracy_val <- sum(diag(tab_cm_VAL$table))/sum(tab_cm_VAL$table)
sensibilidad_val <- mean(tab_cm_VAL$byClass[,1])
especificidad_val <- mean(tab_cm_VAL$byClass[,2])
F1_score_val <- mean(tab_cm_VAL$byClass[,7])
kappa_val <- 0.9231
precision_val <- mean(tab_cm_VAL$byClass[,5])
NPV_val <- mean(tab_cm_VAL$byClass[,4])


######################################
setwd(tablesDirectory)
zscore_log_FPKM_MET500 <- read.csv(file="zscore_log_fpkm_MET500_clean.csv")
phenodata_MET500_clean <- read.csv(file="phenodata_MET500_clean.csv")

rownames(zscore_log_FPKM_MET500) <- zscore_log_FPKM_MET500$X
zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, colnames(zscore_log_FPKM_MET500)!="X"]
colnames(zscore_log_FPKM_MET500) <- str_replace_all(colnames(zscore_log_FPKM_MET500), "[.]", "-")
table(colnames(zscore_log_FPKM_MET500) %in% phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, match(phenodata_MET500_clean$Sample_id, colnames(zscore_log_FPKM_MET500))]
table(colnames(zscore_log_FPKM_MET500)==phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

phenodata_MET500_clean$enfermedad_primaria <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$tissue_set, df_tissue_biopsy_tissue$set)]
phenodata_MET500_clean$enfermedad_primaria <- as.factor(phenodata_MET500_clean$enfermedad_primaria)

zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[rownames(zscore_log_FPKM_MET500) %in% colnames(zscore_log_fpkm.var),]
zscore_log_FPKM_MET500 <- as.data.frame(t(zscore_log_FPKM_MET500))

table(colnames(zscore_log_FPKM_MET500)==colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"])
zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, match(colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"], colnames(zscore_log_FPKM_MET500))]


predRF_MET500 <- predict(model_SVM, zscore_log_FPKM_MET500)
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)

#Accuracy : 0.3747 
#Kappa : 0.3146

tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)
tejidos_MET500 <- paste("Class:",tejidos_MET500)

accuracy_MET500 <- sum(diag(tab_cm_MET500$table))/sum(tab_cm_MET500$table)
sensibilidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,1])
especificidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,2])
F1_score_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,7])
kappa_MET500 <-0.3146 
precision_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,5])
NPV_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,4])

######################################
setwd(resultsDirectory)
saveRDS(tab_cm_MET500, file="tab_cm_MET500_SVM_top1000_zscore_log_FPKM.rds")


######################################
#### STUDY MODEL
######################################
importance <- model_rf$importance
setwd(resultsDirectory)
write.csv(importance, "importance_RF.csv")

######################################
#### STUDY BY CLASS
######################################

metrics_MET500 <- tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass) %in% tejidos_MET500,]

metrics_training <- tab_cm$byClass

metrics_VAL <- tab_cm_VAL$byClass

df_byClass_training <- round(tab_cm$byClass[,c(1,2,3,4,5,7)],4)
rownames(df_byClass_training) <- sub("Class:", "", rownames(df_byClass_training))
colnames(df_byClass_training) <- c("Sensibilidad", "Especificidad", "PPV", "NPV", "Precisión", "F1-score")
df_byClass_VAL <- round(tab_cm_VAL$byClass[,c(1,2,3,4,5,7)],4)
rownames(df_byClass_VAL) <- sub("Class:", "", rownames(df_byClass_VAL))
colnames(df_byClass_VAL) <- c("Sensibilidad", "Especificidad", "PPV", "NPV", "Precisión", "F1-score")
df_byClass_MET500 <- round(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass) %in% tejidos_MET500,c(1,2,3,4,5,7)],4)
rownames(df_byClass_MET500) <- sub("Class:", "", rownames(df_byClass_MET500))
colnames(df_byClass_MET500) <- c("Sensibilidad", "Especificidad", "PPV", "NPV", "Precisión", "F1-score")

setwd(resultsDirectory)
write.csv(df_byClass_training, file="df_byClass_training_SVM.csv")
write.csv(df_byClass_VAL, file="df_byClass_VAL_SVM.csv")
write.csv(df_byClass_MET500, file="df_byClass_MET500_SVM.csv")
