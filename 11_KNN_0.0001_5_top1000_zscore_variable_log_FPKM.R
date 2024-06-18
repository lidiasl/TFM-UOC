#####################################################################
# Genes of pvalue adj < 0.0001 & log2FC>3 & top 1000 in KNN zscore by gene
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
zscore_variable_log_fpkm <- read.csv(file="zscore_variable_log_fpkm.csv")

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_variable_log_fpkm) <- zscore_variable_log_fpkm$X
zscore_variable_log_fpkm <- zscore_variable_log_fpkm[, colnames(zscore_variable_log_fpkm)!="X"]

setwd(resultsDirectory)
genes_top1000 <- readRDS(file="genes_top1000.rds")

setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file = "df_tissue_biopsy_tissue.csv", header = TRUE)
######################################

zscore_variable_log_fpkm.f <- zscore_variable_log_fpkm[rownames(zscore_variable_log_fpkm) %in% genes_top1000,]
zscore_variable_log_fpkm.f <- zscore_variable_log_fpkm.f[,colnames(zscore_variable_log_fpkm.f) %in% phenoData_training_final$sample,]

# Top variables
zscore_variable_log_fpkm.var <- zscore_variable_log_fpkm.f

zscore_variable_log_fpkm.var <- t(zscore_variable_log_fpkm.var)
zscore_variable_log_fpkm.var <- as.data.frame(zscore_variable_log_fpkm.var)

table(rownames(zscore_variable_log_fpkm.var)==phenoData_training_final$sample)
zscore_variable_log_fpkm.var$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_variable_log_fpkm.var), phenoData_training_final$sample)]

######################################
#### KNN
######################################

zscore_variable_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_variable_log_fpkm.var$enfermedad_primaria)

####
tunegrid <- expand.grid(k=seq(1,10,1))

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)



library(doParallel)
registerDoParallel(cores = 16)
set.seed(42) 
toc <- Sys.time()
model_KNN <- caret::train(enfermedad_primaria ~ .,
                          data = zscore_variable_log_fpkm.var,
                          method = "knn", 
                          metric = "Accuracy",
                          tuneGrid = tunegrid,
                          trControl = control)

tic <- Sys.time()

tic-toc #Time difference of 25.03921 secs

bestTune_list <- model_KNN$results
str(bestTune_list)
#k=3

######################################
predKNN <- predict(model_KNN)
tab_cm <- confusionMatrix(predKNN,zscore_variable_log_fpkm.var$enfermedad_primaria)
######################################
setwd(resultsDirectory)
write.csv(bestTune_list, file="bestTune_KNN_0.0001_5_top1000_zscore_variable_log_FPKM.csv")


#### 

setwd(resultsDirectory)
save(model_KNN, file = "model_KNN_variable.RData")
#Accuracy : 0.9763  
#Kappa : 0.9754      

accuracy_train <- sum(diag(tab_cm$table))/sum(tab_cm$table)
sensibilidad_train <- mean(tab_cm$byClass[,1])
especificidad_train <- mean(tab_cm$byClass[,2])
F1_score_train <- mean(tab_cm$byClass[,7])
kappa_train <-0.9754      
precision_train <- mean(tab_cm$byClass[,5])
NPV_train <- mean(tab_cm$byClass[,4])



######################################
########## KNN VALIDATION
######################################

setwd(tablesDirectory)
zscore_variable_log_FPKM_test <- read.csv(file="zscore_variable_log_FPKM_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_variable_log_FPKM_test) <- zscore_variable_log_FPKM_test$X
zscore_variable_log_FPKM_test <- zscore_variable_log_FPKM_test[, colnames(zscore_variable_log_FPKM_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]
phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)

zscore_variable_log_FPKM_test <- zscore_variable_log_FPKM_test[rownames(zscore_variable_log_FPKM_test) %in% colnames(zscore_variable_log_fpkm.var),]
zscore_variable_log_FPKM_test <- as.data.frame(t(zscore_variable_log_FPKM_test))

table(colnames(zscore_variable_log_FPKM_test)==colnames(zscore_variable_log_fpkm.var)[colnames(zscore_variable_log_fpkm.var)!="enfermedad_primaria"])
zscore_variable_log_FPKM_test <- zscore_variable_log_FPKM_test[, match(colnames(zscore_variable_log_fpkm.var)[colnames(zscore_variable_log_fpkm.var)!="enfermedad_primaria"], colnames(zscore_variable_log_FPKM_test))]

pred_val <- predict(model_KNN, zscore_variable_log_FPKM_test)
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)

#Accuracy : 0.8783      
#Kappa : 0.8736  

accuracy_val <- sum(diag(tab_cm_VAL$table))/sum(tab_cm_VAL$table)
sensibilidad_val <- mean(tab_cm_VAL$byClass[,1])
especificidad_val <- mean(tab_cm_VAL$byClass[,2])
F1_score_val <- mean(tab_cm_VAL$byClass[,7])
kappa_val <- 0.8736  
precision_val <- mean(tab_cm_VAL$byClass[,5])
NPV_val <- mean(tab_cm_VAL$byClass[,4])


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

zscore_variable_log_FPKM_MET500 <- zscore_variable_log_FPKM_MET500[rownames(zscore_variable_log_FPKM_MET500) %in% colnames(zscore_variable_log_fpkm.var),]
zscore_variable_log_FPKM_MET500 <- as.data.frame(t(zscore_variable_log_FPKM_MET500))

table(colnames(zscore_variable_log_FPKM_MET500)==colnames(zscore_variable_log_fpkm.var)[colnames(zscore_variable_log_fpkm.var)!="enfermedad_primaria"])
zscore_variable_log_FPKM_MET500 <- zscore_variable_log_FPKM_MET500[, match(colnames(zscore_variable_log_fpkm.var)[colnames(zscore_variable_log_fpkm.var)!="enfermedad_primaria"], colnames(zscore_variable_log_FPKM_MET500))]


predRF_MET500 <- predict(model_KNN, zscore_variable_log_FPKM_MET500)
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)

#Accuracy : 0.4729 
#Kappa :0.4129 

tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)
tejidos_MET500 <- paste("Class:",tejidos_MET500)

accuracy_MET500 <- sum(diag(tab_cm_MET500$table))/sum(tab_cm_MET500$table)
sensibilidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,1])
especificidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,2])
F1_score_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,7])
kappa_MET500 <- 0.4129 
precision_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,5])
NPV_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,4])

######################################
setwd(resultsDirectory)
saveRDS(tab_cm_MET500, file="tab_cm_MET500_KNN_top1000_zscore_variable_log_FPKM.rds")


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
write.csv(df_byClass_training, file="df_byClass_training_kNN_variable.csv")
write.csv(df_byClass_VAL, file="df_byClass_VAL_kNN_variable.csv")
write.csv(df_byClass_MET500, file="df_byClass_MET500_kNN_variable.csv")


