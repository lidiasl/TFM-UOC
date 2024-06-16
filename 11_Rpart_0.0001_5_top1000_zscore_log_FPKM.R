#####################################################################
# Genes of pvalue adj < 0.0001 & log2FC>3 & top 1000 in Decision tree
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
#### DECISION TREE
######################################

zscore_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_log_fpkm.var$enfermedad_primaria)

####
tunegrid <- expand.grid(cp=seq(0,1,0.05))
maxdepth_list <- seq(5,30,5)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)

bestTune_list<- list()

library(doParallel)
registerDoParallel(cores = 16)
set.seed(42) 
toc <- Sys.time()
for (i in 1:length(maxdepth_list)){
  model_rf <- caret::train(enfermedad_primaria ~ .,
                           data = zscore_log_fpkm.var,
                           method = "rpart", # this will use the randomForest::randomForest function
                           metric = "Accuracy", # which metric should be optimized for classification
                           # options to be passed to randomForest
                           tuneGrid = tunegrid,
                           trControl = control,
                           maxdepth = maxdepth_list[i])
  bestTune_list[[i]] <- model_rf[["results"]]
  print(i)
}


tic <- Sys.time()

tic-toc #Time difference of 1.68827 mins


for (i in 1:length(maxdepth_list)){
  bestTune_list[[i]]$maxdepth <- maxdepth_list[i]
}

bestTune_df <- Reduce(rbind, bestTune_list)
bestTune_df$maxdepth <- as.factor(bestTune_df$maxdepth)
ggplot(bestTune_df, aes(x=cp, y=Kappa, col=maxdepth)) + geom_line()
ggplot(bestTune_df, aes(x=cp, y=Accuracy, col=maxdepth)) + geom_line()
######################################
setwd(resultsDirectory)
write.csv(bestTune_df, file="bestTune_RPART_0.0001_5_top1000_zscore_log_FPKM.csv")


#### 

zscore_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_log_fpkm.var$enfermedad_primaria)

library(rpart)
set.seed(42)
model_rpart <- rpart(enfermedad_primaria ~., data = zscore_log_fpkm.var,
                         cp=0,
                         maxdepth=20)


predrpart <- predict(model_rpart, type="class")
tab_cm <- confusionMatrix(predrpart,zscore_log_fpkm.var$enfermedad_primaria)

setwd(resultsDirectory)
save(model_rpart,file = "model_RPART.RData")

#Accuracy : 0.8066  
#Kappa : 0.7992   

accuracy_train <- sum(diag(tab_cm$table))/sum(tab_cm$table)
sensibilidad_train <- mean(tab_cm$byClass[,1])
especificidad_train <- mean(tab_cm$byClass[,2])
F1_score_train <- mean(tab_cm$byClass[,7])
kappa_train <- 0.7992   
precision_train <- mean(tab_cm$byClass[,5])
NPV_train <- mean(tab_cm$byClass[,4])


######################################
########## PROBAR RPART SOBRE VALIDATION
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

pred_val <- predict(model_rpart, zscore_log_FPKM_test, type="class")
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)

#Accuracy : 0.6931
#Kappa : 0.6813    

accuracy_val <- sum(diag(tab_cm_VAL$table))/sum(tab_cm_VAL$table)
sensibilidad_val <- mean(tab_cm_VAL$byClass[,1])
especificidad_val <- mean(tab_cm_VAL$byClass[,2])
F1_score_val <- mean(tab_cm_VAL$byClass[,7])
kappa_val <- 0.6813
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


predRF_MET500 <- predict(model_rpart, zscore_log_FPKM_MET500, type="class")
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)

#Accuracy : 0.3514 
#Kappa : 0.297 

tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)
tejidos_MET500 <- paste("Class:",tejidos_MET500)

accuracy_MET500 <- sum(diag(tab_cm_MET500$table))/sum(tab_cm_MET500$table)
sensibilidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,1])
especificidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,2])
F1_score_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,7])
kappa_MET500 <- 0.297   
precision_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,5])
NPV_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,4])

######################################
setwd(resultsDirectory)
saveRDS(tab_cm_MET500, file="tab_cm_MET500_RPART_top1000_zscore_log_FPKM.rds")

######################################
#### STUDY MODEL
######################################
importance <- varImp(model_rpart)
setwd(resultsDirectory)
write.csv(importance, "importance_RPART.csv")

model_rpart$variable.importance

library(rpart.plot)
rpart.plot(model_rpart)


pdf("arbol_decision.pdf", height=10, width=8)
prp(model_rpart)
dev.off()

importance_values <- as.data.frame(model_rpart$variable.importance)
importance_values$Variable <- rownames(importance_values)
colnames(importance_values) <- c("Importance", "Variable")
top_30_importance <- importance_values[order(importance_values$Importance, decreasing = TRUE), ][1:30, ]

ggplot(top_30_importance, aes(x = Importance, y = reorder(Variable, Importance))) +
  geom_point(size = 3, color = "black") +
  labs(title = "Importancia de las variables (Árbol de decisión)",
       x = "Importancia",
       y = "Variables") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

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
write.csv(df_byClass_training, file="df_byClass_training_rpart.csv")
write.csv(df_byClass_VAL, file="df_byClass_VAL_rpart.csv")
write.csv(df_byClass_MET500, file="df_byClass_MET500_rpart.csv")
