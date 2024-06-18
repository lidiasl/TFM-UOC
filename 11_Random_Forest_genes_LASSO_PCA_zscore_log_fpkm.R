#####################################################################
# RANDOM FOREST genes LASSO PCA zscore log FPKM
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Random Forest Prueba"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"


######################################
# Read data
######################################

setwd(tablesDirectory)
zscore_log_fpkm <- read.csv(file="zscore_log_fpkm.csv")
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")

genes_LASSO <- read.csv(file="nonzero_genes_LASSO_zscore_log_FPKM.csv")


setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
df_tissue_biopsy_tissue <- read.csv(file = "df_tissue_biopsy_tissue.csv", header = TRUE)
######################################

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_log_fpkm) <- zscore_log_fpkm$X
zscore_log_fpkm <- zscore_log_fpkm[, colnames(zscore_log_fpkm)!="X"]
######################################

genes_LASSO <- genes_LASSO$Variable
genes_LASSO <- unique(genes_LASSO) #308


######################################
zscore_log_fpkm.f <- zscore_log_fpkm[rownames(zscore_log_fpkm) %in% genes_LASSO,]
zscore_log_fpkm.f <- zscore_log_fpkm.f[,colnames(zscore_log_fpkm.f) %in% phenoData_training_final$sample,]
zscore_log_fpkm.f <- t(zscore_log_fpkm.f)
zscore_log_fpkm.f <- as.data.frame(zscore_log_fpkm.f)

table(rownames(zscore_log_fpkm.f)==phenoData_training_final$sample)
zscore_log_fpkm.f$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_log_fpkm.f), phenoData_training_final$sample)]


######################################
####### PCA
######################################

set.seed(42)
PCA <- prcomp(as.matrix(zscore_log_fpkm.f[,colnames(zscore_log_fpkm.f)!="enfermedad_primaria"]), 
              scale. = FALSE, 
              center = FALSE)

library(factoextra)
fviz_eig(PCA, addlabels = TRUE, labelsize=0.1,ncp=30)+ 
  geom_hline(yintercept = 0.1, color = "red") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
transformacion_PCA <- PCA$rotation

PCA_training <- PCA$x
table(rownames(PCA_training) == phenoData_training_final$sample)

PCA_training <- as.data.frame(PCA_training)

PCA_training$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(PCA_training), phenoData_training_final$sample)]
PCA_training$enfermedad_primaria <- as.factor(PCA_training$enfermedad_primaria)

####
tunegrid <- expand.grid(mtry=seq(10,50,5))
ntree_list <- seq(500,2000,100)

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

tic-toc #Time difference of 1.616388 hours

for (i in 1:length(ntree_list)){
  bestTune_list[[i]]$ntree <- ntree_list[i]
}

bestTune_df <- Reduce(rbind, bestTune_list)
bestTune_df$ntree <- as.factor(bestTune_df$ntree)
ggplot(bestTune_df, aes(x=mtry, y=Kappa, col=ntree)) + geom_line()
ggplot(bestTune_df, aes(x=mtry, y=Accuracy, col=ntree)) + geom_line()
######################################
setwd("~/Desktop/LIDIA/TCGA_classification/TFM/Results/Random Forest Prueba")
write.csv(bestTune_df, file="bestTune_df_LASSO_PCA_zscore_log_fpkm.csv")

#### 
library(randomForest)
set.seed(42)
model_rf <- randomForest(x = as.matrix(PCA_training[, colnames(PCA_training)!="enfermedad_primaria"]),
                         y=PCA_training$enfermedad_primaria,
                         mtry=10,
                         ntree=1500,
                         keep.forest = TRUE,
                         importance = TRUE)
predRF <- predict(model_rf)
tab_cm <- confusionMatrix(predRF,PCA_training$enfermedad_primaria)
#Accuracy : 0.9711
#Kappa : 0.9699

accuracy_train <- sum(diag(tab_cm$table))/sum(tab_cm$table)
sensibilidad_train <- mean(tab_cm$byClass[,1])
especificidad_train <- mean(tab_cm$byClass[,2])
F1_score_train <- mean(tab_cm$byClass[,7])
kappa_train <- 0.9699    
precision_train <- mean(tab_cm$byClass[,5])


######################################
########## RF VALIDATION
######################################

setwd(tablesDirectory)
zscore_log_FPKM_test <- read.csv(file="zscore_log_fpkm_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_log_FPKM_test) <- zscore_log_FPKM_test$X
zscore_log_FPKM_test <- zscore_log_FPKM_test[, colnames(zscore_log_FPKM_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]


zscore_log_FPKM_test <- zscore_log_FPKM_test[rownames(zscore_log_FPKM_test) %in% genes_LASSO,]
zscore_log_FPKM_test <- as.data.frame(t(zscore_log_FPKM_test))

table(colnames(zscore_log_FPKM_test)==rownames(transformacion_PCA))


zscore_log_FPKM_test_PCA <- as.matrix(zscore_log_FPKM_test) %*% transformacion_PCA
zscore_log_FPKM_test_PCA <- as.data.frame(zscore_log_FPKM_test_PCA)

phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)


pred_val <- predict(model_rf, zscore_log_FPKM_test_PCA)
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)
#Accuracy : 0.9259
#Kappa : 0.9231


accuracy_val <- sum(diag(tab_cm_VAL$table))/sum(tab_cm_VAL$table)
sensibilidad_val <- mean(tab_cm_VAL$byClass[,1])
especificidad_val <- mean(tab_cm_VAL$byClass[,2])
F1_score_val <- mean(tab_cm_VAL$byClass[,7])
kappa_val <- 0.9231      
precision_val <- mean(tab_cm_VAL$byClass[,5])

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

zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[rownames(zscore_log_FPKM_MET500) %in% genes_LASSO,]
zscore_log_FPKM_MET500 <- as.data.frame(t(zscore_log_FPKM_MET500))

table(colnames(zscore_log_FPKM_MET500)==rownames(transformacion_PCA))
table(colnames(zscore_log_FPKM_MET500) %in% rownames(transformacion_PCA))
zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, match(rownames(transformacion_PCA), colnames(zscore_log_FPKM_MET500))]

zscore_log_FPKM_MET500_PCA <- as.matrix(zscore_log_FPKM_MET500) %*% transformacion_PCA
zscore_log_FPKM_MET500_PCA <- as.data.frame(zscore_log_FPKM_MET500_PCA)


predRF_MET500 <- predict(model_rf, zscore_log_FPKM_MET500_PCA)
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)
#Accuracy : 0.4574
#Kappa : 0.406

tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)
tejidos_MET500 <- paste("Class:",tejidos_MET500)

accuracy_MET500 <- sum(diag(tab_cm_MET500$table))/sum(tab_cm_MET500$table)
sensibilidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,1])
especificidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,2])
F1_score_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,7])
kappa_MET500 <- 0.406 
precision_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,5])
