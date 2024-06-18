#####################################################################
# Genes of pvalue adj < 0.0001 & log2FC>5 & top 50 in each tumor
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Random Forest Prueba"
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

genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")
setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file = "df_tissue_biopsy_tissue.csv", header = TRUE)

######################################


genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

names_enfermedad_primaria <- names(genes_list_0.0001_5)

list_top50 <- list()

zscore_log_fpkm <- zscore_log_fpkm[,colnames(zscore_log_fpkm) %in% phenoData_training_final$sample,]

for (i in 1:length(names_enfermedad_primaria)){
  zscore_log_fpkm.f <- zscore_log_fpkm[rownames(zscore_log_fpkm) %in% genes_list_0.0001_5[[i]],]
  zscore_log_fpkm.var <- zscore_log_fpkm.f[order(rowVars(as.matrix(zscore_log_fpkm.f)), decreasing = TRUE)[1:50],]
  list_top50[[i]] <- rownames(zscore_log_fpkm.var)
}

names(list_top50) <- names_enfermedad_primaria

list_top50 <- Reduce(c,list_top50)
list_top50 <- unique(list_top50) #139


zscore_log_fpkm.var <- zscore_log_fpkm[rownames(zscore_log_fpkm)%in%list_top50,]

zscore_log_fpkm.var <- t(zscore_log_fpkm.var)
zscore_log_fpkm.var <- as.data.frame(zscore_log_fpkm.var)

table(rownames(zscore_log_fpkm.var)==phenoData_training_final$sample)
zscore_log_fpkm.var$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_log_fpkm.var), phenoData_training_final$sample)]


######################################
#### RANDOM FOREST
######################################
zscore_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_log_fpkm.var$enfermedad_primaria)

####
tunegrid <- expand.grid(mtry=seq(5,50,5))
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
                           data = zscore_log_fpkm.var,
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

tic-toc #Time difference of 36.03792 mins

for (i in 1:length(ntree_list)){
  bestTune_list[[i]]$ntree <- ntree_list[i]
}

bestTune_df <- Reduce(rbind, bestTune_list)
bestTune_df$ntree <- as.factor(bestTune_df$ntree)
ggplot(bestTune_df, aes(x=mtry, y=Kappa, col=ntree)) + geom_line()
ggplot(bestTune_df, aes(x=mtry, y=Accuracy, col=ntree)) + geom_line()

######################################
setwd(resultsDirectory)
write.csv(bestTune_df, file="bestTune_df_0.0001_5_top50_zscore_log_FPKM.csv")


#### 

zscore_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_log_fpkm.var$enfermedad_primaria)

library(randomForest)
set.seed(42)
model_rf <- randomForest(x = as.matrix(zscore_log_fpkm.var[, colnames(zscore_log_fpkm.var)!="enfermedad_primaria"]),
                         y=zscore_log_fpkm.var$enfermedad_primaria,
                         mtry=10,
                         ntree=1100,
                         keep.forest = TRUE,
                         importance = TRUE)
predRF <- predict(model_rf)
tab_cm <- confusionMatrix(predRF,zscore_log_fpkm.var$enfermedad_primaria)
#Accuracy : 0.9368 
#KKappa : 0.9344      

accuracy_train <- sum(diag(tab_cm$table))/sum(tab_cm$table)
sensibilidad_train <- mean(tab_cm$byClass[,1])
especificidad_train <- mean(tab_cm$byClass[,2])
F1_score_train <- mean(tab_cm$byClass[,7])
kappa_train <- 0.9344
precision_train <- mean(tab_cm$byClass[,5])



######################################
########## RF VALIDATION COHORT
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

pred_val <- predict(model_rf, zscore_log_FPKM_test)
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)

#Accuracy : 0.8995
#Kappa : 0.8956

accuracy_val <- sum(diag(tab_cm_VAL$table))/sum(tab_cm_VAL$table)
sensibilidad_val <- mean(tab_cm_VAL$byClass[,1])
especificidad_val <- mean(tab_cm_VAL$byClass[,2])
F1_score_val <- mean(tab_cm_VAL$byClass[,7])
kappa_val <- 0.8956 
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

zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[rownames(zscore_log_FPKM_MET500) %in% colnames(zscore_log_fpkm.var),]
zscore_log_FPKM_MET500 <- as.data.frame(t(zscore_log_FPKM_MET500))

table(colnames(zscore_log_FPKM_MET500)==colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"])
zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, match(colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"], colnames(zscore_log_FPKM_MET500))]


predRF_MET500 <- predict(model_rf, zscore_log_FPKM_MET500)
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)

#Accuracy : 0.4677
#Kappa : 0.3971 

tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)
tejidos_MET500 <- paste("Class:",tejidos_MET500)

accuracy_MET500 <- sum(diag(tab_cm_MET500$table))/sum(tab_cm_MET500$table)
sensibilidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,1])
especificidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,2])
F1_score_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,7])
kappa_MET500 <- 0.3971    
precision_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,5])

