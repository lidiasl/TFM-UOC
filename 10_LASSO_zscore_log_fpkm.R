#####################################################################
# Select genes applying Lasso ZSCORE LOG FPKM
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/LASSO"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Read training data
######################################
setwd(tablesDirectory)
# Filtered
zscore_log_fpkm <- read.csv(file="zscore_log_fpkm.csv")
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")

genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
df_tissue_biopsy_tissue <- read.csv(file="df_tissue_biopsy_tissue.csv", header=TRUE)
######################################

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_log_fpkm) <- zscore_log_fpkm$X
zscore_log_fpkm <- zscore_log_fpkm[, colnames(zscore_log_fpkm)!="X"]

######################################

genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

genes_list_0.0001_5 <- Reduce(c,genes_list_0.0001_5) 
genes_list_0.0001_5 <- unique(genes_list_0.0001_5) #7339
length(genes_list_0.0001_5) #7339

zscore_log_fpkm <- zscore_log_fpkm[rownames(zscore_log_fpkm) %in% genes_list_0.0001_5,]
zscore_log_fpkm <- t(zscore_log_fpkm)
zscore_log_fpkm <- as.data.frame(zscore_log_fpkm)

table(rownames(zscore_log_fpkm) == phenoData_training_final$sample)

zscore_log_fpkm$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_log_fpkm), phenoData_training_final$sample)]
zscore_log_fpkm$enfermedad_primaria <- as.factor(zscore_log_fpkm$enfermedad_primaria)

######################################
####### LASSO REGRESSION
######################################

library(glmnet)
set.seed(42)
cv_lasso_fit <- cv.glmnet(
  as.matrix(zscore_log_fpkm[, colnames(zscore_log_fpkm) != "enfermedad_primaria"]), 
  zscore_log_fpkm$enfermedad_primaria,
  family = "multinomial",
  alpha = 1
)
plot(cv_lasso_fit)

lambda_min <- cv_lasso_fit$lambda.1se #0.002868027
coef_matrix <- coef(cv_lasso_fit, s = lambda_min)

#List coefs!=0
nonzero_coefs <- lapply(coef_matrix, function(x) {
  coef_df <- as.data.frame(as.matrix(x))
  nonzero <- coef_df[coef_df[, 1]!=0, , drop = FALSE]
  return(nonzero)
})

# Dataframe to add coefs!=0
nonzero_coefs_df <- data.frame()

# Loop to add to coef_dataframe coefs!=0
for (i in 1:length(nonzero_coefs)) {
  nonzero_df <- as.data.frame(nonzero_coefs[[i]])
  nonzero_df$Variable <- rownames(nonzero_df)
  nonzero_coefs_df <- rbind(nonzero_coefs_df, nonzero_df)
}

nonzero_coefs_df <- nonzero_coefs_df[nonzero_coefs_df$Variable!="(Intercept)",]

nrow(nonzero_coefs_df) #319

######################################
#### Save genes LASSO
######################################
setwd(tablesDirectory)
write.csv(nonzero_coefs_df, file="nonzero_genes_LASSO_zscore_log_FPKM.csv")
saveRDS(coef_matrix, file="coef_matrix_LASSO_zscore_log_FPKM.rds")


######################################
########## LASSO VALIDATION COHORT
######################################
genes_LASSO <- nonzero_coefs_df$Variable
genes_LASSO <- unique(genes_LASSO) #308

setwd(tablesDirectory)
zscore_log_FPKM_test <- read.csv(file="zscore_log_fpkm_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_log_FPKM_test) <- zscore_log_FPKM_test$X
zscore_log_FPKM_test <- zscore_log_FPKM_test[, colnames(zscore_log_FPKM_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]
phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)

zscore_log_FPKM_test <- zscore_log_FPKM_test[rownames(zscore_log_FPKM_test) %in% genes_list_0.0001_5,]
zscore_log_FPKM_test <- as.data.frame(t(zscore_log_FPKM_test))

pred_LASSO_val <- predict(cv_lasso_fit, newx = as.matrix(zscore_log_FPKM_test), type = "class", s=lambda_min)
tab_cm_VAL <- confusionMatrix(as.factor(as.vector(pred_LASSO_val)), phenodata_test$enfermedad_primaria)

#Accuracy : 0.9312
#Kappa : 0.9286
######################################
setwd(resultsDirectory)
saveRDS(tab_cm_VAL, file="tab_cm_VAL_LASSO_zscore_log_FPKM.rds")


######################################
########## LASSO MET500 COHORT
######################################

setwd(tablesDirectory)
zscore_log_fpkm_MET500_clean <- read.csv(file="zscore_log_fpkm_MET500_clean.csv")
phenodata_MET500_clean <- read.csv(file="phenodata_MET500_clean.csv")

rownames(zscore_log_fpkm_MET500_clean) <- zscore_log_fpkm_MET500_clean$X
zscore_log_fpkm_MET500_clean <- zscore_log_fpkm_MET500_clean[, colnames(zscore_log_fpkm_MET500_clean)!="X"]
colnames(zscore_log_fpkm_MET500_clean) <- str_replace_all(colnames(zscore_log_fpkm_MET500_clean), "[.]", "-")
table(colnames(zscore_log_fpkm_MET500_clean) %in% phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

zscore_log_fpkm_MET500_clean <- zscore_log_fpkm_MET500_clean[, match(phenodata_MET500_clean$Sample_id, colnames(zscore_log_fpkm_MET500_clean))]
table(colnames(zscore_log_fpkm_MET500_clean)==phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

phenodata_MET500_clean$enfermedad_primaria <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$tissue_set, df_tissue_biopsy_tissue$set)]
phenodata_MET500_clean$enfermedad_primaria <- as.factor(phenodata_MET500_clean$enfermedad_primaria)

zscore_log_fpkm_MET500_clean <- zscore_log_fpkm_MET500_clean[rownames(zscore_log_fpkm_MET500_clean) %in% genes_list_0.0001_5,]
zscore_log_fpkm_MET500_clean <- as.data.frame(t(zscore_log_fpkm_MET500_clean))



pred_LASSO <- predict(cv_lasso_fit, newx = as.matrix(zscore_log_fpkm_MET500_clean), type = "class", s=lambda_min)
pred_LASSO <- as.factor(as.vector(pred_LASSO))

pred_LASSO <- factor(pred_LASSO, levels=union(levels(phenodata_MET500_clean$enfermedad_primaria), levels(pred_LASSO)))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, 
                                                     levels=union(levels(phenodata_MET500_clean$enfermedad_primaria), levels(pred_LASSO)))

tab_cm_MET500 <- confusionMatrix(pred_LASSO, phenodata_MET500_clean$enfermedad_primaria)
#Accuracy : 0.0258
#Kappa : 0.0034
######################################
setwd(resultsDirectory)
saveRDS(tab_cm_MET500, file="tab_cm_MET500_LASSO_zscore_log_FPKM.rds")


