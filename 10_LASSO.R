#####################################################################
# Select genes applying Lasso
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
zscore_variable_log_fpkm <- read.csv(file="zscore_variable_log_fpkm.csv")
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")

genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

######################################

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_variable_log_fpkm) <- zscore_variable_log_fpkm$X
zscore_variable_log_fpkm <- zscore_variable_log_fpkm[, colnames(zscore_variable_log_fpkm)!="X"]

######################################


genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

genes_list_0.0001_5 <- Reduce(c,genes_list_0.0001_5) 
genes_list_0.0001_5 <- unique(genes_list_0.0001_5) #7339
length(genes_list_0.0001_5) #7339

zscore_variable_log_fpkm <- zscore_variable_log_fpkm[rownames(zscore_variable_log_fpkm) %in% genes_list_0.0001_5,]
zscore_variable_log_fpkm <- t(zscore_variable_log_fpkm)
zscore_variable_log_fpkm <- as.data.frame(zscore_variable_log_fpkm)

table(rownames(zscore_variable_log_fpkm) == phenoData_training_final$sample)

zscore_variable_log_fpkm$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_variable_log_fpkm), phenoData_training_final$sample)]
zscore_variable_log_fpkm$enfermedad_primaria <- as.factor(zscore_variable_log_fpkm$enfermedad_primaria)

######################################
####### LASSO REGRESSION
######################################

library(glmnet)
set.seed(42)
cv_lasso_fit <- cv.glmnet(
  as.matrix(zscore_variable_log_fpkm[, colnames(zscore_variable_log_fpkm) != "enfermedad_primaria"]), 
  zscore_variable_log_fpkm$enfermedad_primaria,
  family = "multinomial",
  alpha = 1
)
plot(cv_lasso_fit)

lambda_min <- cv_lasso_fit$lambda.1se #0.002783792
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

nrow(nonzero_coefs_df) #316

######################################
#### Save genes LASSO
######################################
setwd(tablesDirectory)
write.csv(nonzero_coefs_df, file="nonzero_genes_LASSO.csv")
saveRDS(coef_matrix, file="coef_matrix_LASSO.rds")


######################################
########## LASSO VALIDATION COHORT
######################################
genes_LASSO <- nonzero_coefs_df$Variable
genes_LASSO <- unique(genes_LASSO) #308

setwd(tablesDirectory)
zscore_variable_log_FPKM_test <- read.csv(file="zscore_variable_log_FPKM_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_variable_log_FPKM_test) <- zscore_variable_log_FPKM_test$X
zscore_variable_log_FPKM_test <- zscore_variable_log_FPKM_test[, colnames(zscore_variable_log_FPKM_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]


zscore_variable_log_FPKM_test <- zscore_variable_log_FPKM_test[rownames(zscore_variable_log_FPKM_test) %in% genes_list_0.0001_5,]
zscore_variable_log_FPKM_test <- as.data.frame(t(zscore_variable_log_FPKM_test))

pred_LASSO_val <- predict(cv_lasso_fit, newx = as.matrix(zscore_variable_log_FPKM_test), type = "class", s=lambda_min)
table_cont_val <- table(phenodata_test$enfermedad_primaria,pred_LASSO_val)
tab_cont_perc_val <-100*(prop.table(table_cont_val,1))

library(caret)

phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)
pred_LASSO_val <- factor(pred_LASSO_val, levels = levels(phenodata_test$enfermedad_primaria))

confusionMatrix(phenodata_test$enfermedad_primaria,pred_LASSO_val)
#Accuracy : 0.9259
#Kappa : 0.9231
library(MLmetrics)
F1_Score(phenodata_test$enfermedad_primaria,pred_LASSO_val)
#0.8571429
######################################
setwd(resultsDirectory)
write.csv(tab_cont_perc_val, file="tab_cont_perc_VAL_LASSO.csv")

######################################
########## LASSO MET500 COHORT
######################################


setwd(tablesDirectory)
zscore_variable_log_FPKM_MET500 <- read.csv(file="zscore_variable_log_FPKM_MET500.csv")
phenodata_MET500_clean <- read.csv(file="phenodata_MET500_clean.csv")
setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file="df_tissue_biopsy_tissue.csv")

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

zscore_variable_log_FPKM_MET500 <- zscore_variable_log_FPKM_MET500[rownames(zscore_variable_log_FPKM_MET500) %in% genes_list_0.0001_5,]
zscore_variable_log_FPKM_MET500 <- as.data.frame(t(zscore_variable_log_FPKM_MET500))



pred_LASSO <- predict(cv_lasso_fit, newx = as.matrix(zscore_variable_log_FPKM_MET500), type = "class", s=lambda_min)
table_cont <- table(phenodata_MET500_clean$enfermedad_primaria,pred_LASSO)
tab_cont_perc <-100*(prop.table(table_cont,1))

#pred_LASSO
#esófago sarcoma  timoma 
#23     363       1 
phenodata_MET500_clean$enfermedad_primaria <- as.factor(phenodata_MET500_clean$enfermedad_primaria)
pred_LASSO <- factor(pred_LASSO, levels = levels(phenodata_MET500_clean$enfermedad_primaria))

confusionMatrix(phenodata_MET500_clean$enfermedad_primaria,pred_LASSO,mode = "everything")
#Accuracy : 0.1344
#Kappa : -0.0103

######################################
setwd(resultsDirectory)
write.csv(tab_cont_perc, file="tab_cont_perc_MET500_LASSO.csv")
