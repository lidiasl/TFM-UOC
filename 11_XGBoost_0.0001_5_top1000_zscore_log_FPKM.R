#####################################################################
# Genes of pvalue adj < 0.0001 & log2FC>5 & top 1000 in XGBoost
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
#### XGBoost
######################################

zscore_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_log_fpkm.var$enfermedad_primaria)

####
tunegrid <- expand.grid(
  nrounds = seq(500,1500,100),
  eta = c(0.001,0.01,0.05,0.1),
  max_depth = seq(5,30,5),
  gamma = 0,
  colsample_bytree = seq(0.1,0.5,0.1),
  min_child_weight = 1,
  subsample = 1
)


control <- trainControl(method='cv', 
                        number=3, 
                        allowParallel = TRUE)



library(doParallel)
registerDoParallel(cores = 16)
set.seed(42) 
toc <- Sys.time()
model_xgboost <- caret::train(enfermedad_primaria ~ .,
                          data = zscore_log_fpkm.var,
                          method = "xgbTree", # this will use the svmRadial function
                          metric = "Accuracy",
                          trControl = control)

tic <- Sys.time()

tic-toc #
