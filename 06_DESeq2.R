#####################################################################
# DESeq2
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/DESeq2"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Read training data
######################################
setwd(tablesDirectory)
# Filtered
rnaseq_counts_training.f_final <- read.csv(file = "rnaseq_counts_training_f_final.csv")
log_rnaseq_training.f_final <- read.csv(file = "log_rnaseq_training_f_final.csv")
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(rnaseq_counts_training.f_final) <- rnaseq_counts_training.f_final$X
rnaseq_counts_training.f_final <- rnaseq_counts_training.f_final[, colnames(rnaseq_counts_training.f_final)!="X"]

rownames(log_rnaseq_training.f_final) <- log_rnaseq_training.f_final$X
log_rnaseq_training.f_final <- log_rnaseq_training.f_final[, colnames(log_rnaseq_training.f_final)!="X"]

######################################
#CHECK:
table(colnames(rnaseq_counts_training.f_final)==phenoData_training_final$sample)
table(colnames(log_rnaseq_training.f_final)==phenoData_training_final$sample)

phenoData_training_final$primary_disease_set <- as.factor(phenoData_training_final$primary_disease_set)

dds.train <- DESeqDataSetFromMatrix(countData = rnaseq_counts_training.f_final, 
                                    colData = phenoData_training_final, 
                                    design = ~ 0 + primary_disease_set)

#"2024-05-09 09:04:47 CEST"
print(Sys.time())
dds.train <- DESeq(dds.train)
print(Sys.time())
#"2024-05-09 09:30:19 CEST"

#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 11 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing

######################################
setwd(resultsDirectory)
saveRDS(dds.train, file="dds.train.rds")
######################################

plotDispEsts(dds.train)

resultsNames(dds.train)
#[1] "primary_disease_setadrenocortical"     "primary_disease_setB_cell_lymphoma"    "primary_disease_setbladder"           
#[4] "primary_disease_setbrain"              "primary_disease_setbreast"             "primary_disease_setcervical"          
#[7] "primary_disease_setcholangiocarcinoma" "primary_disease_setcolon"              "primary_disease_setesophageal"        
#[10] "primary_disease_sethead_neck"          "primary_disease_setkidney"             "primary_disease_setleukemia"          
#[13] "primary_disease_setliver"              "primary_disease_setlung"               "primary_disease_setmesothelioma"      
#[16] "primary_disease_setovarian"            "primary_disease_setpancreas"           "primary_disease_setparaganglioma"     
#[19] "primary_disease_setprostate"           "primary_disease_setsarcoma"            "primary_disease_setskin"              
#[22] "primary_disease_setstomach"            "primary_disease_settesticles"          "primary_disease_setthymoma"           
#[25] "primary_disease_setthyroid"            "primary_disease_setuterine"            "primary_disease_setuveal"     



#logarithmic fold change log2(treated/untreated)
#res <- results(dds, contrast=c("condition","treated","untreated"))

names_primary_disease_set <- unique(phenoData_training_final$primary_disease_set)

genes_list <- vector("list", length = length(names_primary_disease_set))
names(genes_list) <- names_primary_disease_set

print(Sys.time()) # "2024-05-09 10:29:46 CEST"
for (i in 1:length(names_primary_disease_set)){
  for (j in 1:length(names_primary_disease_set)){
    if (names_primary_disease_set[i]!=names_primary_disease_set[j]){
      #results(dds, contrast=c("condition","treated","untreated"))
      res <- results(dds.train, 
                     contrast = c("primary_disease_set", 
                                  as.character(names_primary_disease_set[i]),
                                  as.character(names_primary_disease_set[j])),
                     pAdjustMethod = c("bonferroni"))
      res <- as.data.frame(res)
      genes_list[[i]][[j]] <- res
    }
  }
  names(genes_list[[i]]) <- names_primary_disease_set
  print(i)
}
print(Sys.time()) # Aproximadamente 13 horas

genes_list[[27]][27] <- list(NULL)
names(genes_list[[27]]) <- names_primary_disease_set



for(i in 1:length(genes_list)){
  genes_list[[i]] <- genes_list[[i]][-i]
}

######################################
setwd(tablesDirectory)
saveRDS(genes_list, file="genes_list.rds")

###
genes_list <- readRDS(file="genes_list.rds")
######################################


############################################################################
########################### Tablas padj y logFC ############################
############################################################################

# Take the union of 26 comparisons for each primary tumor

padj_list <- c(0.001, 0.005, 0.01, 0.05)
log2FC_list <- c(0,1,2,3,4,5)

padj_log2FC <- cbind(rep(padj_list,each=6), rep(log2FC_list,4))
padj_log2FC <- as.data.frame(padj_log2FC)
colnames(padj_log2FC) <- c("padj", "log2FC")

######################################
setwd(tablesDirectory)
write.csv(padj_log2FC, file="padj_log2FC.csv")
######################################

genes_list_combinations <- list()
number_genes_list_combinations <- list()

for (j in 1:nrow(padj_log2FC)){
  genes_df <- lapply(genes_list, function(x){
    lapply(x, function(y){
      rownames(y[y$padj < padj_log2FC$padj[j] & y$log2FoldChange > padj_log2FC$log2FC[j],])
    })
  })
  genes_list_combinations[[j]] <- genes_df
  
  genes_number <- lapply(genes_df, function(x){
    lapply(x, function(y){
      length(y)
    })
  })
  number_genes_list_combinations[[j]] <- genes_number
}

######################################
setwd(tablesDirectory)
saveRDS(genes_list_combinations, file="genes_list_all_combinations_padj_log2FC.rds")
######################################

##### Number of genes

tables_list <- list()

for (i in 1:length(names_primary_disease_set)){
  df <- data.frame()
  for (j in 1:nrow(padj_log2FC)){
    fila <- unlist(number_genes_list_combinations[[j]][[i]])
    df <- rbind(df, fila)
  }
  names <- paste(names_primary_disease_set[i], names_primary_disease_set[-i],sep ="_VS_")
  colnames(df) <- names
  tables_list[[i]] <- cbind(padj_log2FC,df)
}

names(tables_list) <- names_primary_disease_set



### Intersect
for (i in 1:length(names_primary_disease_set)){
  sum_vector <- c()
  for (j in 1:nrow(padj_log2FC)){
    sum_vector[j] <- length(Reduce(intersect, genes_list_df[[j]][[i]]))
  }
  tables_list[[i]]$Intersect <- sum_vector
}

### Intersect >= 24

for (i in 1:length(names_primary_disease_set)){
  sum_vector <- c()
  for (j in 1:nrow(padj_log2FC)){
    all <- unlist(genes_list_df[[j]][[i]])
    table_all <- table(all)
    sum_vector[j] <- length(names(table_all[table_all>=24]))
  }
  tables_list[[i]]$Intersect_24 <- sum_vector
}

### UNION
for (i in 1:length(names_primary_disease_set)){
  sum_vector <- c()
  for (j in 1:nrow(padj_log2FC)){
    sum_vector[j] <-  length(Reduce(unique, genes_list_df[[j]][[i]]))
  }
  tables_list[[i]]$Union <- sum_vector
}

######################################
##### SAVE TABLES
setwd(resultsDirectory)
for (i in 1:length(tables_list)){
  write.csv(tables_list[[i]], file=paste0(names_primary_disease_set[i], "_padj_log2FC.csv"))
}
######################################

