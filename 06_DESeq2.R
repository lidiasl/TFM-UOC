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

genes_list_df <- list()
genes_union_list <- list()

for (j in 1:nrow(padj_log2FC)){
  genes_df <- lapply(genes_list, function(x){
    lapply(x, function(y){
      rownames(y[y$padj <padj_log2FC$padj[j] & y$log2FoldChange > padj_log2FC$log2FC[j],])
    })
  })
  genes_list_df[[j]] <- genes_df
  
  genes_union <- lapply(genes_df, function(x){
    lapply(x, function(y){
      length(y)
    })
  })
  genes_union_list[[j]] <- genes_union
}



##### Numer of genes

tables_list <- list()

for (i in 1:length(names_primary_disease_set)){
  df <- data.frame()
  for (j in 1:nrow(padj_log2FC)){
    fila <- unlist(genes_union_list[[j]][[i]])
    df <- rbind(df, fila)
  }
  names <- paste(names_primary_disease_set[i], names_primary_disease_set[-i],sep ="_VS_")
  colnames(df) <- names
  tables_list[[i]] <- cbind(padj_log2FC,df)
}

names(tables_list) <- names_primary_disease_set


##### UPSET
library(UpSetR)

for (i in 1:length(names_primary_disease_set)){
  sum_vector <- c()
  for (j in 1:nrow(padj_log2FC)){
    upset <- upset(fromList(genes_list_df[[j]][[i]]), nsets=26)
    upset$New_data$sum <- rowSums(upset$New_data)
    sum_vector[j] <-  sum(upset$New_data$sum==26)
  }
  tables_list[[i]]$intersect <- sum_vector
}

# 24

tables_list_new <- tables_list

for (i in 1:length(names_primary_disease_set)){
  sum_vector <- c()
  for (j in 1:nrow(padj_log2FC)){
    upset <- upset(fromList(genes_list_df[[j]][[i]]), nsets=26)
    upset$New_data$sum <- rowSums(upset$New_data)
    sum_vector[j] <-  sum(upset$New_data$sum>=24)
  }
  tables_list_new[[i]]$intersect <- sum_vector
}

#### UNION
for (i in 1:length(names_primary_disease_set)){
  sum_vector <- c()
  for (j in 1:nrow(padj_log2FC)){
    union_genes <- Reduce(unique, genes_list_df[[j]][[i]])
    sum_vector[j] <-  length(union_genes)
  }
  tables_list[[i]]$union <- sum_vector
}

##### SAVE TABLES
setwd(resultsDirectory)
for (i in 1:length(tables_list)){
  write.csv(tables_list[[i]], file=paste0(names_primary_disease_set[i], "_padj_log2FC.csv"))
}


##### SAVE TABLE INTERSECT

table_intersect <- as.data.frame(tables_list[[1]]$intersect)

for (i in 2:length(tables_list)){
  table_intersect <- cbind(table_intersect, tables_list[[i]]$intersect)
}

colnames(table_intersect) <- names_primary_disease_set

table_intersect <- cbind(padj_log2FC,table_intersect)

write.csv(table_intersect, file="table_intersect.csv")

#new
table_intersect_new <- as.data.frame(tables_list_new[[1]]$intersect)

for (i in 2:length(tables_list_new)){
  table_intersect_new <- cbind(table_intersect_new, tables_list_new[[i]]$intersect)
}

colnames(table_intersect_new) <- names_primary_disease_set

table_intersect_new <- cbind(padj_log2FC,table_intersect_new)

write.csv(table_intersect_new, file="table_intersect_new.csv")


##### SAVE TABLE UNION

table_union <- as.data.frame(tables_list[[1]]$union)

for (i in 2:length(tables_list)){
  table_union <- cbind(table_union, tables_list[[i]]$union)
}

colnames(table_union) <- names_primary_disease_set

table_union <- cbind(padj_log2FC,table_union)

write.csv(table_union, file="table_union.csv")




############################################################################
# Genes in each comparison

tables_each_gene_list <- list()

for (i in 1:length(names_primary_disease_set)){
  small_list <- ()
  for (j in 1:nrow(padj_log2FC)){
    small_list <- unlist(genes_list_df[[j]][[i]])
    df <- cbind(df, df$genes %in% columna)
  }
  names <- paste(names_primary_disease_set[i], names_primary_disease_set[-i],sep ="_VS_")
  colnames(df) <- names
  tables_each_gene_list[[i]] <- cbind(padj_log2FC,df)
}

names(tables_each_gene_list) <- names_primary_disease_set


############################################################################

# Study genes intersect in each class
# padj <0.001 & log2FoldChange > 0 in each comparison
genes_list_df <- lapply(genes_list, function(x){
  lapply(x, function(y){
    rownames(y[y$padj <0.01 & y$log2FoldChange > 0,])
  })
})

# upset plot
library(UpSetR)
upset(fromList(genes_list_df$adrenocortical), nsets=26)

genes_list_df2 <- lapply(genes_list, function(x){
  lapply(x, function(y){
    rownames(y[y$padj <0.05 & y$log2FoldChange > 0,])
  })
})

genes_intersect_in_class <- lapply(genes_list_df, function(x){
  Reduce(intersect, x)
})

genes_intersect_in_class2 <- lapply(genes_list_df2, function(x){
  Reduce(intersect, x)
})

genes_union_in_class2 <- lapply(genes_list_df2, function(x){
  Reduce(unique, x)
})

genes_intersect_in_class_plot <- sapply(genes_intersect_in_class, function(x){
  list(length(x))
})
genes_intersect_in_class_plot <- as.data.frame(t(as.data.frame(genes_intersect_in_class_plot)))
colnames(genes_intersect_in_class_plot) <- "number"
genes_intersect_in_class_plot$class <- rownames(genes_intersect_in_class_plot)
genes_intersect_in_class_plot$clase <- set_subset$spanish[match(genes_intersect_in_class_plot$class,set_subset$set)]

ggplot(data=genes_intersect_in_class_plot, aes(x=clase, y=number)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=number), vjust=-0.3, size=3.5)+
  theme_minimal()+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 12, color="black"), 
        axis.text.x = element_text(size = 10, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 10, color="black"))+
  ggtitle("Número de genes por clase (padj <0.01 & log2FoldChange > 0)")+
  xlab("Tumor primario") + ylab("Número de genes")


genes_intersect_in_class_plot2 <- sapply(genes_intersect_in_class2, function(x){
  list(length(x))
})
genes_intersect_in_class_plot2 <- as.data.frame(t(as.data.frame(genes_intersect_in_class_plot2)))
colnames(genes_intersect_in_class_plot2) <- "number"
genes_intersect_in_class_plot2$class <- rownames(genes_intersect_in_class_plot2)
genes_intersect_in_class_plot2$clase <- set_subset$spanish[match(genes_intersect_in_class_plot2$class,set_subset$set)]

ggplot(data=genes_intersect_in_class_plot2, aes(x=clase, y=number)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=number), vjust=-0.3, size=3.5)+
  theme_minimal()+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 12, color="black"), 
        axis.text.x = element_text(size = 10, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 10, color="black"))+
  ggtitle("Número de genes por clase (padj <0.05 & log2FoldChange > 0)")+
  xlab("Tumor primario") + ylab("Número de genes")



genes_union_in_class_plot2 <- sapply(genes_union_in_class2, function(x){
  list(length(x))
})
genes_union_in_class_plot2 <- as.data.frame(t(as.data.frame(genes_union_in_class_plot2)))
colnames(genes_union_in_class_plot2) <- "number"
genes_union_in_class_plot2$class <- rownames(genes_union_in_class_plot2)
genes_union_in_class_plot2$clase <- set_subset$spanish[match(genes_union_in_class_plot2$class,set_subset$set)]

ggplot(data=genes_union_in_class_plot2, aes(x=clase, y=number)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=number), vjust=-0.3, size=2.5)+
  theme_minimal()+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 12, color="black"), 
        axis.text.x = element_text(size = 10, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 10, color="black"))+
  ggtitle("Número de genes por clase (padj <0.05 & log2FoldChange > 0)")+
  xlab("Tumor primario") + ylab("Número de genes")

######################################
setwd(tablesDirectory)
saveRDS(genes_union_in_class2, file="genes_union_list.rds")
######################################