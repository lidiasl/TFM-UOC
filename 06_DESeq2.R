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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Read training data
######################################
setwd(tablesDirectory)
rnaseq_counts_training <- read.csv(file = "rnaseq_counts_training.csv")
# Filtered
rnaseq_counts_training.f <- read.csv(file = "rnaseq_counts_training_f.csv")
# phenoData
phenoData_training <- read.csv(file = "phenoData_training.csv")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(rnaseq_counts_training) <- rnaseq_counts_training$X
rnaseq_counts_training <- rnaseq_counts_training[, colnames(rnaseq_counts_training)!="X"]

rownames(rnaseq_counts_training.f) <- rnaseq_counts_training.f$X
rnaseq_counts_training.f <- rnaseq_counts_training.f[, colnames(rnaseq_counts_training.f)!="X"]

######################################
#CHECK:
table(colnames(rnaseq_counts_training.f)==phenoData_training$sample)

phenoData_training$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)

dds.train <- DESeqDataSetFromMatrix(countData = rnaseq_counts_training.f, 
                                    colData = phenoData_training, 
                                    design = ~ 0 + primary_disease_set)

#"2024-04-18 16:31:58 CEST"
print(Sys.time())
dds.train <- DESeq(dds.train)
print(Sys.time())
#"2024-04-18 17:02:19 CEST"
#replacing outliers and refitting for 11 genes

######################################
setwd(tablesDirectory)
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


names_primary_disease_set <- unique(phenoData_training$primary_disease_set)

genes_list <- vector("list", length = length(names_primary_disease_set))
names(genes_list) <- names_primary_disease_set

print(Sys.time()) # "2024-04-19 09:46:29 CEST"
for (i in 1:length(names_primary_disease_set)){
  for (j in 1:length(names_primary_disease_set)){
    if (names_primary_disease_set[i]!=names_primary_disease_set[j]){
      #results(dds, contrast=c("condition","treated","untreated"))
      res <- results(dds.train, contrast = c("primary_disease_set", 
                                             as.character(names_primary_disease_set[i]), 
                                             as.character(names_primary_disease_set[j])))
      res <- as.data.frame(res)
      genes_list[[i]][[j]] <- res
    }
  }
  names(genes_list[[i]]) <- names_primary_disease_set
  print(i)
}
print(Sys.time()) # Aprox 13 horas y media

genes_list[[27]][27] <- list(NULL)
names(genes_list[[27]]) <- names_primary_disease_set

######################################
setwd(tablesDirectory)
saveRDS(genes_list, file="genes_list.rds")
######################################

# Study genes intersect in each class

for(i in 1:length(genes_list)){
  genes_list[[i]] <- genes_list[[i]][-i]
}

genes_list_df <- lapply(genes_list, function(x){
  lapply(x, function(y){
    rownames(y[y$padj <0.001 & y$log2FoldChange > 0,])
  })
})

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
  ggtitle("Número de genes por clase (padj <0.001 & log2FoldChange > 0)")+
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

######################################
# Volcano plot de ejemplo

with(genes_list$adrenocortical$brain, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(genes_list$adrenocortical$brain, padj<.001 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(genes_list$adrenocortical$brain, padj<.001 & log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))