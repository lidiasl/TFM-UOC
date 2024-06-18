#####################################################################
# Transform data validation ML
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/VAL Preparation"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

#####################################################################
setwd(tablesDirectory)
rnaseq_counts_training.f_final <- read.csv(file = "rnaseq_counts_training_f_final.csv")
rnaseq_counts_test <- read.csv(file = "rnaseq_counts_test.csv")
phenoData_test <- read.csv(file = "phenoData_test.csv")
fpkm_training <- read.csv(file="fpkm_matrix.csv")

rownames(rnaseq_counts_training.f_final) <- rnaseq_counts_training.f_final$X
rnaseq_counts_training.f_final <- rnaseq_counts_training.f_final[, colnames(rnaseq_counts_training.f_final)!="X"]

rownames(rnaseq_counts_test) <- rnaseq_counts_test$X
rnaseq_counts_test <- rnaseq_counts_test[, colnames(rnaseq_counts_test)!="X"]

rownames(fpkm_training) <- fpkm_training$X
fpkm_training <- fpkm_training[, colnames(fpkm_training)!="X"]

rownames(phenoData_test) <- phenoData_test$X
phenoData_test <- phenoData_test[, colnames(phenoData_test)!="X"]

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

#####################################################################


######################################

setwd(tablesDirectory)
gtf <- read.table('gencode.v23.chr_patch_hapl_scaff.basic.annotation.gtf', header = FALSE, sep = '\t')

gtf$gene_id_mod <- sapply(gtf$V9, function(x){
  sub(".*gene_id ([^.]+).*", "\\1", x)
})

gtf.f <- gtf[gtf$V3 == "gene",]

gtf.f <- gtf.f[gtf.f$gene_id_mod %in% rownames(rnaseq_counts_test),]

featureLength <- (gtf.f$V5 - gtf.f$V4)+1

gtf.f$featureLength <- featureLength
gtf.f <- gtf.f[match(rownames(rnaseq_counts_test), gtf$gene_id_mod),]

#install.packages("DGEobj.utils")
library(DGEobj.utils)

fpkm_test <- convertCounts(as.matrix(rnaseq_counts_test),
                             unit       = "fpkm",
                             geneLength = featureLength,
                             log        = FALSE,
                             normalize  = "none")

# Me quedo con los genes de entrenamiento
fpkm_test <- fpkm_test[rownames(fpkm_test) %in% rownames(fpkm_training),]

table(apply(fpkm_test,1,function(x){sd(x)!=0}))
#FALSE  TRUE 
#48 18942


sum_fpkm_per_sample <- colSums(fpkm_test)
scaling_factors <- sum_fpkm_per_sample / 1e6
tpm_test <- t(t(fpkm_test) / scaling_factors)

log_fpkm_test <- log2(fpkm_test+0.001)
log_tpm_test <- log2(tpm_test+0.001)

#####################################################################

setwd(tablesDirectory)
medias_genes_log_fpkm <- readRDS(file="medias_genes_log_fpkm.rds")
sd_genes_log_fpkm <- readRDS(file="sd_genes_log_fpkm.rds")
medias_genes_log_tpm <- readRDS(file = "medias_genes_log_tpm.rds")
sd_genes_log_tpm <- readRDS(file = "sd_genes_log_tpm.rds")
#####################################################################

# zscore variable log

zscore_variable_log_FPKM_test <- apply(log_fpkm_test,2,function(x){
  (x-medias_genes_log_fpkm)/sd_genes_log_fpkm
})

zscore_variable_log_tpm_test <- apply(log_tpm_test,2,function(x){
  (x-medias_genes_log_tpm)/sd_genes_log_tpm
})


# zscore

zscore_fpkm_test <- scale(fpkm_test)
zscore_tpm_test <- scale(tpm_test)


# zscore log

zscore_log_fpkm_test <- scale(log_fpkm_test)
zscore_log_tpm_test <- scale(log_tpm_test)


######################################
setwd(tablesDirectory)
write.csv(fpkm_test, file="fpkm_test.csv")
write.csv(log_fpkm_test, file="log_fpkm_test.csv")
write.csv(zscore_log_fpkm_test, file = "zscore_log_fpkm_test.csv")
write.csv(zscore_variable_log_FPKM_test, file = "zscore_variable_log_FPKM_test.csv")
write.csv(zscore_fpkm_test, file="zscore_fpkm_test.csv")

write.csv(tpm_test, file="tpm_test.csv")
write.csv(log_tpm_test, file="log_tpm_test.csv")
write.csv(zscore_log_tpm_test, file = "zscore_log_tpm_test.csv")
write.csv(zscore_variable_log_tpm_test, file = "zscore_variable_log_tpm_test.csv")
write.csv(zscore_tpm_test, file="zscore_tpm_test.csv")
######################################


######################################
# Filter genes union final 
setwd(tablesDirectory)
genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")

genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

genes_list_0.0001_5 <- Reduce(c,genes_list_0.0001_5) 
genes_union <- unique(genes_list_0.0001_5) #7339



fpkm_test.genes_final <- fpkm_test[rownames(fpkm_test) %in% genes_union,
                                                   match(phenoData_test$sample, colnames(fpkm_test))]
log_fpkm_test.genes_final <- log_fpkm_test[rownames(log_fpkm_test) %in% genes_union,
                                                           match(phenoData_test$sample, colnames(log_fpkm_test))]
zscore_log_fpkm_test.genes_final <- zscore_log_fpkm_test[rownames(zscore_log_fpkm_test) %in% genes_union,
                                                                         match(phenoData_test$sample, colnames(zscore_log_fpkm_test))]
zscore_variable_log_FPKM_test.genes_final <- zscore_variable_log_FPKM_test[rownames(zscore_variable_log_FPKM_test) %in% genes_union,
                                                                                           match(phenoData_test$sample, colnames(zscore_variable_log_FPKM_test))]
zscore_fpkm_test.genes_final <- zscore_fpkm_test[rownames(zscore_fpkm_test) %in% genes_union,
                                                                 match(phenoData_test$sample, colnames(zscore_fpkm_test))]

tpm_test.genes_final <- tpm_test[rownames(tpm_test) %in% genes_union,
                                                 match(phenoData_test$sample, colnames(tpm_test))]
log_tpm_test.genes_final <- log_tpm_test[rownames(log_tpm_test) %in% genes_union,
                                                         match(phenoData_test$sample, colnames(log_tpm_test))]
zscore_log_tpm_test.genes_final <- zscore_log_tpm_test[rownames(zscore_log_tpm_test) %in% genes_union,
                                                                       match(phenoData_test$sample, colnames(zscore_log_tpm_test))]
zscore_variable_log_tpm_test.genes_final <- zscore_variable_log_tpm_test[rownames(zscore_variable_log_tpm_test) %in% genes_union,
                                                                                         match(phenoData_test$sample, colnames(zscore_variable_log_tpm_test))]
zscore_tpm_test.genes_final <- zscore_tpm_test[rownames(zscore_tpm_test) %in% genes_union,
                                                               match(phenoData_test$sample, colnames(zscore_tpm_test))]

######################################
############     PCA      ############
######################################

list_data <- list(fpkm_test.genes_final,
                  log_fpkm_test.genes_final,
                  zscore_log_fpkm_test.genes_final,
                  zscore_variable_log_FPKM_test.genes_final,
                  zscore_fpkm_test.genes_final,
                  tpm_test.genes_final,
                  log_tpm_test.genes_final,
                  zscore_log_tpm_test.genes_final,
                  zscore_variable_log_tpm_test.genes_final,
                  zscore_tpm_test.genes_final)

names(list_data) <- list("fpkm",
                         "log_fpkm",
                         "zscore_log_fpkm",
                         "zscore_variable_log_fpkm",
                         "zscore_fpkm",
                         "tpm",
                         "log_tpm",
                         "zscore_log_tpm",
                         "zscore_variable_log_tpm",
                         "zscore_tpm")
# Cargar la paleta de colores turbo
library(viridisLite)
colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

set.seed(42)
PCA_list <- list()

for (i in 1:length(list_data)){
  PCA <- prcomp(t(as.matrix(list_data[[i]])), 
                scale. = FALSE, 
                center = FALSE)
  autplot_PCA <- autoplot(PCA, data = phenoData_test, colour = "enfermedad_primaria",
                          main = as.character(names(list_data)[i]), size=0.4)+
    labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")+
    scale_color_manual(values = colors_blind)
  PCA_list[[i]] <- autplot_PCA
}


setwd(resultsDirectory)
pdf(file = "PCA_data_ML_test_TCGA.pdf", width = 14, height = 6)
grid.arrange(PCA_list[[1]], PCA_list[[2]], PCA_list[[3]], PCA_list[[4]],PCA_list[[5]],
             PCA_list[[6]], PCA_list[[7]], PCA_list[[8]], PCA_list[[9]],PCA_list[[10]],
             ncol=5)
dev.off()


######################################
############     UMAP     ############
######################################
umap.defaults
#n_neighbors: 15
#min_dist: 0.1

# Bucle listas
n_neighbors_list <- seq(5,100,10)
min_dist_list <- seq(0.1,0.5,0.1)

custom.settings = umap.defaults
custom.settings$n_neighbors = n_neighbors_list[3]
custom.settings$min_dist=min_dist_list[4]

set.seed(42)
UMAP_list <- list()

for (i in 1:length(list_data)){
  umap_i <- umap(t(list_data[[i]]))
  
  df.umap <- data.frame(x = umap_i$layout[,1],
                        y = umap_i$layout[,2],
                        Group = phenoData_test$enfermedad_primaria)
  
  UMAP_list[[i]] <- ggplot(df.umap, aes(x, y, colour = Group)) +
    geom_point(size=0.3) + ggtitle(as.character(names(list_data)[i]))+
    scale_color_manual(values = colors_blind)+
    labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")
}


setwd(resultsDirectory)
pdf(file = "UMAP_data_ML_test_TCGA.pdf", width = 14, height = 6)
grid.arrange(UMAP_list[[1]], UMAP_list[[2]], UMAP_list[[3]], UMAP_list[[4]], UMAP_list[[5]],
             UMAP_list[[6]], UMAP_list[[7]], UMAP_list[[8]], UMAP_list[[9]], UMAP_list[[10]],
             ncol=5)
dev.off()

