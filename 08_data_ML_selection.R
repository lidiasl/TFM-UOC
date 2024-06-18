#####################################################################
# Data selection ML
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Data ML selection"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Read training data
######################################

setwd(tablesDirectory)
fpkm_matrix <- read.csv(file="fpkm_matrix.csv")
log_fpkm_matrix <- read.csv(file="log_fpkm_matrix.csv")
zscore_log_fpkm <- read.csv(file = "zscore_log_fpkm.csv")
zscore_variable_log_fpkm <- read.csv(file = "zscore_variable_log_fpkm.csv")
zscore_fpkm <- read.csv(file="zscore_fpkm.csv")

tpm_matrix <- read.csv(file="tpm_matrix.csv")
log_tpm_matrix <- read.csv(file="log_tpm_matrix.csv")
zscore_log_tpm <- read.csv(file = "zscore_log_tpm.csv")
zscore_variable_log_tpm <- read.csv(file = "zscore_variable_log_tpm.csv")
zscore_tpm <- read.csv(file="zscore_tpm.csv")

rnaseq_counts_training.f_final <- read.csv(file = "rnaseq_counts_training_f_final.csv")
######################################
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")

genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
######################################


rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(rnaseq_counts_training.f_final) <- rnaseq_counts_training.f_final$X
rnaseq_counts_training.f_final <- rnaseq_counts_training.f_final[, colnames(rnaseq_counts_training.f_final)!="X"]

rownames(fpkm_matrix) <- fpkm_matrix$X
fpkm_matrix <- fpkm_matrix[, colnames(fpkm_matrix)!="X"]

rownames(log_fpkm_matrix) <- log_fpkm_matrix$X
log_fpkm_matrix <- log_fpkm_matrix[, colnames(log_fpkm_matrix)!="X"]

rownames(zscore_log_fpkm) <- zscore_log_fpkm$X
zscore_log_fpkm <- zscore_log_fpkm[, colnames(zscore_log_fpkm)!="X"]

rownames(zscore_variable_log_fpkm) <- zscore_variable_log_fpkm$X
zscore_variable_log_fpkm <- zscore_variable_log_fpkm[, colnames(zscore_variable_log_fpkm)!="X"]

rownames(zscore_fpkm) <- zscore_fpkm$X
zscore_fpkm <- zscore_fpkm[, colnames(zscore_fpkm)!="X"]

rownames(tpm_matrix) <- tpm_matrix$X
tpm_matrix <- tpm_matrix[, colnames(tpm_matrix)!="X"]

rownames(log_tpm_matrix) <- log_tpm_matrix$X
log_tpm_matrix <- log_tpm_matrix[, colnames(log_tpm_matrix)!="X"]

rownames(zscore_log_tpm) <- zscore_log_tpm$X
zscore_log_tpm <- zscore_log_tpm[, colnames(zscore_log_tpm)!="X"]

rownames(zscore_variable_log_tpm) <- zscore_variable_log_tpm$X
zscore_variable_log_tpm <- zscore_variable_log_tpm[, colnames(zscore_variable_log_tpm)!="X"]

rownames(zscore_tpm) <- zscore_tpm$X
zscore_tpm <- zscore_tpm[, colnames(zscore_tpm)!="X"]

######################################
#CHECK:
table(colnames(rnaseq_counts_training.f_final)==phenoData_training_final$sample)
######################################

######################################
# Filter genes union final 

genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

genes_list_0.0001_5 <- Reduce(c,genes_list_0.0001_5) 
genes_union <- unique(genes_list_0.0001_5) #7339



fpkm_matrix.genes_final <- fpkm_matrix[rownames(fpkm_matrix) %in% genes_union,
                                               match(phenoData_training_final$sample, colnames(fpkm_matrix))]
log_fpkm_matrix.genes_final <- log_fpkm_matrix[rownames(log_fpkm_matrix) %in% genes_union,
                                               match(phenoData_training_final$sample, colnames(log_fpkm_matrix))]
zscore_log_fpkm.genes_final <- zscore_log_fpkm[rownames(zscore_log_fpkm) %in% genes_union,
                                               match(phenoData_training_final$sample, colnames(zscore_log_fpkm))]
zscore_variable_log_fpkm.genes_final <- zscore_variable_log_fpkm[rownames(zscore_variable_log_fpkm) %in% genes_union,
                                                                 match(phenoData_training_final$sample, colnames(zscore_variable_log_fpkm))]
zscore_fpkm.genes_final <- zscore_fpkm[rownames(zscore_fpkm) %in% genes_union,
                                       match(phenoData_training_final$sample, colnames(zscore_fpkm))]

tpm_matrix.genes_final <- tpm_matrix[rownames(tpm_matrix) %in% genes_union,
                                             match(phenoData_training_final$sample, colnames(tpm_matrix))]
log_tpm_matrix.genes_final <- log_tpm_matrix[rownames(log_tpm_matrix) %in% genes_union,
                                                     match(phenoData_training_final$sample, colnames(log_tpm_matrix))]
zscore_log_tpm.genes_final <- zscore_log_tpm[rownames(zscore_log_tpm) %in% genes_union,
                                                     match(phenoData_training_final$sample, colnames(zscore_log_tpm))]
zscore_variable_log_tpm.genes_final <- zscore_variable_log_tpm[rownames(zscore_variable_log_tpm) %in% genes_union,
                                                                       match(phenoData_training_final$sample, colnames(zscore_variable_log_tpm))]
zscore_tpm.genes_final <- zscore_tpm[rownames(zscore_tpm) %in% genes_union,
                                             match(phenoData_training_final$sample, colnames(zscore_tpm))]


######################################
############     PCA      ############
######################################

list_data <- list(fpkm_matrix.genes_final,
                  log_fpkm_matrix.genes_final,
                  zscore_log_fpkm.genes_final,
                  zscore_variable_log_fpkm.genes_final,
                  zscore_fpkm.genes_final,
                  tpm_matrix.genes_final,
                  log_tpm_matrix.genes_final,
                  zscore_log_tpm.genes_final,
                  zscore_variable_log_tpm.genes_final,
                  zscore_tpm.genes_final)

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
  autplot_PCA <- autoplot(PCA, data = phenoData_training_final, colour = "enfermedad_primaria",
                          main = as.character(names(list_data)[i]), size=0.4)+
    labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")+
    scale_color_manual(values = colors_blind)
  PCA_list[[i]] <- autplot_PCA
}

library("factoextra")
PCA_3 <- prcomp(t(as.matrix(list_data[[3]])), 
              scale. = FALSE, 
              center = FALSE)
fviz_eig(PCA_3, addlabels = TRUE, ylim = c(0, 50))
library(rgl)
phenoData_training_final_foo <- phenoData_training_final
df_color <- as.data.frame(cbind(colors_blind,unique(phenoData_training_final_foo$enfermedad_primaria)))
phenoData_training_final_foo$color <- df_color$colors_blind[match(phenoData_training_final_foo$enfermedad_primaria, df_color$V2)]
plot3d(PCA_3$x[,1:3], col = phenoData_training_final_foo$color)



PCA_4 <- prcomp(t(as.matrix(list_data[[4]])), 
              scale. = FALSE, 
              center = FALSE)
fviz_eig(PCA_4, addlabels = TRUE, ylim = c(0, 50))
plot3d(PCA_4$x[,1:3], col = phenoData_training_final_foo$color)



setwd(resultsDirectory)
pdf(file = "PCA_data_ML_selection.pdf", width = 14, height = 6)
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
                        Group = phenoData_training_final$enfermedad_primaria)
  
  UMAP_list[[i]] <- ggplot(df.umap, aes(x, y, colour = Group)) +
    geom_point(size=0.3) + ggtitle(as.character(names(list_data)[i]))+
    scale_color_manual(values = colors_blind)+
    labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")
}


setwd(resultsDirectory)
pdf(file = "UMAP_data_ML_selection.pdf", width = 14, height = 6)
grid.arrange(UMAP_list[[1]], UMAP_list[[2]], UMAP_list[[3]], UMAP_list[[4]], UMAP_list[[5]],
             UMAP_list[[6]], UMAP_list[[7]], UMAP_list[[8]], UMAP_list[[9]], UMAP_list[[10]],
             ncol=5)
dev.off()



umap_i <- umap(t(list_data[[2]]))
df.umap <- data.frame(x = umap_i$layout[,1],
                      y = umap_i$layout[,2],
                      Group = phenoData_training_final$enfermedad_primaria)
ggplot(df.umap, aes(x, y, colour = Group)) +
  geom_point(size=1) + ggtitle(as.character(names(list_data)[2]))+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()


######################################

