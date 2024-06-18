#####################################################################
# Transform data MET500 ML
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/MET500 Preparation"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

#####################################################################
setwd(tablesDirectory)
fpkm_matrix_training <- read.csv(file="fpkm_matrix.csv")

rownames(fpkm_matrix_training) <- fpkm_matrix_training$X
fpkm_matrix_training <- fpkm_matrix_training[, colnames(fpkm_matrix_training)!="X"]


setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

# Load FPKM for MET500
FPKM_MET500 <- fread("~/Desktop/LIDIA/TCGA_classification/Data/M.mx.txt", sep = "\t", stringsAsFactors = F)
genenames_MET500 <- separate(data = FPKM_MET500, col = sample, into = c("GeneName", "version"), remove = T)
FPKM_MET500 <- FPKM_MET500[,-1]
FPKM_MET500 <- as.data.frame(FPKM_MET500)
rownames(FPKM_MET500) <- genenames_MET500$GeneName

#####################################################################

table(rownames(FPKM_MET500) %in% rownames(fpkm_matrix_training))
#FALSE  TRUE 
#1916 19063 

FPKM_MET500.f <- FPKM_MET500[rownames(FPKM_MET500) %in% rownames(fpkm_matrix_training), ]

#####################################################################
# Samples
setwd(dataDirectory)
phenodata_MET500_1 <- fread(file="MET500_geneExpression_M.meta.plus.txt")
phenodata_MET500_2 <- fread(file="phenodata_complete_MET500.csv")


table(phenodata_MET500_1$Sample_id == colnames(FPKM_MET500.f))
# TRUE: 868

# Muestras:
muestras_polyA <- phenodata_MET500_1$Sample_id[grep("poly", phenodata_MET500_1$run.id)]
muestras_capt <- phenodata_MET500_1$Sample_id[grep("capt", phenodata_MET500_1$run.id)]

muestras_source_polyA <- phenodata_MET500_1$sample_source[grep("poly", phenodata_MET500_1$run.id)]
muestras_source_capt <- phenodata_MET500_1$sample_source[grep("capt", phenodata_MET500_1$run.id)]

table(muestras_source_polyA %in% muestras_source_capt)
#FALSE  TRUE 
#12   372 

muestras_capt_outofpolyA <- muestras_source_capt[!muestras_source_capt %in% muestras_source_polyA]
muestras_capt_outofpolyA <- phenodata_MET500_1$Sample_id[phenodata_MET500_1$sample_source %in% muestras_capt_outofpolyA]

muestras_totales <- c(muestras_polyA, muestras_capt_outofpolyA)


# Filtering
FPKM_MET500.f <- FPKM_MET500.f[, colnames(FPKM_MET500.f) %in% muestras_totales]
phenodata_MET500.f <- phenodata_MET500_1[match(colnames(FPKM_MET500.f), phenodata_MET500_1$run.id),]

## Prepare tissue
table(is.na(phenodata_MET500.f$tissue)) #FALSE
table(is.na(phenodata_MET500.f$biopsy_tissue)) #FALSE

table(phenodata_MET500.f$tissue=="")
#FALSE  TRUE 
#463    33 
table(phenodata_MET500.f$biopsy_tissue=="")
#FALSE 
#496

table(phenodata_MET500.f$sample_type) #tumor

phenodata_MET500.f <- phenodata_MET500.f[phenodata_MET500.f$tissue!="",]
phenodata_MET500.f <- as.data.frame(phenodata_MET500.f)

table(phenodata_MET500.f$tissue==phenodata_MET500.f$biopsy_tissue)
#FALSE  TRUE 
#434    29 

#####################################################################

setwd(tablesDirectory)
write.csv(phenodata_MET500.f, file="phenodata_MET500.f.csv")
write.csv(FPKM_MET500.f, file="FPKM_MET500.f.csv")

#####################################################################

phenodata_MET500_clean <- phenodata_MET500.f

phenodata_MET500_clean <- phenodata_MET500_clean[phenodata_MET500_clean$tissue!=phenodata_MET500_clean$biopsy_tissue,] #434 samples

table(phenodata_MET500_clean$tissue %in% c("gall_bladder", "parotid", "soft_tissue", "lymph_node", "bone_marrow", "other"))
#FALSE  TRUE 
#387    47 

phenodata_MET500_clean <- phenodata_MET500_clean[!phenodata_MET500_clean$tissue %in% c("gall_bladder", "parotid", "soft_tissue", "lymph_node", "bone_marrow", "other"),] #387



FPKM_MET500_clean <- FPKM_MET500.f[,colnames(FPKM_MET500.f) %in% phenodata_MET500_clean$Sample_id]
phenodata_MET500_clean <- phenodata_MET500_clean[match(colnames(FPKM_MET500_clean), phenodata_MET500_clean$run.id),]


setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file="df_tissue_biopsy_tissue.csv", sep=",")
df_tissue_biopsy_tissue <- df_tissue_biopsy_tissue[,colnames(df_tissue_biopsy_tissue)!="X"]

phenodata_MET500_clean$tissue_set <- df_tissue_biopsy_tissue$set[match(phenodata_MET500_clean$tissue, df_tissue_biopsy_tissue$tissue)]
phenodata_MET500_clean$tejido <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$tissue, df_tissue_biopsy_tissue$tissue)]

phenodata_MET500_clean$biopsy_tissue_set <- df_tissue_biopsy_tissue$set[match(phenodata_MET500_clean$biopsy_tissue, df_tissue_biopsy_tissue$tissue)]
phenodata_MET500_clean$biopsia_tejido <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$biopsy_tissue, df_tissue_biopsy_tissue$tissue)]

#####################################################################
setwd(tablesDirectory)
write.csv(phenodata_MET500_clean, file="phenodata_MET500_clean.csv")
write.csv(FPKM_MET500_clean, file="FPKM_MET500_clean.csv")
#####################################################################

setwd(tablesDirectory)
medias_genes_log_fpkm <- readRDS(file="medias_genes_log_fpkm.rds")
sd_genes_log_fpkm <- readRDS(file="sd_genes_log_fpkm.rds")
medias_genes_log_tpm <- readRDS(file = "medias_genes_log_tpm.rds")
sd_genes_log_tpm <- readRDS(file = "sd_genes_log_tpm.rds")
#####################################################################

# TPM
sum_fpkm_per_sample <- colSums(FPKM_MET500_clean)
scaling_factors <- sum_fpkm_per_sample / 1e6
tpm_MET500_clean <- t(t(FPKM_MET500_clean) / scaling_factors)

# logs

log_fpkm_MET500_clean <- log2(FPKM_MET500_clean+0.001)
log_tpm_MET500_clean <- log2(tpm_MET500_clean+0.001)

# zscore log

zscore_log_fpkm_MET500_clean <- scale(log_fpkm_MET500_clean)
zscore_log_tpm_MET500_clean <- scale(log_tpm_MET500_clean)

# zscore variable log

zscore_variable_log_FPKM_MET500_clean <- apply(log_fpkm_MET500_clean,2,function(x){
  (x-medias_genes_log_fpkm)/sd_genes_log_fpkm
})

zscore_variable_log_tpm_MET500_clean <- apply(log_tpm_MET500_clean,2,function(x){
  (x-medias_genes_log_tpm)/sd_genes_log_tpm
})

# zscore

zscore_fpkm_MET500_clean <- scale(FPKM_MET500_clean)
zscore_tpm_MET500_clean <- scale(tpm_MET500_clean)


######################################
setwd(tablesDirectory)
write.csv(log_fpkm_MET500_clean, file="log_fpkm_MET500_clean.csv")
write.csv(zscore_log_fpkm_MET500_clean, file = "zscore_log_fpkm_MET500_clean.csv")
write.csv(zscore_variable_log_FPKM_MET500_clean, file = "zscore_variable_log_FPKM_MET500_clean.csv")
write.csv(zscore_fpkm_MET500_clean, file="zscore_fpkm_MET500_clean.csv")

write.csv(tpm_MET500_clean, file="tpm_MET500_clean.csv")
write.csv(log_tpm_MET500_clean, file="log_tpm_MET500_clean.csv")
write.csv(zscore_log_tpm_MET500_clean, file = "zscore_log_tpm_MET500_clean.csv")
write.csv(zscore_variable_log_tpm_MET500_clean, file = "zscore_variable_log_tpm_MET500_clean.csv")
write.csv(zscore_tpm_MET500_clean, file="zscore_tpm_MET500_clean.csv")
######################################



#####################################################################
############# GRAPHS
#####################################################################

library(viridisLite)
colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

ggplot(data=phenodata_MET500_clean, aes(x=tejido, fill=biopsia_tejido))+
  geom_bar()+ theme(axis.text.x = element_text(size=10, angle=90)) + ggtitle("Muestras MET500 (limpias)") +
  xlab("Tejido") +
  scale_fill_manual(name = "Tejido biopsia",values=c(colors_blind[c(1,2,3,6)],"black", colors_blind[c(9,12)], "grey",colors_blind[c(15,17,19)], "red", colors_blind[27]))


######################################
# Read data
######################################
setwd(tablesDirectory)
genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")
phenodata_MET500_clean <- read.csv(file="phenodata_MET500_clean.csv")
FPKM_MET500_clean <- read.csv(file="FPKM_MET500_clean.csv")
log_fpkm_MET500_clean <- read.csv(file="log_fpkm_MET500_clean.csv")
zscore_log_fpkm_MET500_clean <- read.csv(file = "zscore_log_fpkm_MET500_clean.csv")
zscore_variable_log_FPKM_MET500_clean <- read.csv(file = "zscore_variable_log_FPKM_MET500_clean.csv")
zscore_fpkm_MET500_clean <- read.csv(file="zscore_fpkm_MET500_clean.csv")

tpm_MET500_clean <- read.csv(file="tpm_MET500_clean.csv")
log_tpm_MET500_clean <- read.csv(file="log_tpm_MET500_clean.csv")
zscore_log_tpm_MET500_clean <- read.csv(file = "zscore_log_tpm_MET500_clean.csv")
zscore_variable_log_tpm_MET500_clean <- read.csv(file = "zscore_variable_log_tpm_MET500_clean.csv")
zscore_tpm_MET500_clean <- read.csv(file="zscore_tpm_MET500_clean.csv")
######################################

rownames(phenodata_MET500_clean) <- phenodata_MET500_clean$X
phenodata_MET500_clean <- phenodata_MET500_clean[, colnames(phenodata_MET500_clean)!="X"]

rownames(FPKM_MET500_clean) <- FPKM_MET500_clean$X
FPKM_MET500_clean <- FPKM_MET500_clean[, colnames(FPKM_MET500_clean)!="X"]
colnames(FPKM_MET500_clean) <- str_replace_all(colnames(FPKM_MET500_clean), "[.]", "-")

rownames(log_fpkm_MET500_clean) <- log_fpkm_MET500_clean$X
log_fpkm_MET500_clean <- log_fpkm_MET500_clean[, colnames(log_fpkm_MET500_clean)!="X"]
colnames(log_fpkm_MET500_clean) <- str_replace_all(colnames(log_fpkm_MET500_clean), "[.]", "-")

rownames(zscore_log_fpkm_MET500_clean) <- zscore_log_fpkm_MET500_clean$X
zscore_log_fpkm_MET500_clean <- zscore_log_fpkm_MET500_clean[, colnames(zscore_log_fpkm_MET500_clean)!="X"]
colnames(zscore_log_fpkm_MET500_clean) <- str_replace_all(colnames(zscore_log_fpkm_MET500_clean), "[.]", "-")

rownames(zscore_variable_log_FPKM_MET500_clean) <- zscore_variable_log_FPKM_MET500_clean$X
zscore_variable_log_FPKM_MET500_clean <- zscore_variable_log_FPKM_MET500_clean[, colnames(zscore_variable_log_FPKM_MET500_clean)!="X"]
colnames(zscore_variable_log_FPKM_MET500_clean) <- str_replace_all(colnames(zscore_variable_log_FPKM_MET500_clean), "[.]", "-")

rownames(zscore_fpkm_MET500_clean) <- zscore_fpkm_MET500_clean$X
zscore_fpkm_MET500_clean <- zscore_fpkm_MET500_clean[, colnames(zscore_fpkm_MET500_clean)!="X"]
colnames(zscore_fpkm_MET500_clean) <- str_replace_all(colnames(zscore_fpkm_MET500_clean), "[.]", "-")

rownames(tpm_MET500_clean) <- tpm_MET500_clean$X
tpm_MET500_clean <- tpm_MET500_clean[, colnames(tpm_MET500_clean)!="X"]
colnames(tpm_MET500_clean) <- str_replace_all(colnames(tpm_MET500_clean), "[.]", "-")

rownames(log_tpm_MET500_clean) <- log_tpm_MET500_clean$X
log_tpm_MET500_clean <- log_tpm_MET500_clean[, colnames(log_tpm_MET500_clean)!="X"]
colnames(log_tpm_MET500_clean) <- str_replace_all(colnames(log_tpm_MET500_clean), "[.]", "-")

rownames(zscore_log_tpm_MET500_clean) <- zscore_log_tpm_MET500_clean$X
zscore_log_tpm_MET500_clean <- zscore_log_tpm_MET500_clean[, colnames(zscore_log_tpm_MET500_clean)!="X"]
colnames(zscore_log_tpm_MET500_clean) <- str_replace_all(colnames(zscore_log_tpm_MET500_clean), "[.]", "-")

rownames(zscore_variable_log_tpm_MET500_clean) <- zscore_variable_log_tpm_MET500_clean$X
zscore_variable_log_tpm_MET500_clean <- zscore_variable_log_tpm_MET500_clean[, colnames(zscore_variable_log_tpm_MET500_clean)!="X"]
colnames(zscore_variable_log_tpm_MET500_clean) <- str_replace_all(colnames(zscore_variable_log_tpm_MET500_clean), "[.]", "-")

rownames(zscore_tpm_MET500_clean) <- zscore_tpm_MET500_clean$X
zscore_tpm_MET500_clean <- zscore_tpm_MET500_clean[, colnames(zscore_tpm_MET500_clean)!="X"]
colnames(zscore_tpm_MET500_clean) <- str_replace_all(colnames(zscore_tpm_MET500_clean), "[.]", "-")
######################################
#CHECK:
table(colnames(FPKM_MET500_clean)==phenodata_MET500_clean$Sample_id)
######################################


######################################
# Filter genes union final 

genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

genes_list_0.0001_5 <- Reduce(c,genes_list_0.0001_5) 
genes_union <- unique(genes_list_0.0001_5) #7339



FPKM_MET500_clean.genes_final <- FPKM_MET500_clean[rownames(FPKM_MET500_clean) %in% genes_union,
                                       match(phenodata_MET500_clean$Sample_id, colnames(FPKM_MET500_clean))]
log_fpkm_MET500_clean.genes_final <- log_fpkm_MET500_clean[rownames(log_fpkm_MET500_clean) %in% genes_union,
                                               match(phenodata_MET500_clean$Sample_id, colnames(log_fpkm_MET500_clean))]
zscore_log_fpkm_MET500_clean.genes_final <- zscore_log_fpkm_MET500_clean[rownames(zscore_log_fpkm_MET500_clean) %in% genes_union,
                                               match(phenodata_MET500_clean$Sample_id, colnames(zscore_log_fpkm_MET500_clean))]
zscore_variable_log_FPKM_MET500_clean.genes_final <- zscore_variable_log_FPKM_MET500_clean[rownames(zscore_variable_log_FPKM_MET500_clean) %in% genes_union,
                                                                 match(phenodata_MET500_clean$Sample_id, colnames(zscore_variable_log_FPKM_MET500_clean))]
zscore_fpkm_MET500_clean.genes_final <- zscore_fpkm_MET500_clean[rownames(zscore_fpkm_MET500_clean) %in% genes_union,
                                       match(phenodata_MET500_clean$Sample_id, colnames(zscore_fpkm_MET500_clean))]

tpm_MET500_clean.genes_final <- tpm_MET500_clean[rownames(tpm_MET500_clean) %in% genes_union,
                                     match(phenodata_MET500_clean$Sample_id, colnames(tpm_MET500_clean))]
log_tpm_MET500_clean.genes_final <- log_tpm_MET500_clean[rownames(log_tpm_MET500_clean) %in% genes_union,
                                             match(phenodata_MET500_clean$Sample_id, colnames(log_tpm_MET500_clean))]
zscore_log_tpm_MET500_clean.genes_final <- zscore_log_tpm_MET500_clean[rownames(zscore_log_tpm_MET500_clean) %in% genes_union,
                                             match(phenodata_MET500_clean$Sample_id, colnames(zscore_log_tpm_MET500_clean))]
zscore_variable_log_tpm_MET500_clean.genes_final <- zscore_variable_log_tpm_MET500_clean[rownames(zscore_variable_log_tpm_MET500_clean) %in% genes_union,
                                                               match(phenodata_MET500_clean$Sample_id, colnames(zscore_variable_log_tpm_MET500_clean))]
zscore_tpm_MET500_clean.genes_final <- zscore_tpm_MET500_clean[rownames(zscore_tpm_MET500_clean) %in% genes_union,
                                     match(phenodata_MET500_clean$Sample_id, colnames(zscore_tpm_MET500_clean))]



######################################
############     PCA      ############
######################################

list_data <- list(FPKM_MET500_clean.genes_final,
                  log_fpkm_MET500_clean.genes_final,
                  zscore_log_fpkm_MET500_clean.genes_final,
                  zscore_variable_log_FPKM_MET500_clean.genes_final,
                  zscore_fpkm_MET500_clean.genes_final,
                  tpm_MET500_clean.genes_final,
                  log_tpm_MET500_clean.genes_final,
                  zscore_log_tpm_MET500_clean.genes_final,
                  zscore_variable_log_tpm_MET500_clean.genes_final,
                  zscore_tpm_MET500_clean.genes_final)

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

colors_blind <- c("#467EF4FF", "#21C6E0FF", "#FE8F28FF", "#FEAA33FF", "#424AB3FF", "#CD2C04FF", "#30123BFF", "#ED5510FF", "#E9D539FF", "#F8721CFF", "#D4E735FF", "#1FE9AEFF",
                  "#38F491FF","#18DAC7FF", "#A2FC3CFF", "#83FF52FF", "#E03F08FF")

unique(phenodata_MET500_clean$tejido)
#"mama"           "colon"          "sarcoma"        "próstata"       "vejiga"         "timoma"         "adrenocortical" "estómago"       "páncreas"       "piel"           "ovario"        
# "cabeza cuello"  "riñón"          "esófago"        "pulmón"         "hígado"         "testículos"    
set.seed(42)
PCA_list <- list()

for (i in 1:length(list_data)){
  PCA <- prcomp(t(as.matrix(list_data[[i]])), 
                scale. = FALSE, 
                center = FALSE)
  autplot_PCA <- autoplot(PCA, data = phenodata_MET500_clean, colour = "tejido",
                          main = as.character(names(list_data)[i]), size=0.4)+
    labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")+
    scale_color_manual(values = colors_blind)
  PCA_list[[i]] <- autplot_PCA
}


setwd(resultsDirectory)
pdf(file = "PCA_data_ML_MET500.pdf", width = 14, height = 6)
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
                        Group = phenodata_MET500_clean$tejido)
  
  UMAP_list[[i]] <- ggplot(df.umap, aes(x, y, colour = Group)) +
    geom_point(size=0.3) + ggtitle(as.character(names(list_data)[i]))+
    scale_color_manual(values = colors_blind)+ theme_light() + theme(legend.position = "none")
}


setwd(resultsDirectory)
pdf(file = "UMAP_data_ML_selection.pdf", width = 14, height = 6)
grid.arrange(UMAP_list[[1]], UMAP_list[[2]], UMAP_list[[3]], UMAP_list[[4]], UMAP_list[[5]],
             UMAP_list[[6]], UMAP_list[[7]], UMAP_list[[8]], UMAP_list[[9]], UMAP_list[[10]],
             ncol=5)
dev.off()



