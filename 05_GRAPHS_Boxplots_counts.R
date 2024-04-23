#####################################################################
# GRAPHS Boxplots
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
log_rnaseq_training <- read.csv(file = "log_rnaseq_training.csv")
counts.vsd <- read.csv(file = "counts_vsd.csv")
counts.norm <- read.csv(file = "counts_norm.csv")
# Filtered
rnaseq_counts_training.f <- read.csv(file = "rnaseq_counts_training_f.csv")
log_rnaseq_training.f <- read.csv(file = "log_rnaseq_training_f.csv")
counts.vsd.f <- read.csv(file = "counts_vsd_f_after.csv")
counts.norm.f <- read.csv(file = "counts_norm_f_after.csv")
# phenoData
phenoData_training <- read.csv(file = "phenoData_training.csv")
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(rnaseq_counts_training) <- rnaseq_counts_training$X
rnaseq_counts_training <- rnaseq_counts_training[, colnames(rnaseq_counts_training)!="X"]

rownames(rnaseq_counts_training.f) <- rnaseq_counts_training.f$X
rnaseq_counts_training.f <- rnaseq_counts_training.f[, colnames(rnaseq_counts_training.f)!="X"]

rownames(log_rnaseq_training) <- log_rnaseq_training$X
log_rnaseq_training <- log_rnaseq_training[, colnames(log_rnaseq_training)!="X"]

rownames(log_rnaseq_training.f) <- log_rnaseq_training.f$X
log_rnaseq_training.f <- log_rnaseq_training.f[, colnames(log_rnaseq_training.f)!="X"]

rownames(counts.vsd) <- counts.vsd$X
counts.vsd <- counts.vsd[, colnames(counts.vsd)!="X"]

rownames(counts.vsd.f) <- counts.vsd.f$X
counts.vsd.f <- counts.vsd.f[, colnames(counts.vsd.f)!="X"]

rownames(counts.norm) <- counts.norm$X
counts.norm <- counts.norm[, colnames(counts.norm)!="X"]

rownames(counts.norm.f) <- counts.norm.f$X
counts.norm.f <- counts.norm.f[, colnames(counts.norm.f)!="X"]

######################################



############################################################################
##################              DATOS CRUDOS             ###################
##################              sin filtrar              ###################
##################                  SUMA                 ###################
############################################################################

rnaseq_counts_training.sum <- colSums(rnaseq_counts_training)
rnaseq_counts_training.sum <- as.data.frame(rnaseq_counts_training.sum)
colnames(rnaseq_counts_training.sum) <- "counts"
table(rownames(rnaseq_counts_training.sum)==phenoData_training$sample)

rnaseq_counts_training.sum$enfermedad_primaria <- as.factor(phenoData_training$enfermedad_primaria)
rnaseq_counts_training.sum$raza <- as.factor(phenoData_training$raza)
rnaseq_counts_training.sum$sexo <- as.factor(phenoData_training$sexo)
rnaseq_counts_training.sum$edad <- phenoData_training$edad


str(rnaseq_counts_training.sum)

######
# Sex
######
Boxplot_rnaseq_sex <- ggplot(rnaseq_counts_training.sum, aes(y = counts, x=sexo, color = sexo)) +
  geom_boxplot() +
  scale_color_manual(values = c("deepskyblue3", "coral1")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot counts (Suma)") + xlab("Sexo") + ylab("Suma de los counts (por muestra)")

######
# Race
######

Boxplot_rnaseq_race <- ggplot(rnaseq_counts_training.sum, aes(y = counts, x=raza, color = raza)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot counts (Suma)") + xlab("Raza") + ylab("Suma de los counts (por muestra)")


# Primary disease
# Cargar la paleta de colores turbo
#library(viridisLite)
#colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

Boxplot_rnaseq_primary_disease <- ggplot(rnaseq_counts_training.sum, 
                                         aes(y = counts, x=enfermedad_primaria)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 4),
        legend.position = "none") +
  ggtitle("Boxplot counts (Suma)") + xlab("Tumor primario") + ylab("Suma de los counts (por muestra)")


# Primary disease and sex

Boxplot_rnaseq_primary_disease_sex <- ggplot(rnaseq_counts_training.sum, 
                                                                  aes(y = counts, x=enfermedad_primaria, 
                                                                                                  color = sexo)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("deepskyblue3", "coral1")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 4)) +
  ggtitle("Boxplot counts (Suma)") + xlab("Tumor primario") + ylab("Suma de los counts (por muestra)")


############################################################################
##################              DATOS CRUDOS             ###################
##################              sin filtrar              ###################
##################                  MEDIA                 ###################
############################################################################



rnaseq_counts_training.mean <- apply(rnaseq_counts_training, 2, mean)
rnaseq_counts_training.mean <- as.data.frame(rnaseq_counts_training.mean)
colnames(rnaseq_counts_training.mean) <- "counts"
table(rownames(rnaseq_counts_training.mean)==phenoData_training$sample)

rnaseq_counts_training.mean$enfermedad_primaria <- as.factor(phenoData_training$enfermedad_primaria)
rnaseq_counts_training.mean$raza <- as.factor(phenoData_training$raza)
rnaseq_counts_training.mean$sexo <- as.factor(phenoData_training$sexo)
rnaseq_counts_training.mean$edad <- phenoData_training$edad


str(rnaseq_counts_training.mean)


######
# Sex
######
Boxplot_rnaseq_mean_sex <- ggplot(rnaseq_counts_training.mean, aes(y = counts, x=sexo, color = sexo)) +
  geom_boxplot() +
  scale_color_manual(values = c("deepskyblue3", "coral1")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot counts (Media)") + xlab("Sexo") + ylab("Media de los counts (por muestra)")



######
# Race
######

Boxplot_rnaseq_mean_race <- ggplot(rnaseq_counts_training.mean, aes(y = counts, x=raza, color = raza)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot counts (Media)") + xlab("Raza") + ylab("Media de los counts (por muestra)")


# Primary disease
# Cargar la paleta de colores turbo
#library(viridisLite)
#colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

Boxplot_rnaseq_mean_primary_disease <- ggplot(rnaseq_counts_training.mean, 
                                         aes(y = counts, x=enfermedad_primaria)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 4),
        legend.position = "none") +
  ggtitle("Boxplot counts (Media)") + xlab("Tumor primario") + ylab("Media de los counts (por muestra)")


# Primary disease and sex

Boxplot_rnaseq_mean_primary_disease_sex <- ggplot(rnaseq_counts_training.mean, 
                                             aes(y = counts, x=enfermedad_primaria, 
                                                 color = sexo)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("deepskyblue3", "coral1")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 4)) +
  ggtitle("Boxplot counts (Media)") + xlab("Tumor primario") + ylab("Media de los counts (por muestra)")



################################################################################
################################################################################
#######################              TOTAL            ##########################
################################################################################
################################################################################

library(gridExtra)


setwd(resultsDirectory)
pdf(file = "Boxplot_comparison.pdf", width = 8, height = 10)
grid.arrange(Boxplot_rnaseq_sex, Boxplot_rnaseq_mean_sex,
             Boxplot_rnaseq_race, Boxplot_rnaseq_mean_race,
             Boxplot_rnaseq_primary_disease, Boxplot_rnaseq_mean_primary_disease,
             Boxplot_rnaseq_primary_disease_sex, Boxplot_rnaseq_mean_primary_disease_sex,
             ncol=2)
dev.off()
