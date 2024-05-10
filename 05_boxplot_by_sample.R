#####################################################################
# GRAPHS Boxplots counts norm each sample
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
library(dplyr)


tablesDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Tables"
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Script 05"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"


######################################
# Read training data
######################################
setwd(tablesDirectory)

counts.vsd <- read.csv(file = "counts_vsd.csv")
counts.norm <- read.csv(file = "counts_norm.csv")
# Filtered

counts.vsd.f <- read.csv(file = "counts_vsd_f_after.csv")
counts.norm.f <- read.csv(file = "counts_norm_f_after.csv")
# phenoData
phenoData_training <- read.csv(file = "phenoData_training.csv")
######################################

rownames(phenoData_training) <- phenoData_training$X
phenoData_training <- phenoData_training[, colnames(phenoData_training)!="X"]

rownames(counts.vsd) <- counts.vsd$X
counts.vsd <- counts.vsd[, colnames(counts.vsd)!="X"]

rownames(counts.vsd.f) <- counts.vsd.f$X
counts.vsd.f <- counts.vsd.f[, colnames(counts.vsd.f)!="X"]

rownames(counts.norm) <- counts.norm$X
counts.norm <- counts.norm[, colnames(counts.norm)!="X"]

rownames(counts.norm.f) <- counts.norm.f$X
counts.norm.f <- counts.norm.f[, colnames(counts.norm.f)!="X"]

######################################
################ VSD  ################
######################################

# Primary disease
# Cargar la paleta de colores turbo
library(viridisLite)
colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

counts.vsd.gather <- gather(counts.vsd)
colnames(counts.vsd.gather) <- c("sample", "value")
counts.vsd.gather$enfermedad_primaria <- phenoData_training$enfermedad_primaria[match(counts.vsd.gather$sample, phenoData_training$sample)]


str(counts.vsd.gather)
counts.vsd.gather$sample <- as.factor(counts.vsd.gather$sample)
counts.vsd.gather$enfermedad_primaria <- as.factor(counts.vsd.gather$enfermedad_primaria)

# Median by primary disease
medianas_por_categoria <- counts.vsd.gather %>%
  group_by(enfermedad_primaria) %>%
  summarise(mediana = median(value))


counts.vsd.gather$mediana <- medianas_por_categoria$mediana[match(counts.vsd.gather$enfermedad_primaria, medianas_por_categoria$enfermedad_primaria)]


Boxplot_vsd_by_sample <- ggplot(counts.vsd.gather, 
                                         aes(value, sample)) +
  geom_boxplot(size=0.4) +
  geom_vline(data = counts.vsd.gather, aes(xintercept = mediana), col="red") +
  geom_vline(xintercept = median(counts.vsd.gather$value),col="blue") +
  facet_grid(enfermedad_primaria~., scales = "free", space = "free")+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 12, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot vsd") + 
  xlab("VSD counts") + ylab("Muestra")


setwd(resultsDirectory)
pdf(file = "BOXPLOT POR MUESTRA VSD.pdf", width = 10, height = 80)
print(Boxplot_vsd_by_sample)
dev.off()


######################################
################ NORM ################
######################################

counts.norm.gather <- gather(counts.norm)
colnames(counts.norm.gather) <- c("sample", "value")
counts.norm.gather$enfermedad_primaria <- phenoData_training$enfermedad_primaria[match(counts.norm.gather$sample, phenoData_training$sample)]


str(counts.norm.gather)
counts.norm.gather$sample <- as.factor(counts.norm.gather$sample)
counts.norm.gather$enfermedad_primaria <- as.factor(counts.norm.gather$enfermedad_primaria)

library(dplyr)
medianas_por_categoria_norm <- counts.norm.gather %>%
  group_by(enfermedad_primaria) %>%
  summarise(mediana = median(value))

counts.norm.gather$mediana <- medianas_por_categoria_norm$mediana[match(counts.norm.gather$enfermedad_primaria, medianas_por_categoria_norm$enfermedad_primaria)]


Boxplot_norm_by_sample <- ggplot(counts.norm.gather, 
                                aes(value, sample)) +
  geom_boxplot(size=0.4) +
  geom_vline(data = counts.norm.gather, aes(xintercept = mediana), col="red") +
  geom_vline(xintercept = median(counts.norm.gather$value),col="blue") +
  facet_grid(enfermedad_primaria~., scales = "free", space = "free")+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 12, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot norm") + 
  xlab("Norm counts") + ylab("Muestra")

setwd(resultsDirectory)
pdf(file = "BOXPLOT POR MUESTRA NORM.pdf", width = 10, height = 80)
print(Boxplot_norm_by_sample)
dev.off()



Boxplot_norm_by_sample_zoom <- ggplot(counts.norm.gather, 
                                 aes(value, sample)) +
  geom_boxplot(size=0.4) +
  geom_vline(data = counts.norm.gather, aes(xintercept = mediana), col="red") +
  geom_vline(xintercept = median(counts.norm.gather$value),col="blue") +
  facet_grid(enfermedad_primaria~., scales = "free", space = "free")+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 12, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") + scale_x_continuous(limits=c(0,1500000)) +
  ggtitle("Boxplot norm") + 
  xlab("Norm counts") + ylab("Muestra")


setwd(resultsDirectory)
pdf(file = "BOXPLOT POR MUESTRA NORM ZOOM.pdf", width = 10, height = 80)
print(Boxplot_norm_by_sample_zoom)
dev.off()


######################################
###############  MEDIAN  #############
######################################
medianas_por_muestra <- apply(counts.vsd, 2, median)
names(medianas_por_muestra) <- colnames(counts.vsd)
hist(medianas_por_muestra, breaks = 100)
boxplot(medianas_por_muestra)
3*median(medianas_por_muestra)
min(medianas_por_muestra)
sort(medianas_por_muestra, decreasing = FALSE)[3]

medianas_por_muestra_norm <- apply(counts.norm, 2, median)

library(outliers)
# VSD
grubbs.test(medianas_por_muestra, type=10)
grubbs.test(medianas_por_muestra[names(medianas_por_muestra)!="TCGA_V4_A9E9_01"], type=10)

#NORM
grubbs.test(medianas_por_muestra_norm, type=10)

library(MASS)
cov_matrix <- cov(counts.vsd)
mean_vector <- colMeans(counts.vsd)
counts.vsd$mahalanobis <- mahalanobis(counts.vsd, center = mean_vector, cov = cov_matrix)
threshold <- qchisq(0.975, df = ncol(counts.vsd) - 1)
counts.vsd$outlier <- counts.vsd$mahalanobis > threshold
table(counts.vsd$outlier)
