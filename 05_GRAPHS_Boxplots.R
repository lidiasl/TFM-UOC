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
############################################################################

rnaseq_counts_training.sum <- colSums(rnaseq_counts_training)
rnaseq_counts_training.sum <- as.data.frame(rnaseq_counts_training.sum)
colnames(rnaseq_counts_training.sum) <- "counts"
table(rownames(rnaseq_counts_training.sum)==phenoData_training$sample)

rnaseq_counts_training.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
rnaseq_counts_training.sum$race <- phenoData_training$race
rnaseq_counts_training.sum$race[rnaseq_counts_training.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
rnaseq_counts_training.sum$race <- as.factor(rnaseq_counts_training.sum$race)
rnaseq_counts_training.sum$sex <- as.factor(phenoData_training$sex)
rnaseq_counts_training.sum$age <- phenoData_training$age


str(rnaseq_counts_training.sum)

######
# Sex
######
Boxplot_rnaseq_sex <- ggplot(rnaseq_counts_training.sum, aes(y = counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot rnaseq") + xlab("Sex") + ylab("Counts")

######
# Race
######

Boxplot_rnaseq_race <- ggplot(rnaseq_counts_training.sum, aes(y = counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot rnaseq") + xlab("Race") + ylab("Counts")


# Primary disease
# Cargar la paleta de colores turbo
#library(viridisLite)
#colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

Boxplot_rnaseq_primary_disease <- ggplot(rnaseq_counts_training.sum, 
                                         aes(y = counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot rnaseq") + xlab("Primary disease") + ylab("Counts")


# Primary disease and sex

Boxplot_rnaseq_primary_disease_sex <- ggplot(rnaseq_counts_training.sum, 
                                                                  aes(y = counts, x=primary_disease_set, 
                                                                                                  color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot rnaseq") + xlab("Primary disease") + ylab("Counts")




############################################################################
##################              DATOS CRUDOS             ###################
##################              Filtrados              ###################
############################################################################


rnaseq_counts_training.f.sum <- colSums(rnaseq_counts_training.f)
rnaseq_counts_training.f.sum <- as.data.frame(rnaseq_counts_training.f.sum)
colnames(rnaseq_counts_training.f.sum) <- "counts"
table(rownames(rnaseq_counts_training.f.sum)==phenoData_training$sample)

rnaseq_counts_training.f.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
rnaseq_counts_training.f.sum$race <- phenoData_training$race
rnaseq_counts_training.f.sum$race[rnaseq_counts_training.f.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
rnaseq_counts_training.f.sum$race <- as.factor(rnaseq_counts_training.f.sum$race)
rnaseq_counts_training.f.sum$sex <- as.factor(phenoData_training$sex)
rnaseq_counts_training.f.sum$age <- phenoData_training$age

str(rnaseq_counts_training.f.sum)


######
# Sex
######
Boxplot_F_rnaseq_sex <- ggplot(rnaseq_counts_training.f.sum, aes(y = counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot rnaseq filtered") + xlab("Sex") + ylab("Counts")


######
# Race
######

Boxplot_F_rnaseq_race <- ggplot(rnaseq_counts_training.f.sum, aes(y = counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot rnaseq filtered") + xlab("Race") + ylab("Counts")


# Primary disease
# Cargar la paleta de colores turbo

Boxplot_F_rnaseq_primary_disease <- ggplot(rnaseq_counts_training.f.sum, 
                                         aes(y = counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot rnaseq filtered") + xlab("Primary disease") + ylab("Counts")


# Primary disease and sex

Boxplot_F_rnaseq_primary_disease_sex <- ggplot(rnaseq_counts_training.f.sum, 
                                             aes(y = counts, x=primary_disease_set, 
                                                 color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot rnaseq filtered") + xlab("Primary disease") + ylab("Counts")



############################################################################
##################               DATOS LOG               ###################
##################              sin filtrar              ###################
############################################################################


log_rnaseq_training.sum <- colSums(log_rnaseq_training)
log_rnaseq_training.sum <- as.data.frame(log_rnaseq_training.sum)
colnames(log_rnaseq_training.sum) <- "log_counts"
table(rownames(log_rnaseq_training.sum)==phenoData_training$sample)

log_rnaseq_training.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
log_rnaseq_training.sum$race <- phenoData_training$race
log_rnaseq_training.sum$race[log_rnaseq_training.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
log_rnaseq_training.sum$race <- as.factor(log_rnaseq_training.sum$race)
log_rnaseq_training.sum$sex <- as.factor(phenoData_training$sex)
log_rnaseq_training.sum$age <- phenoData_training$age

str(log_rnaseq_training.sum)


######
# Sex
######
Boxplot_log_rnaseq_sex <- ggplot(log_rnaseq_training.sum, aes(y = log_counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot log rnaseq") + xlab("Sex") + ylab("Log counts")

######
# Race
######

Boxplot_log_rnaseq_race <- ggplot(log_rnaseq_training.sum, aes(y = log_counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot log rnaseq") + xlab("Race") + ylab("Log counts")



# Primary disease

Boxplot_log_rnaseq_primary_disease <- ggplot(log_rnaseq_training.sum, 
                                             aes(y = log_counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot log rnaseq") + xlab("Primary disease") + ylab("Log counts")


# Primary disease and sex

Boxplot_log_rnaseq_primary_disease_sex <- ggplot(log_rnaseq_training.sum, 
                                             aes(y = log_counts, x=primary_disease_set, 
                                                 color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot log rnaseq") + xlab("Primary disease") + ylab("Log counts")


############################################################################
##################               DATOS LOG               ###################
##################               Filtrados               ###################
############################################################################


log_rnaseq_training.f.sum <- colSums(log_rnaseq_training.f)
log_rnaseq_training.f.sum <- as.data.frame(log_rnaseq_training.f.sum)
colnames(log_rnaseq_training.f.sum) <- "log_counts"
table(rownames(log_rnaseq_training.f.sum)==phenoData_training$sample)

log_rnaseq_training.f.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
log_rnaseq_training.f.sum$race <- phenoData_training$race
log_rnaseq_training.f.sum$race[log_rnaseq_training.f.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
log_rnaseq_training.f.sum$race <- as.factor(log_rnaseq_training.f.sum$race)
log_rnaseq_training.f.sum$sex <- as.factor(phenoData_training$sex)
log_rnaseq_training.f.sum$age <- phenoData_training$age

str(log_rnaseq_training.f.sum)

######
# Sex
######
Boxplot_F_log_rnaseq_sex <- ggplot(log_rnaseq_training.f.sum, aes(y = log_counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot log rnaseq filtered") + xlab("Sex") + ylab("Log counts")


######
# Race
######

Boxplot_F_log_rnaseq_race <- ggplot(log_rnaseq_training.f.sum, aes(y = log_counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot log rnaseq filtered") + xlab("Race") + ylab("Log counts")

# Primary disease

Boxplot_F_log_rnaseq_primary_disease <- ggplot(log_rnaseq_training.f.sum, 
                                             aes(y = log_counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot log rnaseq filtered") + xlab("Primary disease") + ylab("Log counts")


# Primary disease and sex

Boxplot_F_log_rnaseq_primary_disease_sex <- ggplot(log_rnaseq_training.f.sum, 
                                                 aes(y = log_counts, x=primary_disease_set, 
                                                     color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot log rnaseq filtered") + xlab("Primary disease") + ylab("Log counts")


############################################################################
##################               DATOS VSD               ###################
##################              sin filtrar              ###################
############################################################################


counts.vsd.sum <- colSums(counts.vsd)
counts.vsd.sum <- as.data.frame(counts.vsd.sum)
colnames(counts.vsd.sum) <- "counts"
table(rownames(counts.vsd.sum)==phenoData_training$sample)

counts.vsd.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
counts.vsd.sum$race <- phenoData_training$race
counts.vsd.sum$race[counts.vsd.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
counts.vsd.sum$race <- as.factor(counts.vsd.sum$race)
counts.vsd.sum$sex <- as.factor(phenoData_training$sex)
counts.vsd.sum$age <- phenoData_training$age


str(counts.vsd.sum)


######
# Sex
######
Boxplot_vsd_sex <- ggplot(counts.vsd.sum, aes(y = counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot vsd") + xlab("Sex") + ylab("Vsd counts")



######
# Race
######

Boxplot_vsd_race <- ggplot(counts.vsd.sum, aes(y = counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot vsd") + xlab("Race") + ylab("Vsd counts")


# Primary disease

Boxplot_vsd_primary_disease <- ggplot(counts.vsd.sum, 
                                      aes(y = counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot vsd") + xlab("Primary disease") + ylab("Vsd counts")


# Primary disease and sex

Boxplot_vsd_primary_disease_sex <- ggplot(counts.vsd.sum, 
                                             aes(y = counts, x=primary_disease_set, 
                                                 color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot vsd") + xlab("Primary disease") + ylab("Vsd counts")


############################################################################
##################               DATOS VSD               ###################
##################               Filtrados               ###################
############################################################################


counts.vsd.f.sum <- colSums(counts.vsd.f)
counts.vsd.f.sum <- as.data.frame(counts.vsd.f.sum)
colnames(counts.vsd.f.sum) <- "counts"
table(rownames(counts.vsd.f.sum)==phenoData_training$sample)

counts.vsd.f.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
counts.vsd.f.sum$race <- phenoData_training$race
counts.vsd.f.sum$race[counts.vsd.f.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
counts.vsd.f.sum$race <- as.factor(counts.vsd.f.sum$race)
counts.vsd.f.sum$sex <- as.factor(phenoData_training$sex)
counts.vsd.f.sum$age <- phenoData_training$age

str(counts.vsd.f.sum)


######
# Sex
######
Boxplot_F_vsd_sex <- ggplot(counts.vsd.f.sum, aes(y = counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot vsd filtered") + xlab("Sex") + ylab("Vsd counts")


######
# Race
######

Boxplot_F_vsd_race <- ggplot(counts.vsd.f.sum, aes(y = counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot vsd filtered") + xlab("Race") + ylab("Vsd counts")


# Primary disease

Boxplot_F_vsd_primary_disease <- ggplot(counts.vsd.f.sum, 
                                      aes(y = counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot vsd filtered") + xlab("Primary disease") + ylab("Vsd counts")


# Primary disease and sex

Boxplot_F_vsd_primary_disease_sex <- ggplot(counts.vsd.f.sum, 
                                          aes(y = counts, x=primary_disease_set, 
                                              color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot vsd filtered") + xlab("Primary disease") + ylab("Vsd counts")



############################################################################
##################               DATOS norm               ###################
##################              sin filtrar              ###################
############################################################################


counts.norm.sum <- colSums(counts.norm)
counts.norm.sum <- as.data.frame(counts.norm.sum)
colnames(counts.norm.sum) <- "counts"
table(rownames(counts.norm.sum)==phenoData_training$sample)

counts.norm.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
counts.norm.sum$race <- phenoData_training$race
counts.norm.sum$race[counts.norm.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
counts.norm.sum$race <- as.factor(counts.norm.sum$race)
counts.norm.sum$sex <- as.factor(phenoData_training$sex)
counts.norm.sum$age <- phenoData_training$age


str(counts.norm.sum)


######
# Sex
######
Boxplot_norm_sex <- ggplot(counts.norm.sum, aes(y = counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot norm") + xlab("Sex") + ylab("Norm counts")


######
# Race
######

Boxplot_norm_race <- ggplot(counts.norm.sum, aes(y = counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot norm") + xlab("Race") + ylab("Norm counts")


# Primary disease

Boxplot_norm_primary_disease <- ggplot(counts.norm.sum, 
                                      aes(y = counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot norm") + xlab("Primary disease") + ylab("Norm counts")


# Primary disease and sex

Boxplot_norm_primary_disease_sex <- ggplot(counts.norm.sum, 
                                          aes(y = counts, x=primary_disease_set, 
                                              color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot norm") + xlab("Primary disease") + ylab("Norm counts")



############################################################################
##################               DATOS norm               ##################
##################               Filtrados               ###################
############################################################################


counts.norm.f.sum <- colSums(counts.norm.f)
counts.norm.f.sum <- as.data.frame(counts.norm.f.sum)
colnames(counts.norm.f.sum) <- "counts"
table(rownames(counts.norm.f.sum)==phenoData_training$sample)

counts.norm.f.sum$primary_disease_set <- as.factor(phenoData_training$primary_disease_set)
counts.norm.f.sum$race <- phenoData_training$race
counts.norm.f.sum$race[counts.norm.f.sum$race=="BLACK_OR_AFRICAN_AMERICAN"] <- "BLACK"
counts.norm.f.sum$race <- as.factor(counts.norm.f.sum$race)
counts.norm.f.sum$sex <- as.factor(phenoData_training$sex)
counts.norm.f.sum$age <- phenoData_training$age


str(counts.norm.f.sum)


######
# Sex
######
Boxplot_F_norm_sex <- ggplot(counts.norm.f.sum, aes(y = counts, x=sex, color = sex)) +
  geom_boxplot() +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black"), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot norm filtered") + xlab("Sex") + ylab("Norm counts")

######
# Race
######

Boxplot_F_norm_race <- ggplot(counts.norm.f.sum, aes(y = counts, x=race, color = race)) +
  geom_boxplot() +
  scale_color_manual(values = c("#287D8EFF", "#FDE725FF", "#482677FF")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 8, color="black"), 
        axis.text.x = element_text(size = 6, color="black"), 
        axis.text.y = element_text(size = 6, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.position = "none") +
  ggtitle("Boxplot norm filtered") + xlab("Race") + ylab("Norm counts")


# Primary disease

Boxplot_F_norm_primary_disease <- ggplot(counts.norm.f.sum, 
                                       aes(y = counts, x=primary_disease_set)) +
  geom_boxplot(size=0.4) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  ggtitle("Boxplot norm filtered") + xlab("Primary disease") + ylab("Norm counts")


# Primary disease and sex

Boxplot_F_norm_primary_disease_sex <- ggplot(counts.norm.f.sum, 
                                           aes(y = counts, x=primary_disease_set, 
                                               color = sex)) +
  geom_boxplot(size=0.4) +
  scale_color_manual(values = c("coral1", "deepskyblue3")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 10, color="black"), 
        axis.text.x = element_text(size = 8, color="black", angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 8, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2),
        legend.text = element_text(size = 6)) +
  ggtitle("Boxplot norm filtered") + xlab("Primary disease") + ylab("Norm counts")


################################################################################
################################################################################
#######################          TOTAL BY SEX         ##########################
################################################################################
################################################################################

library(gridExtra)

#Sex
setwd(resultsDirectory)
pdf(file = "Boxplot_comparison_by_sex.pdf", width = 9, height = 4)
grid.arrange(Boxplot_rnaseq_sex, Boxplot_log_rnaseq_sex, Boxplot_vsd_sex, Boxplot_norm_sex,
             Boxplot_F_rnaseq_sex, Boxplot_F_log_rnaseq_sex, Boxplot_F_vsd_sex, Boxplot_F_norm_sex,
             ncol=4)
dev.off()

#Race
setwd(resultsDirectory)
pdf(file = "Boxplot_comparison_by_race.pdf", width = 10, height = 4)
grid.arrange(Boxplot_rnaseq_race, Boxplot_log_rnaseq_race, Boxplot_vsd_race, Boxplot_norm_race,
             Boxplot_F_rnaseq_race, Boxplot_F_log_rnaseq_race, Boxplot_F_vsd_race, Boxplot_F_norm_race,
             ncol=4)
dev.off()

#Primary disease
setwd(resultsDirectory)
pdf(file = "Boxplot_comparison_by_primary_disease.pdf", width = 8, height = 10)
grid.arrange(Boxplot_rnaseq_primary_disease, Boxplot_F_rnaseq_primary_disease,
             Boxplot_log_rnaseq_primary_disease, Boxplot_F_log_rnaseq_primary_disease,
             Boxplot_vsd_primary_disease, Boxplot_F_vsd_primary_disease,
             Boxplot_norm_primary_disease, Boxplot_F_norm_primary_disease,
             ncol=2)
dev.off()

#Primary disease & sex
setwd(resultsDirectory)
pdf(file = "Boxplot_comparison_by_primary_disease_sex.pdf", width = 10, height = 10)
grid.arrange(Boxplot_rnaseq_primary_disease_sex, Boxplot_F_rnaseq_primary_disease_sex,
             Boxplot_log_rnaseq_primary_disease_sex, Boxplot_F_log_rnaseq_primary_disease_sex,
             Boxplot_vsd_primary_disease_sex, Boxplot_F_vsd_primary_disease_sex,
             Boxplot_norm_primary_disease_sex, Boxplot_F_norm_primary_disease_sex,
             ncol=2)
dev.off()

