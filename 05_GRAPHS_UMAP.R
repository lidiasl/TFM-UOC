#####################################################################
# GRAPHS UMAP
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
counts.norm.log <- read.csv(file = "counts_norm_log.csv")
# Filtered
rnaseq_counts_training.f <- read.csv(file = "rnaseq_counts_training_f.csv")
log_rnaseq_training.f <- read.csv(file = "log_rnaseq_training_f.csv")
counts.vsd.f <- read.csv(file = "counts_vsd_f_after.csv")
counts.norm.f <- read.csv(file = "counts_norm_f_after.csv")
counts.norm.f.log <- read.csv(file = "counts_norm_f_after_log.csv")
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

rownames(counts.norm.log) <- counts.norm.log$X
counts.norm.log <- counts.norm.log[, colnames(counts.norm.log)!="X"]

rownames(counts.norm.f.log) <- counts.norm.f.log$X
counts.norm.f.log <- counts.norm.f.log[, colnames(counts.norm.f.log)!="X"]

######################################

# Cargar la paleta de colores turbo
library(viridisLite)
colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)



#######################
#### UMAP
#######################
library(umap)
umap.defaults
#n_neighbors: 15
#min_dist: 0.1

# Bucle listas
n_neighbors_list <- seq(5,100,10)
min_dist_list <- seq(0.1,0.5,0.1)

custom.settings = umap.defaults

############################################################################
##################              DATOS CRUDOS             ###################
##################              sin filtrar              ###################
############################################################################

table(colnames(rnaseq_counts_training)==phenoData_training$sample)

list_UMAP_rnaseq <- list()
k<-1

for(i in 1:length(n_neighbors_list)){
  for(j in 1:length(min_dist_list)){
    custom.settings$n_neighbors = n_neighbors_list[i]
    custom.settings$min_dist=min_dist_list[j]
    
    rnaseq.umap <- umap(t(rnaseq_counts_training), config=custom.settings)
    
    df.umap <- data.frame(x = rnaseq.umap$layout[,1],
                          y = rnaseq.umap$layout[,2],
                          Group = phenoData_training$primary_disease_set)
    
    list_UMAP_rnaseq[[k]] <- ggplot(df.umap, aes(x, y, colour = Group)) +
      geom_point() + ggtitle("UMAP")
    k <- k+1
  }
}

setwd(resultsDirectory)
pdf(file = "UMAP_rnaseq_optimizing.pdf", width = 80, height = 80)
grid.arrange(list_UMAP_rnaseq[[1]],list_UMAP_rnaseq[[2]],list_UMAP_rnaseq[[3]],list_UMAP_rnaseq[[4]],list_UMAP_rnaseq[[5]],
             list_UMAP_rnaseq[[6]],list_UMAP_rnaseq[[7]],list_UMAP_rnaseq[[8]],list_UMAP_rnaseq[[9]], list_UMAP_rnaseq[[10]],
             list_UMAP_rnaseq[[11]],list_UMAP_rnaseq[[12]],list_UMAP_rnaseq[[13]],list_UMAP_rnaseq[[14]],list_UMAP_rnaseq[[15]],
             list_UMAP_rnaseq[[16]],list_UMAP_rnaseq[[17]],list_UMAP_rnaseq[[18]],list_UMAP_rnaseq[[19]],list_UMAP_rnaseq[[20]],
             list_UMAP_rnaseq[[21]],list_UMAP_rnaseq[[22]],list_UMAP_rnaseq[[23]],list_UMAP_rnaseq[[24]],list_UMAP_rnaseq[[25]],
             list_UMAP_rnaseq[[26]],list_UMAP_rnaseq[[27]],list_UMAP_rnaseq[[28]],list_UMAP_rnaseq[[29]], list_UMAP_rnaseq[[30]],
             list_UMAP_rnaseq[[31]],list_UMAP_rnaseq[[32]],list_UMAP_rnaseq[[33]],list_UMAP_rnaseq[[34]],list_UMAP_rnaseq[[35]],
             list_UMAP_rnaseq[[36]],list_UMAP_rnaseq[[37]],list_UMAP_rnaseq[[38]],list_UMAP_rnaseq[[39]], list_UMAP_rnaseq[[40]],
             list_UMAP_rnaseq[[41]],list_UMAP_rnaseq[[42]],list_UMAP_rnaseq[[43]],list_UMAP_rnaseq[[44]],list_UMAP_rnaseq[[45]],
             list_UMAP_rnaseq[[46]],list_UMAP_rnaseq[[47]],list_UMAP_rnaseq[[48]],list_UMAP_rnaseq[[49]], list_UMAP_rnaseq[[50]],
             ncol=5)
dev.off()


custom.settings$n_neighbors = n_neighbors_list[3]
custom.settings$min_dist=min_dist_list[4]

rnaseq.umap <- umap(t(rnaseq_counts_training), config=custom.settings)

df.umap_rnaseq <- data.frame(x = rnaseq.umap$layout[,1],
                      y = rnaseq.umap$layout[,2])

table(rownames(df.umap_rnaseq)==phenoData_training$sample)

df.umap_rnaseq$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.umap_rnaseq$raza <- as.factor(phenoData_training$raza)
df.umap_rnaseq$sexo <- as.factor(phenoData_training$sexo)
df.umap_rnaseq$edad <- phenoData_training$edad


UMAP_rnaseq_sex <- ggplot(df.umap_rnaseq, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Counts") +
  scale_color_manual(values = c("deepskyblue3", "coral1")) + theme_light()+
  labs(colour = "Sexo")

UMAP_rnaseq_age <- ggplot(df.umap_rnaseq, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Counts") + theme_light()+
  labs(colour = "Edad")

UMAP_rnaseq_race <- ggplot(df.umap_rnaseq, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Counts") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF")) + theme_light()+
  labs(colour = "Raza")

# Cargar la paleta de colores turbo
library(viridisLite)
colors_blind <- turbo(27, alpha = 1, begin = 0, end = 1, direction = 1)

UMAP_rnaseq_primary_disease <- ggplot(df.umap_rnaseq, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Counts") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################              DATOS CRUDOS             ###################
##################              Filtrados              ###################
############################################################################

rnaseq.F.umap <- umap(t(rnaseq_counts_training.f), config=custom.settings)

df.F.umap_rnaseq <- data.frame(x = rnaseq.F.umap$layout[,1],
                             y = rnaseq.F.umap$layout[,2])

table(rownames(df.F.umap_rnaseq)==phenoData_training$sample)

df.F.umap_rnaseq$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.F.umap_rnaseq$raza <- as.factor(phenoData_training$raza)
df.F.umap_rnaseq$sexo <- as.factor(phenoData_training$sexo)
df.F.umap_rnaseq$edad <- phenoData_training$edad


UMAP_F_rnaseq_sex <- ggplot(df.F.umap_rnaseq, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Counts filtrados") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_F_rnaseq_age <- ggplot(df.F.umap_rnaseq, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Counts filtrados")+ theme_light()+
  labs(colour = "Edad")

UMAP_F_rnaseq_race <- ggplot(df.F.umap_rnaseq, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Counts filtrados") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_F_rnaseq_primary_disease <- ggplot(df.F.umap_rnaseq, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Counts filtrados") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################               DATOS LOG               ###################
##################              sin filtrar              ###################
############################################################################


log.rnaseq.umap <- umap(t(log_rnaseq_training), config=custom.settings)

df.umap_log_rnaseq <- data.frame(x = log.rnaseq.umap$layout[,1],
                               y = log.rnaseq.umap$layout[,2])

table(rownames(df.umap_log_rnaseq)==phenoData_training$sample)

df.umap_log_rnaseq$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.umap_log_rnaseq$raza <- as.factor(phenoData_training$raza)
df.umap_log_rnaseq$sexo <- as.factor(phenoData_training$sexo)
df.umap_log_rnaseq$edad <- phenoData_training$edad


UMAP_log_rnaseq_sex <- ggplot(df.umap_log_rnaseq, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Log counts") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_log_rnaseq_age <- ggplot(df.umap_log_rnaseq, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Log counts")+ theme_light()+
  labs(colour = "Edad")

UMAP_log_rnaseq_race <- ggplot(df.umap_log_rnaseq, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Log counts") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_log_rnaseq_primary_disease <- ggplot(df.umap_log_rnaseq, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Log counts") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################               DATOS LOG               ###################
##################               Filtrados               ###################
############################################################################


log.F.rnaseq.umap <- umap(t(log_rnaseq_training.f), config=custom.settings)

df.F.umap_log_rnaseq <- data.frame(x = log.F.rnaseq.umap$layout[,1],
                                 y = log.F.rnaseq.umap$layout[,2])


table(rownames(df.F.umap_log_rnaseq)==phenoData_training$sample)

df.F.umap_log_rnaseq$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.F.umap_log_rnaseq$raza <- as.factor(phenoData_training$raza)
df.F.umap_log_rnaseq$sexo <- as.factor(phenoData_training$sexo)
df.F.umap_log_rnaseq$edad <- phenoData_training$edad



UMAP_F_log_rnaseq_sex <- ggplot(df.F.umap_log_rnaseq, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Log counts filtrados") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_F_log_rnaseq_age <- ggplot(df.F.umap_log_rnaseq, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Log counts filtrados")+ theme_light()+
  labs(colour = "Edad")

UMAP_F_log_rnaseq_race <- ggplot(df.F.umap_log_rnaseq, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Log counts filtrados") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_F_log_rnaseq_primary_disease <- ggplot(df.F.umap_log_rnaseq, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Log counts filtrados") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################               DATOS VSD               ###################
##################              sin filtrar              ###################
############################################################################


vsd.umap <- umap(t(counts.vsd), config=custom.settings)

df.umap_vsd <- data.frame(x = vsd.umap$layout[,1], 
                          y = vsd.umap$layout[,2])

table(rownames(df.umap_vsd)==phenoData_training$sample)

df.umap_vsd$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.umap_vsd$raza <- as.factor(phenoData_training$raza)
df.umap_vsd$sexo <- as.factor(phenoData_training$sexo)
df.umap_vsd$edad <- phenoData_training$edad


UMAP_vsd_sex <- ggplot(df.umap_vsd, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Vsd") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_vsd_age <- ggplot(df.umap_vsd, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Vsd")+ theme_light()+
  labs(colour = "Edad")

UMAP_vsd_race <- ggplot(df.umap_vsd, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Vsd") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_vsd_primary_disease <- ggplot(df.umap_vsd, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Vsd") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################               DATOS VSD               ###################
##################               Filtrados               ###################
############################################################################


vsd.F.umap <- umap(t(counts.vsd.f), config=custom.settings)

df.F.umap_vsd <- data.frame(x = vsd.F.umap$layout[,1], 
                            y = vsd.F.umap$layout[,2])

table(rownames(df.F.umap_vsd)==phenoData_training$sample)

df.F.umap_vsd$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.F.umap_vsd$raza <- as.factor(phenoData_training$raza)
df.F.umap_vsd$sexo <- as.factor(phenoData_training$sexo)
df.F.umap_vsd$edad <- phenoData_training$edad


UMAP_F_vsd_sex <- ggplot(df.F.umap_vsd, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Vsd filtrados") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_F_vsd_age <- ggplot(df.F.umap_vsd, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Vsd filtrados")+ theme_light()+
  labs(colour = "Edad")

UMAP_F_vsd_race <- ggplot(df.F.umap_vsd, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Vsd filtrados") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_F_vsd_primary_disease <- ggplot(df.F.umap_vsd, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Vsd filtrados") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################               DATOS norm               ###################
##################              sin filtrar              ###################
############################################################################


norm.umap <- umap(t(counts.norm), config=custom.settings)

df.umap_norm <- data.frame(x = norm.umap$layout[,1], 
                            y = norm.umap$layout[,2])


table(rownames(df.umap_norm)==phenoData_training$sample)

df.umap_norm$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.umap_norm$raza <- as.factor(phenoData_training$raza)
df.umap_norm$sexo <- as.factor(phenoData_training$sexo)
df.umap_norm$edad <- phenoData_training$edad


UMAP_norm_sex <- ggplot(df.umap_norm, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Normalizados") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_norm_age <- ggplot(df.umap_norm, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Normalizados")+ theme_light()+
  labs(colour = "Edad")

UMAP_norm_race <- ggplot(df.umap_norm, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Normalizados") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_norm_primary_disease <- ggplot(df.umap_norm, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Normalizados") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")



############################################################################
##################               DATOS norm               ##################
##################               Filtrados               ###################
############################################################################


norm.F.umap <- umap(t(counts.norm.f), config=custom.settings)

df.F.umap_norm <- data.frame(x = norm.F.umap$layout[,1], 
                           y = norm.F.umap$layout[,2])

table(rownames(df.F.umap_norm)==phenoData_training$sample)

df.F.umap_norm$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.F.umap_norm$raza <- as.factor(phenoData_training$raza)
df.F.umap_norm$sexo <- as.factor(phenoData_training$sexo)
df.F.umap_norm$edad <- phenoData_training$edad


UMAP_F_norm_sex <- ggplot(df.F.umap_norm, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Normalizados filtrados") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_F_norm_age <- ggplot(df.F.umap_norm, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Normalizados filtrados")+ theme_light()+
  labs(colour = "Edad")

UMAP_F_norm_race <- ggplot(df.F.umap_norm, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Normalizados filtrados") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_F_norm_primary_disease <- ggplot(df.F.umap_norm, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Normalizados filtrados") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################            DATOS norm log             ###################
##################              sin filtrar              ###################
############################################################################


norm.log.umap <- umap(t(counts.norm.log), config=custom.settings)

df.umap_norm_log <- data.frame(x = norm.log.umap$layout[,1], 
                             y = norm.log.umap$layout[,2])


table(rownames(df.umap_norm_log)==phenoData_training$sample)

df.umap_norm_log$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.umap_norm_log$raza <- as.factor(phenoData_training$raza)
df.umap_norm_log$sexo <- as.factor(phenoData_training$sexo)
df.umap_norm_log$edad <- phenoData_training$edad

UMAP_norm_log_sex <- ggplot(df.umap_norm_log, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Normalizados log") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_norm_log_age <- ggplot(df.umap_norm_log, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Normalizados log")+ theme_light()+
  labs(colour = "Edad")

UMAP_norm_log_race <- ggplot(df.umap_norm_log, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Normalizados log") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_norm_log_primary_disease <- ggplot(df.umap_norm_log, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Normalizados log") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")


############################################################################
##################            DATOS norm log             ###################
##################              Filtraados               ###################
############################################################################

norm.log.F.umap <- umap(t(counts.norm.f.log), config=custom.settings)

df.F.umap_norm_log <- data.frame(x = norm.log.F.umap$layout[,1], 
                               y = norm.log.F.umap$layout[,2])


table(rownames(df.F.umap_norm_log)==phenoData_training$sample)

df.F.umap_norm_log$enfermedad_primaria <- phenoData_training$enfermedad_primaria
df.F.umap_norm_log$raza <- as.factor(phenoData_training$raza)
df.F.umap_norm_log$sexo <- as.factor(phenoData_training$sexo)
df.F.umap_norm_log$edad <- phenoData_training$edad


UMAP_F_norm_log_sex <- ggplot(df.F.umap_norm_log, aes(x, y, colour = sexo)) +
  geom_point(size=0.3) + ggtitle("Normalizados log filtrados") +
  scale_color_manual(values = c("deepskyblue3", "coral1"))+ theme_light()+
  labs(colour = "Sexo")

UMAP_F_norm_log_age <- ggplot(df.F.umap_norm_log, aes(x, y, colour = edad)) +
  geom_point(size=0.3) + ggtitle("Normalizados log filtrados")+ theme_light()+
  labs(colour = "Edad")

UMAP_F_norm_log_race <- ggplot(df.F.umap_norm_log, aes(x, y, colour = raza)) +
  geom_point(size=0.3) + ggtitle("Normalizados log filtrados") +
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+ theme_light()+
  labs(colour = "Raza")

UMAP_F_norm_log_primary_disease <- ggplot(df.F.umap_norm_log, aes(x, y, colour = enfermedad_primaria)) +
  geom_point(size=0.3) + ggtitle("Normalizados log filtrados") +
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")



################################################################################
################################################################################
#######################              TOTAL            ##########################
################################################################################
################################################################################

library(gridExtra)

#Sex
setwd(resultsDirectory)
pdf(file = "UMAP_comparison_by_sex.pdf", width = 8, height = 10)
grid.arrange(UMAP_rnaseq_sex, UMAP_F_rnaseq_sex,
             UMAP_log_rnaseq_sex, UMAP_F_log_rnaseq_sex,
             UMAP_vsd_sex, UMAP_F_vsd_sex,
             UMAP_norm_sex, UMAP_F_norm_sex,
             UMAP_norm_log_sex, UMAP_F_norm_log_sex,
             ncol=2)
dev.off()

#Race
setwd(resultsDirectory)
pdf(file = "UMAP_comparison_by_race.pdf", width = 8, height = 10)
grid.arrange(UMAP_rnaseq_race, UMAP_F_rnaseq_race,
             UMAP_log_rnaseq_race, UMAP_F_log_rnaseq_race,
             UMAP_vsd_race, UMAP_F_vsd_race,
             UMAP_norm_race, UMAP_F_norm_race,
             UMAP_norm_log_race, UMAP_F_norm_log_race,
             ncol=2)
dev.off()

#Age
setwd(resultsDirectory)
pdf(file = "UMAP_comparison_by_age.pdf", width = 8, height = 10)
grid.arrange(UMAP_rnaseq_age, UMAP_F_rnaseq_age,
             UMAP_log_rnaseq_age, UMAP_F_log_rnaseq_age,
             UMAP_vsd_age, UMAP_F_vsd_age,
             UMAP_norm_age, UMAP_F_norm_age,
             UMAP_norm_log_age, UMAP_F_norm_log_age,
             ncol=2)
dev.off()

#Primary disease
setwd(resultsDirectory)
pdf(file = "UMAP_comparison_by_primary_disease.pdf", width = 8, height = 10)
grid.arrange(UMAP_rnaseq_primary_disease, UMAP_F_rnaseq_primary_disease,
             UMAP_log_rnaseq_primary_disease, UMAP_F_log_rnaseq_primary_disease,
             UMAP_vsd_primary_disease, UMAP_F_vsd_primary_disease,
             UMAP_norm_primary_disease, UMAP_F_norm_primary_disease,
             UMAP_norm_log_primary_disease, UMAP_F_norm_log_primary_disease,
             ncol=2)
dev.off()






