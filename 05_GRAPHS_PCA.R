#####################################################################
# GRAPHS PCA
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


############################################################################
##################              DATOS CRUDOS             ###################
##################              sin filtrar              ###################
############################################################################


PCA_rnaseq <- prcomp(t(as.matrix(rnaseq_counts_training)), 
                                      scale. = FALSE, 
                                      center = FALSE)


PCA_rnaseq_sex <- autoplot(PCA_rnaseq, data = phenoData_training, colour = "sexo",
                                       main = "Counts", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_rnaseq_race <- autoplot(PCA_rnaseq, data = phenoData_training, colour = "raza",
                           main = "Counts", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_rnaseq_age <- autoplot(PCA_rnaseq, data = phenoData_training, colour = "edad",
                            main = "Counts", size=0.4)+
  labs(colour = "Edad") + theme_light()

PCA_rnaseq_primary_disease <- autoplot(PCA_rnaseq, data = phenoData_training, colour = "enfermedad_primaria",
                                       main = "Counts", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")

PCA_rnaseq_primary_disease_noun <- autoplot(PCA_rnaseq, data = phenoData_training, colour = "enfermedad_primaria",
                                       main = "Counts", label=TRUE)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light() + theme(legend.position = "none")
# Muestra de vejiga:TCGA_CF_A3MF_01

############################################################################
##################              DATOS CRUDOS             ###################
##################              Filtrados              ###################
############################################################################


PCA_F_rnaseq <- prcomp(t(as.matrix(rnaseq_counts_training.f)), 
                     scale. = FALSE, 
                     center = FALSE)


PCA_F_rnaseq_sex <- autoplot(PCA_F_rnaseq, data = phenoData_training, colour = "sexo",
                           main = "Counts filtrados", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_F_rnaseq_race <- autoplot(PCA_F_rnaseq, data = phenoData_training, colour = "raza",
                            main = "Counts filtrados", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_F_rnaseq_age <- autoplot(PCA_F_rnaseq, data = phenoData_training, colour = "edad",
                           main = "Counts filtrados", size=0.4)+
  labs(colour = "Edad") + theme_light()

PCA_F_rnaseq_primary_disease <- autoplot(PCA_F_rnaseq, data = phenoData_training, colour = "enfermedad_primaria",
                                       main = "Counts filtrados", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")


############################################################################
##################               DATOS LOG               ###################
##################              sin filtrar              ###################
############################################################################

PCA_log_rnaseq <- prcomp(t(as.matrix(log_rnaseq_training)), 
                     scale. = FALSE, 
                     center = FALSE)


PCA_log_rnaseq_sex <- autoplot(PCA_log_rnaseq, data = phenoData_training, colour = "sexo",
                           main = "Log counts", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_log_rnaseq_race <- autoplot(PCA_log_rnaseq, data = phenoData_training, colour = "raza",
                            main = "Log counts", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_log_rnaseq_age <- autoplot(PCA_log_rnaseq, data = phenoData_training, colour = "edad",
                           main = "Log counts", size=0.4)+
  labs(colour = "Edad") + theme_light()

PCA_log_rnaseq_primary_disease <- autoplot(PCA_log_rnaseq, data = phenoData_training, colour = "enfermedad_primaria",
                                       main = "Log counts", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")
 

############################################################################
##################               DATOS LOG               ###################
##################               Filtrados               ###################
############################################################################


PCA_F_log_rnaseq <- prcomp(t(as.matrix(log_rnaseq_training.f)), 
                         scale. = FALSE, 
                         center = FALSE)


PCA_F_log_rnaseq_sex <- autoplot(PCA_F_log_rnaseq, data = phenoData_training, colour = "sexo",
                               main = "Log counts filtrados", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_F_log_rnaseq_race <- autoplot(PCA_F_log_rnaseq, data = phenoData_training, colour = "raza",
                                main = "Log counts filtrados", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_F_log_rnaseq_age <- autoplot(PCA_F_log_rnaseq, data = phenoData_training, colour = "edad",
                               main = "Log counts filtrados", size=0.4)+
  labs(colour = "Edad") + theme_light()

PCA_F_log_rnaseq_primary_disease <- autoplot(PCA_F_log_rnaseq, data = phenoData_training, colour = "enfermedad_primaria",
                                           main = "Log counts filtrados", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")


############################################################################
##################               DATOS VSD               ###################
##################              sin filtrar              ###################
############################################################################


PCA_vsd <- prcomp(t(as.matrix(counts.vsd)), 
                           scale. = FALSE, 
                           center = FALSE)


PCA_vsd_sex <- autoplot(PCA_vsd, data = phenoData_training, colour = "sexo",
                               main = "Vsd", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_vsd_race <- autoplot(PCA_vsd, data = phenoData_training, colour = "raza",
                                main = "Vsd", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_vsd_age <- autoplot(PCA_vsd, data = phenoData_training, colour = "edad",
                               main = "Vsd", size=0.4)+
  labs(colour = "Edad") + theme_light()


PCA_vsd_primary_disease <- autoplot(PCA_vsd, data = phenoData_training, colour = "enfermedad_primaria",
                                           main = "Vsd", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")



############################################################################
##################               DATOS VSD               ###################
##################               Filtrados               ###################
############################################################################

PCA_F_vsd <- prcomp(t(as.matrix(counts.vsd.f)), 
                  scale. = FALSE, 
                  center = FALSE)


PCA_F_vsd_sex <- autoplot(PCA_F_vsd, data = phenoData_training, colour = "sexo",
                        main = "Vsd filtrados", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_F_vsd_race <- autoplot(PCA_F_vsd, data = phenoData_training, colour = "raza",
                         main = "Vsd filtrados", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_F_vsd_age <- autoplot(PCA_F_vsd, data = phenoData_training, colour = "edad",
                        main = "Vsd filtrados", size=0.4)+
  labs(colour = "Edad") + theme_light()


PCA_F_vsd_primary_disease <- autoplot(PCA_F_vsd, data = phenoData_training, colour = "enfermedad_primaria",
                                    main = "Vsd filtrados", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")


############################################################################
##################               DATOS norm               ###################
##################              sin filtrar              ###################
############################################################################


PCA_norm <- prcomp(t(as.matrix(counts.norm)), 
                  scale. = FALSE, 
                  center = FALSE)


PCA_norm_sex <- autoplot(PCA_norm, data = phenoData_training, colour = "sexo",
                        main = "Normalizados", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_norm_race <- autoplot(PCA_norm, data = phenoData_training, colour = "raza",
                         main = "Normalizados", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_norm_age <- autoplot(PCA_norm, data = phenoData_training, colour = "edad",
                        main = "Normalizados", size=0.4)+
  labs(colour = "Edad") + theme_light()


PCA_norm_primary_disease <- autoplot(PCA_norm, data = phenoData_training, colour = "enfermedad_primaria",
                                    main = "Normalizados", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")


############################################################################
##################               DATOS norm               ##################
##################               Filtrados               ###################
############################################################################


PCA_F_norm <- prcomp(t(as.matrix(counts.norm.f)), 
                   scale. = FALSE, 
                   center = FALSE)


PCA_F_norm_sex <- autoplot(PCA_F_norm, data = phenoData_training, colour = "sexo",
                         main = "Normalizados filtrados", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_F_norm_race <- autoplot(PCA_F_norm, data = phenoData_training, colour = "raza",
                          main = "Normalizados filtrados", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_F_norm_age <- autoplot(PCA_F_norm, data = phenoData_training, colour = "edad",
                         main = "Normalizados filtrados", size=0.4)+
  labs(colour = "Edad") + theme_light()


PCA_F_norm_primary_disease <- autoplot(PCA_F_norm, data = phenoData_training, colour = "enfermedad_primaria",
                                     main = "Normalizados filtrados", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")



############################################################################
##################            DATOS norm log             ###################
##################              sin filtrar              ###################
############################################################################

PCA_norm_log <- prcomp(t(as.matrix(counts.norm.log)), 
                   scale. = FALSE, 
                   center = FALSE)


PCA_norm_log_sex <- autoplot(PCA_norm_log, data = phenoData_training, colour = "sexo",
                         main = "Normalizados log", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_norm_log_race <- autoplot(PCA_norm_log, data = phenoData_training, colour = "raza",
                          main = "Normalizados log", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_norm_log_age <- autoplot(PCA_norm_log, data = phenoData_training, colour = "edad",
                         main = "Normalizados log", size=0.4)+
  labs(colour = "Edad") + theme_light()


PCA_norm_log_primary_disease <- autoplot(PCA_norm_log, data = phenoData_training, colour = "enfermedad_primaria",
                                     main = "Normalizados log", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")



############################################################################
##################            DATOS norm log             ###################
##################              Filtraados               ###################
############################################################################


PCA_F_norm_log <- prcomp(t(as.matrix(counts.norm.f.log)), 
                       scale. = FALSE, 
                       center = FALSE)


PCA_F_norm_log_sex <- autoplot(PCA_F_norm_log, data = phenoData_training, colour = "sexo",
                             main = "Normalizados log filtrados", size=0.4)+
  scale_color_manual(values = c("deepskyblue3", "coral1"))+
  labs(colour = "Sexo") + theme_light()


PCA_F_norm_log_race <- autoplot(PCA_F_norm_log, data = phenoData_training, colour = "raza",
                              main = "Normalizados log filtrados", size=0.4)+
  scale_color_manual(values = c("#287D8EFF", "#482677FF", "#FDE725FF"))+
  labs(colour = "Raza") + theme_light()


PCA_F_norm_log_age <- autoplot(PCA_F_norm_log, data = phenoData_training, colour = "edad",
                             main = "Normalizados log filtrados", size=0.4)+
  labs(colour = "Edad") + theme_light()


PCA_F_norm_log_primary_disease <- autoplot(PCA_F_norm_log, data = phenoData_training, colour = "enfermedad_primaria",
                                         main = "Normalizados log filtrados", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()+ theme(legend.position = "none")



################################################################################
################################################################################
#######################              TOTAL            ##########################
################################################################################
################################################################################

library(gridExtra)

#Sex
setwd(resultsDirectory)
pdf(file = "PCA_comparison_by_sex.pdf", width = 8, height = 10)
grid.arrange(PCA_rnaseq_sex, PCA_F_rnaseq_sex,
             PCA_log_rnaseq_sex, PCA_F_log_rnaseq_sex,
             PCA_vsd_sex, PCA_F_vsd_sex,
             PCA_norm_sex, PCA_F_norm_sex,
             PCA_norm_log_sex, PCA_F_norm_log_sex,
             ncol=2)
dev.off()

#Race
setwd(resultsDirectory)
pdf(file = "PCA_comparison_by_race.pdf", width = 8, height = 10)
grid.arrange(PCA_rnaseq_race, PCA_F_rnaseq_race,
             PCA_log_rnaseq_race, PCA_F_log_rnaseq_race,
             PCA_vsd_race, PCA_F_vsd_race,
             PCA_norm_race, PCA_F_norm_race,
             PCA_norm_log_race, PCA_F_norm_log_race,
             ncol=2)
dev.off()

#Age
setwd(resultsDirectory)
pdf(file = "PCA_comparison_by_age.pdf", width = 8, height = 10)
grid.arrange(PCA_rnaseq_age, PCA_F_rnaseq_age,
             PCA_log_rnaseq_age, PCA_F_log_rnaseq_age,
             PCA_vsd_age, PCA_F_vsd_age,
             PCA_norm_age, PCA_F_norm_age,
             PCA_norm_log_age, PCA_F_norm_log_age,
             ncol=2)
dev.off()

#Primary disease
setwd(resultsDirectory)
pdf(file = "PCA_comparison_by_primary_disease.pdf", width = 8, height = 10)
grid.arrange(PCA_rnaseq_primary_disease, PCA_F_rnaseq_primary_disease,
             PCA_log_rnaseq_primary_disease, PCA_F_log_rnaseq_primary_disease,
             PCA_vsd_primary_disease, PCA_F_vsd_primary_disease,
             PCA_norm_primary_disease, PCA_F_norm_primary_disease,
             PCA_norm_log_primary_disease, PCA_F_norm_log_primary_disease,
             ncol=2)
dev.off()


setwd(resultsDirectory)
pdf(file = "foo_legend.pdf", width = 8, height = 10)
autoplot(PCA_F_norm_log, data = phenoData_training, colour = "enfermedad_primaria",
         main = "Normalizados log filtrados", size=0.5)+
  scale_color_manual(values = colors_blind)+
  labs(colour = "Tumor primario") + theme_light()
dev.off()
