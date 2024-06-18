#####################################################################
# Samples study and  representation
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

tablesDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Tables"
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

#####################################################################
#### LOAD TCGA DATA
#####################################################################

#### Load the phenotypic dataset:
phenoData_TCGA <- fread("~/Desktop/LIDIA/TCGA_classification/Data/TCGA_phenotype_denseDataOnlyDownload.tsv", sep = "\t", stringsAsFactors = F)
colnames(phenoData_TCGA)[4] <- "primary_disease"

phenoData_TCGA$sample <- str_replace_all(phenoData_TCGA$sample, "-", "_")
phenoData_TCGA$sample_type <- str_replace_all(phenoData_TCGA$sample_type, " ", "_")
phenoData_TCGA$sample_type <- str_replace_all(phenoData_TCGA$sample_type, "-", "_")
phenoData_TCGA$sample_type <- str_replace_all(phenoData_TCGA$sample_type, "___", "_")
phenoData_TCGA$primary_disease <- str_replace_all(phenoData_TCGA$primary_disease, " ", "_")
phenoData_TCGA$primary_disease <- str_replace_all(phenoData_TCGA$primary_disease, "-", "_")

# Table of sample_type
table(phenoData_TCGA$sample_type)

# Complete cases
table(complete.cases(phenoData_TCGA)==TRUE)
#FALSE  TRUE 
#72 12732 

phenoData_TCGA$sample_type[complete.cases(phenoData_TCGA)==FALSE] <- "NA"
phenoData_TCGA$sample_type_id[complete.cases(phenoData_TCGA)==FALSE] <- "NA"


#### Phenodata with gender (and other variables):
phenoData_TCGA_gender <- fread("~/Desktop/LIDIA/TCGA_classification/Data/Survival_SupplementalTable_S1_20171025_xena_sp.txt", sep = "\t")
phenoData_TCGA_gender$sample <- str_replace_all(phenoData_TCGA_gender$sample, "-", "_")


table(phenoData_TCGA$sample %in% phenoData_TCGA_gender$sample)
#FALSE  TRUE 
#213 12591 

table(phenoData_TCGA_gender$sample %in% phenoData_TCGA$sample)
#TRUE 
#12591 

### Tissue source site codes

table_tissue_source_site_codes <- fread("~/Desktop/LIDIA/TCGA_classification/Data/tissue_source_site_codes.csv", sep = ",", stringsAsFactors = F)

df_source_site <- as.data.frame(table(table_tissue_source_site_codes$`Source Site`))

write.table(df_source_site, file="~/Desktop/LIDIA/TCGA_classification/Data/tissue_source_site_df_modified.csv")

#### Load the gene expression dataset: (unit log2(expected_count+1))
rnaseq <- fread("~/Desktop/LIDIA/TCGA_classification/Data/tcga_gene_expected_count", sep = "\t", stringsAsFactors = F)
genenames <- separate(data = rnaseq, col = sample, into = c("GeneName", "version"), remove = T)
rnaseq <- rnaseq[,-1]
rnaseq <- as.data.frame(rnaseq)
rownames(rnaseq) <- genenames$GeneName
colnames(rnaseq) <- str_replace_all(colnames(rnaseq), "-", "_")

# Samples phenodata_TCGA in rnaseq
table(phenoData_TCGA$sample %in% colnames(rnaseq))
#FALSE  TRUE 
#2275 10529 

# Samples with gender in rnaseq
table(phenoData_TCGA_gender$sample %in% colnames(rnaseq))
#FALSE  TRUE 
#2100 10491

####################
# Filtering 
####################
phenoData_TCGA <- as.data.frame(phenoData_TCGA)
phenoData_TCGA_in_rnaseq <- phenoData_TCGA[phenoData_TCGA$sample %in% colnames(rnaseq),]
phenoData_TCGA_in_rnaseq_with_gender <- phenoData_TCGA_in_rnaseq[phenoData_TCGA_in_rnaseq$sample %in% phenoData_TCGA_gender$sample,]

######################################
# Write phenoData_TCGA
######################################
setwd(tablesDirectory)
write.csv(phenoData_TCGA, file = "phenoData_TCGA.csv")
write.csv(phenoData_TCGA_in_rnaseq_with_gender, file = "phenoData_TCGA_in_rnaseq_with_gender.csv")

#####################################################################
#### GRAPHICS PRE-ANALYSIS
#####################################################################


# Counts by sample type (raw data)
phenoData_TCGA$sample_type <- factor(phenoData_TCGA$sample_type,
                                         levels = names(sort(table(phenoData_TCGA$sample_type), decreasing = FALSE)))

Counts_by_sample_type_raw_data <- ggplot(phenoData_TCGA, aes(sample_type)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Counts by sample type (raw data)") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=1.2) +
  theme(plot.background = element_rect(fill = "white"), 
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 6, color="black"), 
    axis.text.x = element_text(size = 4, color="black"), 
    axis.text.y = element_text(size = 3.5, color="black")) +
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by sample type raw data.png", plot = Counts_by_sample_type_raw_data, width = 4, height = 2, units = "in")


# Counts by sample type by rnaseq
phenoData_TCGA$rnaseq <- ifelse(phenoData_TCGA$sample %in% colnames(rnaseq), "in", "out")

Counts_by_sample_type_raw_data_in_rnaseq <- ggplot(phenoData_TCGA, aes(sample_type, fill=rnaseq, color=rnaseq)) + geom_bar(alpha=0.4) + 
  ggtitle("Counts by sample type in rnaseq") + ylab("Counts") + xlab("Sample type") +  
  geom_bar(alpha = 0.4, show.legend = FALSE) + 
  theme_minimal() +
  scale_fill_manual(values = c("darkgreen","red")) +
  scale_color_manual(values = c("darkgreen","red")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 6, color="black"), 
        axis.text.x = element_text(size = 4, color="black"), 
        axis.text.y = element_text(size = 4, color="black"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by sample type raw data in rnaseq.png", plot = Counts_by_sample_type_raw_data_in_rnaseq, width = 4, height = 2, units = "in")



# Counts by sample type in phenoData_TCGA_gender
phenoData_TCGA$gender <- ifelse(phenoData_TCGA$sample %in% phenoData_TCGA_gender$sample, "yes", "no")

Counts_by_sample_type_raw_data_with_gender <- ggplot(phenoData_TCGA, aes(sample_type, fill=gender, color=gender)) + geom_bar(alpha=0.4) + 
  ggtitle("Counts by sample type with gender") + ylab("Counts") + xlab("Sample type") +  
  geom_bar(alpha = 0.4, show.legend = FALSE) + 
  theme_minimal() +
  scale_fill_manual(values = c("red","darkgreen")) +
  scale_color_manual(values = c("red","darkgreen")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 6, color="black"), 
        axis.text.x = element_text(size = 4, color="black"), 
        axis.text.y = element_text(size = 4, color="black"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by sample type raw data with gender.png", plot = Counts_by_sample_type_raw_data_with_gender, width = 4, height = 2, units = "in")


# Counts by primary disease

phenoData_TCGA$primary_disease <- factor(phenoData_TCGA$primary_disease,
                                         levels = names(sort(table(phenoData_TCGA$primary_disease), decreasing = FALSE)))

Counts_by_primary_disease_raw_data <- ggplot(phenoData_TCGA, aes(primary_disease)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Counts by primary disease (raw data)") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease raw data.png", plot = Counts_by_primary_disease_raw_data, width = 8, height = 5, units = "in")


# Counts by primary_disease in rnaseq with gender

phenoData_TCGA_in_rnaseq_with_gender$primary_disease <- factor(phenoData_TCGA_in_rnaseq_with_gender$primary_disease,
                                                   levels = names(sort(table(phenoData_TCGA_in_rnaseq_with_gender$primary_disease), 
                                                                       decreasing = FALSE)))

Counts_by_primary_disease_in_rnaseq_with_gender <- ggplot(phenoData_TCGA_in_rnaseq_with_gender, aes(primary_disease)) + 
  geom_bar(color="blue", fill="blue", alpha=0.2) +
  ggtitle("Counts by primary disease (in rnaseq with gender)") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease in rnaseq with gender.png", plot = Counts_by_primary_disease_in_rnaseq_with_gender, width = 8, 
       height = 5, units = "in")


#####################################################################
#### Separate primary and metastatic
#####################################################################

phenoData_TCGA.f <- phenoData_TCGA_in_rnaseq_with_gender

# Filter table by primary tumor sites:
phenoData_TCGA_primary <- phenoData_TCGA.f[which(phenoData_TCGA.f$sample_type == "Primary_Tumor" | phenoData_TCGA.f$sample_type == "Primary_Blood_Derived_Cancer_Peripheral_Blood"),]
phenoData_TCGA_metastatic <- phenoData_TCGA.f[which(phenoData_TCGA.f$sample_type == "Metastatic"),]


# Table of primary_disease TCGA PRIMARY in rnaseq
phenoData_TCGA_primary <- as.data.frame(phenoData_TCGA_primary)
#Check
table(phenoData_TCGA_primary$sample %in% colnames(rnaseq))
#TRUE 
#9332 


################
# Add variables
################
setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

# PRIMARY
phenoData_TCGA_primary$sex <- phenoData_TCGA_gender$gender[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$race <- phenoData_TCGA_gender$race[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$tumor_stage <- phenoData_TCGA_gender$ajcc_pathologic_tumor_stage[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$histological_type <- phenoData_TCGA_gender$histological_type[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$vital_status <- phenoData_TCGA_gender$vital_status[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$age <- phenoData_TCGA_gender$age_at_initial_pathologic_diagnosis[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$primary_disease_set <- set_subset$set[match(phenoData_TCGA_primary$primary_disease, set_subset$subset)]
phenoData_TCGA_primary$TSS <- sapply(phenoData_TCGA_primary$sample, function(x){
  str_split(x, "_")[[1]][2]
})
phenoData_TCGA_primary$BCR <- table_tissue_source_site_codes$BCR[match(phenoData_TCGA_primary$TSS, table_tissue_source_site_codes$`TSS Code`)]
phenoData_TCGA_primary$source_site <- table_tissue_source_site_codes$`Source Site`[match(phenoData_TCGA_primary$TSS, table_tissue_source_site_codes$`TSS Code`)]
phenoData_TCGA_primary$sexo <- ifelse(phenoData_TCGA_primary$sex=="MALE", "HOMBRE", "MUJER")
phenoData_TCGA_primary$race <- as.factor(phenoData_TCGA_primary$race)
phenoData_TCGA_primary$raza <- factor(phenoData_TCGA_primary$race, 
                                      levels=c("", "[Not Evaluated]", "[Unknown]", "AMERICAN INDIAN OR ALASKA NATIVE", "ASIAN", "BLACK OR AFRICAN AMERICAN",
                                               "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "WHITE"),
                                      labels= c("", "[No Evaluado]", "[Desconocido]", 
                                                "INDIO AMERICANO O NATIVO DE ALASKA", "ASIÁTICO", 
                                                "NEGRO O AFROAMERICANO", "HAWAIANO NATIVO U OTRO ISLEÑO DEL PACÍFICO",
                                                "BLANCO"))
phenoData_TCGA_primary$edad <- phenoData_TCGA_primary$age
phenoData_TCGA_primary$raza <- str_replace_all(phenoData_TCGA_primary$raza, " ", "_")
phenoData_TCGA_primary$enfermedad_primaria <- set_subset$spanish[match(phenoData_TCGA_primary$primary_disease_set, set_subset$set)]



# METASTATIC
phenoData_TCGA_metastatic$sex <- phenoData_TCGA_gender$gender[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$race <- phenoData_TCGA_gender$race[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$tumor_stage <- phenoData_TCGA_gender$ajcc_pathologic_tumor_stage[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$histological_type <- phenoData_TCGA_gender$histological_type[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$vital_status <- phenoData_TCGA_gender$vital_status[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$age <- phenoData_TCGA_gender$age_at_initial_pathologic_diagnosis[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$primary_disease_set <- set_subset$set[match(phenoData_TCGA_metastatic$primary_disease, set_subset$subset)]
phenoData_TCGA_metastatic$TSS <- sapply(phenoData_TCGA_metastatic$sample, function(x){
  str_split(x, "_")[[1]][2]
})
phenoData_TCGA_metastatic$BCR <- table_tissue_source_site_codes$BCR[match(phenoData_TCGA_metastatic$TSS, table_tissue_source_site_codes$`TSS Code`)]
phenoData_TCGA_metastatic$source_site <- table_tissue_source_site_codes$`Source Site`[match(phenoData_TCGA_metastatic$TSS, table_tissue_source_site_codes$`TSS Code`)]
phenoData_TCGA_metastatic$sexo <- ifelse(phenoData_TCGA_metastatic$sex=="MALE", "HOMBRE", "MUJER")
phenoData_TCGA_metastatic$race <- as.factor(phenoData_TCGA_metastatic$race)
phenoData_TCGA_metastatic$raza <- factor(phenoData_TCGA_metastatic$race, 
                                      levels=c("", "[Not Evaluated]", "[Unknown]", "ASIAN", "BLACK OR AFRICAN AMERICAN", "WHITE"),
                                      labels= c("", "[No Evaluado]", "[Desconocido]", "ASIÁTICO", "NEGRO O AFROAMERICANO", "BLANCO"))
phenoData_TCGA_metastatic$edad <- phenoData_TCGA_metastatic$age
phenoData_TCGA_metastatic$raza <- str_replace_all(phenoData_TCGA_metastatic$raza, " ", "_")
phenoData_TCGA_metastatic$enfermedad_primaria <- set_subset$spanish[match(phenoData_TCGA_metastatic$primary_disease_set, set_subset$set)]

######################################
# Write phenoData_TCGA_primary and phenoData_TCGA_metastatic
######################################
setwd(tablesDirectory)
write.csv(phenoData_TCGA_primary, file = "phenoData_TCGA_primary.csv")
write.csv(phenoData_TCGA_metastatic, file = "phenoData_TCGA_metastatic.csv")


#####################################################################
# Graphs
#####################################################################

library(lessR)
setwd(resultsDirectory)
phenoData_TCGA_primary_graph <- phenoData_TCGA_primary[phenoData_TCGA_primary$raza %in% races_list,]
pdf(file="PieChart_razas.pdf", width=8, height = 8)
PieChart(raza, hole = 0, data = phenoData_TCGA_primary_graph,
         fill = c("#287D8EFF", "#482677FF", "#808080", "#FF0000","#FDE725FF"), main = "")
dev.off()
# Counts by primary disease (Primary T)
phenoData_TCGA_primary$enfermedad_primaria <- factor(phenoData_TCGA_primary$enfermedad_primaria,
                                                     levels = names(sort(table(phenoData_TCGA_primary$enfermedad_primaria), decreasing = FALSE)))

Counts_by_primary_disease_set_PT <- ggplot(phenoData_TCGA_primary, aes(enfermedad_primaria)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Muestras por tumores primarios (conjunto Tumores primarios)") + ylab("Número de muestras") + xlab("Tumor primario") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 12))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Muestras por tumores primarios conjunto PT.png", plot = Counts_by_primary_disease_set_PT, width = 8, height = 5, units = "in")


setwd(resultsDirectory)
pdf(file = "Muestras por tumores primarios conjunto PT.pdf", width = 6, height = 4)
print(Counts_by_primary_disease_set_PT)
dev.off()

# Counts by primary disease (Metastatic T)
phenoData_TCGA_metastatic$enfermedad_primaria <- factor(phenoData_TCGA_metastatic$enfermedad_primaria,
                                                     levels = names(sort(table(phenoData_TCGA_metastatic$enfermedad_primaria), decreasing = FALSE)))

Counts_by_primary_disease_set_MetT <- ggplot(phenoData_TCGA_metastatic, aes(enfermedad_primaria)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Muestras por tumores primarios (conjunto Metástasis)") + ylab("Número de muestras") + xlab("Tumor primario") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 12))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Muestras por tumores primarios conjunto MetT.png", plot = Counts_by_primary_disease_set_MetT, width = 8, height = 4, units = "in")


setwd(resultsDirectory)
pdf(file = "Muestras por tumores primarios conjunto MetT.pdf", width = 6, height = 4)
print(Counts_by_primary_disease_set_MetT)
dev.off()

######################################
# Filter samples with low and high counts
######################################

# The data is in log2(expected_count+1), we have to change it (we want
# counts without the transformation)
rnaseq_counts <- round(2^(rnaseq)-1, 0)

rnaseq_counts_primary <- rnaseq_counts[,colnames(rnaseq_counts) %in% phenoData_TCGA_primary$sample]
rnaseq_counts_primary <- rnaseq_counts_primary[, match(phenoData_TCGA_primary$sample, colnames(rnaseq_counts_primary))]

log_rnaseq_counts_primary <- rnaseq[,colnames(rnaseq) %in% phenoData_TCGA_primary$sample]
log_rnaseq_counts_primary <- log_rnaseq_counts_primary[, match(phenoData_TCGA_primary$sample, colnames(log_rnaseq_counts_primary))]


######################################
# Write rnaseq and log rnaseq PT
######################################
setwd(tablesDirectory)
write.csv(rnaseq_counts_primary, file = "rnaseq_counts_primary.csv")
write.csv(log_rnaseq_counts_primary, file = "log_rnaseq_counts_primary.csv")

#######################################
# Filtering races
#######################################
phenoData_TCGA_primary$raza <- as.factor(phenoData_TCGA_primary$raza)
phenoData_TCGA_primary$raza <- factor(phenoData_TCGA_primary$raza,
                                      levels=c("BLANCO", 
                                               "NEGRO_O_AFROAMERICANO", 
                                               "ASIÁTICO", 
                                               "INDIO_AMERICANO_O_NATIVO_DE_ALASKA", 
                                               "HAWAIANO_NATIVO_U_OTRO_ISLEÑO_DEL_PACÍFICO",
                                               "", "[No Evaluado]", "[Desconocido]"))
phenoData_TCGA_primary_split <- split(x = phenoData_TCGA_primary, f = phenoData_TCGA_primary$enfermedad_primaria)


races_list <- c("BLANCO", "NEGRO_O_AFROAMERICANO", "ASIÁTICO", "INDIO_AMERICANO_O_NATIVO_DE_ALASKA", "HAWAIANO_NATIVO_U_OTRO_ISLEÑO_DEL_PACÍFICO")

length(names(phenoData_TCGA_primary_split))
# 27

setwd(resultsDirectory)
pdf(file = "muestras_TCGA_primary_raza.pdf", width = 15, height = 15)
par(mfrow=c(5,6))
for (i in 1:length(names(phenoData_TCGA_primary_split))){
  counts <- table(phenoData_TCGA_primary_split[[i]]$raza)
  counts.df <- as.data.frame(cbind(races_list, rep(0, length(races_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$races_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$races_list
  barplot(vector, 
          main = names(phenoData_TCGA_primary_split)[i], 
          ylab = "Nº muestras", 
          col = c("#482677FF", "#FDE725FF", "#287D8EFF", "#FF0000", "#808080"),
          names.arg = c("BLANCO"="BLANCO",
                        "NEGRO_O_AFROAMERICANO"="NEGRO", 
                        "ASIÁTICO"="ASIÁTICO", 
                        "INDIO_AMERICANO_O_NATIVO_DE_ALASKA"="INDIO", 
                        "HAWAIANO_NATIVO_U_OTRO_ISLEÑO_DEL_PACÍFICO"="HAWAIANO"),
          cex.axis = 0.9,
          cex.lab = 1.4,
          cex.main = 1.5,
          cex.names =0.9,
          las = 2)
}

dev.off()

# FILTER BY RACE
phenoData_TCGA_primary <- phenoData_TCGA_primary[phenoData_TCGA_primary$raza %in% c("BLANCO", "NEGRO_O_AFROAMERICANO", "ASIÁTICO"),]
rnaseq_counts_primary <- rnaseq_counts_primary[, match(phenoData_TCGA_primary$sample, colnames(rnaseq_counts_primary))]
log_rnaseq_counts_primary <- log_rnaseq_counts_primary[, match(phenoData_TCGA_primary$sample, colnames(log_rnaseq_counts_primary))]

#######################################
# Identifying outliers by Sum raw counts by col
#######################################
counts_by_sample <- colSums(rnaseq_counts_primary)
counts_by_sample.df <- as.data.frame(counts_by_sample)

min(counts_by_sample)
max(counts_by_sample)

Q1 <- quantile(counts_by_sample, c(0.25)); Q1
Q3 <- quantile(counts_by_sample, c(0.75)); Q3
IQR <- Q3-Q1

boxplot(counts_by_sample.df$counts_by_sample)

hist_counts_by_sample_PT <- ggplot(data = counts_by_sample.df, aes(x=counts_by_sample)) + 
  geom_histogram(fill="black", colour="black", alpha = 0.25, bins = 100) +
  geom_vline(xintercept = Q1-1.5*IQR, col = "red") +
  geom_vline(xintercept = Q3+1.5*IQR, col = "red") +
  ggtitle("Histograma suma de counts por muestra") +
  xlab("Counts totales por muestra") + ylab("Número de muestras")

setwd(resultsDirectory)
ggsave("Histograma counts PT.png", plot = hist_counts_by_sample_PT, width = 8, height = 5, units = "in")


table(counts_by_sample<(Q1-1.5*IQR)) # 9 outliers
table(counts_by_sample>(Q3+1.5*IQR)) # 119 outliers

table(rownames(counts_by_sample.df)==phenoData_TCGA_primary$sample)
counts_by_sample.df$enfermedad_primaria <- phenoData_TCGA_primary$enfermedad_primaria
counts_by_sample.df$TSS <- phenoData_TCGA_primary$TSS
counts_by_sample.df$BCR <- phenoData_TCGA_primary$BCR
counts_by_sample.df$source_site <- phenoData_TCGA_primary$source_site

### Boxplot counts by sample by PT
boxplot_counts_by_sample_by_PT <- ggplot(data = counts_by_sample.df, aes(x=counts_by_sample, y= enfermedad_primaria)) + 
  geom_boxplot() + ylab("Tumor primario") + xlab("Suma de counts por muestra")

setwd(resultsDirectory)
ggsave("Boxplot counts por PT.png", plot = boxplot_counts_by_sample_by_PT, width = 8, height = 5, units = "in")


### Boxplot counts by TSS

setwd(resultsDirectory)
pdf(file = "Boxplot_por_TSS_PT.pdf", width = 15, height = 4)
print(ggplot(data = counts_by_sample.df, aes(x=counts_by_sample, y= TSS)) + 
    geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_por_BCR_PT.pdf", width = 15, height = 4)
print(ggplot(data = counts_by_sample.df, aes(x=counts_by_sample, y= BCR)) + 
        geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_por_source_site_PT.pdf", width = 15, height = 4)
print(ggplot(data = counts_by_sample.df, aes(x=counts_by_sample, y= source_site)) + 
        geom_boxplot()+theme_minimal()+
        theme(plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white"),
              text = element_text(size = 6, color="black"), 
              axis.text.x = element_text(size = 4, color="black", angle = 90, vjust = 0.5), 
              axis.text.y = element_text(size = 3.5, color="black")) + coord_flip())
dev.off()



### With log counts
log_counts_by_sample <- colSums(log_rnaseq_counts_primary)
log_counts_by_sample.df <- as.data.frame(log_counts_by_sample)

table(rownames(log_counts_by_sample.df)==phenoData_TCGA_primary$sample)
log_counts_by_sample.df$enfermedad_primaria <- phenoData_TCGA_primary$enfermedad_primaria
log_counts_by_sample.df$TSS <- phenoData_TCGA_primary$TSS
log_counts_by_sample.df$BCR <- phenoData_TCGA_primary$BCR
log_counts_by_sample.df$source_site <- phenoData_TCGA_primary$source_site


boxplot_log_counts_by_sample_by_PT <- ggplot(data = log_counts_by_sample.df, aes(x=log_counts_by_sample, y= enfermedad_primaria)) + 
  geom_boxplot() + ylab("Tumor primario") + xlab("Suma de log counts por muestra")

setwd(resultsDirectory)
ggsave("Boxplot log counts por PT.png", plot = boxplot_log_counts_by_sample_by_PT, width = 8, height = 5, units = "in")



setwd(resultsDirectory)
pdf(file = "Boxplot_log_counts_por_TSS_PT.pdf", width = 15, height = 4)
print(ggplot(data = log_counts_by_sample.df, aes(x=log_counts_by_sample, y= TSS)) + 
        geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_log_counts_por_BCR_PT.pdf", width = 15, height = 4)
print(ggplot(data = log_counts_by_sample.df, aes(x=log_counts_by_sample, y= BCR)) + 
        geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_log_counts_por_source_site_PT.pdf", width = 15, height = 4)
print(ggplot(data = log_counts_by_sample.df, aes(x=log_counts_by_sample, y= source_site)) + 
        geom_boxplot()+theme_minimal()+
        theme(plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white"),
              text = element_text(size = 6, color="black"), 
              axis.text.x = element_text(size = 4, color="black", angle = 90, vjust = 0.5), 
              axis.text.y = element_text(size = 3.5, color="black")) + coord_flip())
dev.off()

samples_down_outliers <- names(counts_by_sample[counts_by_sample<(Q1-1.5*IQR)])
samples_up_outliers <- names(counts_by_sample[counts_by_sample>(Q3+1.5*IQR)])

######################################
### Filtering
######################################
rnaseq_counts_primary.f <- rnaseq_counts_primary[, !colnames(rnaseq_counts_primary) %in% c(samples_down_outliers,samples_up_outliers)]
log_rnaseq_counts_primary.f <- log_rnaseq_counts_primary[, !colnames(log_rnaseq_counts_primary) %in% c(samples_down_outliers,samples_up_outliers)]
phenoData_TCGA_primary.f <- phenoData_TCGA_primary[!phenoData_TCGA_primary$sample %in% c(samples_down_outliers,samples_up_outliers), ]

rnaseq_counts_primary.f <- rnaseq_counts_primary.f[,match(phenoData_TCGA_primary.f$sample, colnames(rnaseq_counts_primary.f))]
log_rnaseq_counts_primary.f <- log_rnaseq_counts_primary.f[,match(phenoData_TCGA_primary.f$sample, colnames(log_rnaseq_counts_primary.f))]


### Boxplot by sample by PT
counts_by_sample.f <- colSums(rnaseq_counts_primary.f)
counts_by_sample.f.df <- as.data.frame(counts_by_sample.f)

table(rownames(counts_by_sample.f.df)==phenoData_TCGA_primary.f$sample)
counts_by_sample.f.df$enfermedad_primaria <- phenoData_TCGA_primary.f$enfermedad_primaria
counts_by_sample.f.df$TSS <- phenoData_TCGA_primary.f$TSS
counts_by_sample.f.df$BCR <- phenoData_TCGA_primary.f$BCR
counts_by_sample.f.df$source_site <- phenoData_TCGA_primary.f$source_site

boxplot_counts_by_sample_by_PT.f <- ggplot(data = counts_by_sample.f.df, aes(x=counts_by_sample.f, y= enfermedad_primaria)) + 
  geom_boxplot() + ylab("Tumor primario") + xlab("Suma de counts por muestra")

setwd(resultsDirectory)
ggsave("Boxplot counts por PT filtradas.png", plot = boxplot_counts_by_sample_by_PT.f, width = 8, height = 5, units = "in")


### Boxplot by sample by PT log counts
log_counts_by_sample.f <- colSums(log_rnaseq_counts_primary.f)
log_counts_by_sample.f.df <- as.data.frame(log_counts_by_sample.f)

table(rownames(log_counts_by_sample.f.df)==phenoData_TCGA_primary.f$sample)
log_counts_by_sample.f.df$enfermedad_primaria <- phenoData_TCGA_primary.f$enfermedad_primaria
log_counts_by_sample.f.df$TSS <- phenoData_TCGA_primary.f$TSS
log_counts_by_sample.f.df$BCR <- phenoData_TCGA_primary.f$BCR
log_counts_by_sample.f.df$source_site <- phenoData_TCGA_primary.f$source_site

boxplot_log_counts_by_sample_by_PT.f <- ggplot(data = log_counts_by_sample.f.df, aes(x=log_counts_by_sample.f, y= enfermedad_primaria)) + 
  geom_boxplot() + ylab("Tumor primario") + xlab("Suma de log counts por muestra")

setwd(resultsDirectory)
ggsave("Boxplot log counts por PT filtradas.png", plot = boxplot_log_counts_by_sample_by_PT.f, width = 8, height = 5, units = "in")

######################################
# Write phenoData_TCGA_primary filter
######################################
setwd(tablesDirectory)
write.csv(phenoData_TCGA_primary.f, file = "phenoData_TCGA_primary_f.csv")
write.csv(rnaseq_counts_primary.f, file = "rnaseq_counts_primary_f.csv")
write.csv(log_rnaseq_counts_primary.f, file = "log_rnaseq_counts_primary_f.csv")
