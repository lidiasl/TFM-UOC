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


phenoData_TCGA_primary <- phenoData_TCGA_primary[phenoData_TCGA_primary$race %in% 
                                                            c("WHITE", "BLACK OR AFRICAN AMERICAN", "ASIAN"),]
phenoData_TCGA_metastatic <- phenoData_TCGA_metastatic[phenoData_TCGA_metastatic$race %in% 
                                                         c("WHITE", "BLACK OR AFRICAN AMERICAN", "ASIAN"),]

phenoData_TCGA_primary$race <- str_replace_all(phenoData_TCGA_primary$race, " ", "_")
phenoData_TCGA_metastatic$race <- str_replace_all(phenoData_TCGA_metastatic$race, " ", "_")


######################################
# Write phenoData_TCGA_primary and phenoData_TCGA_metastatic
######################################
setwd(tablesDirectory)
write.csv(phenoData_TCGA_primary, file = "phenoData_TCGA_primary.csv")
write.csv(phenoData_TCGA_metastatic, file = "phenoData_TCGA_metastatic.csv")

setwd(tablesDirectory)
phenoData_TCGA_primary <- read.csv(file="phenoData_TCGA_primary.csv")

phenoData_TCGA_primary <- phenoData_TCGA_primary[, colnames(phenoData_TCGA_primary)!="X"]
rownames(phenoData_TCGA_primary) <- phenoData_TCGA_primary$sample

#####################################################################
# Graphs
#####################################################################

# Counts by primary disease (Primary T)
phenoData_TCGA_primary$primary_disease_set <- factor(phenoData_TCGA_primary$primary_disease_set,
                                                     levels = names(sort(table(phenoData_TCGA_primary$primary_disease_set), decreasing = FALSE)))

Counts_by_primary_disease_set_PT <- ggplot(phenoData_TCGA_primary, aes(primary_disease_set)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Counts by primary disease (set Primary Tumor)") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease SET PT.png", plot = Counts_by_primary_disease_set_PT, width = 8, height = 5, units = "in")


setwd(resultsDirectory)
pdf(file = "Counts by primary disease SET PT.pdf", width = 10, height = 8)
print(Counts_by_primary_disease_set_PT)
dev.off()

# Counts by primary disease (Metastatic T)
phenoData_TCGA_metastatic$primary_disease_set <- factor(phenoData_TCGA_metastatic$primary_disease_set,
                                                     levels = names(sort(table(phenoData_TCGA_metastatic$primary_disease_set), decreasing = FALSE)))

Counts_by_primary_disease_set_MetT <- ggplot(phenoData_TCGA_metastatic, aes(primary_disease_set)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Counts by primary disease (set Metastatic Tumor)") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease SET MetT.png", plot = Counts_by_primary_disease_set_MetT, width = 8, height = 5, units = "in")


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
# Identifying outliers by Sum raw counts by col
#######################################
setwd(tablesDirectory)
rnaseq_counts_primary <- read.csv(file = "rnaseq_counts_primary.csv")
rownames(rnaseq_counts_primary) <- rnaseq_counts_primary$X
rnaseq_counts_primary <- rnaseq_counts_primary[, colnames(rnaseq_counts_primary)!="X"]


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
counts_by_sample.df$primary_disease_set <- phenoData_TCGA_primary$primary_disease_set
counts_by_sample.df$TSS <- phenoData_TCGA_primary$TSS
counts_by_sample.df$BCR <- phenoData_TCGA_primary$BCR
counts_by_sample.df$source_site <- phenoData_TCGA_primary$source_site

### Boxplot counts by sample by PT
boxplot_counts_by_sample_by_PT <- ggplot(data = counts_by_sample.df, aes(x=counts_by_sample, y= primary_disease_set)) + 
  geom_boxplot() 

setwd(resultsDirectory)
ggsave("Boxplot counts by PT.png", plot = boxplot_counts_by_sample_by_PT, width = 8, height = 5, units = "in")


### Boxplot counts by TSS

setwd(resultsDirectory)
pdf(file = "Boxplot_by_TSS_PT.pdf", width = 15, height = 4)
print(ggplot(data = counts_by_sample.df, aes(x=counts_by_sample, y= TSS)) + 
    geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_by_BCR_PT.pdf", width = 15, height = 4)
print(ggplot(data = counts_by_sample.df, aes(x=counts_by_sample, y= BCR)) + 
        geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_by_source_site_PT.pdf", width = 15, height = 4)
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
log_counts_by_sample.df$primary_disease_set <- phenoData_TCGA_primary$primary_disease_set
log_counts_by_sample.df$TSS <- phenoData_TCGA_primary$TSS
log_counts_by_sample.df$BCR <- phenoData_TCGA_primary$BCR
log_counts_by_sample.df$source_site <- phenoData_TCGA_primary$source_site


boxplot_log_counts_by_sample_by_PT <- ggplot(data = log_counts_by_sample.df, aes(x=log_counts_by_sample, y= primary_disease_set)) + 
  geom_boxplot() 

setwd(resultsDirectory)
ggsave("Boxplot log counts by PT.png", plot = boxplot_log_counts_by_sample_by_PT, width = 8, height = 5, units = "in")



setwd(resultsDirectory)
pdf(file = "Boxplot_log_counts_by_TSS_PT.pdf", width = 15, height = 4)
print(ggplot(data = log_counts_by_sample.df, aes(x=log_counts_by_sample, y= TSS)) + 
        geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_log_counts_by_BCR_PT.pdf", width = 15, height = 4)
print(ggplot(data = log_counts_by_sample.df, aes(x=log_counts_by_sample, y= BCR)) + 
        geom_boxplot()+theme_minimal()+coord_flip())
dev.off()

setwd(resultsDirectory)
pdf(file = "Boxplot_log_counts_by_source_site_PT.pdf", width = 15, height = 4)
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
counts_by_sample.f.df$primary_disease_set <- phenoData_TCGA_primary.f$primary_disease_set
counts_by_sample.f.df$TSS <- phenoData_TCGA_primary.f$TSS
counts_by_sample.f.df$BCR <- phenoData_TCGA_primary.f$BCR
counts_by_sample.f.df$source_site <- phenoData_TCGA_primary.f$source_site

boxplot_counts_by_sample_by_PT.f <- ggplot(data = counts_by_sample.f.df, aes(x=counts_by_sample.f, y= primary_disease_set)) + 
  geom_boxplot() 

setwd(resultsDirectory)
ggsave("Boxplot counts by PT filter.png", plot = boxplot_counts_by_sample_by_PT.f, width = 8, height = 5, units = "in")



setwd(resultsDirectory)
pdf(file = "Boxplot_by_source_site_PT_filter.pdf", width = 15, height = 4)
print(ggplot(data = counts_by_sample.f.df, aes(x=counts_by_sample.f, y= source_site)) + 
        geom_boxplot()+theme_minimal()+
        theme(plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white"),
              text = element_text(size = 6, color="black"), 
              axis.text.x = element_text(size = 4, color="black", angle = 90, vjust = 0.5), 
              axis.text.y = element_text(size = 3.5, color="black")) + coord_flip())
dev.off()



### Boxplot by sample by PT log counts
log_counts_by_sample.f <- colSums(log_rnaseq_counts_primary.f)
log_counts_by_sample.f.df <- as.data.frame(log_counts_by_sample.f)

table(rownames(log_counts_by_sample.f.df)==phenoData_TCGA_primary.f$sample)
log_counts_by_sample.f.df$primary_disease_set <- phenoData_TCGA_primary.f$primary_disease_set
log_counts_by_sample.f.df$TSS <- phenoData_TCGA_primary.f$TSS
log_counts_by_sample.f.df$BCR <- phenoData_TCGA_primary.f$BCR
log_counts_by_sample.f.df$source_site <- phenoData_TCGA_primary.f$source_site

boxplot_log_counts_by_sample_by_PT.f <- ggplot(data = log_counts_by_sample.f.df, aes(x=log_counts_by_sample.f, y= primary_disease_set)) + 
  geom_boxplot() 

setwd(resultsDirectory)
ggsave("Boxplot log counts by PT filter.png", plot = boxplot_log_counts_by_sample_by_PT.f, width = 8, height = 5, units = "in")


setwd(resultsDirectory)
pdf(file = "Boxplot_log_counts_by_source_site_PT_filter.pdf", width = 15, height = 4)
print(ggplot(data = log_counts_by_sample.f.df, aes(x=log_counts_by_sample.f, y= source_site)) + 
        geom_boxplot()+theme_minimal()+
        theme(plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white"),
              text = element_text(size = 6, color="black"), 
              axis.text.x = element_text(size = 4, color="black", angle = 90, vjust = 0.5), 
              axis.text.y = element_text(size = 3.5, color="black")) + coord_flip())
dev.off()

######################################
# Write phenoData_TCGA_primary filter
######################################
setwd(tablesDirectory)
write.csv(phenoData_TCGA_primary.f, file = "phenoData_TCGA_primary_f.csv")
write.csv(rnaseq_counts_primary.f, file = "rnaseq_counts_primary_f.csv")
write.csv(log_rnaseq_counts_primary.f, file = "log_rnaseq_counts_primary_f.csv")
