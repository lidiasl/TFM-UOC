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


# Filtering 
phenoData_TCGA <- as.data.frame(phenoData_TCGA)
phenoData_TCGA_in_rnaseq <- phenoData_TCGA[phenoData_TCGA$sample %in% colnames(rnaseq),]
phenoData_TCGA_in_rnaseq_with_gender <- phenoData_TCGA_in_rnaseq[phenoData_TCGA_in_rnaseq$sample %in% phenoData_TCGA_gender$sample,]

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


# Table of primary_disease
table(phenoData_TCGA$primary_disease)

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


# Table of primary_disease in rnaseq

phenoData_TCGA_in_rnaseq$primary_disease <- factor(phenoData_TCGA_in_rnaseq$primary_disease,
                                         levels = names(sort(table(phenoData_TCGA_in_rnaseq$primary_disease), decreasing = FALSE)))

Counts_by_primary_disease_in_rnaseq <- ggplot(phenoData_TCGA_in_rnaseq, aes(primary_disease)) + geom_bar(color="blue", fill="blue", alpha=0.2) +
  ggtitle("Counts by primary disease (in rnaseq)") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease in rnaseq.png", plot = Counts_by_primary_disease_in_rnaseq, width = 8, height = 5, units = "in")

# Table of primary_disease in rnaseq

phenoData_TCGA_in_rnaseq_with_gender$primary_disease <- factor(phenoData_TCGA_in_rnaseq_with_gender$primary_disease,
                                                   levels = names(sort(table(phenoData_TCGA_in_rnaseq_with_gender$primary_disease), decreasing = FALSE)))

Counts_by_primary_disease_in_rnaseq_with_gender <- ggplot(phenoData_TCGA_in_rnaseq_with_gender, aes(primary_disease)) + geom_bar(color="blue", fill="blue", alpha=0.2) +
  ggtitle("Counts by primary disease (in rnaseq with gender)") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease in rnaseq with gender.png", plot = Counts_by_primary_disease_in_rnaseq_with_gender, width = 8, height = 5, units = "in")



######################################
# Write phenoData_TCGA
######################################
setwd(tablesDirectory)
write.csv(phenoData_TCGA, file = "phenoData_TCGA.csv")
write.csv(phenoData_TCGA_in_rnaseq_with_gender, file = "phenoData_TCGA_in_rnaseq_with_gender.csv")
######################################

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


#####################################################################
# Study by primary disease
#####################################################################

phenoData_TCGA_primary$sex <- phenoData_TCGA_gender$gender[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$race <- phenoData_TCGA_gender$race[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$tumor_stage <- phenoData_TCGA_gender$ajcc_pathologic_tumor_stage[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$histological_type <- phenoData_TCGA_gender$histological_type[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$vital_status <- phenoData_TCGA_gender$vital_status[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$age <- phenoData_TCGA_gender$age_at_initial_pathologic_diagnosis[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]

phenoData_TCGA_metastatic$sex <- phenoData_TCGA_gender$gender[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$race <- phenoData_TCGA_gender$race[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$tumor_stage <- phenoData_TCGA_gender$ajcc_pathologic_tumor_stage[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$histological_type <- phenoData_TCGA_gender$histological_type[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$vital_status <- phenoData_TCGA_gender$vital_status[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$age <- phenoData_TCGA_gender$age_at_initial_pathologic_diagnosis[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]

phenoData_TCGA_primary <- phenoData_TCGA_primary[phenoData_TCGA_primary$race %in% 
                                                            c("WHITE", "BLACK OR AFRICAN AMERICAN",
                                                              "ASIAN", "AMERICAN INDIAN OR ALASKA NATIVE"),]
phenoData_TCGA_metastatic <- phenoData_TCGA_metastatic[phenoData_TCGA_metastatic$race %in% 
                                                         c("WHITE", "BLACK OR AFRICAN AMERICAN",
                                                           "ASIAN", "AMERICAN INDIAN OR ALASKA NATIVE"),]


phenoData_TCGA_primary$primary_disease <- factor(phenoData_TCGA_primary$primary_disease,
                                                 levels = names(sort(table(phenoData_TCGA_primary$primary_disease), decreasing = FALSE)))

Counts_by_primary_disease_PT <- ggplot(phenoData_TCGA_primary, aes(primary_disease)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Counts by primary disease Primary Tumor") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease PT.png", plot = Counts_by_primary_disease_PT, width = 8, height = 5, units = "in")


#####################################################################
# Group in disease set
#####################################################################

setwd(tablesDirectory)

set_subset <- read.csv(file = "set_subset.csv", header = TRUE)
phenoData_TCGA_primary$primary_disease_set <- set_subset$set[match(phenoData_TCGA_primary$primary_disease, set_subset$subset)]

phenoData_TCGA_primary$primary_disease_set <- factor(phenoData_TCGA_primary$primary_disease_set,
                                                 levels = names(sort(table(phenoData_TCGA_primary$primary_disease_set), decreasing = FALSE)))

Counts_by_primary_disease_PT_set <- ggplot(phenoData_TCGA_primary, aes(primary_disease_set)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Counts by primary disease Primary Tumor SET") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease PT SET.png", plot = Counts_by_primary_disease_PT_set, width = 8, height = 5, units = "in")



phenoData_TCGA_primary_split <- split(x = phenoData_TCGA_primary, f = phenoData_TCGA_primary$primary_disease_set)




races_list <- c("AMERICAN INDIAN OR ALASKA NATIVE", 
                "ASIAN", 
              "BLACK OR AFRICAN AMERICAN", 
              "WHITE")

pdf(file = "samples_PT_race.pdf", width = 15, height = 15)
par(mfrow=c(5,6))
for (i in 1:length(names(phenoData_TCGA_primary_split))){
  counts <- table(phenoData_TCGA_primary_split[[i]]$race)
  counts.df <- as.data.frame(cbind(races_list, rep(0, length(races_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$races_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$races_list
  barplot(vector, 
          main = names(phenoData_TCGA_primary_split)[i], 
          ylab = "Nº muestras", 
          col = c("deepskyblue3", "coral1", "mediumpurple2", "palegreen3"),
          names.arg = c("AMERICAN INDIAN OR ALASKA NATIVE" = "A.I/A.N", 
                        "ASIAN" = "Asia", 
                        "BLACK OR AFRICAN AMERICAN" = "Black", 
                        "WHITE" = "White"),
          cex.axis = 1,
          cex.lab = 1.4,
          cex.main = 1.5,
          cex.names =1.2,
          las = 2)
}

dev.off()



#####################################################################
# Create subset of samples by primary disease set
#####################################################################

names(table(phenoData_TCGA_primary_split$esophageal$race)[table(phenoData_TCGA_primary_split$esophageal$race)<=5])
n <- 21
subset_PT <- list()

for (i in 1:length(names(phenoData_TCGA_primary_split))){
  races_split <- split(phenoData_TCGA_primary_split[[i]], phenoData_TCGA_primary_split[[i]]$race)
  names <- names(table(phenoData_TCGA_primary_split[[i]]$race)[table(phenoData_TCGA_primary_split[[i]]$race)<=5])
  subset_PT[[i]] <- phenoData_TCGA_primary_split[[i]][phenoData_TCGA_primary_split[[i]]$race %in% names,]
  races_split <- races_split[!names(races_split) %in% names]
  while (length(races_split)>1){
    n_samples <- n - nrow(subset_PT[[i]])
    n_races <- length(unique(phenoData_TCGA_primary_split[[i]]$race)) - length(names)
    n_random <- floor(n_samples/n_races)
    test <- any(lapply(races_split, function(df) nrow(df) < n_random))
    if (test==TRUE){
      race_min <- names(which.min(sapply(races_split, nrow)))
      subset_PT[[i]] <- rbind(subset_PT[[i]], phenoData_TCGA_primary_split[[i]][phenoData_TCGA_primary_split[[i]]$race %in% race_min,])
      races_split <- races_split[!names(races_split) %in% race_min]
      names <- c(names, race_min)
    } else{
      for (j in 1:length(races_split)){
        index <- sample(seq_len(nrow(races_split[[j]])), size = n_random)
        subset_PT[[i]] <- rbind(subset_PT[[i]], races_split[[j]][index,])
      }
      if (nrow(subset_PT[[i]]) < n){
        if(n - nrow(subset_PT[[i]])==1){
          names_ult <- names(sort(sapply(races_split,nrow), decreasing = TRUE)[1])
          subset_PT[[i]] <- rbind(subset_PT[[i]], races_split[[j]][index,])
        }else{
          names_ult <- names(sort(sapply(races_split,nrow), decreasing = TRUE)[1:(n - nrow(subset_PT[[i]]))])
        }
        
        
      }
      races_split <- NULL
    }
  }
  if (length(races_split)==1){
    n_samples <- n - nrow(subset_PT[[i]])
    index <- sample(seq_len(nrow(races_split[[1]])), size = n_samples)
    subset_PT[[i]] <- rbind(subset_PT[[i]], races_split[[1]][index,])
  }
}



#####################################################################
# Separate TRAIN and TEST
#####################################################################

# Create training and testing datasets
# For each primary disease, select 80% for train and 20% for test

## set the seed to make your partition reproducible
set.seed(123)

selection <- lapply(phenoData_TCGA_primary.f_split, function(x){
  smp_size <- floor(0.8 * nrow(x)) #  Returns the largest integer that is smaller than or equal to value passed to it as argument 
  #"sample" takes a sample of the specified size from the elements of "seq_len(nrow(x))" without replacement.
  train_ind <- sample(seq_len(nrow(x)), size = smp_size)
  
  train_primary <- x[train_ind, ]
  test_primary <- x[-train_ind, ]
  train_primary$set <- "Train"
  test_primary$set <- "Test"
  rbind(train_primary,test_primary)
})

#Check:
#prop.table(table(selection[["acute_myeloid_leukemia"]]$set))

train_vs_test <- Reduce(rbind, selection)

################################################################################
#Add variables phenodata TCGA
################################################################################

train_vs_test$sex <- phenoData_TCGA_gender$gender[match(train_vs_test$sample, phenoData_TCGA_gender$sample)]
train_vs_test$race <- phenoData_TCGA_gender$race[match(train_vs_test$sample, phenoData_TCGA_gender$sample)]
train_vs_test$tumor_stage <- phenoData_TCGA_gender$ajcc_pathologic_tumor_stage[match(train_vs_test$sample, phenoData_TCGA_gender$sample)]
train_vs_test$histological_type <- phenoData_TCGA_gender$histological_type[match(train_vs_test$sample, phenoData_TCGA_gender$sample)]
train_vs_test$vital_status <- phenoData_TCGA_gender$vital_status[match(train_vs_test$sample, phenoData_TCGA_gender$sample)]


phenodata_train_primary <- train_vs_test[train_vs_test$set == "Train",]
phenodata_test_primary <- train_vs_test[train_vs_test$set == "Test",]

# Get number of samples of the minimal set of primary disease in training set:
min_disease <- names(table(phenodata_train_primary$primary_disease)[which.min(table(phenodata_train_primary$primary_disease))]) #"cholangiocarcinoma"
max_disease <- names(table(phenodata_train_primary$primary_disease)[which.max(table(phenodata_train_primary$primary_disease))]) #"breast_invasive_carcinoma"
min_number <- nrow(phenodata_train_primary[phenodata_train_primary$primary_disease==min_disease,]) #28
max_number <- nrow(phenodata_train_primary[phenodata_train_primary$primary_disease==max_disease,]) #873


#####################################################################
# LOAD DATA MET500
#####################################################################

# Load FPKM for MET500
FPKM_MET500 <- fread("/usr/local/Projects/TCGA_classification/Data/M.mx.txt", sep = "\t", stringsAsFactors = F)
genenames_MET500 <- separate(data = FPKM_MET500, col = sample, into = c("GeneName", "version"), remove = T)
FPKM_MET500 <- FPKM_MET500[,-1]
FPKM_MET500 <- as.data.frame(FPKM_MET500)
rownames(FPKM_MET500) <- genenames_MET500$GeneName

# Filter genes common in both datasets:
table(rownames(rnaseq) %in% rownames(FPKM_MET500))
#FALSE  TRUE 
#41435 19063 

rnaseq <- rnaseq[rownames(rnaseq) %in% rownames(FPKM_MET500),]

#Check:
table(rownames(rnaseq) %in% rownames(FPKM_MET500))
#TRUE 
#19063 


# The data is in log2(expected_count+1), we have to change it (we want
# counts without the transformation)
rnaseq_counts <- round(2^(rnaseq)-1, 0)

# Remove genenames
rm(genenames)
rm(genenames_MET500)

# Split in training, test dataset and metastatic:
rnaseq_counts_primary_train <- rnaseq_counts[,colnames(rnaseq_counts) %in% phenodata_train_primary$sample]
phenodata_train_primary <- phenodata_train_primary[phenodata_train_primary$sample %in% colnames(rnaseq_counts_primary_train),]
phenodata_train_primary <- phenodata_train_primary[match(colnames(rnaseq_counts_primary_train),phenodata_train_primary$sample),]

rnaseq_counts_primary_test <- rnaseq_counts[,colnames(rnaseq_counts) %in% phenodata_test_primary$sample]
phenodata_test_primary <- phenodata_test_primary[phenodata_test_primary$sample %in% colnames(rnaseq_counts_primary_test),]
phenodata_test_primary <- phenodata_test_primary[match(colnames(rnaseq_counts_primary_test),phenodata_test_primary$sample),]


Counts_by_primary_disease_set <- ggplot(train_vs_test, aes(primary_disease,fill=set)) + geom_bar(alpha=0.4) + coord_flip() + 
  ggtitle("Counts by primary disease (Train vs test) / PT rnaseq") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label = after_stat(count), color=set), position = position_stack(vjust = 0.5), size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = c("red1","blue")) +
  scale_color_manual(values = c("red4","darkblue")) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white")) + labs(fill="Set")

setwd(resultsDirectory)
ggsave("Counts by primary disease set.png", plot = Counts_by_primary_disease_set, width = 10, height = 5, units = "in")

table(train_vs_test$gender)
#no  yes 
#22 9332 

table(train_vs_test$primary_disease[train_vs_test$gender=="no"])



table(train_vs_test$sex)
ggplot(train_vs_test, aes(sex,fill=sex)) + geom_bar(alpha=1) + 
  ggtitle("Samples by sex (train and test)") + ylab("Counts") + xlab("Sex") + theme_minimal() + 
  geom_text(stat='count', aes(label = after_stat(count)), position = position_stack(vjust = 1.03), size = 3.5, show.legend = FALSE)

table(train_vs_test$race)
ggplot(train_vs_test, aes(race,fill=race)) + geom_bar(alpha=1) + 
  ggtitle("Samples by race (train and test)") + ylab("Counts") + xlab("Sex") + theme_minimal() + coord_flip() +  
  geom_text(stat='count', aes(label = after_stat(count)), position = position_stack(vjust = 1.03), size = 3.5, show.legend = FALSE)

table(train_vs_test$tumor_stage)
ggplot(train_vs_test, aes(tumor_stage,fill=tumor_stage)) + geom_bar(alpha=1) + 
  ggtitle("Samples by tumor_stage (train and test)") + ylab("Counts") + xlab("Sex") + theme_minimal() + coord_flip() +  
  geom_text(stat='count', aes(label = after_stat(count)), position = position_stack(vjust = 1.03), size = 3.5, show.legend = FALSE)


length(unique(train_vs_test$histological_type)) # 145

################################################################################
################################################################################

#RNA-seq count distribution

ggplot(rnaseq_counts_primary_train, aes(x = TCGA_19_1787_01)) +
  geom_histogram(stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  ggtitle("TCGA_19_1787_01")

#If we zoom in close to zero, we can see a large number of genes with counts of zero:

ggplot(rnaseq_counts_primary_train, aes(x = TCGA_19_1787_01)) +
  geom_histogram(stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  ggtitle("TCGA_19_1787_01") +
  xlim(-5, 500)

# Now the histogram for a specific gene:
rnaseq_counts_primary_train.t <- t(rnaseq_counts_primary_train)
ggplot(rnaseq_counts_primary_train.t, aes(x = ENSG00000167578)) +
  geom_histogram(stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of samples") +
  ggtitle("ENSG00000167578")

# Check if mean < variance (to use NB model) for each subset of replicates
# for each gene

mean_counts <- apply(rnaseq_counts_primary_train[,colnames(rnaseq_counts_primary_train) %in%
                                            phenodata_train_primary$sample[phenodata_train_primary$primary_disease == "acute_myeloid_leukemia"]], 1, mean)

variance_counts <- apply(rnaseq_counts_primary_train[,colnames(rnaseq_counts_primary_train) %in%
                                                phenodata_train_primary$sample[phenodata_train_primary$primary_disease == "acute_myeloid_leukemia"]], 1, var)

df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts), color="red") +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Variance vs mean") +
  xlab("mean") +
  ylab("variance")

################################################################################
#Add variables phenodata TCGA
################################################################################

################################################################################
################################################################################

# Filtering by expression
# Filtering out lowly expressed genes
smallestGroupSize <- min(table(phenodata_train_primary$primary_disease)) #28

nrow(rnaseq_counts_primary_train) #19063

keep.train <- rowSums(rnaseq_counts_primary_train >= 10) >= floor(0.9*smallestGroupSize)

table(keep.train)
# keep 18472 genes
# Remove 591 genes

# Filter by expression

rnaseq_counts_primary_train <- rnaseq_counts_primary_train[keep.train,]


################################################################################
################################################################################

# Identifying outliers by Sum raw counts by col

counts_by_sample <- colSums(rnaseq_counts_primary_train)
counts_by_sample.df <- as.data.frame(counts_by_sample)

min(counts_by_sample)
max(counts_by_sample)

Q1 <- quantile(counts_by_sample, c(0.25)); Q1
Q3 <- quantile(counts_by_sample, c(0.75)); Q3
IQR <- Q3-Q1

boxplot(counts_by_sample.df$counts_by_sample)

ggplot(data = counts_by_sample.df, aes(x=counts_by_sample)) + 
  geom_histogram(fill="black", colour="black", alpha = 0.25, bins = 100) +
  geom_vline(xintercept = Q1-1.5*IQR, col = "red") +
  geom_vline(xintercept = Q3+1.5*IQR, col = "red") +
  ggtitle("Histogram of counts by sample") +
  xlab("Counts") + ylab("Numer of samples")

table(counts_by_sample<(Q1-1.5*IQR)) # 6 outliers
table(counts_by_sample>(Q3+1.5*IQR)) # 108 outliers



samples_down_outliers <- names(counts_by_sample[counts_by_sample<(Q1-1.5*IQR)])
samples_up_outliers <- names(counts_by_sample[counts_by_sample>(Q3+1.5*IQR)])

table_diseases_up_outliers <- table(phenodata_train_primary$primary_disease[phenodata_train_primary$sample %in% samples_up_outliers])
diseases_up_outliers <- names(table(phenodata_train_primary$primary_disease[phenodata_train_primary$sample %in% samples_up_outliers]))
table_diseases <- table(phenodata_train_primary$primary_disease[phenodata_train_primary$primary_disease %in% diseases_up_outliers])

# Outliers by disease in percentage:
table_diseases_up_outliers/table_diseases*100

setwd(QCDirectory)
write.table(table_diseases_up_outliers/table_diseases*100, file = "outliers_by_disease_percentage.csv", 
            sep = " ", dec = ".", row.names = FALSE,
            col.names = TRUE)

rnaseq_primary_train <- rnaseq_primary_train[, !colnames(rnaseq_primary_train) %in% c(samples_down_outliers,samples_up_outliers)]
phenodata_train_primary <- phenodata_train_primary[!phenodata_train_primary$sample %in% c(samples_down_outliers,samples_up_outliers), ]

################################################################################
################################################################################

# Identifying outliers (samples) by correlation depending on counts

dds <- DESeqDataSetFromMatrix(countData = rnaseq_primary_train, colData = phenodata_train_primary, design = ~ primary_disease)
dds <- estimateSizeFactors(dds)
## Transform counts for data visualization
## Do vst to take less time:
rld <- vst(dds, blind=TRUE)
rld_mat <- assay(rld)

correlation_matrix <- cor(rld_mat) 
table(correlation_matrix>=0.7)
# True 51131239
# False 3170922

#Rate:
51131239/(3170922 + 51131239) # 94%

# Como la mayoría de genes no están diferentcialmente expresados, la correlación entre las muestras 
# debería de ser igual o mayor al 0.80.

table(correlation_matrix>=0.6)
# True 54229667
# False 72494
54229667/(72494 + 54229667) # 99.8%

table(correlation_matrix>=0.8)
# True 31799053
# False 22503108
31799053/(31799053+22503108) # 58.6%

################################################################################
# PCA general
################################################################################
prior.general.pca <- prcomp(t(rnaseq_primary_train), 
                    center = TRUE, 
                    scale. = TRUE) 

pca.general.plot <- autoplot(prior.general.pca, 
                     data = phenodata_train_primary, 
                     colour = 'primary_disease') 

pca.general.plot


################################################################################
# Filter test data
################################################################################

rnaseq_primary_test <- rnaseq_primary_test[match(rownames(rnaseq_primary_train), rownames(rnaseq_primary_test)),]

rnaseq_metastatic <- rnaseq_metastatic[match(rownames(rnaseq_primary_train), rownames(rnaseq_metastatic)),]

################################################################################
# Save final data
################################################################################
setwd(tablesDirectory)
write.csv(phenodata_train_primary, file = "phenodata_train_primary.csv")
write.csv(phenodata_test_primary, file = "phenodata_test_primary.csv")
write.csv(phenodata_metastatic, file = "phenodata_metastatic.csv")

write.csv(rnaseq_primary_train, file = "rnaseq_primary_train.csv")
write.csv(rnaseq_primary_test, file = "rnaseq_primary_test.csv")
write.csv(rnaseq_metastatic, file = "rnaseq_metastatic.csv")






################################################################################
# Read final data
################################################################################
setwd(tablesDirectory)
phenodata_train_primary <- read.csv(file = "phenodata_train_primary.csv")
rnaseq_primary_train <- read.csv(file = "rnaseq_primary_train.csv")

# 10 samples of each condition:

library(dplyr)
phenodata_train_primary_divided <- phenodata_train_primary %>% group_by(primary_disease)
phenodata_train_primary_divided <- group_split(phenodata_train_primary_divided)
phenodata_train_primary_divided <- lapply(phenodata_train_primary_divided, function(x){
  x[1:20,]
})

phenodata_train_primary_divided <- Reduce(rbind,phenodata_train_primary_divided)

# PCA general
prior.general.pca.divided <- prcomp(t(rnaseq_primary_train[, colnames(rnaseq_primary_train) %in% phenodata_train_primary_divided$sample]), 
                            center = TRUE, 
                            scale. = TRUE) 

pca.general.divided.plot <- autoplot(prior.general.pca.divided, 
                             data = phenodata_train_primary_divided, 
                             colour = 'primary_disease') 

pca.general.divided.plot
