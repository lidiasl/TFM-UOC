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

#########################################################################

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
setwd(tablesDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

phenoData_TCGA_primary$sex <- phenoData_TCGA_gender$gender[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$race <- phenoData_TCGA_gender$race[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$tumor_stage <- phenoData_TCGA_gender$ajcc_pathologic_tumor_stage[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$histological_type <- phenoData_TCGA_gender$histological_type[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$vital_status <- phenoData_TCGA_gender$vital_status[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$age <- phenoData_TCGA_gender$age_at_initial_pathologic_diagnosis[match(phenoData_TCGA_primary$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_primary$primary_disease_set <- set_subset$set[match(phenoData_TCGA_primary$primary_disease, set_subset$subset)]



phenoData_TCGA_metastatic$sex <- phenoData_TCGA_gender$gender[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$race <- phenoData_TCGA_gender$race[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$tumor_stage <- phenoData_TCGA_gender$ajcc_pathologic_tumor_stage[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$histological_type <- phenoData_TCGA_gender$histological_type[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$vital_status <- phenoData_TCGA_gender$vital_status[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$age <- phenoData_TCGA_gender$age_at_initial_pathologic_diagnosis[match(phenoData_TCGA_metastatic$sample, phenoData_TCGA_gender$sample)]
phenoData_TCGA_metastatic$primary_disease_set <- set_subset$set[match(phenoData_TCGA_metastatic$primary_disease, set_subset$subset)]



phenoData_TCGA_primary <- phenoData_TCGA_primary[phenoData_TCGA_primary$race %in% 
                                                            c("WHITE", "BLACK OR AFRICAN AMERICAN",
                                                              "ASIAN", "AMERICAN INDIAN OR ALASKA NATIVE"),]
phenoData_TCGA_metastatic <- phenoData_TCGA_metastatic[phenoData_TCGA_metastatic$race %in% 
                                                         c("WHITE", "BLACK OR AFRICAN AMERICAN",
                                                           "ASIAN", "AMERICAN INDIAN OR ALASKA NATIVE"),]


#####################################################################
# Group in disease set
#####################################################################

setwd(tablesDirectory)

phenoData_TCGA_primary$primary_disease_set <- factor(phenoData_TCGA_primary$primary_disease_set,
                                                     levels = names(sort(table(phenoData_TCGA_primary$primary_disease_set), decreasing = FALSE)))

Counts_by_primary_disease_set_PT <- ggplot(phenoData_TCGA_primary, aes(primary_disease_set)) + geom_bar(color="blue", fill="blue", alpha=0.2) + 
  ggtitle("Counts by primary disease set Primary Tumor") + ylab("Counts") + xlab("Sample type") + theme_minimal() + 
  geom_text(stat='count', aes(label=after_stat(count)), hjust=-0.1, size=2.5) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))+
  scale_y_continuous(expand = c(0.1, 0.5)) + coord_flip()

setwd(resultsDirectory)
ggsave("Counts by primary disease set PT.png", plot = Counts_by_primary_disease_set_PT, width = 8, height = 5, units = "in")


# BY RACE
phenoData_TCGA_primary_split <- split(x = phenoData_TCGA_primary, f = phenoData_TCGA_primary$primary_disease_set)


races_list <- c("AMERICAN INDIAN OR ALASKA NATIVE", 
                "ASIAN", 
              "BLACK OR AFRICAN AMERICAN", 
              "WHITE")

length(names(phenoData_TCGA_primary_split))
# 27

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
sort(table(phenoData_TCGA_primary$primary_disease_set))
# min: cholangiocarcinoma (36)

n <- 21
subset_PT <- list()

set.seed(123)

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
        subset_PT[[i]] <- rbind(subset_PT[[i]], add_by_sex(n_random, races_split[[j]]))
        races_split[[j]] <- races_split[[j]][!races_split[[j]]$sample %in% subset_PT[[i]]$sample,]
      }
      if (nrow(subset_PT[[i]]) < n){
        if(n - nrow(subset_PT[[i]])==1){
          names_ult <- names(sort(sapply(races_split,nrow), decreasing = TRUE)[1])
          add <- lapply(races_split[names(races_split) %in% names_ult], function(x){
            x[sample(seq_len(nrow(x)), size = 1),]
          })
          subset_PT[[i]] <- rbind(subset_PT[[i]], Reduce(rbind,add))
        }else{
          names_ult <- names(sort(sapply(races_split,nrow), decreasing = TRUE)[1:(n - nrow(subset_PT[[i]]))])
          add <- lapply(races_split[names(races_split) %in% names_ult], function(x){
            x[sample(seq_len(nrow(x)), size = 1),]
          })
          subset_PT[[i]] <- rbind(subset_PT[[i]], Reduce(rbind,add))
        }
      }
      races_split <- NULL
    }
  }
  if (length(races_split)==1){
    n_samples <- n - nrow(subset_PT[[i]])
    subset_PT[[i]] <- rbind(subset_PT[[i]], add_by_sex(n_samples,races_split[[1]]))
  }
}



########### DEF FUNCTION
add_by_sex <- function(m, races.df) {
  subset.df <- data.frame()
  if (length(names(table(races.df$sex)))==1){
    index <- sample(seq_len(nrow(races.df)), size = m)
    subset.df  <- races.df[index,]
  }else{
    name_max <- names(which.max(table(races.df$sex)))[1]
    male_df <- races.df[races.df$sex == "MALE",]
    female_df <- races.df[races.df$sex == "FEMALE",]
    list_sex <- list(male_df, female_df)
    names(list_sex) <- c("MALE", "FEMALE")
    test <- any(lapply(list_sex, function(df) nrow(df) < floor(m/2)))
    if (test == TRUE){
      name_min <- names(which.min(sapply(list_sex, nrow)))
      subset.df <- rbind(subset.df, races.df[races.df$sex %in% name_min,])
      races.new.df <- races.df[!races.df$sex %in% name_min,]
      index <- sample(seq_len(nrow(races.new.df)), size = m-nrow(races.df[races.df$sex %in% name_min,]))
      subset.df <- rbind(subset.df, races.new.df[index,])
    } else {
      for (k in 1:length(list_sex)){
        if (names(list_sex)[k]==name_max){
          index <- sample(seq_len(nrow(list_sex[[k]])), size = ceiling(m/2))
          subset.df <- rbind(subset.df, list_sex[[k]][index,])
        } else{
          index <- sample(seq_len(nrow(list_sex[[k]])), size = floor(m/2))
          subset.df <- rbind(subset.df, list_sex[[k]][index,])
        }
      }
    }
  }
  return(subset.df)
}



#######################################
names(subset_PT) <- names(phenoData_TCGA_primary_split)

setwd(resultsDirectory)
pdf(file = "samples_subset_PT_race.pdf", width = 15, height = 15)
par(mfrow=c(5,6))
for (i in 1:length(names(subset_PT))){
  counts <- table(subset_PT[[i]]$race)
  counts.df <- as.data.frame(cbind(races_list, rep(0, length(races_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$races_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$races_list
  barplot(vector, 
          main = names(subset_PT)[i], 
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


sex_list <- c("MALE", 
                "FEMALE")

setwd(resultsDirectory)
pdf(file = "samples_subset_PT_sex.pdf", width = 15, height = 15)
par(mfrow=c(5,6))
for (i in 1:length(names(subset_PT))){
  counts <- table(subset_PT[[i]]$sex)
  counts.df <- as.data.frame(cbind(sex_list, rep(0, length(sex_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$sex_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$sex_list
  barplot(vector, 
          main = names(subset_PT)[i], 
          ylab = "Nº muestras", 
          col = c("deepskyblue3", "coral1"),
          cex.axis = 1,
          cex.lab = 1.4,
          cex.main = 1.5,
          cex.names =1.2,
          las = 2)
}

dev.off()


train_primary <- Reduce(rbind, subset_PT)

######################################
# Write train primary
######################################
setwd(tablesDirectory)
write.csv(train_primary, file = "train_primary.csv")

train_primary <- read.csv(file = "train_primary.csv")
train_primary <- train_primary[, colnames(train_primary) != "X"]
#####################################################################
# RNAseq train primary
#####################################################################

# Load FPKM for MET500
FPKM_MET500 <- fread("~/Desktop/LIDIA/TCGA_classification/Data/M.mx.txt", sep = "\t", stringsAsFactors = F)
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


# Split in training, test dataset and metastatic:
rnaseq_counts_primary_train <- rnaseq_counts[,colnames(rnaseq_counts) %in% train_primary$sample]
train_primary <- train_primary[train_primary$sample %in% colnames(rnaseq_counts_primary_train),]
train_primary <- train_primary[match(colnames(rnaseq_counts_primary_train),train_primary$sample),]


# Study sex
phenoData_TCGA_primary_split <- split(x = train_primary, f = train_primary$primary_disease_set)
lapply(phenoData_TCGA_primary_split, function(x){
  table(x$sex)
})
