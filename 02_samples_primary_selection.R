#####################################################################
# Samples study Primary Tumor and selection of samples
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Script 02"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"


######################################
# Read data
######################################
setwd(tablesDirectory)
phenoData_TCGA <- read.csv(file = "phenoData_TCGA.csv")
phenoData_TCGA_in_rnaseq_with_gender <- read.csv(file = "phenoData_TCGA_in_rnaseq_with_gender.csv")
phenoData_TCGA_primary_f <- read.csv(file = "phenoData_TCGA_primary_f.csv")
rnaseq_counts_primary_f <- read.csv(file = "rnaseq_counts_primary_f.csv")
log_rnaseq_counts_primary_f <- read.csv(file = "log_rnaseq_counts_primary_f.csv")

rownames(rnaseq_counts_primary_f) <- rnaseq_counts_primary_f$X
rnaseq_counts_primary_f <- rnaseq_counts_primary_f[, colnames(rnaseq_counts_primary_f)!="X"]

rownames(log_rnaseq_counts_primary_f) <- log_rnaseq_counts_primary_f$X 
log_rnaseq_counts_primary_f <- log_rnaseq_counts_primary_f[, colnames(log_rnaseq_counts_primary_f)!="X"]

phenoData_TCGA_primary_f <- phenoData_TCGA_primary_f[, colnames(phenoData_TCGA_primary_f)!="X"]
rownames(phenoData_TCGA_primary_f) <- phenoData_TCGA_primary_f$sample

phenoData_TCGA_in_rnaseq_with_gender <- phenoData_TCGA_in_rnaseq_with_gender[, colnames(phenoData_TCGA_in_rnaseq_with_gender)!="X"]
rownames(phenoData_TCGA_in_rnaseq_with_gender) <- phenoData_TCGA_in_rnaseq_with_gender$sample

setwd(dataDirectory)
set_subset <- read.csv(file = "set_subset.csv", header = TRUE)

# Load FPKM for MET500
FPKM_MET500 <- fread("~/Desktop/LIDIA/TCGA_classification/Data/M.mx.txt", sep = "\t", stringsAsFactors = F)
genenames_MET500 <- separate(data = FPKM_MET500, col = sample, into = c("GeneName", "version"), remove = T)
FPKM_MET500 <- FPKM_MET500[,-1]
FPKM_MET500 <- as.data.frame(FPKM_MET500)
rownames(FPKM_MET500) <- genenames_MET500$GeneName

# Filter genes common in both datasets:
table(rownames(rnaseq_counts_primary_f) %in% rownames(FPKM_MET500))
#FALSE  TRUE 
#41435 19063 

rnaseq_counts_primary_f <- rnaseq_counts_primary_f[rownames(rnaseq_counts_primary_f) %in% rownames(FPKM_MET500),]

#Check:
table(rownames(rnaseq_counts_primary_f) %in% rownames(FPKM_MET500))
#TRUE 
#19063 


#####################################################################
# Filter breast cancer and select only FEMALE
#####################################################################

phenoData_breast <- phenoData_TCGA_primary_split$breast

table(phenoData_breast$sex)
#FEMALE   MALE 
#974     11 

setwd(tablesDirectory)

BRCA_subtypes <- read.csv(file = "BRCA_subtypes.csv")
BRCA_subtypes$sampleID <- str_replace_all(BRCA_subtypes$sampleID, "-", "_")

table(phenoData_breast$sample %in% BRCA_subtypes$sampleID)

counts_by_sample <- colSums(rnaseq_counts_primary_f[,match(phenoData_breast$sample, colnames(rnaseq_counts_primary_f))])
counts_by_sample.df <- as.data.frame(counts_by_sample)

table(rownames(counts_by_sample.df)==phenoData_breast$sample)
phenoData_breast$counts <- counts_by_sample.df$counts_by_sample

phenoData_breast$subtype <- BRCA_subtypes$PR_Status_nature2012[match(phenoData_breast$sample, BRCA_subtypes$sampleID)]

phenoData_breast <- phenoData_breast[phenoData_breast$subtype!="",]

### See if the samples differenciate by sex
#PCA log counts
PCA_log_counts_breast <- prcomp(t(as.matrix(log_rnaseq_counts_primary_f[,match(phenoData_breast$sample, 
                                                                               colnames(log_rnaseq_counts_primary_f))])))

autoplot(PCA_log_counts_breast, 
         label = F, 
         data = as.data.frame(phenoData_breast), colour = "subtype", main = "PCA Subtype Breast Cancer", legend = F)


names_male_breast <- rownames(phenoData_TCGA_primary_split$breast[phenoData_TCGA_primary_split$breast$sex=="MALE",])
table(names_male_breast %in% phenoData_breast$sample)
#FALSE  TRUE 
#3     8 

setwd(tablesDirectory)
write.csv(phenoData_breast, file = "phenoData_breast.csv")


#####################################################################
# Create subset of samples by primary disease set
#####################################################################
phenoData_TCGA_primary_split <- split(x = phenoData_TCGA_primary_f, f = phenoData_TCGA_primary_f$primary_disease_set)

sort(table(phenoData_TCGA_primary_f$primary_disease_set))
# min: cholangiocarcinoma (36)

n <- 36
subset_PT <- list()

set.seed(42)

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



########### DEF FUNCTION SELECTION BY SEX
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

subset_PT_español <- subset_PT
names(subset_PT_español) <- set_subset$spanish[match(names(subset_PT_español), set_subset$set)]

######################################
# Write subset_PT
######################################
setwd(tablesDirectory)
saveRDS(subset_PT, file = "subset_PT.rds")
saveRDS(subset_PT_español, file = "subset_PT_español.rds")

######################################

races_list <- c("BLANCO", "NEGRO_O_AFROAMERICANO", "ASIÁTICO")


setwd(resultsDirectory)
pdf(file = "muestras_selecPT_raza.pdf", width = 8, height = 12)
par(mfrow=c(6,5))
for (i in 1:length(names(subset_PT_español))){
  counts <- table(subset_PT_español[[i]]$raza)
  counts.df <- as.data.frame(cbind(races_list, rep(0, length(races_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$races_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$races_list
  barplot(vector, 
          main = names(subset_PT_español)[i], 
          ylab = "Nº muestras", 
          col = c("#482677FF", "#FDE725FF", "#287D8EFF"),
          names.arg = c("BLANCO"="BLANCO", 
                        "NEGRO_O_AFROAMERICANO"="NEGRO", 
                        "ASIÁTICO"="ASIÁTICO"),
          cex.axis = 1,
          cex.lab = 1,
          cex.main = 1,
          cex.names =0.8,
          las = 2)
}

dev.off()


sex_list <- c("HOMBRE", 
              "MUJER")


setwd(resultsDirectory)
pdf(file = "muestras_selecPT_sexo.pdf", width = 8, height = 12)
par(mfrow=c(5,6))
for (i in 1:length(names(subset_PT_español))){
  counts <- table(subset_PT_español[[i]]$sexo)
  counts.df <- as.data.frame(cbind(sex_list, rep(0, length(sex_list))))
  for (j in (names(counts))){
    counts.df$V2[counts.df$sex_list==j] <- counts[j]
  }
  vector <- as.numeric(counts.df$V2)
  names(vector) <- counts.df$sex_list
  barplot(vector, 
          main = names(subset_PT_español)[i], 
          ylab = "Nº muestras", 
          col = c("deepskyblue3", "coral1"),
          cex.axis = 1,
          cex.lab = 1,
          cex.main = 1,
          cex.names =0.8,
          las = 2)
}

dev.off()


phenoData_select_primary <- Reduce(rbind, subset_PT)

table(phenoData_select_primary$sample[phenoData_select_primary$primary_disease_set=="breast"] %in% BRCA_subtypes$sampleID)
#TRUE 
#36 

BRCA_subtypes$PR_Status_nature2012[BRCA_subtypes$sampleID %in% phenoData_select_primary$sample[phenoData_select_primary$primary_disease_set=="breast"]]


######################################
# Write select primary
######################################
setwd(tablesDirectory)
write.csv(phenoData_select_primary, file = "phenoData_select_primary.csv")

######################################



######################################
# Filter by training samples:
######################################
log_rnaseq_primary_select <- log_rnaseq_counts_primary_f[,colnames(log_rnaseq_counts_primary_f) %in% phenoData_select_primary$sample]
rnaseq_counts_primary_select <- rnaseq_counts_primary_f[,colnames(rnaseq_counts_primary_f) %in% phenoData_select_primary$sample]

log_rnaseq_primary_select <- log_rnaseq_primary_select[, match(phenoData_select_primary$sample, colnames(log_rnaseq_primary_select))]
rnaseq_counts_primary_select <- rnaseq_counts_primary_select[, match(phenoData_select_primary$sample, colnames(rnaseq_counts_primary_select))]

log_rnaseq_primary_select <- log_rnaseq_primary_select[rownames(log_rnaseq_primary_select) %in% rownames(rnaseq_counts_primary_select),]

######################################
# Write rnaseq counts and log counts
######################################
setwd(tablesDirectory)
write.csv(log_rnaseq_primary_select, file = "log_rnaseq_primary_select.csv")
write.csv(rnaseq_counts_primary_select, file = "rnaseq_counts_primary_select.csv")
######################################


################################################################################
# GRAPHS
################################################################################

counts_by_sample_PT <- colSums(rnaseq_counts_primary_train)
counts_by_sample_PT.df <- as.data.frame(counts_by_sample_PT)

table(rownames(counts_by_sample_PT.df)==phenoData_train_primary$sample)
counts_by_sample_PT.df$primary_disease_set <- phenoData_train_primary$primary_disease_set

boxplot_counts_by_sample_by_PT_train <- ggplot(data = counts_by_sample_PT.df, aes(x=counts_by_sample_PT, y= primary_disease_set)) + 
  geom_boxplot() 

setwd(resultsDirectory)
ggsave("Boxplot counts by sample by PT TRAIN.png", plot = boxplot_counts_by_sample_by_PT_train, width = 8, height = 5, units = "in")

######################################

setwd(resultsDirectory)
pdf(file = "hist_selecPT.pdf", width = 10, height = 8)
hist(as.matrix(rnaseq_counts_primary_train), breaks = 1000, main = "Count distribution by sample")
dev.off()


setwd(resultsDirectory)
#RNA-seq count distribution
pdf(file = "hist_selecPT_one_sample.pdf", width = 10, height = 8)
ggplot(rnaseq_counts_primary_train, aes(x = TCGA_OR_A5J6_01)) +
  geom_histogram(stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  ggtitle("TCGA_OR_A5J6_01")
#If we zoom in close to zero, we can see a large number of genes with counts of zero:
ggplot(rnaseq_counts_primary_train, aes(x = TCGA_OR_A5J6_01)) +
  geom_histogram(stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  ggtitle("TCGA_OR_A5J6_01") +
  xlim(-5, 500)
dev.off()

###################
# Now the histogram for a specific gene:
setwd(resultsDirectory)
pdf(file = "hist_selecPT_one_gene.pdf", width = 10, height = 8)
rnaseq_counts_primary_train.t <- t(rnaseq_counts_primary_train)
ggplot(rnaseq_counts_primary_train.t, aes(x = ENSG00000167578)) +
  geom_histogram(stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of samples") +
  ggtitle("ENSG00000167578")
dev.off()

###################
# PCA
###################

PCA_selecPT <- autoplot(prcomp(t(as.matrix(rnaseq_counts_primary_train))), label = F, 
         data = as.data.frame(phenoData_train_primary), colour = "primary_disease_set", main = "PCA", legend = F)

setwd(resultsDirectory)
ggsave("PCA_selecPT.png", plot = PCA_selecPT, width = 8, height = 6, units = "in")


PCA_selecPT.log <- autoplot(prcomp(t(as.matrix(log_rnaseq_primary_train))), label = F, 
         data = as.data.frame(phenoData_train_primary), colour = "primary_disease_set", main = "PCA", legend = F)

setwd(resultsDirectory)
ggsave("PCA_selecPT_log.png", plot = PCA_selecPT.log, width = 8, height = 6, units = "in")



###################
# UMAP log
###################

umap.defaults
#n_neighbors: 15
#min_dist: 0.1

# Bucle listas
n_neighbors_list <- seq(5,100,10)
min_dist_list <- seq(0.1,0.5,0.1)

custom.settings = umap.defaults

list_UMAP <- list()
k<-1

for(i in 1:length(n_neighbors_list)){
  for(j in 1:length(min_dist_list)){
    custom.settings$n_neighbors = n_neighbors_list[i]
    custom.settings$min_dist=min_dist_list[j]
    
    log.rnaseq.umap <- umap(t(log_rnaseq_primary_train), config=custom.settings)
    
    df.log.umap <- data.frame(x = log.rnaseq.umap$layout[,1],
                                 y = log.rnaseq.umap$layout[,2],
                                 Group = phenoData_train_primary$primary_disease_set)
    
    list_UMAP[[k]] <- ggplot(df.log.umap, aes(x, y, colour = Group)) +
      geom_point() + ggtitle("UMAP log rnaseq")
    k <- k+1
  }
}


setwd(resultsDirectory)
pdf(file = "UMAP_log_rnaseq.pdf", width = 100, height = 100)
grid.arrange(list_UMAP[[1]],list_UMAP[[2]],list_UMAP[[3]],list_UMAP[[4]],list_UMAP[[5]],list_UMAP[[6]],list_UMAP[[7]],list_UMAP[[8]],list_UMAP[[9]], list_UMAP[[10]],
             list_UMAP[[11]],list_UMAP[[12]],list_UMAP[[13]],list_UMAP[[14]],list_UMAP[[15]],list_UMAP[[16]],list_UMAP[[17]],list_UMAP[[18]],list_UMAP[[19]],list_UMAP[[20]],
             list_UMAP[[21]],list_UMAP[[22]],list_UMAP[[23]],list_UMAP[[24]],list_UMAP[[25]],list_UMAP[[26]],list_UMAP[[27]],list_UMAP[[28]],list_UMAP[[29]], list_UMAP[[30]],
             list_UMAP[[31]],list_UMAP[[32]],list_UMAP[[33]],list_UMAP[[34]],list_UMAP[[35]],list_UMAP[[36]],list_UMAP[[37]],list_UMAP[[38]],list_UMAP[[39]], list_UMAP[[40]],
             list_UMAP[[41]],list_UMAP[[42]],list_UMAP[[43]],list_UMAP[[44]],list_UMAP[[45]],list_UMAP[[46]],list_UMAP[[47]],list_UMAP[[48]],list_UMAP[[49]], list_UMAP[[50]],
             ncol=5)
dev.off()

#########################

setwd(resultsDirectory)
ggsave("UMAP_log_rnaseq5.png", plot = list_UMAP[[5]], width = 10, height = 6, units = "in")

ggsave("UMAP_log_rnaseq10.png", plot = list_UMAP[[10]], width = 10, height = 6, units = "in")

#########################

custom.settings = umap.defaults

list_UMAP_counts <- list()
k<-1

for(i in 1:length(n_neighbors_list)){
  for(j in 1:length(min_dist_list)){
    custom.settings$n_neighbors = n_neighbors_list[i]
    custom.settings$min_dist=min_dist_list[j]
    
    rnaseq.umap <- umap(t(rnaseq_counts_primary_train), config=custom.settings)
    
    df.umap <- data.frame(x = rnaseq.umap$layout[,1],
                              y = rnaseq.umap$layout[,2],
                              Group = phenoData_train_primary$primary_disease_set)
    
    list_UMAP_counts[[k]] <- ggplot(df.umap, aes(x, y, colour = Group)) +
      geom_point() + ggtitle("UMAP rnaseq")
    k <- k+1
  }
}

setwd(resultsDirectory)
pdf(file = "UMAP_rnaseq.pdf", width = 100, height = 100)
grid.arrange(list_UMAP_counts[[1]],list_UMAP_counts[[2]],list_UMAP_counts[[3]],list_UMAP_counts[[4]],list_UMAP_counts[[5]],
             list_UMAP_counts[[6]],list_UMAP_counts[[7]],list_UMAP_counts[[8]],list_UMAP_counts[[9]], list_UMAP_counts[[10]],
             list_UMAP_counts[[11]],list_UMAP_counts[[12]],list_UMAP_counts[[13]],list_UMAP_counts[[14]],list_UMAP_counts[[15]],
             list_UMAP_counts[[16]],list_UMAP_counts[[17]],list_UMAP_counts[[18]],list_UMAP_counts[[19]],list_UMAP_counts[[20]],
             list_UMAP_counts[[21]],list_UMAP_counts[[22]],list_UMAP_counts[[23]],list_UMAP_counts[[24]],list_UMAP_counts[[25]],
             list_UMAP_counts[[26]],list_UMAP_counts[[27]],list_UMAP_counts[[28]],list_UMAP_counts[[29]], list_UMAP_counts[[30]],
             list_UMAP_counts[[31]],list_UMAP_counts[[32]],list_UMAP_counts[[33]],list_UMAP_counts[[34]],list_UMAP_counts[[35]],
             list_UMAP_counts[[36]],list_UMAP_counts[[37]],list_UMAP_counts[[38]],list_UMAP_counts[[39]], list_UMAP_counts[[40]],
             list_UMAP_counts[[41]],list_UMAP_counts[[42]],list_UMAP_counts[[43]],list_UMAP_counts[[44]],list_UMAP_counts[[45]],
             list_UMAP_counts[[46]],list_UMAP_counts[[47]],list_UMAP_counts[[48]],list_UMAP_counts[[49]], list_UMAP_counts[[50]],
             ncol=5)
dev.off()

setwd(resultsDirectory)
ggsave("UMAP_rnaseq4.png", plot = list_UMAP_counts[[4]], width = 10, height = 6, units = "in")

ggsave("UMAP_rnaseq18.png", plot = list_UMAP_counts[[18]], width = 10, height = 6, units = "in")


###################
# Check if mean < variance (to use NB model) for each subset of replicates
# for each gene

mean_counts <- apply(rnaseq_counts_primary_train[,colnames(rnaseq_counts_primary_train) %in%
                                                   phenoData_train_primary$sample[phenoData_train_primary$primary_disease_set == "lung"]], 1, mean)

variance_counts <- apply(rnaseq_counts_primary_train[,colnames(rnaseq_counts_primary_train) %in%
                                                       phenoData_train_primary$sample[phenoData_train_primary$primary_disease_set == "lung"]], 1, var)

df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts), color="red") +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Variance vs mean") +
  xlab("mean") +
  ylab("variance")


