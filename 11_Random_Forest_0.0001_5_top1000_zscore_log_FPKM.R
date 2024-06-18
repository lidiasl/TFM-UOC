#####################################################################
# Genes of pvalue adj < 0.0001 & log2FC>5 & top 1000 in Random Forest
# zscore_log_FPKM
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
library(fgsea)


# Directories
tablesDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Tables"
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/FINAL MODELS"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
setwd(tablesDirectory)
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")
zscore_log_fpkm <- read.csv(file="zscore_log_fpkm.csv")

rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

rownames(zscore_log_fpkm) <- zscore_log_fpkm$X
zscore_log_fpkm <- zscore_log_fpkm[, colnames(zscore_log_fpkm)!="X"]

genes_list_0.0001_5 <- readRDS(file="genes_list_0.0001_5.rds")
setwd(dataDirectory)
df_tissue_biopsy_tissue <- read.csv(file = "df_tissue_biopsy_tissue.csv", header = TRUE)

setwd("~/Desktop/LIDIA/TCGA_classification/TFM/Results/Human Protein")
genes_union_enhanced <- readRDS(file="genes_union_enhanced.rds")
######################################


genes_list_0.0001_5 <- lapply(genes_list_0.0001_5, function(x){
  Reduce(c,x)
})

genes_list_0.0001_5 <- Reduce(c,genes_list_0.0001_5) 
genes_list_0.0001_5 <- unique(genes_list_0.0001_5) #7339
length(genes_list_0.0001_5) #7339


table(genes_union_enhanced %in% genes_list_0.0001_5)
#FALSE  TRUE 
#1663  2474 

genes_intersect <- genes_union_enhanced[genes_union_enhanced %in% genes_list_0.0001_5]


zscore_log_fpkm.f <- zscore_log_fpkm[rownames(zscore_log_fpkm) %in% genes_intersect,]
zscore_log_fpkm.f <- zscore_log_fpkm.f[,colnames(zscore_log_fpkm.f) %in% phenoData_training_final$sample,]

# Top variables
zscore_log_fpkm.var <- zscore_log_fpkm.f[order(rowVars(as.matrix(zscore_log_fpkm.f)), decreasing = TRUE)[1:1000],]

zscore_log_fpkm.var <- t(zscore_log_fpkm.var)
zscore_log_fpkm.var <- as.data.frame(zscore_log_fpkm.var)

table(rownames(zscore_log_fpkm.var)==phenoData_training_final$sample)
zscore_log_fpkm.var$enfermedad_primaria <- phenoData_training_final$enfermedad_primaria[match(rownames(zscore_log_fpkm.var), phenoData_training_final$sample)]

genes_top1000 <- colnames(zscore_log_fpkm.var)
genes_top1000 <- genes_top1000[genes_top1000 != "enfermedad_primaria"]

setwd(resultsDirectory)
saveRDS(genes_top1000, file="genes_top1000.rds")
######################################
#### RANDOM FOREST
######################################
zscore_log_fpkm.var$enfermedad_primaria <- as.factor(zscore_log_fpkm.var$enfermedad_primaria)

####
tunegrid <- expand.grid(mtry=seq(10,285,25))
ntree_list <- seq(500,1500,100)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)

bestTune_list<- list()

library(doParallel)
registerDoParallel(cores = 16)
set.seed(42) 
toc <- Sys.time()
for (i in 1:length(ntree_list)){
  model_rf <- caret::train(enfermedad_primaria ~ .,
                           data = zscore_log_fpkm.var,
                           method = "rf", 
                           metric = "Accuracy", 
                           tuneGrid = tunegrid,
                           trControl = control,
                           ntree = ntree_list[i],
                           keep.forest = TRUE,
                           importance = TRUE)
  bestTune_list[[i]] <- model_rf[["results"]]
  print(i)
}


tic <- Sys.time()

tic-toc #Time difference of 4.097734 hours


for (i in 1:length(ntree_list)){
  bestTune_list[[i]]$ntree <- ntree_list[i]
}

bestTune_df <- Reduce(rbind, bestTune_list)
bestTune_df$ntree <- as.factor(bestTune_df$ntree)
ggplot(bestTune_df, aes(x=mtry, y=Kappa, col=ntree)) + geom_line()
ggplot(bestTune_df, aes(x=mtry, y=Accuracy, col=ntree)) + geom_line()

######################################
setwd(resultsDirectory)
write.csv(bestTune_df, file="bestTune_df_0.0001_5_top1000_zscore_log_FPKM.csv")


#### 

library(randomForest)
set.seed(42)
model_rf <- randomForest(x = as.matrix(zscore_log_fpkm.var[, colnames(zscore_log_fpkm.var)!="enfermedad_primaria"]),
                         y=zscore_log_fpkm.var$enfermedad_primaria,
                         mtry=10,
                         ntree=1400,
                         keep.forest = TRUE,
                         importance = TRUE)
predRF <- predict(model_rf)
tab_cm <- confusionMatrix(predRF,zscore_log_fpkm.var$enfermedad_primaria)

setwd(resultsDirectory)
save(model_rf,file = "model_rf.RData")
#Accuracy : 0.9421  
#Kappa : 0.9399

accuracy_train <- sum(diag(tab_cm$table))/sum(tab_cm$table)
sensibilidad_train <- mean(tab_cm$byClass[,1])
especificidad_train <- mean(tab_cm$byClass[,2])
F1_score_train <- mean(tab_cm$byClass[,7])
kappa_train <- 0.9399    
precision_train <- mean(tab_cm$byClass[,5])
NPV_train <- mean(tab_cm$byClass[,4])

######################################
########## RF VALIDATION
######################################

setwd(tablesDirectory)
zscore_log_FPKM_test <- read.csv(file="zscore_log_fpkm_test.csv")
phenodata_test <- read.csv(file="phenoData_test.csv")

rownames(zscore_log_FPKM_test) <- zscore_log_FPKM_test$X
zscore_log_FPKM_test <- zscore_log_FPKM_test[, colnames(zscore_log_FPKM_test)!="X"]

rownames(phenodata_test) <- phenodata_test$X
phenodata_test <- phenodata_test[, colnames(phenodata_test)!="X"]
phenodata_test$enfermedad_primaria <- as.factor(phenodata_test$enfermedad_primaria)

zscore_log_FPKM_test <- zscore_log_FPKM_test[rownames(zscore_log_FPKM_test) %in% colnames(zscore_log_fpkm.var),]
zscore_log_FPKM_test <- as.data.frame(t(zscore_log_FPKM_test))

table(colnames(zscore_log_FPKM_test)==colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"])
zscore_log_FPKM_test <- zscore_log_FPKM_test[, match(colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"], colnames(zscore_log_FPKM_test))]

pred_val <- predict(model_rf, zscore_log_FPKM_test)
tab_cm_VAL <- confusionMatrix(pred_val, phenodata_test$enfermedad_primaria)

#Accuracy : 0.8995
#Kappa : 0.8956 

accuracy_val <- sum(diag(tab_cm_VAL$table))/sum(tab_cm_VAL$table)
sensibilidad_val <- mean(tab_cm_VAL$byClass[,1])
especificidad_val <- mean(tab_cm_VAL$byClass[,2])
F1_score_val <- mean(tab_cm_VAL$byClass[,7])
kappa_val <- 0.8956 
precision_val <- mean(tab_cm_VAL$byClass[,5])
NPV_val <- mean(tab_cm_VAL$byClass[,4])

######################################
setwd(tablesDirectory)
zscore_log_FPKM_MET500 <- read.csv(file="zscore_log_fpkm_MET500_clean.csv")
phenodata_MET500_clean <- read.csv(file="phenodata_MET500_clean.csv")

rownames(zscore_log_FPKM_MET500) <- zscore_log_FPKM_MET500$X
zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, colnames(zscore_log_FPKM_MET500)!="X"]
colnames(zscore_log_FPKM_MET500) <- str_replace_all(colnames(zscore_log_FPKM_MET500), "[.]", "-")
table(colnames(zscore_log_FPKM_MET500) %in% phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, match(phenodata_MET500_clean$Sample_id, colnames(zscore_log_FPKM_MET500))]
table(colnames(zscore_log_FPKM_MET500)==phenodata_MET500_clean$Sample_id)
#TRUE 
#387 

phenodata_MET500_clean$enfermedad_primaria <- df_tissue_biopsy_tissue$spanish[match(phenodata_MET500_clean$tissue_set, df_tissue_biopsy_tissue$set)]
phenodata_MET500_clean$enfermedad_primaria <- as.factor(phenodata_MET500_clean$enfermedad_primaria)

zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[rownames(zscore_log_FPKM_MET500) %in% colnames(zscore_log_fpkm.var),]
zscore_log_FPKM_MET500 <- as.data.frame(t(zscore_log_FPKM_MET500))

table(colnames(zscore_log_FPKM_MET500)==colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"])
zscore_log_FPKM_MET500 <- zscore_log_FPKM_MET500[, match(colnames(zscore_log_fpkm.var)[colnames(zscore_log_fpkm.var)!="enfermedad_primaria"], colnames(zscore_log_FPKM_MET500))]


predRF_MET500 <- predict(model_rf, zscore_log_FPKM_MET500)
table(levels(phenodata_MET500_clean$enfermedad_primaria) %in% levels(predRF_MET500))
phenodata_MET500_clean$enfermedad_primaria <- factor(phenodata_MET500_clean$enfermedad_primaria, levels=levels(predRF_MET500))
tab_cm_MET500 <- confusionMatrix(predRF_MET500, phenodata_MET500_clean$enfermedad_primaria)

#Accuracy : 0.6305
#Kappa : 0.5789   

tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)
tejidos_MET500 <- paste("Class:",tejidos_MET500)

accuracy_MET500 <- sum(diag(tab_cm_MET500$table))/sum(tab_cm_MET500$table)
sensibilidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,1])
especificidad_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,2])
F1_score_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,7])
kappa_MET500 <- 0.5789 
precision_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,5])
NPV_MET500 <- mean(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass)%in%tejidos_MET500,4])

######################################
setwd(resultsDirectory)
saveRDS(tab_cm_MET500, file="tab_cm_MET500_RF_top1000_zscore_log_FPKM.rds")


######################################
#### STUDY MODEL
######################################
setwd(resultsDirectory)
load(file = "model_rf.RData")

importance <- model_rf$importance
setwd(resultsDirectory)
write.csv(importance, "importance_RF.csv")

varImpPlot(model_rf,main = "Importancia de las variables (Random Forest)")


######################################
#### STUDY BY CLASS
######################################
setwd(resultsDirectory)
tab_cm_MET500 <- readRDS(file="tab_cm_MET500_RF_top1000_zscore_log_FPKM.rds")
tejidos_MET500 <- unique(phenodata_MET500_clean$tejido)

tab_cm_MET500$table <- tab_cm_MET500$table[,colnames(tab_cm_MET500$table) %in% tejidos_MET500]
table_cm <- round(prop.table(tab_cm_MET500$table, margin=2),4)*100
table_cm <- matrix(table_cm, nrow=nrow(table_cm), ncol=ncol(table_cm))
rownames(table_cm) <- rownames(prop.table(tab_cm_MET500$table, margin=2))
colnames(table_cm) <- colnames(prop.table(tab_cm_MET500$table, margin=2))


library(plot.matrix)
colors <- colorRampPalette(c("white","red"))(20)
pdf(file="tabla_cm_RF.pdf")
par(mar=c(8,8,8,8),cex.axis=0.7)
plot(table_cm, las=2, digits=0, text.cell=list(cex=0.5), fmt.key="%1f", col=colors)
dev.off()

table_numbers <- matrix(tab_cm_MET500$table, nrow=nrow(tab_cm_MET500$table), ncol=ncol(tab_cm_MET500$table))
rownames(table_numbers) <- rownames(tab_cm_MET500$table)
colnames(table_numbers) <- colnames(tab_cm_MET500$table)
library(plot.matrix)
colors <- colorRampPalette(c("white","red"))(3)
pdf(file="tabla_cm_RF_numbers.pdf")
par(mar=c(8,8,8,8),cex.axis=0.7)
plot(table_numbers, las=2, digits=0, text.cell=list(cex=0.5), fmt.key="%1f", col="white")
dev.off()


metrics_MET500 <- tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass) %in% tejidos_MET500,]

metrics_training <- tab_cm$byClass

metrics_VAL <- tab_cm_VAL$byClass


# Adrenocortical
metrics_training_adrenocortical <- metrics_training[1,c(1,2,3,4,5,7)]
metrics_VAL_adrenocortical <- metrics_VAL[1,c(1,2,3,4,5,7)]
metrics_MET500_adrenocortical <- metrics_MET500[1,c(1,2,3,4,5,7)]
df_adrenocortical <- rbind(metrics_training_adrenocortical, metrics_VAL_adrenocortical,metrics_VAL_adrenocortical)
df_adrenocortical <- as.data.frame(df_adrenocortical)

df_adrenocortical <- rbind(rep(0,6),df_adrenocortical)
df_adrenocortical <- rbind(rep(1,6),df_adrenocortical)
rownames(df_adrenocortical)[1:2] <- c("Max","Min")
library(fmsb)
par(mar=c(0,0,0,0))
radarchart(df_adrenocortical, 
           axistype = 1,
           pcol = c("red", "blue", "green"),   
           plwd = 2,                
           cglcol = "grey",         
           cglty = 1,             
           axislabcol = "grey",   
           caxislabels = c(0.00, 0.25, 0.50, 0.75,1.00), 
           vlcex = 1.2,         
           cex.main = 1.5,       
           cex.lab = 1.2)     
legend(x = 1, y = 1, legend = c("Entrenamiento","Test", "MET500"), 
       col = c("red", "blue", "green"), lty = 1, lwd = 2)

library(fmsb)

for (i in 1:length(tejidos_MET500)){
  
  metrics_training_i <- metrics_training[rownames(metrics_training)==tejidos_MET500[i],c(1,2,3,4,5,7)]
  metrics_VAL_i <- metrics_VAL[rownames(metrics_VAL)==tejidos_MET500[i],c(1,2,3,4,5,7)]
  metrics_MET500_i <- metrics_MET500[rownames(metrics_MET500)==tejidos_MET500[i],c(1,2,3,4,5,7)]
  df <- rbind(metrics_training_i, metrics_VAL_i,metrics_MET500_i)
  df <- as.data.frame(df)
  
  df <- rbind(rep(0,6),df)
  df <- rbind(rep(1,6),df)
  rownames(df)[1:2] <- c("Max","Min")
  colnames(df) <- c("Sensibilidad", "Especificidad", "PPV", "NPV", "Precisión", "F1")
  pdf(file=paste("radar plots RF", tejidos_MET500[i],".pdf" ), width=8)
  radarchart(df, 
             axistype = 1,
             pcol = c("#984EA3", "#FFD92F", "#E41A1C"),   
             plwd = 5,              
             cglcol = "grey",        
             cglty = 1,           
             axislabcol = "grey",  
             caxislabels = c(0.00, 0.25, 0.50, 0.75,1.00),
             vlcex =1.5,              
             cex.main = 1.5,           
             cex.lab = 1.2,
             cex=4)
  title(main = sub("Class:", "", tejidos_MET500[i]), cex.main = 2)
  dev.off()
}

radarchart(df, 
           axistype = 1,
           pcol = c("#984EA3", "#FFD92F", "#E41A1C"),  
           plwd = 5,                 
           cglcol = "grey",         
           cglty = 1,                
           axislabcol = "grey",    
           caxislabels = c(0.00, 0.25, 0.50, 0.75,1.00), 
           vlcex =1.5,             
           cex.main = 1.5,           
           cex.lab = 1.2,
           cex=4)
title(main = sub("Class:", "", tejidos_MET500[i]), cex.main = 2)
legend(x = 1, y = 1, legend = c("Entrenamiento","Test", "MET500"), 
       col = c("#984EA3", "#FFD92F", "#E41A1C"), lty = c(1,2,3), lwd = 2)


df_byClass_training <- round(tab_cm$byClass[,c(1,2,3,4,5,7)],4)
rownames(df_byClass_training) <- sub("Class:", "", rownames(df_byClass_training))
colnames(df_byClass_training) <- c("Sensibilidad", "Especificidad", "PPV", "NPV", "Precisión", "F1-score")
df_byClass_VAL <- round(tab_cm_VAL$byClass[,c(1,2,3,4,5,7)],4)
rownames(df_byClass_VAL) <- sub("Class:", "", rownames(df_byClass_VAL))
colnames(df_byClass_VAL) <- c("Sensibilidad", "Especificidad", "PPV", "NPV", "Precisión", "F1-score")
df_byClass_MET500 <- round(tab_cm_MET500$byClass[rownames(tab_cm_MET500$byClass) %in% tejidos_MET500,c(1,2,3,4,5,7)],4)
rownames(df_byClass_MET500) <- sub("Class:", "", rownames(df_byClass_MET500))
colnames(df_byClass_MET500) <- c("Sensibilidad", "Especificidad", "PPV", "NPV", "Precisión", "F1-score")

setwd(resultsDirectory)
write.csv(df_byClass_training, file="df_byClass_training_RF.csv")
write.csv(df_byClass_VAL, file="df_byClass_VAL_RF.csv")
write.csv(df_byClass_MET500, file="df_byClass_MET500_RF.csv")
######################################
#### STUDY SITE OF BIOPSY
######################################
df_biopsia <- as.data.frame(predRF_MET500)
df_biopsia <- cbind(df_biopsia, phenodata_MET500_clean$enfermedad_primaria, phenodata_MET500_clean$biopsia_tejido)
colnames(df_biopsia) <- c("Prediccion", "Enfermedad_primaria", "Tejido_biopsia")
table(df_biopsia$Prediccion==df_biopsia$Enfermedad_primaria)
#FALSE  TRUE 
#143   244 

df_biopsia <- df_biopsia[df_biopsia$Prediccion!=df_biopsia$Enfermedad_primaria,]
table(df_biopsia$Prediccion==df_biopsia$Tejido_biopsia)
#FALSE  TRUE 
#114    29 

predRF_MET500_prob <- predict(model_rf, zscore_log_FPKM_MET500, type="prob")
predRF_MET500_prob <- as.data.frame(predRF_MET500_prob)
predRF_MET500_prob$prediccion <- df_biopsia$predRF_MET500
predRF_MET500_prob$enfermedad_primaria <- phenodata_MET500_clean$enfermedad_primaria
