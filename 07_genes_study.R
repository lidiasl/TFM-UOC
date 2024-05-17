#####################################################################
# Study what important genes are in my results
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/Genes study"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# Genes list
setwd(tablesDirectory)
genes_list <- readRDS(file="genes_list.rds")
# phenoData
phenoData_training_final <- read.csv(file = "phenoData_training_final.csv")


rownames(phenoData_training_final) <- phenoData_training_final$X
phenoData_training_final <- phenoData_training_final[, colnames(phenoData_training_final)!="X"]

######################################

#### Adrenocortical
# https://www.proteinatlas.org/humanproteome/tissue/adrenal+gland#protein_expression_of_genes_elevated_in_adrenal_gland
#CYP11B2, CYP11B1, CYP17A1, HSD3B2, GML, CYP21A2, CCN3, MC2R, CYP11A1, ADGRV1, STAR, AKR1B1
adrenocortical_genes_reference <- c("CYP11B2", "CYP11B1", "CYP17A1", "HSD3B2", "GML", "CYP21A2", "CCN3", "MC2R", "CYP11A1", "ADGRV1", "STAR", "AKR1B1")

#### B_cell_lymphoma
# https://www.proteinatlas.org/humanproteome/tissue/lymphoid+tissue
#CD1B, PSMB11, RAG1, TRBV12-5, TRAJ61, CD1E, TRBV15, PRSS16, SH2D1A, TRBV6-1, PTCRA, TRBV12-4
#B_cell_lymphoma_genes_reference <- c("CD1B", "PSMB11", "RAG1", "TRBV12-5", "TRAJ61", "CD1E", "TRBV15", "PRSS16", "SH2D1A", "TRBV6-1", "PTCRA", "TRBV12-4")


#### bladder
#https://www.proteinatlas.org/humanproteome/tissue/urinary+bladder
#IGHV4-34, UPK2, IGHV5-51, MMP13, TNFAIP6, UPK1A, IL23A
bladder_genes_reference <- c("IGHV4-34", "UPK2", "IGHV5-51", "MMP13", "TNFAIP6", "UPK1A", "IL23A")

#### brain
#https://www.proteinatlas.org/humanproteome/brain/human+brain
#HCRT, AVP, PMCH, GRM4, BARHL1, GABRA6, MOG, OXT, FGF3, NEUROD6, SLC6A3, HAPLN2
brain_genes_reference <- c("HCRT", "AVP", "PMCH", "GRM4", "BARHL1", "GABRA6", "MOG", "OXT", "FGF3", "NEUROD6", "SLC6A3", "HAPLN2")

#### breast
#https://www.proteinatlas.org/humanproteome/tissue/breast
#LALBA, CSN2, CSN1S1, SULT1C3, BTN1A1, ANKRD30A, UGT2B11, ACSM1, UGT2B28, ABCC11, SERHL2, CST9
breast_genes_reference <- c("LALBA", "CSN2", "CSN1S1", "SULT1C3", "BTN1A1", "ANKRD30A", "UGT2B11", "ACSM1", "UGT2B28", "ABCC11", "SERHL2", "CST9")

#### cervical
#https://www.proteinatlas.org/humanproteome/tissue/cervix#protein_expression_of_genes_elevated_in_cervix
#MUC16, SLPI, WFDC2, TMPRSS11D, IVL, ESR1, ESR1, COL1A2, COL3A1, HOXA11, ASRGL1, MSX1, BPIFB1
cervical_genes_reference <- c("MUC16", "SLPI", "WFDC2", "TMPRSS11D", "IVL", "ESR1", "ESR1", "COL1A2", "COL3A1", "HOXA11", "ASRGL1", "MSX1", "BPIFB1")

#### cholangiocarcinoma
#

#### colon
#https://www.proteinatlas.org/humanproteome/tissue/intestine
#TMPRSS15, DEFA6, DEFA5, MLN, RBP2, GIP, S100G, LCT, FABP2, SI, ALPI, INSL5
colon_genes_reference <- c("TMPRSS15", "DEFA6", "DEFA5", "MLN", "RBP2", "GIP", "S100G", "LCT", "FABP2", "SI", "ALPI", "INSL5")

#### esophageal
#https://www.proteinatlas.org/humanproteome/tissue/esophagus
#CAPN14, UGT1A7, DYNAP, KRT4, CSTB, TGM3, LYPD2, MUC22, PADI1, TGM1, GBP6, ADH7
esophageal_genes_reference <- c("CAPN14", "UGT1A7", "DYNAP", "KRT4", "CSTB", "TGM3", "LYPD2", "MUC22", "PADI1", "TGM1", "GBP6", "ADH7")

#### head_neck
#

#### kidney
#https://www.proteinatlas.org/humanproteome/tissue/kidney
#SLC12A1, TMEM174, UMOD, SLC7A13, SLC22A12, NPHS2, SLC12A3, TMEM207, SLC22A2, SLC6A18, MCCD1, SLC34A1
kidney_genes_reference <- c("SLC12A1", "TMEM174", "UMOD", "SLC7A13", "SLC22A12", "NPHS2", "SLC12A3", "TMEM207", "SLC22A2", "SLC6A18", "MCCD1", "SLC34A1")

#### leukemia

#### liver
#https://www.proteinatlas.org/humanproteome/tissue/liver
#SPP2, AHSG, CFHR2, MBL2, F9, CFHR5, A1BG, SERPINC1, APOA2, F2, SLC10A1, CFHR3
liver_genes_reference <- c("SPP2", "AHSG", "CFHR2", "MBL2", "F9", "CFHR5", "A1BG", "SERPINC1", "APOA2", "F2", "SLC10A1", "CFHR3")

#### lung
#https://www.proteinatlas.org/humanproteome/tissue/lung
#SFTPA1, SFTPC, SFTPA2, SCGB3A2, SFTPB, SFTPD, AGER, NAPSA, MS4A15, RTKN2, SCGB1A1, SLC34A2
lung_genes_reference <- c("SFTPA1", "SFTPC", "SFTPA2", "SCGB3A2", "SFTPB", "SFTPD", "AGER", "NAPSA", "MS4A15", "RTKN2", "SCGB1A1", "SLC34A2")

#### mesothelioma

#### ovarian
#https://www.proteinatlas.org/humanproteome/tissue/ovary
#PADI6, ZP4, KLHDC8A, LHX9, GREB1
ovarian_genes_reference <- c("PADI6", "ZP4", "KLHDC8A", "LHX9", "GREB1")

#### pancreas
#https://www.proteinatlas.org/humanproteome/tissue/pancreas
#CELA2A, CPA1, CTRB1, CELA3A, PNLIP, AMY2A, PRSS1, CLPS, CTRC, PNLIPRP1, CELA2B, CTRB2
pancreas_genes_reference <- c("CELA2A", "CPA1", "CTRB1", "CELA3A", "PNLIP", "AMY2A", "PRSS1", "CLPS", "CTRC", "PNLIPRP1", "CELA2B", "CTRB2")

#### paraganglioma

#### prostate
#https://www.proteinatlas.org/humanproteome/tissue/prostate
#TGM4, KLK3, KLK2, ACP3, RLN1, KLK4, MSMB, SLC45A3, SP8, CHRNA2, STEAP2, NKX3-1
prostate_genes_reference <- c("TGM4", "KLK3", "KLK2", "ACP3", "RLN1", "KLK4", "MSMB", "SLC45A3", "SP8", "CHRNA2", "STEAP2", "NKX3-1")

#### sarcoma

#### skin
#https://www.proteinatlas.org/humanproteome/tissue/skin
#KRTAP4-7, KRTAP9-3, KRTAP4-9, KRTAP2-1, KRTAP9-4, KRTAP9-8, KRTAP2-2, KRTAP4-6, KRTAP1-1, 
# KRTAP3-1, KRTAP1-3, KRTAP4-12
skin_genes_reference <- c("KRTAP4-7"," KRTAP9-3", "KRTAP4-9", "KRTAP2-1", "KRTAP9-4", "KRTAP9-8", "KRTAP2-2", "KRTAP4-6", "KRTAP1-1", 
                           "KRTAP3-1", "KRTAP1-3", "KRTAP4-12")

#### stomach
#https://www.proteinatlas.org/humanproteome/tissue/stomach
#PGA5, CBLIF, GKN1, PGA4, PGA3, GAST, ATP4B, LIPF, MUC5AC, ATP4A, GKN2, TFF1
stomach_genes_reference <- c("PGA5", "CBLIF", "GKN1", "PGA4", "PGA3", "GAST", "ATP4B", "LIPF", "MUC5AC", "ATP4A", "GKN2", "TFF1")

#### testicles
#https://www.proteinatlas.org/humanproteome/tissue/testis
#HMGB4, CXorf51B, LELP1, TNP1, TPD52L3, SHCBP1L, ACTRT2, PRR30, PRM1, ACTL7B, CETN1, ACTL9
testicles_genes_reference <- c("HMGB4", "CXorf51B", "LELP1", "TNP1", "TPD52L3", "SHCBP1L", "ACTRT2", "PRR30", "PRM1", "ACTL7B", "CETN1", "ACTL9")

#### thymoma

#### thyroid
#https://www.proteinatlas.org/humanproteome/tissue/thyroid+gland
#TG, TPO, SLC26A4, BMP8A, IYD, SLC26A7, TSHR, FOXE1, SCUBE3, IGFBPL1, CRYGN, PDE8B
thyroid_genes_reference <- c("TG", "TPO", "SLC26A4", "BMP8A", "IYD", "SLC26A7", "TSHR", "FOXE1", "SCUBE3", "IGFBPL1", "CRYGN", "PDE8B")

#### uterine

#### uveal

###################################################
setwd(tablesDirectory)
gtf <- read.table('gencode.v23.chr_patch_hapl_scaff.basic.annotation.gtf', header = FALSE, sep = '\t')

gtf$gene_id_mod <- sapply(gtf$V9, function(x){
  sub(".*gene_id ([^.]+).*", "\\1", x)
})

gtf$gene_symbol <- sapply(gtf$V9, function(x){
  sub(".*gene_name ([^.]+).*", "\\1", x)
})

gtf$gene_symbol <- sapply(gtf$gene_symbol, function(x){
  str_split(x, ";")[[1]][1]
})

gtf.f <- gtf[gtf$V3 == "gene",]


adrenocortical_genes_reference <- adrenocortical_genes_reference[adrenocortical_genes_reference %in% gtf.f$gene_symbol] # One false CCN3
adrenocortical_genes_reference_ID <- gtf.f$gene_id_mod[match(adrenocortical_genes_reference, gtf.f$gene_symbol)]

B_cell_lymphoma_genes_reference <- B_cell_lymphoma_genes_reference[B_cell_lymphoma_genes_reference %in% gtf.f$gene_symbol]
B_cell_lymphoma_genes_reference_ID <- gtf.f$gene_id_mod[match(B_cell_lymphoma_genes_reference, gtf.f$gene_symbol)]

bladder_genes_reference <- bladder_genes_reference[bladder_genes_reference %in% gtf.f$gene_symbol]
bladder_genes_reference_ID <- gtf.f$gene_id_mod[match(bladder_genes_reference, gtf.f$gene_symbol)]

brain_genes_reference <- brain_genes_reference[brain_genes_reference %in% gtf.f$gene_symbol]
brain_genes_reference_ID <- gtf.f$gene_id_mod[match(brain_genes_reference, gtf.f$gene_symbol)]

breast_genes_reference <- breast_genes_reference[breast_genes_reference %in% gtf.f$gene_symbol]
breast_genes_reference_ID <- gtf.f$gene_id_mod[match(breast_genes_reference, gtf.f$gene_symbol)]

cervical_genes_reference <- cervical_genes_reference[cervical_genes_reference %in% gtf.f$gene_symbol]
cervical_genes_reference_ID <- gtf.f$gene_id_mod[match(cervical_genes_reference, gtf.f$gene_symbol)]

colon_genes_reference <- colon_genes_reference[colon_genes_reference %in% gtf.f$gene_symbol]
colon_genes_reference_ID <- gtf.f$gene_id_mod[match(colon_genes_reference, gtf.f$gene_symbol)]

esophageal_genes_reference <- esophageal_genes_reference[esophageal_genes_reference %in% gtf.f$gene_symbol]
esophageal_genes_reference_ID <- gtf.f$gene_id_mod[match(esophageal_genes_reference, gtf.f$gene_symbol)]

kidney_genes_reference <- kidney_genes_reference[kidney_genes_reference %in% gtf.f$gene_symbol]
kidney_genes_reference_ID <- gtf.f$gene_id_mod[match(kidney_genes_reference, gtf.f$gene_symbol)]

liver_genes_reference <- liver_genes_reference[liver_genes_reference %in% gtf.f$gene_symbol]
liver_genes_reference_ID <- gtf.f$gene_id_mod[match(liver_genes_reference, gtf.f$gene_symbol)]

lung_genes_reference <- lung_genes_reference[lung_genes_reference %in% gtf.f$gene_symbol]
lung_genes_reference_ID <- gtf.f$gene_id_mod[match(lung_genes_reference, gtf.f$gene_symbol)]

ovarian_genes_reference <- ovarian_genes_reference[ovarian_genes_reference %in% gtf.f$gene_symbol]
ovarian_genes_reference_ID <- gtf.f$gene_id_mod[match(ovarian_genes_reference, gtf.f$gene_symbol)]

pancreas_genes_reference <- pancreas_genes_reference[pancreas_genes_reference %in% gtf.f$gene_symbol]
pancreas_genes_reference_ID <- gtf.f$gene_id_mod[match(pancreas_genes_reference, gtf.f$gene_symbol)]

prostate_genes_reference <- prostate_genes_reference[prostate_genes_reference %in% gtf.f$gene_symbol] # One false ACP3
prostate_genes_reference_ID <- gtf.f$gene_id_mod[match(prostate_genes_reference, gtf.f$gene_symbol)]

skin_genes_reference <- skin_genes_reference[skin_genes_reference %in% gtf.f$gene_symbol] # One false KRTAP9-3
skin_genes_reference_ID <- gtf.f$gene_id_mod[match(skin_genes_reference, gtf.f$gene_symbol)]

stomach_genes_reference <- stomach_genes_reference[stomach_genes_reference %in% gtf.f$gene_symbol] # One false CBLIF
stomach_genes_reference_ID <- gtf.f$gene_id_mod[match(stomach_genes_reference, gtf.f$gene_symbol)]

testicles_genes_reference <- testicles_genes_reference[testicles_genes_reference %in% gtf.f$gene_symbol]
testicles_genes_reference_ID <- gtf.f$gene_id_mod[match(testicles_genes_reference, gtf.f$gene_symbol)]

thyroid_genes_reference <- thyroid_genes_reference[thyroid_genes_reference %in% gtf.f$gene_symbol]
thyroid_genes_reference_ID <- gtf.f$gene_id_mod[match(thyroid_genes_reference, gtf.f$gene_symbol)]

###################################################




list_genes_reference_ID <- list(adrenocortical_genes_reference_ID,
                                bladder_genes_reference_ID,
                                brain_genes_reference_ID,
                                breast_genes_reference_ID,
                                cervical_genes_reference_ID,
                                colon_genes_reference_ID,
                                esophageal_genes_reference_ID,
                                kidney_genes_reference_ID,
                                liver_genes_reference_ID,
                                lung_genes_reference_ID,
                                ovarian_genes_reference_ID,
                                pancreas_genes_reference_ID,
                                prostate_genes_reference_ID,
                                skin_genes_reference_ID,
                                stomach_genes_reference_ID,
                                testicles_genes_reference_ID,
                                thyroid_genes_reference_ID)
names(list_genes_reference_ID) <- c("adrenocortical","bladder","brain","breast","cervical",
                                    "colon","esophageal","kidney","liver","lung","ovarian","pancreas","prostate",
                                    "skin", "stomach","testicles","thyroid")


setwd(tablesDirectory)
padj_log2FC <- read.csv(file="padj_log2FC.csv")
padj_log2FC <- padj_log2FC[,colnames(padj_log2FC)!="X"]


setwd(tablesDirectory)
genes_list_combinations <- readRDS(file="genes_list_all_combinations_padj_log2FC.rds")


#### UNION
genes_included_union <- data.frame()

for (i in 1:nrow(padj_log2FC)){
  genes_union <- lapply(genes_list_combinations[[i]], function(x){
    Reduce(unique,x)
  })
  
  vector <- c()
  
  for (j in 1:length(names(genes_union))){
    if (names(genes_union)[j] %in% names(list_genes_reference_ID)){
      value <- sum(list_genes_reference_ID[[names(genes_union)[j]]] %in% genes_union[[j]]==TRUE)/length(list_genes_reference_ID[[names(genes_union)[j]]])*100
      vector <- c(vector, value)
    }
  }
  
  genes_included_union <- rbind(genes_included_union, vector)
}

colnames(genes_included_union) <- names(genes_union)[names(genes_union) %in% names(list_genes_reference_ID)]
genes_included_union <- cbind(padj_log2FC,genes_included_union)

setwd(resultsDirectory)
write.csv(genes_included_union, file="genes_included_union.csv")

#### INTERSECT
genes_included_intersect <- data.frame()

for (i in 1:nrow(padj_log2FC)){
  genes_intersect <- lapply(genes_list_combinations[[i]], function(x){
    Reduce(intersect,x)
  })
  
  vector <- c()
  
  for (j in 1:length(names(genes_intersect))){
    if (names(genes_intersect)[j] %in% names(list_genes_reference_ID)){
      value <- sum(list_genes_reference_ID[[names(genes_intersect)[j]]] %in% genes_intersect[[j]]==TRUE)/length(list_genes_reference_ID[[names(genes_intersect)[j]]])*100
      vector <- c(vector, value)
    }
  }
  
  genes_included_intersect <- rbind(genes_included_intersect, vector)
}


colnames(genes_included_intersect) <- names(genes_intersect)[names(genes_intersect) %in% names(list_genes_reference_ID)]
genes_included_intersect <- cbind(padj_log2FC,genes_included_intersect)

setwd(resultsDirectory)
write.csv(genes_included_intersect, file="genes_included_intersect.csv")


#INTERSECT NEW
genes_included_intersect_new <- data.frame()

for (i in 1:nrow(padj_log2FC)){
  genes_intersect_new <- lapply(genes_list_combinations[[i]], function(x){
    all <- unlist(x)
    table_all <- table(all)
    names(table_all[table_all>=24])
  })
  
  vector <- c()
  
  for (j in 1:length(names(genes_intersect_new))){
    if (names(genes_intersect)[j] %in% names(list_genes_reference_ID)){
      value <- sum(list_genes_reference_ID[[names(genes_intersect_new)[j]]] %in% genes_intersect_new[[j]]==TRUE)/length(list_genes_reference_ID[[names(genes_intersect_new)[j]]])*100
      vector <- c(vector, value)
    }
  }
  
  genes_included_intersect_new <- rbind(genes_included_intersect_new, vector)
}


colnames(genes_included_intersect_new) <- names(genes_intersect_new)[names(genes_intersect_new) %in% names(list_genes_reference_ID)]
genes_included_intersect_new <- cbind(padj_log2FC,genes_included_intersect_new)

setwd(resultsDirectory)
write.csv(genes_included_intersect_new, file="genes_included_intersect_24.csv")

######################################
### Select combination
######################################

max_number_genes_disease <- apply(genes_included_union[,3:ncol(genes_included_union)],2,max)

combinations <- data.frame()

for (i in 1:nrow(genes_included_union)){
  if ((sum(genes_included_union[i,3:ncol(genes_included_union)]==max_number_genes_disease))==(ncol(genes_included_union)-2)){
    combinations <- rbind(combinations, genes_included_union[i,])
  }
}

setwd(resultsDirectory)
write.csv(combinations, file="combinations_opt.csv")


######################################
### Save genes for this combination
######################################

genes_list_final <- lapply(genes_list, function(x){
  lapply(x, function(y){
    rownames(y[y$padj < 0.01 & y$log2FoldChange > 1,])
  })
})


setwd(tablesDirectory)
saveRDS(genes_list_final, file="genes_list_final.rds")
