#####################################################################
# fgsea after DESeq2
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
resultsDirectory <- "~/Desktop/LIDIA/TCGA_classification/TFM/Results/fgsea"
dataDirectory <- "~/Desktop/LIDIA/TCGA_classification/Data"

######################################
# C8
setwd(dataDirectory)

data(examplePathways)
data("exampleRanks")

c8.all.v2023.2.Hs.symbols <- readLines("c8.all.v2023.2.Hs.symbols.gmt")

gene_sets_C8 <- lapply(c8.all.v2023.2.Hs.symbols, function(linea) {
  elementos <- strsplit(linea, "\t")[[1]]  
  list(
    ID = elementos[1],
    WEB = elementos[2],
    GENES = elementos[-c(1,2)]  
  )
})

names(gene_sets_C8) <- sapply(gene_sets_C8, function(x){x$ID})
gene_sets_C8 <- lapply(gene_sets_C8, function(x){
  x <- x$GENES
  })


c5.all.v2023.2.Hs.symbols <- readLines("c5.all.v2023.2.Hs.symbols.gmt")
gene_sets_C5 <- lapply(c5.all.v2023.2.Hs.symbols, function(linea) {
  elementos <- strsplit(linea, "\t")[[1]]  
  list(
    ID = elementos[1],
    WEB = elementos[2],
    GENES = elementos[-c(1,2)]  
  )
})
names(gene_sets_C5) <- sapply(gene_sets_C5, function(x){x$ID})
gene_sets_C5 <- lapply(gene_sets_C5, function(x){
  x <- x$GENES
})
######################################
# Genes list
setwd(tablesDirectory)
genes_list <- readRDS(file="genes_list.rds")
######################################


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

######################################

# First comparison
# cervix vs uterine
genes_list_cervix_uterine <- genes_list$cervical$uterine
genes_list_cervix_uterine <- genes_list_cervix_uterine[order(genes_list_cervix_uterine$stat, decreasing = TRUE),]
genes_list_cervix_uterine$symbol <- gtf.f$gene_symbol[match(rownames(genes_list_cervix_uterine), gtf.f$gene_id_mod)]

library(tidyverse)
ranks <- genes_list_cervix_uterine[, c("symbol", "stat")]
ranks <- deframe(ranks)
head(ranks, 20)
# fgsea
fgseaRes <- fgsea(pathways=gene_sets_C5, stats=ranks)


# bladder vs uterine
genes_list_bladder_uterine <- genes_list$bladder$uterine
genes_list_bladder_uterine <- genes_list_bladder_uterine[order(genes_list_bladder_uterine$stat, decreasing = TRUE),]
genes_list_bladder_uterine$symbol <- gtf.f$gene_symbol[match(rownames(genes_list_bladder_uterine), gtf.f$gene_id_mod)]

library(tidyverse)
ranks <- genes_list_bladder_uterine[, c("symbol", "stat")]
ranks <- deframe(ranks)
length(unique(names(ranks)))
ranks <- ranks[unique(names(ranks))]
head(ranks, 20)
# fgsea
fgseaRes <- fgsea(pathways=gene_sets_C8, stats=ranks)



# adrenocortical vs liver
genes_list_adrenocortical_liver <- genes_list$adrenocortical$liver
genes_list_adrenocortical_liver <- genes_list_adrenocortical_liver[order(genes_list_adrenocortical_liver$stat, decreasing = TRUE),]
genes_list_adrenocortical_liver$symbol <- gtf.f$gene_symbol[match(rownames(genes_list_adrenocortical_liver), gtf.f$gene_id_mod)]

library(tidyverse)
ranks <- genes_list_adrenocortical_liver[, c("symbol", "stat")]
ranks <- deframe(ranks)
length(unique(names(ranks)))
ranks <- ranks[unique(names(ranks))]
head(ranks, 20)
# fgsea
fgseaRes <- fgsea(pathways=gene_sets_C8, stats=ranks)
