setwd("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/")
library(reshape2)
library(ggplot2)
library(plotly)
library(tidyverse)
library(dplyr)
library(plyr)
library(gridExtra)
library(edgeR)
library(DESeq2)
library(ggfortify)

# Merge data frames of same samples
myFilesRP1 <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/RP", pattern = "*.32_33.count", full.names=TRUE)
myFilesRP2 <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/RP", pattern = "*.26_31.count", full.names=TRUE)

for (I in myFilesRP1){
  S <- (strsplit(I, split='\\.')[[1]][1])
  for (L in myFilesRP2){
    D <- (strsplit(L, split='\\.')[[1]][1])
    if (S == D){
      count_data1 <- read.table(I)
      count_data2 <- read.table(L)
      count_data_merged <- merge(count_data1, count_data2, by = 'V1', all=TRUE)
      count_data_merged[is.na(count_data_merged)] <- 0
      count_data_merged$counts <- count_data_merged$V2.x + count_data_merged$V2.y
      #count_data_merged <- count_data_merged[count_data_merged$counts >= 1, ]
      count_data_merged$V2.x = NULL
      count_data_merged$V2.y = NULL
      Z <- (strsplit(S, split='\\//')[[1]][2])
      colnames(count_data_merged) <- c('Gene_ID', Z)
      csvFile <- paste(D, '_merged', '.csv', sep='')
      write.csv(count_data_merged, csvFile, row.names = FALSE)
    }
  }
}

myFilesRPa <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/RP", pattern = "*merged.csv", full.names=TRUE)
import.list <- llply(myFilesRPa, read.csv)
RP_merged <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Gene_ID"), import.list, accumulate=F)
RP_merged[is.na(RP_merged)] <- 0
RP_merged <- RP_merged[rowSums(RP_merged[,c(2:16)]) >= 20, ]
write.csv(RP_merged, 'RP_counts.csv', row.names = FALSE)

#Merge new and old sequencing
RP_counts_new <- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/RP_counts.csv")
RP_counts_old <- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/RP_counts.csv")
RP_count_merged <- merge(x = RP_counts_new, y = RP_counts_old, by = 'Gene_ID', all = TRUE)
RP_counts_merged <- matrix(nrow = NROW(RP_count_merged), ncol = 16)
RP_counts_merged[,1] <- RP_count_merged$Gene_ID
for (i in 2:16){
RP_counts_merged[,i] <- RP_count_merged[,i] + RP_count_merged[,i+15]
}
RP_counts_merged[is.na(RP_counts_merged)] <- 0
colnames(RP_counts_merged) <- c("Gene_ID","RP113","RP114","RP115","RP116","RP117","RP118","RP126","RP127","RP128","RP129","RP139","RP140","RP141","RP142","RP143")
write.csv(RP_counts_merged, 'RP_counts_merged.csv', row.names = FALSE)

#TR
myFilesTR <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/TR/", full.names=TRUE)
myFilesTR
for (i in myFilesTR){
  S <- (strsplit(i, split='\\//')[[1]][2])
  D <- (strsplit(S, split='\\.')[[1]][1])
  count_data <- read.table(i)
  #count_data <- count_data[count_data[,2] >= 10, ]
  colnames(count_data) <- c('Gene_ID', D)
  csvFile <- paste(S, '.csv', sep='')
  write.csv(count_data, csvFile, row.names = FALSE)
}

myFilesTRa <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/", pattern = "*.count.csv", full.names=TRUE)
import.list <- llply(myFilesTRa, read.csv)
TR_merged <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Gene_ID"), import.list, accumulate=F)
TR_merged[is.na(TR_merged)] <- 0
TR_merged <- TR_merged[rowSums(TR_merged[,c(2:16)]) >= 20, ] #Remove gene for which there are less than 10 reads for all samples together
write.csv(TR_merged, 'TR_counts.csv', row.names = FALSE)

#Merge new and old sequencing
TR_counts_new <- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/TR_counts.csv")
TR_counts_old <- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/TR_count/TR_counts.csv")
TR_count_merged <- merge(x = TR_counts_new, y = TR_counts_old, by = 'Gene_ID', all = TRUE)
TR_counts_merged <- matrix(nrow = NROW(TR_count_merged), ncol = 16)
TR_counts_merged[,1] <- TR_count_merged$Gene_ID
for (i in 2:16){
  TR_counts_merged[,i] <- TR_count_merged[,i] + TR_count_merged[,i+15]
}
TR_counts_merged[is.na(TR_counts_merged)] <- 0
colnames(TR_counts_merged) <- c("Gene_ID","TR113","TR114","TR115","TR116","TR117","TR118","TR126","TR127","TR128","TR129","TR139","TR140","TR141","TR142","TR143")
write.csv(TR_counts_merged, 'TR_counts_merged.csv', row.names = FALSE)

###Normalization read counts for TR samples

TR_merged <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/TR_counts_merged.csv')
design <- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/TR_design.csv", header = TRUE)
TR_merged$Gene_ID <- NULL
TR_raw_counts <- as.matrix(TR_merged)

#EdgeR TMM
dge2 <- DGEList(TR_raw_counts)
dge2
dge2 <- calcNormFactors(dge2, method = "TMM")
TR_pseudo_TMM <- log2(cpm(dge2) + 1)
TR_TMM_norm <- data.frame(TR_pseudo_TMM)
TR_counts_merged <- read.csv('TR_counts_merged.csv')
TR_TMM_norm$Gene_ID <- TR_counts_merged[,1]
TR_TMM_norm <- TR_TMM_norm[, c(16, 1:15)]
write.csv(TR_TMM_norm, "TR_TMM_norm.csv", row.names = FALSE)

#PCA analysis
TR_pseudo_TMM <- read.csv('TR_TMM_norm.csv')
TR_pseudo_TMM$Gene_ID <- NULL
TR_pseudo_TMM <- as.matrix(TR_pseudo_TMM)
pca <- prcomp(t(TR_pseudo_TMM ))
summary(pca)
pdf('TR_TMM_PCA_PC1_PC2.pdf')
autoplot(prcomp(t(TR_pseudo_TMM)), data = design, colour = 'groups', size = 2, shape = 'batches') + theme_bw()
dev.off()
pdf('TR_TMM_PCA_PC1_PC3.pdf')
autoplot(prcomp(t(TR_pseudo_TMM)), data = design, colour = 'groups', size = 2, x=2, y=3, shape = 'batches') + theme_bw()
dev.off()

###Normalization read counts for RP samples

RP_merged <- read.csv('RP_counts_merged.csv')
design <- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/RP_design.csv", header = TRUE)
RP_merged$Gene_ID <- NULL
RP_raw_counts <- as.matrix(RP_merged)

#EdgeR TMM
dge2 <- DGEList(RP_raw_counts)
dge2
dge2 <- calcNormFactors(dge2, method = "TMM")
RP_pseudo_TMM <- log2(cpm(dge2) + 1)
write.csv(RP_pseudo_TMM, "RP_TMM_norm.csv", row.names = FALSE)
RP_pseudo_TMM <- read.csv('RP_TMM_norm.csv')
RP_counts_merged <- read.csv('RP_counts_merged.csv')
RP_pseudo_TMM$Gene_ID <- RP_counts_merged[,1]
RP_pseudo_TMM <- RP_pseudo_TMM[, c(16, 1:15)]
write.csv(RP_pseudo_TMM, "RP_TMM_norm.csv", row.names = FALSE)

#PCA
RP_pseudo_TMM <- read.csv('RP_TMM_norm.csv')
RP_pseudo_TMM$Gene_ID <- NULL
pdf('RP_TMM_PCA_PC1_PC2.pdf')
autoplot(prcomp(t(RP_pseudo_TMM)), data = design, colour = 'groups', size = 2, shape = 'batches') + theme_bw()
dev.off()
pdf('RP_TMM_PCA_PC1_PC3.pdf')
autoplot(prcomp(t(RP_pseudo_TMM)), data = design, colour = 'groups', size = 2, x=1, y=3, shape = 'batches') + theme_bw()
dev.off()
pdf('RP_TMM_PCA_PC2_PC3.pdf')
autoplot(prcomp(t(RP_pseudo_TMM)), data = design, colour = 'groups', size = 2, x=2, y=3, shape = 'batches') + theme_bw()
dev.off()
pdf('RP_TMM_PCA_PC1_PC4.pdf')
autoplot(prcomp(t(RP_pseudo_TMM)), data = design, colour = 'groups', size = 2, x=1, y=4, shape = 'batches') + theme_bw()
dev.off()

### TPM calculation for RP

RP_TMM_norm <- read.csv("RP_TMM_norm.csv")
gene_length <- read.table("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/annotations/Mus_musculus.GRCm38.100.geneCDSlength")
colnames(gene_length) <- c("Gene_ID", "gene_length")
gene_length <- gene_length[gene_length$Gene_ID %in% RP_TMM_norm$Gene_ID, ] #To select rows that are in both tables
gene_length["gene_length"]=gene_length["gene_length"]/1000 #Divide length genes by 1000 to have length in kb
RP_TMM_norm_c <- RP_TMM_norm[RP_TMM_norm$Gene_ID %in% gene_length$Gene_ID, ]

RP_TPM <- matrix(ncol=ncol(RP_TMM_norm_c), nrow=nrow(RP_TMM_norm_c))
RP_TPM <- data.frame(RP_TPM)
for(i in 2:ncol(RP_TMM_norm_c)) {
  col <- RP_TMM_norm_c[,i]
  col_RPK = col/gene_length["gene_length"]
  per_million_scaling_factor = (sum(col_RPK))/1000000
  col_TPM = col/per_million_scaling_factor
  RP_TPM[,i] <- col_TPM
}
RP_TPM[,1] <- gene_length[,1]
write.csv(RP_TPM, "RP_TPM.csv", row.names = FALSE)

### TPM calculation for TR

TR_TMM_norm <- read.csv("TR_TMM_norm.csv")
gene_length <- read.table("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/annotations/Mus_musculus.GRCm38.100.geneCDSlength")
colnames(gene_length) <- c("Gene_ID", "gene_length")
gene_length <- gene_length[gene_length$Gene_ID %in% TR_TMM_norm$Gene_ID, ] #To select rows that are in both tables
gene_length["gene_length"]=gene_length["gene_length"]/1000 #Divide length genes by 1000 to have length in kb
TR_TMM_norm_c <- TR_TMM_norm[TR_TMM_norm$Gene_ID %in% gene_length$Gene_ID, ]

TR_TPM <- matrix(ncol=ncol(TR_TMM_norm_c), nrow=nrow(TR_TMM_norm_c))
TR_TPM <- data.frame(TR_TPM)
for(i in 2:ncol(TR_TMM_norm_c)) {
  col <- TR_TMM_norm_c[,i]
  col_RPK = col/gene_length["gene_length"]
  per_million_scaling_factor = (sum(col_RPK))/1000000
  col_TPM = col/per_million_scaling_factor
  TR_TPM[,i] <- col_TPM
}
TR_TPM[,1] <- gene_length[,1]
colnames(TR_TPM) <- colnames(TR_TMM_norm)
write.csv(TR_TPM, "TR_TPM.csv", row.names = FALSE)

### Translation Efficiency
setwd("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/")

#Merge TR and RP data
RP_TMM_norm <- read.csv('RP_TPM.csv')
TR_TMM_norm <- read.csv('TR_TPM.csv')
RP_TR_merged <- merge(RP_TMM_norm, TR_TMM_norm, all=TRUE, by="Gene_ID")
RP_TR_merged[is.na(RP_TR_merged)] <- 0
write.csv(RP_TR_merged, "RP_TR_merged.csv", row.names = FALSE)
RP_TR_merged <- read.csv('RP_TR_merged.csv')

#Calculate TE (RP/TR)
RP_TR_ratio <- matrix(ncol = 16, nrow = NROW(RP_TR_merged))
for (i in 2:16){
  RP_TR_ratio[,i] <- (RP_TR_merged[,i]/RP_TR_merged[,i+15])
}
RP_TR_ratio <- apply(RP_TR_ratio, 2, as.numeric)
RP_TR_ratio <- data.frame(RP_TR_ratio)
colnames(RP_TR_ratio) <- colnames(RP_TR_merged[,1:16])
RP_TR_ratio$Gene_ID <- RP_TR_merged[,1]
write.csv(RP_TR_ratio, 'RP_TR_ratio.csv', row.names = FALSE)

### relative TE (knock-down/control RPF)
ctrl_kd_ratio <- matrix(ncol = 16, nrow = NROW(RP_TR_ratio))
ctrl_kd_ratio[,2] <- RP_TR_ratio[,2]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,3] <- RP_TR_ratio[,5]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,4] <- RP_TR_ratio[,8]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,5] <- RP_TR_ratio[,4]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,6] <- RP_TR_ratio[,7]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,7] <- RP_TR_ratio[,10]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,8] <- RP_TR_ratio[,3]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,9] <- RP_TR_ratio[,6]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,10] <- RP_TR_ratio[,9]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,11] <- RP_TR_ratio[,11]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,12] <- RP_TR_ratio[,12]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,13] <- RP_TR_ratio[,13]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,14] <- RP_TR_ratio[,14]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,15] <- RP_TR_ratio[,15]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio[,16] <- RP_TR_ratio[,16]/rowMeans(RP_TR_ratio[,c(2,4,5,7,8,10)])
ctrl_kd_ratio <- apply(ctrl_kd_ratio, 2, as.numeric)
ctrl_kd_ratio <- data.frame(ctrl_kd_ratio)
colnames(ctrl_kd_ratio) <- c('Gene_ID','Scr1','Scr2','Scr3','Eif2d_nf1','Eif2d_nf2','Eif2d_nf3','Denr_shRNA2_1','Denr_shRNA2_2','Denr_shRNA2_3','Eif2d_sh3_1','Eif2d_sh3_2','Eif2d_sh3_3','Eif2d_sh4_1','Eif2d_sh4_2','Eif2d_sh4_3')
ctrl_kd_ratio$Gene_ID <- RP_TR_merged[,1]
write.csv(ctrl_kd_ratio, 'relative_TE_ctrl_kd_TPM_TMM.csv', row.names = FALSE)
