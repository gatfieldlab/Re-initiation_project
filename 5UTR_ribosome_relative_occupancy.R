setwd('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/')
library(plyr)
library(edgeR)
library(plotly)
library(ggplot2)
library(ggfortify)
library(tidyr)
library(dplyr)
library(gridExtra)

myFilesRP1 <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/", pattern = "*.32_33.5UTR.count", full.names=TRUE)
myFilesRP2 <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/", pattern = "*.26_31.5UTR.count", full.names=TRUE)

for (I in myFilesRP1){
  S <- (strsplit(I, split='\\.')[[1]][1])
  for (L in myFilesRP2){
    D <- (strsplit(L, split='\\.')[[1]][1])
    if (S == D){
      count_data1 <- read.table(I)
      count_data2 <- read.table(L)
      count_data_merged <- merge(count_data1, count_data2,by = 'V1', all=TRUE)
      count_data_merged[is.na(count_data_merged)] <- 0
      count_data_merged$counts <- count_data_merged$V2.x + count_data_merged$V2.y
      #count_data_merged <- count_data_merged[count_data_merged$counts >= 1, ]
      count_data_merged$V2.x = NULL
      count_data_merged$V2.y = NULL
      Z <- (strsplit(S, split='\\//')[[1]][2])
      colnames(count_data_merged) <- c('Gene_ID', Z)
      csvFile <- paste(D, '5UTR_merged', '.csv', sep='')
      write.csv(count_data_merged, csvFile, row.names = FALSE)
    }
  }
}

#Merge RP117 files
RP117_1 <- read_csv('RP117_346_2_0015UTR_merged.csv')
RP117_2 <- read_csv('RP117_346_3_0015UTR_merged.csv')
RP117 <- merge(RP117_1, RP117_2, by = 'Gene_ID', all = TRUE)
RP117$total <- RP117[,2] + RP117[,3]
RP117$X2.x <- NULL
RP117$X2.y <- NULL
write_csv(RP117, 'RP117_346_0015UTR_merged.csv')

#Merge old and new sequencing
myFilesRPa <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/", pattern = "*5UTR_merged.csv", full.names=TRUE)
myFilesRPb <- list.files("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/5UTR_count/RP", pattern = "*5UTR_merged.csv", full.names=TRUE)
import.lista <- llply(myFilesRPa, read.csv)
import.listb <- llply(myFilesRPb, read.csv)
RP_5UTR_a <- Reduce(function(x, y) merge(x, y, all=TRUE, by='Gene_ID'), import.lista, accumulate=F)
RP_5UTR_b <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Gene_ID"), import.listb, accumulate=F)
RP_5UTR_merged <- merge(RP_5UTR_a, RP_5UTR_b, by = 'Gene_ID' , all=TRUE)
RP_5UTR <- matrix(nrow=NROW(RP_5UTR_merged), ncol = 16)
RP_5UTR[,1] <- RP_5UTR_merged[,1]
for (i in 2:16){
  RP_5UTR[,i] <- RP_5UTR_merged[,i] + RP_5UTR_merged[,i+15]
}
colnames(RP_5UTR) <- c('Gene_ID','Scr_1','Denrsh2_1','Eif2dsh1_1','Scr_2','Denrsh2_2','Eif2dsh1_2','Scr_3','Denrsh2_3','Eif2dsh1_3','Eif2dsh3_1','Eif2dsh3_2','Eif2dsh3_3','Eif2dsh4_1','Eif2dsh4_2','Eif2dsh4_3')
RP_5UTR <- data.frame(RP_5UTR)
write.csv(RP_5UTR, 'RP_5UTR_counts.csv', row.names = FALSE)

RP_5UTR_merged<- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/RP_5UTR_counts.csv", header = TRUE)
RP_5UTR_merged <- RP_5UTR_merged[rowSums(RP_5UTR_merged >= 50) >= 10,]
design <- read.csv("~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/RP_design.csv", header = TRUE)
RP_5UTR_merged$Gene_ID <- NULL
colnames(RP_5UTR_merged) <- design[,1]
RP_5UTR_merged[is.na(RP_5UTR_merged)] <- 0
RP_raw_counts <- as.matrix(RP_5UTR_merged)

for (i in 1:15){
RP_raw_counts[,i] <- as.numeric(as.character(RP_raw_counts[,i])) 
}

#Calculate ratio 5UTR count/CDS count (raw counts)

RP_CDS_raw_count <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/RP_counts_merged.csv')
RP_5UTR_raw_count <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/RP_5UTR_counts.csv')

CDS_5UTR_raw_merged <- merge(RP_CDS_raw_count, RP_5UTR_raw_count, all=TRUE, by="Gene_ID")
CDS_5UTR_raw_merged <- CDS_5UTR_raw_merged %>% drop_na()
write.csv(CDS_5UTR_raw_merged, "CDS_5UTR_raw_merged.csv", row.names = FALSE)
CDS_5UTR_raw_ratio <- matrix(ncol = 16, nrow = nrow(CDS_5UTR_raw_merged))
CDS_5UTR_raw_ratio[,1] <- CDS_5UTR_raw_merged[,1]
for (i in 2:16){
  CDS_5UTR_raw_ratio[,i] <- CDS_5UTR_raw_merged[,i+15]/CDS_5UTR_raw_merged[,i]
}
CDS_5UTR_raw_ratio <- data.frame(CDS_5UTR_raw_ratio)
colnames(CDS_5UTR_raw_ratio) <- colnames(RP_5UTR_raw_count[,1:16])
write.csv(CDS_5UTR_raw_ratio, '5UTR_CDS_raw_ratio.csv', row.names = FALSE)

CDS_5UTR_raw_ratio <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/5UTR_CDS_raw_ratio.csv')
Gene_name <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/Xtail/Gene_ID_name.csv') 
colnames(Gene_name) <- c('Gene_ID', 'Gene_name')
CDS_5UTR_raw_ratio <- merge(CDS_5UTR_raw_ratio, Gene_name, by = 'Gene_ID', all = FALSE)
CDS_5UTR_raw_ratio <- CDS_5UTR_raw_ratio[,c(1,17,2:16)]
CDS_5UTR_raw_ratio <- CDS_5UTR_raw_ratio %>% filter_all(all_vars(!is.infinite(.)))
CDS_5UTR_raw_ratio$mean_Control <- rowMeans(CDS_5UTR_raw_ratio[,c(3,5,6,8,9,11)])
CDS_5UTR_raw_ratio$mean_Denrsh2 <- rowMeans(CDS_5UTR_raw_ratio[,c(4,7,10)])
CDS_5UTR_raw_ratio$mean_Eif2dsh <- rowMeans(CDS_5UTR_raw_ratio[,c(12:17)])
CDS_5UTR_raw_ratio_means <- CDS_5UTR_raw_ratio[,c(1,2,18:22)]
write.csv(CDS_5UTR_raw_ratio_means, '5UTR_CDS_raw_ratio_mean.csv', row.names = FALSE)

#Calculate log2 ((5UTR/CDS)knock-down/(5UTR/CDS)ctrl)
CDS_5UTR_ratio <- read.csv('5UTR_CDS_raw_ratio_mean.csv')
CDS_5UTR_ratio$ratio_Denrsh2_Ctrl <- log2(CDS_5UTR_ratio[,4]/CDS_5UTR_ratio[,3])
CDS_5UTR_ratio$ratio_Eif2dsh_Ctrl <- log2(CDS_5UTR_ratio[,5]/CDS_5UTR_ratio[,3])
Gene_ID <- data.frame(CDS_5UTR_ratio$Gene_ID)
ratio_Denrsh2_Ctrl <- data.frame(CDS_5UTR_ratio$ratio_Denrsh2_Ctrl)
ratio_Eif2dsh_Ctrl <- data.frame(CDS_5UTR_ratio$ratio_Eif2dsh_Ctrl)
relative_5UTR_translation <- data.frame(Gene_ID,ratio_Denrsh2_Ctrl,ratio_Eif2dsh_Ctrl,ratio_Eif2dsh3_Ctrl,ratio_Eif2dsh4_Ctrl)
colnames(relative_5UTR_translation) <- c('Gene_ID','ratio_Denrsh2_Ctrl','ratio_Eif2dsh_Ctrl', 'ratio_Eif2dsh3_Ctrl','ratio_Eif2dsh4_Ctrl')
relative_5UTR_translation <- merge(relative_5UTR_translation, Gene_name, by = 'Gene_ID', all = FALSE)
relative_5UTR_translation <- relative_5UTR_translation[,c(1,6,2:5)]
write.csv(relative_5UTR_translation,'5UTR_CDS_raw_log2_ratio_vs_ctrl.csv' , row.names = FALSE)

#plot
pdf('relative_5UTR_translation_plots.pdf',  width=10, height=20)
Denrsh2_plot <- ggplot(ratio_Denrsh2_Ctrl, aes(x = relative_5UTR_translation$ratio_Denrsh2_Ctrl)) + geom_histogram(breaks=seq(-2, 2.5, by = 0.1), fill='#00AF6696') + geom_vline(xintercept = 0, linetype="dotted", lwd=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("")  + geom_vline(xintercept = 0.398957688410847, linetype="dotted", color = 'red3', lwd=1) + scale_y_continuous(limits = c(0,1000))
Eif2dsh_plot <- ggplot(ratio_Eif2dsh_Ctrl, aes(x = relative_5UTR_translation$ratio_Eif2dsh_Ctrl)) + geom_histogram(breaks=seq(-2, 2.5, by = 0.1), fill='#00B2EE96') + geom_vline(xintercept = 0, linetype="dotted", lwd=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("") + geom_vline(xintercept = 0.112624108690127, linetype="dotted", color = 'red3', lwd=1) + scale_y_continuous(limits = c(0,1000))
x <- grid.arrange(Denrsh2_plot, Eif2dsh_plot, nrow = 2)
dev.off()

#Number of genes with increased number of reads on 5_UTR
ratio_values <- matrix(nrow = 2, ncol = 4)
colnames(ratio_values) <- c('< 0', '> 0', '= 0', 'median')
rownames(ratio_values) <- c('Denrsh2/Ctrl','Eif2dsh/Ctrl')
ratio_values[1,1] <- NROW(relative_5UTR_translation[relative_5UTR_translation$ratio_Denrsh2_Ctrl < 0,])
ratio_values[1,2] <- NROW(relative_5UTR_translation[relative_5UTR_translation$ratio_Denrsh2_Ctrl > 0,])
ratio_values[1,3] <- NROW(relative_5UTR_translation[relative_5UTR_translation$ratio_Denrsh2_Ctrl == 0,])
m <- relative_5UTR_translation$ratio_Denrsh2_Ctrl 
ratio_values[1,4] <- median(m, na.rm = TRUE)
ratio_values[2,1] <- NROW(relative_5UTR_translation[relative_5UTR_translation$ratio_Eif2dsh_Ctrl  < 0,])
ratio_values[2,2] <- NROW(relative_5UTR_translation[relative_5UTR_translation$ratio_Eif2dsh_Ctrl  > 0,])
n <- relative_5UTR_translation$ratio_Eif2dsh_Ctrl
ratio_values[2,3] <- NROW(relative_5UTR_translation[relative_5UTR_translation$ratio_Eif2dsh_Ctrl  == 0,])
ratio_values[2,4] <- median(n, na.rm = TRUE)
write.csv(ratio_values, '5UTR_relative_translation_median.csv')

#Plot Denr sh against Eif2d sh
UTR_coverage <- read.csv('5UTR_CDS_raw_log2_ratio_vs_ctrl.csv')
UTR_coverage[sapply(UTR_coverage, is.infinite)] <- NA
UTR_coverage <- UTR_coverage %>% drop_na()
corr <- cor(UTR_coverage$ratio_Denrsh2_Ctrl, UTR_coverage$ratio_Eif2dsh_Ctrl, method = 'spearman')
pdf('5UTR_coverage_Denrsh2_vs_Eif2dsh.pdf')
plot(UTR_coverage$ratio_Eif2dsh_Ctrl, UTR_coverage$ratio_Denrsh2_Ctrl, type = 'p', xlab = 'Eif2d shRNAs', ylab = 'Denr shRNA 2', col='gray35', main = corr, xlim = c(-2,3))
model <- lm(UTR_coverage$ratio_Denrsh2_Ctrl ~ UTR_coverage$ratio_Eif2dsh_Ctrl, data = UTR_coverage)
abline(model)
dev.off()
summary(model)


#Wilcoxon test
log2_ratio <- read.csv('5UTR_CDS_raw_ratio_mean.csv')
wilcox.test(log2_ratio$mean_Control, log2_ratio$mean_Denrsh2, paired=TRUE)
wilcox.test(log2_ratio$mean_Control, log2_ratio$mean_Eif2dsh, paired=TRUE)


# 5' UTR TE all transcripts with 5' UTR coverage (TMM normalization)

RP_5UTR_count <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/count_data/new_sequencing/5_UTR/RP_5UTR_counts.csv')
RP_5UTR_count$ribo_ctrl_counts <- rowSums(RP_5UTR_count[,c(2,4,5,7,8,10)])
RP_5UTR_count$ribo_Denrsh2_counts <- rowSums(RP_5UTR_count[,c(3,6,9)])
RP_5UTR_count_only <- RP_5UTR_count[,c(17,18)]
RP_5UTR_count_matrix <- as.matrix(RP_5UTR_count_only)

dge2 <- DGEList(RP_5UTR_count_matrix)
dge2
dge2 <- calcNormFactors(dge2, method = "TMM")
Ctrl_Denrsh2_RP_5UTR_counts_pseudo_TMM <- log2(cpm(dge2) + 1)
Ctrl_Denrsh2_RP_5UTR_counts_TMM_norm <- data.frame(Ctrl_Denrsh2_RP_5UTR_counts_pseudo_TMM)
Ctrl_Denrsh2_RP_5UTR_counts_TMM_norm$Gene_ID <- RP_5UTR_count$Gene_ID
write.csv(Ctrl_Denrsh2_RP_5UTR_counts_TMM_norm, "Ctrl_Denrsh2_RP_5UTR_counts_TMM_norm.csv", row.names = FALSE)

Ctrl_Denrsh2_RNA_counts_TMM_norm <- read.csv('/home/rmeurs/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/uORF_annotations/Ctrl_Denrsh2_RNA_counts_TMM_norm.csv')

UTR_RNA_counts_TMM <- merge(Ctrl_Denrsh2_RP_5UTR_counts_TMM_norm, Ctrl_Denrsh2_RNA_counts_TMM_norm, by = 'Gene_ID')
UTR_RNA_counts_TMM$TE_ctrl <- UTR_RNA_counts_TMM$ribo_ctrl_counts/UTR_RNA_counts_TMM$RNA_ctrl_counts
UTR_RNA_counts_TMM$TE_Denrsh2 <- UTR_RNA_counts_TMM$ribo_Denrsh2_counts/UTR_RNA_counts_TMM$RNA_Denrsh2_counts
UTR_RNA_counts_TMM$log2_TE_FC <- log2(UTR_RNA_counts_TMM$TE_Denrsh2/UTR_RNA_counts_TMM$TE_ctrl)
UTR_RNA_counts_TMM <- na.omit(UTR_RNA_counts_TMM)
Gene_name <- read.csv('/home/rmeurs/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Gene_ID_name.csv')
UTR_RNA_counts_TMM <- merge(UTR_RNA_counts_TMM, Gene_name, by = 'Gene_ID')
write.csv(UTR_RNA_counts_TMM, '5UTR_RP_RNA_TMM_counts_TE_Denrsh2.csv', row.names = FALSE)

median(UTR_RNA_counts_TMM$TE_Denrsh2)
median(UTR_RNA_counts_TMM$TE_ctrl)
median(UTR_RNA_counts_TMM$log2_TE_FC)
wilcox.test(UTR_RNA_counts_TMM$TE_ctrl, UTR_RNA_counts_TMM$TE_Denrsh2, paired=FALSE) 

pdf('log2_TE_FC_5UTR_TMM_ctrl_Denrsh2_histogram.pdf')
ggplot(UTR_RNA_counts_TMM, aes(x = log2_TE_FC)) + 
  geom_histogram(breaks=seq(-5, 5, by = 0.1), fill='#00AF6696') + 
  geom_vline(xintercept = 0, linetype="dotted", lwd=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.back  = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab("log2 TE FC") + geom_vline(xintercept = -0.0101016, linetype="dotted", color = 'red3', lwd=1) + 
  scale_y_continuous(limits = c(0,3000)) 
dev.off() 

# 5' UTR TE transcripts with translated uORF TMM

translated_uORFs <- read.table('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/uORF_annotations/RP_ctrl.translated_uORF.longest.3.Gene_ID.txt')
colnames(translated_uORFs) <- 'Gene_ID'
translated_uORFs <- unique(translated_uORFs)
UTR_RNA_counts_translated_uORFs <- merge(UTR_RNA_counts_TMM,translated_uORFs, by = 'Gene_ID')

median(UTR_RNA_counts_translated_uORFs$TE_Denrsh2)
median(UTR_RNA_counts_translated_uORFs$TE_ctrl)
median(UTR_RNA_counts_translated_uORFs$log2_TE_FC)
wilcox.test(UTR_RNA_counts_translated_uORFs$TE_ctrl, UTR_RNA_counts_translated_uORFs$TE_Denrsh2, paired=FALSE) 

pdf('log2_TE_FC_5UTR_TMM_translated_uORF_transcripts_ctrl_Denrsh2_histogram_2.pdf')
ggplot(UTR_RNA_counts_translated_uORFs, aes(x = log2_TE_FC)) + 
  geom_histogram(breaks=seq(-2.5, 2.5, by = 0.05), fill='#c883c8d1') + 
  geom_vline(xintercept = 0, linetype="dotted", lwd=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.back  = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab("log2 TE FC") + geom_vline(xintercept = -0.01255897, linetype="dotted", color = '#c883c8d1', lwd=1) + 
  scale_y_continuous(limits = c(0,1000)) 
dev.off()  

# 5' UTR TE DENR-responsive transcripts with translated uORF TMM
DENR_responsive <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv')
UTR_RNA_counts_DENR_responsive <- merge(DENR_responsive,UTR_RNA_counts_TMM, by = 'Gene_ID')
median(UTR_RNA_counts_DENR_responsive$TE_Denrsh2)
median(UTR_RNA_counts_DENR_responsive$TE_ctrl)
median(UTR_RNA_counts_DENR_responsive$log2_TE_FC)
wilcox.test(UTR_RNA_counts_DENR_responsive$TE_ctrl, UTR_RNA_counts_DENR_responsive$TE_Denrsh2, paired=FALSE) 

pdf('log2_TE_FC_5UTR_TMM_DENR_responsive_transcripts_ctrl_Denrsh2_histogram_2.pdf')
ggplot(UTR_RNA_counts_DENR_responsive, aes(x = log2_TE_FC)) + 
  geom_histogram(breaks=seq(-2.5, 2.5, by = 0.05), fill='#00AF6696') + 
  geom_vline(xintercept = 0, linetype="dotted", lwd=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.back  = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab("log2 TE FC") + geom_vline(xintercept = -0.06102209, linetype="dotted", color = 'red3', lwd=1) + 
  scale_y_continuous(limits = c(0,1000)) 
dev.off()  

DENR_responsive <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv')
UTR_RNA_counts_DENR_responsive_uORFs <- merge(DENR_responsive,UTR_RNA_counts_translated_uORFs, by = 'Gene_ID')

median(UTR_RNA_counts_DENR_responsive_uORFs$TE_Denrsh2)
median(UTR_RNA_counts_DENR_responsive_uORFs$TE_ctrl)
median(UTR_RNA_counts_DENR_responsive_uORFs$log2_TE_FC)
wilcox.test(UTR_RNA_counts_DENR_responsive_uORFs$TE_ctrl, UTR_RNA_counts_DENR_responsive_uORFs$TE_Denrsh2, paired=FALSE) 

pdf('log2_TE_FC_5UTR_TMM_DENR_responsive_transcripts_translated_uORFs_ctrl_Denrsh2_histogram_2.pdf')
ggplot(UTR_RNA_counts_DENR_responsive_uORFs, aes(x = log2_TE_FC)) + 
  geom_histogram(breaks=seq(-2.5, 2.5, by = 0.05), fill='#00AF6696') + 
  geom_vline(xintercept = 0, linetype="dotted", lwd=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.back  = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab("log2 TE FC") + geom_vline(xintercept = -0.05533796, linetype="dotted", color = 'red3', lwd=1) + 
  scale_y_continuous(limits = c(0,1000)) 
dev.off()  

wilcox.test(UTR_RNA_counts_DENR_responsive_uORFs$log2FoldChange, UTR_RNA_counts_translated_uORFs$log2_TE_FC) 
wilcox.test(UTR_RNA_counts_DENR_responsive$log2FoldChange, UTR_RNA_counts_translated_uORFs$log2_TE_FC)
wilcox.test(UTR_RNA_counts_TMM$log2_TE_FC, UTR_RNA_counts_translated_uORFs$log2_TE_FC) 
wilcox.test(UTR_RNA_counts_DENR_responsive_uORFs$log2FoldChange, UTR_RNA_counts_TMM$log2_TE_FC, paired=FALSE) 
wilcox.test(UTR_RNA_counts_DENR_responsive$log2FoldChange, UTR_RNA_counts_TMM$log2_TE_FC, paired=FALSE) 
