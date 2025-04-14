setwd('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/')
library(DESeq2)
library(apeglm)
library(ashr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggExtra)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)

###Denr shRNA2
ribo_counts = read.table('ribo_counts_Denr.txt', row.names = 1)
rna_counts = read.table('rna_counts_Denr.txt', row.names = 1)
sample_info = read.delim('sample_info_Denr.txt')
sample_info$SeqType <- relevel(as.factor(sample_info$SeqType), "RNA")

#Create DESeq2 object for the combined dataset of Ribo-seq and RNA-seq counts. The interaction term should be included in the linear model design as follows:
ddsMat = DESeqDataSetFromMatrix(
  countData=cbind(ribo_counts,rna_counts),
  colData=sample_info,
  design=~Condition+SeqType+Condition:SeqType)

#Run DESeq2:
ddsMat = DESeq(ddsMat)

#Obtain fold changes for TE:
resultsNames(ddsMat)
res_Denr = results(ddsMat, name='ConditionDenrsh2.SeqTypeRIBO')
write.table(res_Denr[which(res_Denr$padj<0.05), ],
            'DTEGs_Denrsh2_Ctrl.txt', quote=F)
write.table(res_Denr, 'DTEGs_Denrsh2_Ctrl_all.csv', quote=F)
res_Denr <- read.table('DTEGs_Denrsh2_Ctrl.txt', header = TRUE, row.names = NULL)
colnames(res_Denr) <- c('Gene_ID','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
gene_name <- read.csv('Gene_ID_name.csv')
res_gene_name <- merge(res_Denr ,gene_name, by='Gene_ID')
res_gene_name <- res_gene_name[,c(1,8,2:7)]
write.csv(res_gene_name[which(res_gene_name$padj<0.05), ],
          'DTEGs_Denrsh2_vs_Ctrl.csv', quote=F, row.names = FALSE)

#Run DESeq2 for mRNA counts in order to obtain DTGs:
ddsMat_rna = DESeqDataSetFromMatrix(countData=rna_counts, colData=sample_info[which(sample_info$SeqType == 'RNA'),], design=~Condition)
ddsMat_rna = DESeq(ddsMat_rna)
resultsNames(ddsMat_rna)
res_rna_Denr = results(ddsMat_rna, name="Condition_Denrsh2_vs_Ctrl")
res_rna_Denr = lfcShrink(ddsMat_rna, res=res_rna_Denr, coef = 2)
write.csv(res_rna_Denr[which(res_rna_Denr$padj<0.05), ],
          'DTGs_Denrsh2_vs_Ctrl.csv', quote=F)
write.csv(res_rna_Denr,
          'DTGs_Denrsh2_vs_Ctrl_all.csv', quote=F)

#Run DESeq2 for RPFs (Ribo-seq counts):
ddsMat_ribo = DESeqDataSetFromMatrix(countData=ribo_counts, colData=sample_info[which(sample_info$SeqType == 'RIBO'),], design=~Condition)
ddsMat_ribo = DESeq(ddsMat_ribo)
resultsNames(ddsMat_ribo)
res_ribo_Denr = results(ddsMat_ribo,name="Condition_Denrsh2_vs_Ctrl")
res_ribo_Denr =lfcShrink(ddsMat_ribo, res=res_ribo_Denr, coef=2)
write.csv(res_ribo_Denr[which(res_ribo_Denr$padj<0.05), ],
          'Ribo_Denrsh2_vs_Ctrl.csv', quote=F)
write.csv(res_ribo_Denr,
          'Ribo_Denr_vs_Ctrl_all.csv', quote=F)

###Eif2d shRNAs
ribo_counts = read.table('ribo_counts_Eif2dsh.txt', row.names = 1)
rna_counts = read.table('rna_counts_Eif2dsh.txt', row.names = 1)
sample_info = read.delim('sample_info_Eif2dsh.txt')
sample_info$SeqType <- relevel(as.factor(sample_info$SeqType), "RNA")

#Create DESeq2 object for the combined dataset of Ribo-seq and RNA-seq counts. The interaction term should be included in the linear model design as follows:
ddsMat = DESeqDataSetFromMatrix(
  countData=cbind(ribo_counts,rna_counts),
  colData=sample_info,
  design=~Condition+SeqType+Condition:SeqType)

#Run DESeq2:
ddsMat = DESeq(ddsMat)

#Obtain fold changes for TE:
resultsNames(ddsMat)
res_Eif2dsh = results(ddsMat, name='ConditionEif2dsh.SeqTypeRIBO')
write.table(res_Eif2dsh[which(res_Eif2dsh$padj<0.05), ],
            'DTEGs_Eif2dsh_Ctrl.txt', quote=F)
write.table(res_Eif2dsh, 'DTEGs_Eif2dsh_Ctrl_all.csv', quote=F)
res_Eif2dsh <- read.table('DTEGs_Eif2dsh_Ctrl.txt', header = TRUE, row.names = NULL)
colnames(res_Eif2dsh) <- c('Gene_ID','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
gene_name <- read.csv('Gene_ID_name.csv')
res_gene_name <- merge(res_Eif2dsh,gene_name, by='Gene_ID', all.x = TRUE)
res_gene_name <- res_gene_name[,c(1,8,2:7)]
write.csv(res_gene_name[which(res_gene_name$padj<0.05), ],
          'DTEGs_Eif2dsh_vs_Ctrl.csv', quote=F, row.names = FALSE)

#Run DESeq2 for mRNA counts in order to obtain DTGs:
ddsMat_rna = DESeqDataSetFromMatrix(countData=rna_counts, colData=sample_info[which(sample_info$SeqType == 'RNA'),], design=~Condition)
ddsMat_rna = DESeq(ddsMat_rna)
resultsNames(ddsMat_rna)
res_rna_Eif2dsh = results(ddsMat_rna, name="Condition_Eif2dsh_vs_Ctrl")
res_rna_Eif2dsh = lfcShrink(ddsMat_rna, res=res_rna_Eif2dsh, coef = 2)
write.csv(res_rna_Eif2dsh[which(res_rna_Eif2dsh$padj<0.05), ],
          'DTGs_Eif2dsh_vs_Ctrl.csv', quote=F)
write.csv(res_rna_Eif2dsh,
          'DTGs_Eif2dsh_vs_Ctrl_all.csv', quote=F)

#Run DESeq2 for RPFs (Ribo-seq counts):
ddsMat_ribo = DESeqDataSetFromMatrix(countData=ribo_counts, colData=sample_info[which(sample_info$SeqType == 'RIBO'),], design=~Condition)
ddsMat_ribo = DESeq(ddsMat_ribo)
resultsNames(ddsMat_ribo)
res_ribo_Eif2dsh = results(ddsMat_ribo,name="Condition_Eif2dsh_vs_Ctrl")
res_ribo_Eif2dsh =lfcShrink(ddsMat_ribo, res=res_ribo_Eif2dsh, coef=2)
write.csv(res_ribo_Eif2dsh[which(res_ribo_Eif2dsh$padj<0.05), ],
          'Ribo_Eif2dsh_vs_Ctrl.csv', quote=F)
write.csv(res_ribo_Eif2dsh,
          'Ribo_Eif2dsh_vs_Ctrl_all.csv', quote=F)

#Obtain genes for each regulation class described in Figure 1D, E

#Denr shRNA2
res_Denr <- read.csv('DTEGs_Denrsh2_Ctrl_all.csv')
res_ribo_Denr <- read.csv('Ribo_Denr_vs_Ctrl_all.csv', row.names = 1)
res_rna_Denr <- read.csv('DTGs_Denrsh2_vs_Ctrl_all.csv', row.names = 1)
forwarded_Denr = rownames(res_Denr)[which(res_Denr$padj > 0.1
                                          & res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj < 0.1)]
exclusive_Denr = rownames(res_Denr)[which(res_Denr$padj < 0.1
                                          & res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj > 0.1)]
both_Denr = rownames(res_Denr)[which(res_Denr$padj < 0.1 &
                                       res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj < 0.1)]
intensified_Denr = rownames(res_Denr[both_Denr[which(res_Denr[both_Denr,2]*res_rna_Denr[both_Denr,2] > 0)],])
buffered_Denr = rownames(res_Denr[both_Denr[which(res_Denr[both_Denr,2]*res_rna_Denr[both_Denr,2] < 0)],])
buffered_Denr = c(rownames(res_Denr)[which(res_Denr$padj < 0.1 & res_ribo_Denr$padj > 0.1 & res_rna_Denr$padj < 0.1)], buffered_Denr)
max_val_Denr = max(res_ribo_Denr[,2],res_rna_Denr[,2],na.rm = T)

res_Denr <- read.table('DTEGs_Denrsh2_Ctrl_all.csv', row.names = NULL)
colnames(res_Denr) <- c('Gene_ID','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
write.csv(forwarded_Denr, 'Denr_forwarded_0.1.csv', row.names = FALSE)
forwarded_Denr <- read.csv('Denr_forwarded_0.1.csv')
colnames(forwarded_Denr) <- 'Gene_ID'
forwarded_Denr <- merge(forwarded_Denr, gene_name, by = 'Gene_ID')
forwarded_Denr <- merge(forwarded_Denr, res_Denr, by = 'Gene_ID')
write.csv(forwarded_Denr, 'Denr_forwarded_0.1.csv', row.names = FALSE)
write.csv(buffered_Denr, 'Denr_buffered_0.1.csv', row.names = FALSE)
buffered_Denr <- read.csv('Denr_buffered_0.1.csv')
colnames(buffered_Denr) <- 'Gene_ID'
buffered_Denr <- merge(buffered_Denr, gene_name, by = 'Gene_ID')
buffered_Denr <- merge(buffered_Denr, res_Denr, by = 'Gene_ID')
write.csv(buffered_Denr, 'Denr_buffered_0.1.csv', row.names = FALSE)
write.csv(both_Denr, 'Denr_both_0.1.csv', row.names = FALSE)
both_Denr <- read.csv('Denr_both_0.1.csv')
colnames(both_Denr) <- 'Gene_ID'
both_Denr <- merge(both_Denr, gene_name, by = 'Gene_ID')
both_Denr <- merge(both_Denr, res_Denr, by = 'Gene_ID')
write.csv(both_Denr, 'Denr_both_0.1.csv', row.names = FALSE)
write.csv(exclusive_Denr, 'Denr_exclusive_0.1.csv', row.names = FALSE)
exclusive_Denr <- read.csv('Denr_exclusive_0.1.csv')
colnames(exclusive_Denr) <- 'Gene_ID'
exclusive_Denr <- merge(exclusive_Denr, gene_name, by = 'Gene_ID')
exclusive_Denr <- merge(exclusive_Denr, res_Denr, by = 'Gene_ID')
write.csv(exclusive_Denr, 'Denr_exclusive_0.1.csv', row.names = FALSE)
write.csv(intensified_Denr, 'Denr_intensified_0.1.csv', row.names = FALSE)
intensified_Denr <- read.csv('Denr_intensified_0.1.csv')
colnames(intensified_Denr) <- 'Gene_ID'
intensified_Denr <- merge(intensified_Denr, gene_name, by = 'Gene_ID')
intensified_Denr <- merge(intensified_Denr, res_Denr, by = 'Gene_ID')
write.csv(intensified_Denr, 'Denr_intensified_0.1.csv', row.names = FALSE)

exclu <- read.csv('Denr_exclusive_0.1.csv', row.names = NULL)
intens <- read.csv('Denr_intensified_0.1.csv', row.names = NULL)
buffer <- read.csv('Denr_buffered_0.1.csv', row.names = NULL)
all <- merge(exclu,intens, all = TRUE)
all <- merge(all, buffer, all = TRUE)
all_pos <- all[all$log2FoldChange > 0,]
write.csv(all_pos, 'Denrsh2_all_pos_TE_FC_FDR_0.1.csv', row.names = FALSE)
all_neg <- all[all$log2FoldChange < 0,]
write.csv(all_neg, 'Denrsh2_all_neg_TE_FC_FDR_0.1.csv', row.names = FALSE)

#Eif2d shRNAs
res_Eif2dsh <- read.csv('DTEGs_Eif2dsh_Ctrl_all.csv')
res_ribo_Eif2dsh <- read.csv('Ribo_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
res_rna_Eif2dsh <- read.csv('DTGs_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
forwarded_Eif2dsh = rownames(res_Eif2dsh)[which(res_Eif2dsh$padj > 0.1
                                                & res_ribo_Eif2dsh$padj < 0.1 & res_rna_Eif2dsh$padj < 0.1)]
exclusive_Eif2dsh = rownames(res_Eif2dsh)[which(res_Eif2dsh$padj < 0.1
                                                & res_ribo_Eif2dsh$padj < 0.1 & res_rna_Eif2dsh$padj > 0.1)]
both_Eif2dsh = rownames(res_Eif2dsh)[which(res_Eif2dsh$padj < 0.1 &
                                             res_ribo_Eif2dsh$padj < 0.1 & res_rna_Eif2dsh$padj < 0.1)]
intensified_Eif2dsh = rownames(res_Eif2dsh[both_Eif2dsh[which(res_Eif2dsh[both_Eif2dsh,2]*res_rna_Eif2dsh[both_Eif2dsh,2] > 0)],])
buffered_Eif2dsh = rownames(res_Eif2dsh[both_Eif2dsh[which(res_Eif2dsh[both_Eif2dsh,2]
                                                           *res_rna_Eif2dsh[both_Eif2dsh,2] < 0)],])
buffered_Eif2dsh = c(rownames(res_Eif2dsh)[which(res_Eif2dsh$padj < 0.1
                                                 & res_ribo_Eif2dsh$padj > 0.1 & res_rna_Eif2dsh$padj < 0.1)],
                     buffered_Eif2dsh)
max_val_Eif2dsh = max(res_ribo_Eif2dsh[,3],res_rna_Eif2dsh[,3],na.rm = T)

res_Eif2dsh <- read.csv('DTEGs_Eif2dsh_Ctrl_all.csv', row.names = NULL)
colnames(res_Eif2dsh) <- c('Gene_ID','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
write.csv(forwarded_Eif2dsh, 'Eif2dsh_forwarded_0.1.csv', row.names = FALSE)
forwarded_Eif2dsh <- read.csv('Eif2dsh_forwarded_0.1.csv')
colnames(forwarded_Eif2dsh) <- 'Gene_ID'
forwarded_Eif2dsh <- merge(forwarded_Eif2dsh, gene_name, by = 'Gene_ID')
forwarded_Eif2dsh <- merge(forwarded_Eif2dsh, res_Eif2dsh, by = 'Gene_ID')
write.csv(forwarded_Eif2dsh, 'Eif2dsh_forwarded_0.1.csv', row.names = FALSE)
write.csv(buffered_Eif2dsh, 'Eif2dsh_buffered_0.1.csv', row.names = FALSE)
buffered_Eif2dsh <- read.csv('Eif2dsh_buffered_0.1.csv')
colnames(buffered_Eif2dsh) <- 'Gene_ID'
buffered_Eif2dsh <- merge(buffered_Eif2dsh, gene_name, by = 'Gene_ID')
buffered_Eif2dsh <- merge(buffered_Eif2dsh, res_Eif2dsh, by = 'Gene_ID')
write.csv(buffered_Eif2dsh, 'Eif2dsh_buffered_0.1.csv', row.names = FALSE)
write.csv(both_Eif2dsh, 'Eif2dsh_both_0.1.csv', row.names = FALSE)
both_Eif2dsh <- read.csv('Eif2dsh_both_0.1.csv')
colnames(both_Eif2dsh) <- 'Gene_ID'
both_Eif2dsh <- merge(both_Eif2dsh, gene_name, by = 'Gene_ID')
both_Eif2dsh <- merge(both_Eif2dsh, res_Eif2dsh, by = 'Gene_ID')
write.csv(both_Eif2dsh, 'Eif2dsh_both_0.1.csv', row.names = FALSE)
write.csv(exclusive_Eif2dsh, 'Eif2dsh_exclusive_0.1.csv', row.names = FALSE)
exclusive_Eif2dsh <- read.csv('Eif2dsh_exclusive_0.1.csv')
colnames(exclusive_Eif2dsh) <- 'Gene_ID'
exclusive_Eif2dsh <- merge(exclusive_Eif2dsh, gene_name, by = 'Gene_ID')
exclusive_Eif2dsh <- merge(exclusive_Eif2dsh, res_Eif2dsh, by = 'Gene_ID')
write.csv(exclusive_Eif2dsh, 'Eif2dsh_exclusive_0.1.csv', row.names = FALSE)
write.csv(intensified_Eif2dsh, 'Eif2dsh_intensified_0.1.csv', row.names = FALSE)
intensified_Eif2dsh <- read.csv('Eif2dsh_intensified_0.1.csv')
colnames(intensified_Eif2dsh) <- 'Gene_ID'
intensified_Eif2dsh <- merge(intensified_Eif2dsh, gene_name, by = 'Gene_ID')
intensified_Eif2dsh <- merge(intensified_Eif2dsh, res_Eif2dsh, by = 'Gene_ID')
write.csv(intensified_Eif2dsh, 'Eif2dsh_intensified_0.1.csv', row.names = FALSE)

exclu <- read.csv('Eif2dsh_exclusive_0.1.csv', row.names = NULL)
intens <- read.csv('Eif2dsh_intensified_0.1.csv', row.names = NULL)
buffer <- read.csv('Eif2dsh_buffered_0.1.csv', row.names = NULL)
all <- merge(exclu,intens, all = TRUE)
all <- merge(all, buffer, all = TRUE)
all_pos <- all[all$log2FoldChange > 0,]
write.csv(all_pos, 'Eif2dsh_all_pos_TE_FC_FDR_0.1.csv', row.names = FALSE)
all_neg <- all[all$log2FoldChange < 0,]
write.csv(all_neg, 'Eif2dsh_all_neg_TE_FC_FDR_0.1.csv', row.names = FALSE)

#Visualize the global translational and transcriptional regulation as in Figure 1E
res_rna_Denr <- read.csv('DTGs_Denrsh2_vs_Ctrl_all.csv', row.names = 1)
res_ribo_Denr <- read.csv('Ribo_Denr_vs_Ctrl_all.csv', row.names = 1)
pdf('DTG_DTEG_Denr_plot_correct.pdf')
plot(y=res_ribo_Denr[,2],x=res_rna_Denr[,2],
     xlab="RNA-seq log2 fold change",
     ylab = "Ribo-seq log2 fold change", asp=1,
     pch=16,
     col=rgb(128/255,128/255,128/255,0.1), ylim=
       c(-max_val_Denr,max_val_Denr), xlim=c(-max_val_Denr,max_val_Denr),
     cex=0.4)
abline(a=0,b=1,col="gray")
abline(h=0,v=0,col="gray")
points(y=res_ribo_Denr[intensified_Denr,2], x=res_rna_Denr
       [intensified_Denr,2], cex = 0.25,
       pch=16,col='skyblue')
points(y=res_ribo_Denr[buffered_Denr,2], x=res_rna_Denr
       [buffered_Denr,2], cex = 0.25,
       pch=16,col='purple')
points(y=res_ribo_Denr[forwarded_Denr,2], x=res_rna_Denr
       [forwarded_Denr,2], cex = 0.25,
       pch=16,col='orange')
points(y=res_ribo_Denr[exclusive_Denr,2], x=res_rna_Denr
       [exclusive_Denr,2], cex = 0.25,
       pch=16,col='red')
text(abs_losses, percent_losses, labels=namebank, cex= 0.7, pos=3)
dev.off()

res_rna_Eif2dsh <- read.csv('DTGs_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
res_ribo_Eif2dsh <- read.csv('Ribo_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
pdf('DTG_DTEG_Eif2dsh_plot_correct.pdf')
plot(y=res_ribo_Eif2dsh[,2],x=res_rna_Eif2dsh[,2],
     xlab="RNA-seq log2 fold change",
     ylab = "Ribo-seq log2 fold change", asp=1,
     pch=16,
     col=rgb(128/255,128/255,128/255,0.1), ylim=
       c(-6,6), xlim=c(-6,6),
     cex=0.4)
abline(a=0,b=1,col="gray")
abline(h=0,v=0,col="gray")
points(y=res_ribo_Eif2dsh[intensified_Eif2dsh,2], x=res_rna_Eif2dsh
       [intensified_Eif2dsh,2], cex = 0.25,
       pch=16,col='skyblue')
points(y=res_ribo_Eif2dsh[buffered_Eif2dsh,2], x=res_rna_Eif2dsh
       [buffered_Eif2dsh,2], cex = 0.25,
       pch=16,col='purple')
points(y=res_ribo_Eif2dsh[forwarded_Eif2dsh,2], x=res_rna_Eif2dsh
       [forwarded_Eif2dsh,2], cex = 0.25,
       pch=16,col='orange')
points(y=res_ribo_Eif2dsh[exclusive_Eif2dsh,2], x=res_rna_Eif2dsh
       [exclusive_Eif2dsh,2], cex = 0.25,
       pch=16,col='red')
dev.off()


# GO enrichment analysis all negative Denr shRNA2
all_neg <- read.csv('Denrsh2_all_negative_TE_FC_FDR_0.1.csv', row.names = NULL)
all_neg_genes <- all_neg[,1]
ego <- enrichGO(gene = all_neg_genes, keyType = "ENSEMBL",OrgDb='org.Mm.eg.db', ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, 'GO_analysis_BP_Denrsh2_all_negative.csv')
pdf('GO_analysis_plots_BP_Denrsh2_all_negative_FDR_0.1.pdf',height=12,width=13)
#category netplot
signif_res_lFC <- res_Denr$log2FoldChange # To color genes by log2 fold changes
cnetplot(ego, categorySize="pvalue", showCategory = 5, foldChange= signif_res_lFC, vertex.label.font=6)
#barplot
barplot(ego, showCategory=20)
dev.off()

# GO enrichment analysis all negative Eif2d shRNAs
all_neg <- read.csv('Eif2dsh_all_negative_TE_FC_FDR_0.1.csv', row.names = NULL)
exclusive_neg_genes <- exclusive_neg[,1]
ego <- enrichGO(gene = exclusive_neg_genes, keyType = "ENSEMBL",OrgDb='org.Mm.eg.db', ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary,  'GO_analysis_BP_Eif2dsh_all_negative_FDR_0.1.csv')
pdf('GO_analysis_plots_BP_Eif2dsh_all_negative_FDR_0.1.pdf',height=12,width=13)
#category netplot
signif_res_lFC <- res_Eif2dsh$log2FoldChange # To color genes by log2 fold changes
cnetplot(ego, categorySize="pvalue", showCategory = 5, foldChange= signif_res_lFC, vertex.label.font=6)
#barplot
barplot(ego, showCategory=20)
dev.off()

#Visualize the global TE and RNA abundance

res_Denr <- read.csv('DTEGs_Denrsh2_Ctrl_all.csv', row.names = 1)
res_ribo_Denr <- read.csv('Ribo_Denr_vs_Ctrl_all.csv', row.names = 1)
gene_name <- read.csv('Gene_ID_name.csv')
res_rna_Denr <- read.csv('DTGs_Denrsh2_vs_Ctrl_all.csv', row.names = 1)
forwarded_Denr_0.1 = rownames(res_Denr)[which(res_Denr$padj > 0.1
                                              & res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj < 0.1)]
exclusive_Denr_0.1 = rownames(res_Denr)[which(res_Denr$padj < 0.1
                                              & res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj > 0.1)]
both_Denr_0.1 = rownames(res_Denr)[which(res_Denr$padj < 0.1 &
                                           res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj < 0.1)]
intensified_Denr_0.1 = rownames(res_Denr[both_Denr_0.1[which(res_Denr[both_Denr_0.1,
                                                                      2]*res_rna_Denr[both_Denr_0.1,2] > 0)],])
buffered_Denr_0.1 = rownames(res_Denr[both_Denr_0.1[which(res_Denr[both_Denr_0.1,2]
                                                          *res_rna_Denr[both_Denr_0.1,2] < 0)],])
buffered_Denr_0.1 = c(rownames(res_Denr)[which(res_Denr$padj < 0.1
                                               & res_ribo_Denr$padj > 0.1 & res_rna_Denr$padj < 0.1)],
                      buffered_Denr_0.1)
max_val_Denr = max(res_ribo_Denr[,2],res_rna_Denr[,2],na.rm = T)

klhdc8a_1 <- 'ENSMUSG00000042115'
Asb8_1 <- 'ENSMUSG00000048175' 
klhdc8a_1 <- 'ENSMUSG00000042115'
Hoxa3_1 <- 'ENSMUSG00000079560' 
Med23_1 <- 'ENSMUSG00000019984'
Rps20_1 <- 'ENSMUSG00000028234'
Ndc80_1 <- 'ENSMUSG00000024056'
Cenpa_1 <- 'ENSMUSG00000029177'

pdf('DTG_TE_Eif2d_scatter_density_plot_colors.pdf')
df_Eif2d <- data.frame(x = res_rna_Eif2d[,2], y = res_Eif2d[,2])
df_no_na <- na.omit(df_Eif2d) 
intensified <- data.frame(y=res_Eif2d[intensified_Denr_0.1,2], x=res_rna_Denr[intensified_Denr_0.1,2])
buffered <- data.frame(y=res_Eif2d[buffered_Eif2d_0.1,2], x=res_rna_Denr[buffered_Denr_0.1,2])
forwarded <- data.frame(y=res_Denr[forwarded_Denr_0.1,2], x=res_rna_Denr[forwarded_Denr_0.1,2])
exclusive <- data.frame(y=res_Denr[exclusive_Denr_0.1,2], x=res_rna_Denr[exclusive_Denr_0.1,2])
klhdc8a <- data.frame(y=res_Denr[klhdc8a_1,2], x=res_rna_Denr[klhdc8a_1,2])
Asb8 <- data.frame(y=res_Denr[Asb8_1,2], x=res_rna_Denr[Asb8_1,2])
DTEGs <- bind_rows(intensified,buffered,exclusive)
p <- ggplot(df_Denr, mapping = aes(x, y)) + geom_point(color = rgb(128/255,128/255,128/255,0.1), pch = 16, cex = 0.5)+ ylim(-3,3) + xlim(-3,3) + theme_classic() + geom_hline(yintercept=0, color = "gray", size=0.2) + geom_vline(xintercept=0, color = "gray", size=0.2)  + labs(x = "RNA-seq log2 fold change", y = 'TE log2 fold change') + 
  geom_point(buffered, mapping = aes(x,y), cex = 1.5, pch=16,color='#4dbe56') +
  geom_point(exclusive, mapping = aes(x,y), cex = 1.5, pch=16,color='#4dbe56') +
  geom_point(forwarded, mapping = aes(x,y), cex = 1, pch=16,color='#bc72ab') +
  geom_point(klhdc8a, mapping = aes(), cex = 3, pch=16,color='black') +
  geom_point(Asb8, mapping = aes(), cex = 3, pch=16,color='black') +
  geom_point(Ssbp3, mapping = aes(), cex = 3, pch=16,color='black')
xplot <- ggdensity(forwarded, "x", color = '#bc72ab', xlim=c(-3,3))
yplot <- ggdensity(DTEGs, "y", color = '#4dbe56') + coord_flip(xlim=c(-3,3))
p <- p + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
plot_grid(xplot, NULL, p, yplot, ncol = 2, align = "hv",
          rel_widths = c(2, 1), rel_heights = c(1, 2))
dev.off()

res_Eif2d <- read.csv('DTEGs_Eif2dsh_Ctrl_all.csv', row.names = 1)
res_ribo_Eif2d <- read.csv('Ribo_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
gene_name <- read.csv('Gene_ID_name.csv')
res_rna_Eif2d <- read.csv('DTGs_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
forwarded_Eif2d_0.1 = rownames(res_Eif2d)[which(res_Eif2d$padj > 0.1
                                                & res_ribo_Eif2d$padj < 0.1 & res_rna_Eif2d$padj < 0.1)]
exclusive_Eif2d_0.1 = rownames(res_Eif2d)[which(res_Eif2d$padj < 0.1
                                                & res_ribo_Eif2d$padj < 0.1 & res_rna_Eif2d$padj > 0.1)]
both_Eif2d_0.1 = rownames(res_Eif2d)[which(res_Eif2d$padj < 0.1 &
                                             res_ribo_Eif2d$padj < 0.1 & res_rna_Eif2d$padj < 0.1)]
intensified_Eif2d_0.1 = rownames(res_Eif2d[both_Eif2d_0.1[which(res_Eif2d[both_Eif2d_0.1,
                                                                          2]*res_rna_Eif2d[both_Eif2d_0.1,2] > 0)],])
buffered_Eif2d_0.1 = rownames(res_Eif2d[both_Eif2d_0.1[which(res_Eif2d[both_Eif2d_0.1,2]
                                                             *res_rna_Eif2d[both_Eif2d_0.1,2] < 0)],])
buffered_Eif2d_0.1 = c(rownames(res_Eif2d)[which(res_Eif2d$padj < 0.1
                                                 & res_ribo_Eif2d$padj > 0.1 & res_rna_Eif2d$padj < 0.1)],
                       buffered_Eif2d_0.1)
max_val_Eif2d = max(res_ribo_Eif2d[,2],res_rna_Eif2d[,2],na.rm = T)

pdf('DTG_TE_Eif2d_scatter_density_plot_Klhdc8a_Asb8_Med23_Rps20_Ndc80_Hoxa3_Cenpa.pdf')
df_Eif2d <- data.frame(x = res_rna_Eif2d[,2], y = res_Eif2d[,2])
df_no_na <- na.omit(df_Eif2d) 
intensified <- data.frame(y=res_Eif2d[intensified_Eif2d_0.1,2], x=res_rna_Eif2d[intensified_Eif2d_0.1,2])
buffered <- data.frame(y=res_Eif2d[buffered_Eif2d_0.1,2], x=res_rna_Eif2d[buffered_Eif2d_0.1,2])
forwarded <- data.frame(y=res_Eif2d[forwarded_Eif2d_0.1,2], x=res_rna_Eif2d[forwarded_Eif2d_0.1,2])
exclusive <- data.frame(y=res_Eif2d[exclusive_Eif2d_0.1,2], x=res_rna_Eif2d[exclusive_Eif2d_0.1,2])
klhdc8a <- data.frame(y=res_Eif2d[klhdc8a_1,2], x=res_rna_Eif2d[klhdc8a_1,2])
Asb8 <- data.frame(y=res_Eif2d[Asb8_1,2], x=res_rna_Eif2d[Asb8_1,2])
Med23 <- data.frame(y=res_Eif2d[Med23_1,2], x=res_rna_Eif2d[Med23_1,2])
Cenpa <- data.frame(y=res_Eif2d[Cenpa_1,2], x=res_rna_Eif2d[Cenpa_1,2])
Rps20 <- data.frame(y=res_Eif2d[Rps20_1,2], x=res_rna_Eif2d[Rps20_1,2])
Hoxa3 <- data.frame(y=res_Eif2d[Hoxa3_1,2], x=res_rna_Eif2d[Hoxa3_1,2])
Ndc80 <- data.frame(y=res_Eif2d[Ndc80_1,2], x=res_rna_Eif2d[Ndc80_1,2])
DTEGs <- bind_rows(intensified,buffered,exclusive)
ggplot(df_Eif2d, mapping = aes(x, y)) + geom_point(color = rgb(128/255,128/255,128/255,0.1), pch = 16, cex = 0.5)+ ylim(-3,3) + xlim(-3,3) + theme_classic() + geom_hline(yintercept=0, color = "gray", size=0.2) + geom_vline(xintercept=0, color = "gray", size=0.2)  + labs(x = "RNA-seq log2 fold change", y = 'TE log2 fold change') + 
  geom_point(buffered, mapping = aes(x,y), cex = 2, pch=16,color='#2b9cdeff') +
  geom_point(exclusive, mapping = aes(x,y), cex = 2, pch=16,color='#2b9cdeff') +
  geom_point(forwarded, mapping = aes(x,y), cex = 1.5, pch=16,color='#bc72abff') +
  geom_point(klhdc8a, mapping = aes(), cex = 2, pch=16,color='black') +
  geom_point(Asb8, mapping = aes(), cex = 2, pch=16,color='black') +
  geom_point(Med23, mapping = aes(), cex = 2, pch=16,color='black') +
  geom_point(Cenpa, mapping = aes(), cex = 2, pch=16,color='black') +
  geom_point(Rps20, mapping = aes(), cex = 2, pch=16,color='black') +
  geom_point(Ndc80, mapping = aes(), cex = 2, pch=16,color='black') 
dev.off()

# With density plots

res_Denr <- read.csv('DTEGs_Denrsh2_Ctrl_all.csv', row.names = 1)
res_ribo_Denr <- read.csv('Ribo_Denr_vs_Ctrl_all.csv', row.names = 1)
gene_name <- read.csv('Gene_ID_name.csv')
res_rna_Denr <- read.csv('DTGs_Denrsh2_vs_Ctrl_all.csv', row.names = 1)
forwarded_Denr_0.1 = rownames(res_Denr)[which(res_Denr$padj > 0.1
                                              & res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj < 0.1)]
exclusive_Denr_0.1 = rownames(res_Denr)[which(res_Denr$padj < 0.1
                                              & res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj > 0.1)]
both_Denr_0.1 = rownames(res_Denr)[which(res_Denr$padj < 0.1 &
                                           res_ribo_Denr$padj < 0.1 & res_rna_Denr$padj < 0.1)]
intensified_Denr_0.1 = rownames(res_Denr[both_Denr_0.1[which(res_Denr[both_Denr_0.1,
                                                                      2]*res_rna_Denr[both_Denr_0.1,2] > 0)],])
buffered_Denr_0.1 = rownames(res_Denr[both_Denr_0.1[which(res_Denr[both_Denr_0.1,2]
                                                          *res_rna_Denr[both_Denr_0.1,2] < 0)],])
buffered_Denr_0.1 = c(rownames(res_Denr)[which(res_Denr$padj < 0.1
                                               & res_ribo_Denr$padj > 0.1 & res_rna_Denr$padj < 0.1)],
                      buffered_Denr_0.1)
max_val_Denr = max(res_ribo_Denr[,2],res_rna_Denr[,2],na.rm = T)


pdf('DTG_TE_Denr_scatter_density_plot_colors_Klhdc8a_Asb8_Ssbp3.pdf')
df_Denr <- data.frame(x = res_rna_Denr[,2], y = res_Denr[,2])
df_no_na <- na.omit(df_Denr) 
intensified <- data.frame(y=res_Denr[intensified_Denr_0.1,2], x=res_rna_Denr[intensified_Denr_0.1,2])
buffered <- data.frame(y=res_Denr[buffered_Denr_0.1,2], x=res_rna_Denr[buffered_Denr_0.1,2])
forwarded <- data.frame(y=res_Denr[forwarded_Denr_0.1,2], x=res_rna_Denr[forwarded_Denr_0.1,2])
exclusive <- data.frame(y=res_Denr[exclusive_Denr_0.1,2], x=res_rna_Denr[exclusive_Denr_0.1,2])
klhdc8a <- data.frame(y=res_Denr[klhdc8a_1,2], x=res_rna_Denr[klhdc8a_1,2])
Asb8 <- data.frame(y=res_Denr[Asb8_1,2], x=res_rna_Denr[Asb8_1,2])
DTEGs <- bind_rows(intensified,buffered,exclusive)
p <- ggplot(df_Denr, mapping = aes(x, y)) + geom_point(color = rgb(128/255,128/255,128/255,0.1), pch = 16, cex = 0.5)+ ylim(-3,3) + xlim(-3,3) + theme_classic() + geom_hline(yintercept=0, color = "gray", size=0.2) + geom_vline(xintercept=0, color = "gray", size=0.2)  + labs(x = "RNA-seq log2 fold change", y = 'TE log2 fold change') + 
  geom_point(buffered, mapping = aes(x,y), cex = 1.5, pch=16,color='#4dbe56') +
  geom_point(exclusive, mapping = aes(x,y), cex = 1.5, pch=16,color='#4dbe56') +
  geom_point(forwarded, mapping = aes(x,y), cex = 1, pch=16,color='#bc72ab') +
  geom_point(klhdc8a, mapping = aes(), cex = 3, pch=16,color='black') +
  geom_point(Asb8, mapping = aes(), cex = 3, pch=16,color='black') 
xplot <- ggdensity(forwarded, "x", color = '#bc72ab', xlim=c(-3,3))
yplot <- ggdensity(DTEGs, "y", color = '#4dbe56') + coord_flip(xlim=c(-3,3))
p <- p + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
plot_grid(xplot, NULL, p, yplot, ncol = 2, align = "hv",
          rel_widths = c(2, 1), rel_heights = c(1, 2))
dev.off()


res_Eif2d <- read.csv('DTEGs_Eif2dsh_Ctrl_all.csv', row.names = 1)
res_ribo_Eif2d <- read.csv('Ribo_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
gene_name <- read.csv('Gene_ID_name.csv')
res_rna_Eif2d <- read.csv('DTGs_Eif2dsh_vs_Ctrl_all.csv', row.names = 1)
forwarded_Eif2d_0.1 = rownames(res_Eif2d)[which(res_Eif2d$padj > 0.1
                                                & res_ribo_Eif2d$padj < 0.1 & res_rna_Eif2d$padj < 0.1)]
exclusive_Eif2d_0.1 = rownames(res_Eif2d)[which(res_Eif2d$padj < 0.1
                                                & res_ribo_Eif2d$padj < 0.1 & res_rna_Eif2d$padj > 0.1)]
both_Eif2d_0.1 = rownames(res_Eif2d)[which(res_Eif2d$padj < 0.1 &
                                             res_ribo_Eif2d$padj < 0.1 & res_rna_Eif2d$padj < 0.1)]
intensified_Eif2d_0.1 = rownames(res_Eif2d[both_Eif2d_0.1[which(res_Eif2d[both_Eif2d_0.1,
                                                                          2]*res_rna_Eif2d[both_Eif2d_0.1,2] > 0)],])
buffered_Eif2d_0.1 = rownames(res_Eif2d[both_Eif2d_0.1[which(res_Eif2d[both_Eif2d_0.1,2]
                                                             *res_rna_Eif2d[both_Eif2d_0.1,2] < 0)],])
buffered_Eif2d_0.1 = c(rownames(res_Eif2d)[which(res_Eif2d$padj < 0.1
                                                 & res_ribo_Eif2d$padj > 0.1 & res_rna_Eif2d$padj < 0.1)],
                       buffered_Eif2d_0.1)
max_val_Eif2d = max(res_ribo_Eif2d[,2],res_rna_Eif2d[,2],na.rm = T)


pdf('DTG_TE_Eif2d_scatter_density_plot_colors_Klhdc8a_Asb8_Ssbp3.pdf')
df_Eif2d <- data.frame(x = res_rna_Eif2d[,2], y = res_Eif2d[,2])
df_no_na <- na.omit(df_Eif2d) 
intensified <- data.frame(y=res_Eif2d[intensified_Eif2d_0.1,2], x=res_rna_Eif2d[intensified_Eif2d_0.1,2])
buffered <- data.frame(y=res_Eif2d[buffered_Eif2d_0.1,2], x=res_rna_Eif2d[buffered_Eif2d_0.1,2])
forwarded <- data.frame(y=res_Eif2d[forwarded_Eif2d_0.1,2], x=res_rna_Eif2d[forwarded_Eif2d_0.1,2])
exclusive <- data.frame(y=res_Eif2d[exclusive_Eif2d_0.1,2], x=res_rna_Eif2d[exclusive_Eif2d_0.1,2])
klhdc8a <- data.frame(y=res_Eif2d[klhdc8a_1,2], x=res_rna_Eif2d[klhdc8a_1,2])
Asb8 <- data.frame(y=res_Eif2d[Asb8_1,2], x=res_rna_Eif2d[Asb8_1,2])
DTEGs <- bind_rows(intensified,buffered,exclusive)
p <- ggplot(df_Eif2d, mapping = aes(x, y)) + geom_point(color = rgb(128/255,128/255,128/255,0.1), pch = 16, cex = 0.5)+ ylim(-3,3) + xlim(-3,3) + theme_classic() + geom_hline(yintercept=0, color = "gray", size=0.2) + geom_vline(xintercept=0, color = "gray", size=0.2)  + labs(x = "RNA-seq log2 fold change", y = 'TE log2 fold change') + 
  geom_point(buffered, mapping = aes(x,y), cex = 1.5, pch=16,color='#2b9cdeff') +
  geom_point(exclusive, mapping = aes(x,y), cex = 1.5, pch=16,color='#2b9cdeff') +
  geom_point(forwarded, mapping = aes(x,y), cex = 1, pch=16,color='#bc72ab') +
  geom_point(klhdc8a, mapping = aes(), cex = 3, pch=16,color='black') +
  geom_point(Asb8, mapping = aes(), cex = 3, pch=16,color='black') 
xplot <- ggdensity(forwarded, "x", color = '#bc72ab', xlim=c(-3,3))
yplot <- ggdensity(DTEGs, "y", color = '#2b9cdeff') + coord_flip(xlim=c(-3,3))
p <- p + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
plot_grid(xplot, NULL, p, yplot, ncol = 2, align = "hv",
          rel_widths = c(2, 1), rel_heights = c(1, 2))
dev.off()
