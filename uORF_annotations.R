setwd('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/uORF_annotations/')
require(dplyr)
library(stringr)
library(ggplot2)
library(vioplot)
library(reshape2)

# Select gene ID and transcript ID of genes that have only one transcript isoform expressed or select the most representative one (based on prepare_cds.sh from Bulak script)
prepared_cds <-read.table('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/prepare_cds/Ctrl/Eif2d_prepared_cds.tsv')
single_proteins <- prepared_cds[prepared_cds$V2 == 'SINGLE_PROTEIN',]
single_proteins_list <- single_proteins[single_proteins$V5 != 'composite', ]
single_proteins_list <- single_proteins_list[,c(1,5)]
colnames(single_proteins_list) <- c('Gene_ID','Transcript_ID')
write.csv(single_proteins_list,'GeneID_transcriptID_single_isoform.csv')
non_single_longest <- read.table('GeneID_transcriptID_non_single_isoform_sameCDS_longest.txt')
colnames(non_single_longest) <- c('Gene_ID','Transcript_ID')
single_isoforms_non_single_isoforms_longest <- rbind(single_proteins_list, non_single_longest)
write.table(single_isoforms_non_single_isoforms_longest,'GeneID_transcriptID_single_isoform_sameCDS_longest.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
seq_5UTR <- read.table('../annotations/Mus_musculus.GRCm38.100.5UTR.seq.PC.csv')
colnames(seq_5UTR) <- c('Gene_ID','Transcript_ID','seq')
single_isoforms_non_single_isoforms_longest_5UTR <- merge(single_isoforms_non_single_isoforms_longest, seq_5UTR, by = c('Gene_ID', 'Transcript_ID'))
write.table(single_isoforms_non_single_isoforms_longest_5UTR , 'GeneID_transcriptID_5UTRseq_single_isoform_sameCDS_longest.txt')

#Annotate potential uORFs with potential_uORF_annotation.py (overlapping uORF not included)
#count number of reads on uorf. For those that are overlapping, count/uORF_length and take the one with highest score. normalize to the average CDS coverage for checking...
#Bed intersect for counting reads mapping to uORFs : bedtools intersect -b bamfile.bam -a GeneID_transcriptID_5UTR_single_isoform_sameCDS_longest_uORFs_16_PC.bed -bed -c -f 0.5 > output.count.bed
#GeneID_transcriptID_5UTR_single_isoform_sameCDS_longest_uORFs_16_PC.bed has a 16nt A site shift

count <- read.table('RP_ctrl.mouse_cDNA.uORF_16.longest.count.bed')
colnames(count) <- c('IDs', 'start_uORF', 'stop_uORF', 'protein_coding', 'counts')
count <- count[order(count$IDs, count$stop_uORF), ]
best_uORF <- as.data.frame(c())
Gene_IDs <- unique(as.character(count$IDs))
for (i in Gene_IDs) {
  List = list()
  Name = list()
  for (n in 1:nrow(count)) {
    IDs <- as.character(count[n,1])
    if (IDs == i) {
      uORF_stop <- as.integer(count[n,3])
      next_uORF_stop <- as.integer(count[n+1,3])
      uORF_start <- as.integer(count[n,2])
      count_reads <- as.integer(count[n,5])
      length_uORF <- uORF_stop - uORF_start
      ratio <- (count_reads/length_uORF)
      if (uORF_stop == next_uORF_stop){
        List <- c(List, ratio)
        Name <- c(Name, uORF_start)
      }else{
        List <- c(List, ratio)
        Name <- c(Name, uORF_start)
        names(List) <- Name
#        List <- c(List, 'uORF_start' = ratio)
        max_ratio <- max(unlist(List))
        p = match(max_ratio,List)
        best_uORF_start = names(List)[p]
        new <- c(IDs, best_uORF_start, uORF_stop)
        best_uORF <- rbind(best_uORF, new)
        List = list()
        Name = list()
      }
    }
  }
}
write.table(best_uORF, 'best_uORFs.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
best_uORF <- read.table('best_uORFs.txt')
colnames(best_uORF) <- c('IDs', 'start_uORF', 'stop_uORF')
uORFs <- merge(best_uORF,count, by = c('IDs','start_uORF','stop_uORF'))
write.csv(uORFs, 'Ctrl_uORF_counts.csv', row.names = FALSE)
uORFs <- uORFs[,c(1:3,5)]
uORFs$length_uORF <- as.integer(uORFs$stop_uORF) - as.integer(uORFs$start_uORF)
uORFs$ratio <- uORFs$count/uORFs$length_uORF
#threshold : ratio min 0.333 (at least 1 read every 3nt)
translated_uORFs <- uORFs[uORFs$ratio >= (1/3),]
write.csv(translated_uORFs,'transled_uORFs_count.csv')
translated_uORFs$Gene_ID <- substring(translated_uORFs[,1], 1, 18)
translated_uORFs$Transcript_ID <- substring(translated_uORFs[,1], 20, 37)
Gene_ID_selected_transcripts <- read.table('GeneID_transcriptID_single_isoform_sameCDS_longest.txt')
colnames(Gene_ID_selected_transcripts) <- c('Gene_ID', 'Transcript_ID')
non_translated_uORF_transcripts <- anti_join(Gene_ID_selected_transcripts, translated_uORFs, by = "Gene_ID") # no translated uORF on 5' UTR
untranslated_uORFs <- uORFs[uORFs$ratio < (1/3),]
colnames(untranslated_uORFs) <- c('Gene_ID', 'start_uORF', 'stop_uORF','count','length','ratio count/length')
no_uORF_transcripts <- anti_join(non_translated_uORF_transcripts, untranslated_uORFs, by = "Gene_ID") # no uORF on 5' UTR
write.table(untranslated_uORFs, 'RP_ctrl.untranslated_uORF.longest.count3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE) # non translated uORFs
write.table(translated_uORFs, 'RP_ctrl.translated_uORF.longest.count3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE) # translated uORFs
write.table(non_translated_uORF_transcripts,'Transcripts_no_translated_uORF.txt', quote = FALSE, row.names = FALSE, col.names = FALSE) # no translated uORF on 5' UTR
write.table(no_uORF_transcripts,'Transcripts_no_uORF.txt', quote = FALSE, row.names = FALSE, col.names = FALSE) # no uORF on 5' UTR

#translated uORFs vs untranslated uORFs
translated_uORFs <- read.table('RP_ctrl.translated_uORF.longest.count3.txt')
Gene_ID_translated <- substring(translated_uORFs[,1], 1, 18)
write.table(Gene_ID_translated, 'RP_ctrl.translated_uORF.longest.3.Gene_ID.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
untranslated_uORFs <- read.table('RP_ctrl.untranslated_uORF.longest.count3.txt')
Gene_ID_untranslated <- substring(untranslated_uORFs[,1], 1, 18)
write.table(Gene_ID_untranslated, 'RP_ctrl.untranslated_uORF.longest.3.Gene_ID.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
untranslated_uORF <- read.table('RP_ctrl.untranslated_uORF.longest.3.Gene_ID.txt')
colnames(untranslated_uORF) <- 'Gene_ID'
translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.3.Gene_ID.txt')
colnames(translated_uORF) <- 'Gene_ID'
transcripts_no_uORF <- read.table('Transcripts_no_translated_uORF.txt')
colnames(transcripts_no_uORF) <- 'Gene_ID'
Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv')
Denrsh2_neg_translated_uORF <- merge(translated_uORF,Denrsh2_neg, by = 'Gene_ID')
Denrsh2_neg_untranslated_uORF <- merge(untranslated_uORF,Denrsh2_neg, by = 'Gene_ID')
Denrsh2_neg_transcripts_no_uORF <- merge(transcripts_no_uORF,Denrsh2_neg, by = 'Gene_ID')
Denrsh2_neg_number_uORFs <- table(Denrsh2_neg_translated_uORF$Gene_ID)
more_than_one <- which(Denrsh2_neg_number_uORFs > 1)
one <- which(Denrsh2_neg_number_uORFs == 1)
fisher <- read.csv('fisher_test_Denrsh2.csv')
fisher.test(fisher)

#Transcript with only untranslated uORF (transcript with translated uORF if only 1 of the uORF is translated)
translated_uORFs <- read.csv('RP_ctrl.translated_uORF.longest.3.Gene_ID.txt')
translated_uORFs <- unique(translated_uORFs)
write.table(translated_uORFs, 'RP_ctrl.translated_uORF.longest.3.Gene_ID.unique.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
colnames(translated_uORFs) <- 'Gene_ID'
untranslated_uORF <- read.csv('RP_ctrl.untranslated_uORF.longest.3.Gene_ID.txt')
colnames(untranslated_uORF) <- 'Gene_ID'
transcripts_without_translated_uORF <- unique(untranslated_uORF)
write.table(transcripts_without_translated_uORF, 'RP_ctrl.untranslated_uORF.longest.3.Gene_ID.unique.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
transcripts_without_translated_uORF <- anti_join(transcripts_without_translated_uORF, translated_uORFs, by = "Gene_ID")
write.table(transcripts_without_translated_uORF, 'RP_ctrl.untranslated_uORF_transcripts.longest.3.Gene_ID.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

#Transcript without translated uORF vs transcript with translated uORF
transcripts_without_translated_uORF <- read.table('RP_ctrl.untranslated_uORF_transcripts.longest.3.Gene_ID.txt')
colnames(transcripts_without_translated_uORF) <- 'Gene_ID'
translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.3.Gene_ID.txt')
transcripts_with_translated_uORF <- unique(translated_uORF)
colnames(transcripts_with_translated_uORF) <- 'Gene_ID'
write.table(transcripts_with_translated_uORF,'RP_ctrl.translated_uORF_transcripts.longest.3.Gene_ID.csv')
transcripts_no_uORF <- read.table('Transcripts_no_translated_uORF.txt') #no uORF on 5' UTR
colnames(transcripts_no_uORF) <- c('Gene_ID','Transcript_ID')
Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Eif2dsh_all_positive_TE_FC_FDR_0.1.csv')
Denrsh2_neg_transcripts_uORF <- merge(transcripts_with_translated_uORF,Denrsh2_neg, by = 'Gene_ID')
Denrsh2_neg_unstranslated_uORF <- merge(transcripts_without_translated_uORF,Denrsh2_neg, by = 'Gene_ID')
Denrsh2_neg_transcripts_no_uORF <- merge(transcripts_no_uORF,Denrsh2_neg, by = 'Gene_ID')

#All expressed transcripts
#Transcript without translated uORF vs transcript with translated uORF
all <- read.table('GeneID_transcriptID_single_isoform_sameCDS_longest.txt')
all <- all[,1]
write.table(all, 'GeneID_single_isoform_sameCDS_longest.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
all <- read.table('GeneID_single_isoform_sameCDS_longest.txt')
colnames(all) <- 'Gene_ID'
all_transcripts_uORF <- merge(transcripts_with_translated_uORF,all, by = 'Gene_ID')
all_transcripts_no_uORF <- merge(transcripts_no_uORF,all, by = 'Gene_ID')
all_transcripts_untranslated_uORF <- merge(transcripts_without_translated_uORF,all, by = 'Gene_ID')

#fisher test
fisher <- read.csv('fisher_test_Denrsh2.csv')
fisher.test(fisher)

#translated uORFs vs untranslated uORFs
all <- read.table('GeneID_single_isoform_sameCDS_longest.txt')
colnames(all) <- 'Gene_ID'
untranslated_uORF <- read.table('RP_ctrl.untranslated_uORF.longest.Gene_ID.txt')
colnames(untranslated_uORF) <- 'Gene_ID'
translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.min5.Gene_ID.txt')
colnames(translated_uORF) <- 'Gene_ID'
all_translated_uORF <- merge(translated_uORF,all, by = 'Gene_ID')
all_untranslated_uORF <- merge(untranslated_uORF,all, by = 'Gene_ID')

#transcript decreased TE
transcripts_without_translated_uORF <- read.table('RP_ctrl.untranslated_uORF_transcripts.longest.Gene_ID.txt')
colnames(transcripts_without_translated_uORF) <- 'Gene_ID'
translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.min5.Gene_ID.txt')
transcripts_with_translated_uORF <- unique(translated_uORF)
colnames(transcripts_with_translated_uORF) <- 'Gene_ID'
Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Eif2dsh3_all_neg_TE_FC.csv')
Denrsh2_neg_lowest_TE_FC <- Denrsh2_neg[Denrsh2_neg$log2FoldChange < -0.5,]
Denrsh2_neg_transcripts_uORF <- merge(transcripts_with_translated_uORF,Denrsh2_neg_lowest_TE_FC, by = 'Gene_ID')
Denrsh2_neg_transcripts_no_uORF <- merge(transcripts_without_translated_uORF,Denrsh2_neg_lowest_TE_FC, by = 'Gene_ID')
#fisher test
fisher <- read.csv('fisher_test_Denrsh2_lowest_TE_FC.csv')
fisher.test(fisher)

##uORFs Start codon

translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.count3.txt')
translated_uORF$Gene_ID <- substring(translated_uORF[,1], 1, 18) # take GeneID + transcript ID
translated_uORF <- translated_uORF[,c(7,1:4)]
colnames(translated_uORF) <- c('Gene_ID','IDs','start_uORF','stop_uORF', 'reads_counts')
untranslated_uORF <- read.table('RP_ctrl.untranslated_uORF.longest.count3.txt')
untranslated_uORF$Gene_ID <- substring(untranslated_uORF[,1], 1, 18)
untranslated_uORF <- untranslated_uORF[,c(7,1:4)]
colnames(untranslated_uORF) <- c('Gene_ID','IDs','start_uORF','stop_uORF', 'reads_counts')
Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv')
Denrsh2_neg_translated_uORFs <- merge(Denrsh2_neg,translated_uORF, by = 'Gene_ID')
Denrsh2_neg_translated_uORFs <- Denrsh2_neg_translated_uORFs[,c(9:12)]
write.table(Denrsh2_neg_translated_uORFs,'Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
Denrsh2_neg_untranslated_uORFs <- merge(untranslated_uORF, Denrsh2_neg, by = 'Gene_ID')
Denrsh2_neg_untranslated_uORFs <- Denrsh2_neg_untranslated_uORFs[,c(2:5)]
write.table(Denrsh2_neg_untranslated_uORFs,'Denrsh2_neg_TE_FC_FDR_0.1_untranslated_uORFs.3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
Eif2dsh_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Eif2dsh_all_negative_TE_FC_FDR_0.1.csv')
Eif2dsh_neg_translated_uORFs <- merge(Eif2dsh_neg,translated_uORF, by = 'Gene_ID')
Eif2dsh_neg_translated_uORFs <- Eif2dsh_neg_translated_uORFs[,c(9:12)]
write.table(Eif2dsh_neg_translated_uORFs,'Eif2dsh_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
Eif2dsh_neg_untranslated_uORFs <- merge(untranslated_uORF,Eif2dsh_neg, by = 'Gene_ID')
Eif2dsh_neg_untranslated_uORFs <- Eif2dsh_neg_untranslated_uORFs[,c(2:5)]
write.table(Eif2dsh_neg_untranslated_uORFs,'Eif2dsh_neg_TE_FC_FDR_0.1_untranslated_uORFs.3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
all <- read.table('GeneID_single_isoform_sameCDS_longest.txt')
colnames(all) <-'Gene_ID'
all_translated_uORFs <- merge(all,translated_uORF, by = 'Gene_ID')
all_translated_uORFs <- all_translated_uORFs[,c(2:5)]
write.table(all_translated_uORFs,'all_translated_uORFs.3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
all_untranslated_uORFs <- merge(untranslated_uORF,all, by = 'Gene_ID')
all_untranslated_uORFs <- all_untranslated_uORFs[,c(1,3,4)]
write.table(all_untranslated_uORFs,'all_untranslated_uORFs.3.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

#Add uORF sequence with merge_sequence.py then count proportion of each start codon with start_site_count.py
#Fisher test
fisher <- read.csv('fisher_test_Denrsh2_start_codon.csv')
fisher.test(fisher)

##number of uORFs

all <- read.table('GeneID_single_isoform_sameCDS_longest.txt')
colnames(all) <- 'Gene_ID'
translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.3.Gene_ID.txt')
colnames(translated_uORF) <- 'Gene_ID'
all_translated_uORF <- merge(translated_uORF,all, by = 'Gene_ID')
all_number_uORFs_gene <- table(all_translated_uORF$Gene_ID)
all_number_uORFs <- data.frame(c(1:28),c(sum(all_number_uORFs_gene == 1),sum(all_number_uORFs_gene == 2),sum(all_number_uORFs_gene == 3),sum(all_number_uORFs_gene == 4),sum(all_number_uORFs_gene == 5),sum(all_number_uORFs_gene == 6),sum(all_number_uORFs_gene == 7),sum(all_number_uORFs_gene == 8),sum(all_number_uORFs_gene == 9),sum(all_number_uORFs_gene == 10),sum(all_number_uORFs_gene == 11),sum(all_number_uORFs_gene == 12),sum(all_number_uORFs_gene == 13),sum(all_number_uORFs_gene == 14),sum(all_number_uORFs_gene == 15),sum(all_number_uORFs_gene == 16),sum(all_number_uORFs_gene == 17),sum(all_number_uORFs_gene == 18),sum(all_number_uORFs_gene == 19),sum(all_number_uORFs_gene == 20),sum(all_number_uORFs_gene == 21),sum(all_number_uORFs_gene == 22),sum(all_number_uORFs_gene == 23),sum(all_number_uORFs_gene == 24),sum(all_number_uORFs_gene == 25),sum(all_number_uORFs_gene == 26),sum(all_number_uORFs_gene == 27),sum(all_number_uORFs_gene == 28)))
write.csv(all_number_uORFs, 'Nber_uORFs_per_gene_translated_uORFs.csv', row.names = FALSE)
plot(all_number_uORFs$c.1.28., all_number_uORFs$c.sum.all_number_uORFs_gene....1...sum.all_number_uORFs_gene....)
mean(all_number_uORFs_gene)
median(all_number_uORFs_gene)

Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv')
Denrsh2_neg_translated_uORF <- merge(all_translated_uORF,Denrsh2_neg, by = 'Gene_ID')
Denr_number_uORFs_gene <- table(Denrsh2_neg_translated_uORF$Gene_ID)
Denr_number_uORFs <- data.frame(c(1:28),c(sum(Denr_number_uORFs_gene == 1),sum(Denr_number_uORFs_gene == 2),sum(Denr_number_uORFs_gene == 3),sum(Denr_number_uORFs_gene == 4),sum(Denr_number_uORFs_gene == 5),sum(Denr_number_uORFs_gene == 6),sum(Denr_number_uORFs_gene == 7),sum(Denr_number_uORFs_gene == 8),sum(Denr_number_uORFs_gene == 9),sum(Denr_number_uORFs_gene == 10),sum(Denr_number_uORFs_gene == 11),sum(Denr_number_uORFs_gene == 12),sum(Denr_number_uORFs_gene == 13),sum(Denr_number_uORFs_gene == 14),sum(Denr_number_uORFs_gene == 15),sum(Denr_number_uORFs_gene == 16),sum(Denr_number_uORFs_gene == 17),sum(Denr_number_uORFs_gene == 18),sum(Denr_number_uORFs_gene == 19),sum(Denr_number_uORFs_gene == 20),sum(Denr_number_uORFs_gene == 21),sum(Denr_number_uORFs_gene == 22),sum(Denr_number_uORFs_gene == 23),sum(Denr_number_uORFs_gene == 24),sum(Denr_number_uORFs_gene == 25),sum(Denr_number_uORFs_gene == 26),sum(Denr_number_uORFs_gene == 27),sum(Denr_number_uORFs_gene == 28)))
write.csv(Denr_number_uORFs, 'Nber_uORFs_per_DENR_responsive_gene_translated_uORFs.csv', row.names = FALSE)
plot(Denr_number_uORFs$c.1.28., Denr_number_uORFs$c.sum.Denr_number_uORFs_gene....1...sum.Denr_number_uORFs_gene....)
mean(Denr_number_uORFs_gene)
median(Denr_number_uORFs_gene)

##length of uORFs

translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.count3.txt')
colnames(translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads','uORF_length','ratio')
translated_uORF$Gene_ID <- substring(translated_uORF[,1], 1, 18)
Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv')
Denrsh2_neg_translated_uORF <- merge(translated_uORF,Denrsh2_neg, by = 'Gene_ID')
write.csv(Denrsh2_neg_translated_uORF, 'Denrsh2_neg_uORF_length.csv', row.names = FALSE)
median(Denrsh2_neg_translated_uORF$uORF_length)
mean(Denrsh2_neg_translated_uORF$uORF_length)
Denrsh2_neg_length_uORFs <- table(Denrsh2_neg_translated_uORF$uORF_length)
length_uORFs <- as.data.frame(Denrsh2_neg_length_uORFs)
length_uORFs$Group <- 'Denr_sh2' 
write.table(length_uORFs,'length_uORFs_Denrsh2_all_negative_TE_FC_FDR_0.1.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

Eif2dsh_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Eif2dsh_all_negative_TE_FC_FDR_0.1.csv')
Eif2dsh_neg_translated_uORF <- merge(translated_uORF,Eif2dsh_neg, by = 'Gene_ID')
write.csv(Eif2dsh_neg_translated_uORF, 'Eif2dsh_neg_uORF_length.csv', row.names = FALSE)
median(Eif2dsh_neg_translated_uORF$uORF_length)
mean(Eif2dsh_neg_translated_uORF$uORF_length)
Eif2dsh_neg_length_uORFs <- table(Eif2dsh_neg_translated_uORF$uORF_length)
length_uORFs <- as.data.frame(Eif2dsh_neg_length_uORFs)
write.table(length_uORFs,'length_uORFs_Eif2dsh_all_negative_TE_FC_FDR_0.1.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

all <- read.table('GeneID_single_isoform_sameCDS_longest.txt')
colnames(all) <- 'Gene_ID'
translated_uORF <- read.table('RP_ctrl.translated_uORF.longest.count3.txt')
colnames(translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads','uORF_length','ratio')
translated_uORF$Gene_ID <- substring(translated_uORF[,1], 1, 18)
all_translated_uORF <- merge(translated_uORF,all, by = 'Gene_ID')
write.csv(all_translated_uORF, 'all_translated_uORF_length.csv', row.names = FALSE)
median(all_translated_uORF$uORF_length)
mean(all_translated_uORF$uORF_length)
all_length_uORFs <- table(all_translated_uORF$uORF_length)
length_uORFs_all <- as.data.frame(all_length_uORFs)
length_uORFs_all$Group <- 'all'
write.table(length_uORFs_all,'length_uORFs_all.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
length_uORFs_all$Group <- 'all'
length_uORFs_all_Dsh2 <- rbind(length_uORFs,length_uORFs_all)
write.table(length_uORFs_all_Dsh2,'length_uORFs_all_Denrsh2.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

length_uORFs <- read.csv('Wilcoxon_test_Denrsh2_length_uORFs.csv')
wilcox.test(length_uORFs$Denrsh2_exclusive, length_uORFs$all_transcripts, paired=TRUE) 

data <- read.table('length_uORFs_all_Denrsh2.txt')
colnames(data) <- c('Length','Value','Group')
data$Length_aa <- data$Length/3
pdf('uORF_length_distribution_all_vs_Dsh2.pdf', height=5,width=12)
data %>%
  ggplot(aes(x=Length_aa, y=Value)) +
  geom_bar(aes(fill=Group),stat="identity", position=position_dodge(), color = 'grey32') + 
  scale_fill_manual(values=c( 'grey42', '#00AF6696')) +
  theme_classic() +
  ylab('Frequency') +
  xlab('uORF length (aa)') + 
  theme(legend.text = element_text(size=16), legend.key.width= unit(1, 'cm'),legend.title=element_blank(),axis.line.x = element_line(size = 1),axis.line.y = element_line(size = 1),axis.text.x = element_text(size = 22),axis.text.y = element_text(size = 22),axis.title.y = element_text(size = 20),,axis.title.x = element_text(size = 20))
dev.off()

library(KSgeneral)
disc_ks_test(length_uORFs$Denrsh2_exclusive,length_uORFs$all_transcripts)
ks.test(all_translated_uORF$uORF_length, Eif2dsh_neg_translated_uORF$uORF_length)
ks.test(all_translated_uORF$uORF_length, Denrsh2_neg_translated_uORF$uORF_length)

##5UTR length

UTR <- read.table('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/annotations/Mus_musculus.GRCm38.100.5UTR.PC.bed')
UTR <- UTR[order(UTR$V1, UTR$V3),]
UTR_longest <- as.data.frame(c())
for (n in 1:(nrow(UTR)-1)) {
  geneID <- substring(UTR[n,1], 1, 18)
  transcriptID <- substring(UTR[n,1], 20, 37)
  lengthUTR <- as.integer(UTR[n,3])
  NextGeneID <- substring(UTR[(n+1),1], 1, 18)
  if (geneID != NextGeneID){
    new <- c(geneID,transcriptID,lengthUTR)
    UTR_longest <- rbind(UTR_longest, new)
  }
}
write.table(UTR_longest,'Mus_musculus.GRCm38.100.5UTR.longest.PC.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
UTR_longest <- read.table('Mus_musculus.GRCm38.100.5UTR.longest.PC.txt')
UTR_length <- UTR_longest[,c(1,3)]
colnames(UTR_length) <- c('Gene_ID','UTR_length')
Eif2dsh_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Eif2dsh_all_negative_TE_FC_FDR_0.1.csv', row.names = NULL)
Eif2dsh_UTR_length <- merge(Eif2dsh_neg, UTR_length, by = 'Gene_ID')
Eif2dsh_UTR_length$UTR_length <- as.integer(Eif2dsh_UTR_length$UTR_length)
mean(Eif2dsh_UTR_length$UTR_length)
median(Eif2dsh_UTR_length$UTR_length)
Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv', row.names = NULL)
Denrsh2_UTR_length <- merge(Denrsh2_neg, UTR_length, by = 'Gene_ID', all.x = TRUE)
Denrsh2_UTR_length$UTR_length <- as.integer(Denrsh2_UTR_length$UTR_length)
Denrsh2_UTR_length <- na.omit(Denrsh2_UTR_length)
mean(Denrsh2_UTR_length$UTR_length)
median(Denrsh2_UTR_length$UTR_length)

all <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/uORF_annotations/GeneID_single_isoform_sameCDS_longest.txt', row.names = NULL, col.names = 'Gene_ID')
all_length <- merge(all, UTR_length, by = 'Gene_ID')
all_length$UTR_length <- as.integer(all_length$UTR_length)
mean(all_length$UTR_length)
median(all_length$UTR_length)

#Kolmogorov-Smirnov test
ks.test(Eif2dsh_UTR_length$UTR_length, all_length$UTR_length)
ks.test(Denrsh2_UTR_length$UTR_length, all_length$UTR_length)

translated_uORFs <- read.table('all_translated_uORFs.3.txt')
translated_uORFs$Gene_ID <- substring(translated_uORFs[,1], 1, 18)
translated_uORFs <- unique(translated_uORFs[,c(1,5)])
translated_uORFs_5UTR_length <- merge(translated_uORFs, all_length, by= "Gene_ID")
translated_uORFs_5UTR_length$UTR_length <- as.integer(translated_uORFs_5UTR_length$UTR_length)
mean(translated_uORFs_5UTR_length$UTR_length)
median(translated_uORFs_5UTR_length$UTR_length)

Denrsh2_neg_translated_uORF_5UTR_length <- merge(Denrsh2_neg, translated_uORFs_5UTR_length, by = 'Gene_ID')
Denrsh2_neg_translated_uORF_5UTR_length$UTR_length <- as.integer(Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)
mean(Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)
median(Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)

non_Denrsh2_neg_translated_uORF_5UTR_length <- anti_join(translated_uORFs_5UTR_length, Denrsh2_neg_translated_uORF_5UTR_length,by = 'Gene_ID')
mean(non_Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)
median(non_Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)

#Kolmogorov-Smirnov test
ks.test(all_length$UTR_length, non_Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)
ks.test(all_length$UTR_length, Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)
ks.test(all_length$UTR_length, Denrsh2_UTR_length$UTR_length)
ks.test(non_Denrsh2_neg_translated_uORF_5UTR_length$UTR_length, Denrsh2_UTR_length$UTR_length)
ks.test(non_Denrsh2_neg_translated_uORF_5UTR_length$UTR_length, Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)
ks.test(Denrsh2_UTR_length$UTR_length, non_Denrsh2_neg_translated_uORF_5UTR_length$UTR_length)

#plots
UTR_length$Type <- 'all_uORFs'
Denrsh2_UTR_length$Type <- 'Denrsh2'
Denrsh2_UTR_length <- Denrsh2_UTR_length[,c(1,9,10)]
Eif2dsh_UTR_length$Type <- 'Eif2dsh'
Eif2dsh_UTR_length <- Eif2dsh_UTR_length[,c(1,9,10)]
all_uORFs <- rbind(UTR_length, Denrsh2_UTR_length, Eif2dsh_UTR_length)

pdf('5UTR_length_all_negative_FDR_0.1_boxplot_dotplot.pdf', height=6,width=7)
ggplot(all_uORFs,aes(x=Type, y=UTR_length)) +
  geom_jitter(aes(x=Type, y=UTR_length, fill = factor(Type)),color = 'grey32', size=1, alpha=0.6, width = 0.2, stroke = 1, shape = 20) +
  geom_boxplot(aes( fill = factor(Type)),outlier.shape=20, outlier.size=1) +
  scale_y_continuous(limits = c(0,2500)) +
  scale_fill_manual(values = c('grey32','#00AF6696','#00B2EE96')) +
  theme_classic()
dev.off()

pdf('5UTR_length_all_negative_FDR_0.1_violinplot_dotplot.pdf', height=5,width=8)
ggplot(all_uORFs,aes(x=Type, y=UTR_length)) +
  geom_jitter(aes(x=Type, y=UTR_length, fill = factor(Type)), color = 'grey32', size=0.5, alpha=0.6, width = 0.3, stroke = 1, shape = 20) +
  geom_violin(aes( fill = factor(Type)),outlier.shape=20, outlier.size=1,alpha = 0.8) +
  scale_y_continuous(limits = c(0,2500)) +
  scale_fill_manual(values = c('grey32','#00AF6696','#00B2EE96')) +
  theme_classic()
dev.off()

#DENR-responsive vs non responsive translated uORF containing transcripts

all_length$Type <- 'all expressed transcripts'
all_length <- all_length[,c(1,2,4)]
Denrsh2_neg_translated_uORF_5UTR_length$Type <- 'Denr-responsive'
Denrsh2_neg_translated_uORF_5UTR_length <- Denrsh2_neg_translated_uORF_5UTR_length[,c(1,10,11)]
non_Denrsh2_neg_translated_uORF_5UTR_length$Type <- 'other translated uORFs'
non_Denrsh2_neg_translated_uORF_5UTR_length <- non_Denrsh2_neg_translated_uORF_5UTR_length[,c(1,3,4)]
all_uORFs <- rbind(all_length, non_Denrsh2_neg_translated_uORF_5UTR_length, Denrsh2_neg_translated_uORF_5UTR_length)

pdf('5UTR_length_DENR_translated_uORFs_FDR_0.1_violinplot_dotplot.pdf', height=5,width=8)
ggplot(all_uORFs,aes(x=factor(Type, levels = c('all expressed transcripts','other translated uORFs','Denr-responsive')), y=UTR_length)) +
  geom_jitter(aes(x=Type, y=UTR_length, fill = factor(Type)), color = 'grey32', size=0.5, alpha=0.6, width = 0.3, stroke = 1, shape = 20) +
  geom_violin(aes( fill = factor(Type)),outlier.shape=20, outlier.size=1,alpha = 0.8) +
  scale_y_continuous(limits = c(0,2500)) +
  scale_fill_manual(values = c('grey32','#00AF6696','#c883c8ff')) +
  theme_classic()
dev.off()

Denrsh2_UTR_length <- Denrsh2_UTR_length[,c(1:3)]
all_uORFs <- rbind(all_length, non_Denrsh2_neg_translated_uORF_5UTR_length, Denrsh2_UTR_length)

pdf('5UTR_length_DENR_&_translated_uORFs_FDR_0.1_violinplot_dotplot.pdf', height=5,width=8)
ggplot(all_uORFs,aes(x=factor(Type, levels = c('all expressed transcripts','other translated uORFs','Denrsh2')), y=UTR_length)) +
  geom_jitter(aes(x=Type, y=UTR_length, fill = factor(Type)), color = 'grey32', size=0.5, alpha=0.6, width = 0.3, stroke = 1, shape = 20) +
  geom_violin(aes( fill = factor(Type)),outlier.shape=20, outlier.size=1,alpha = 0.8) +
  scale_y_continuous(limits = c(0,2500)) +
  scale_fill_manual(values = c('grey32','#00AF6696','#c883c8ff')) +
  theme_classic()
dev.off()

##translated uORFs GC content

seq <- read.table('GeneID_transcriptID_5UTR_single_isoform_sameCDS_longest_uORFs_seq.txt')
colnames(seq) <- c('Gene_ID','uORF_start','uORF_stop','sequence')
seq$sequence <- substring(seq$sequence, 1, nchar(seq$sequence)-3)

Denrsh2_neg_translated_uORF <- read.table('Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Denrsh2_neg_translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads')
Denrsh2_neg_translated_uORF$Gene_ID <- substring(Denrsh2_neg_translated_uORF[,1], 1, 18)
Denrsh2_neg_translated_uORF$uORF_start <- Denrsh2_neg_translated_uORF$uORF_start + 16
Denrsh2_neg_translated_uORF$uORF_stop <- Denrsh2_neg_translated_uORF$uORF_stop + 16
Denrsh2_neg_translated_uORF_seq <- merge(seq,Denrsh2_neg_translated_uORF, by = c('Gene_ID','uORF_start','uORF_stop'))
Denrsh2_neg_translated_uORF_seq$G_content <- str_count(Denrsh2_neg_translated_uORF_seq$sequence, "G")
Denrsh2_neg_translated_uORF_seq$C_content <- str_count(Denrsh2_neg_translated_uORF_seq$sequence, "C")
Denrsh2_neg_translated_uORF_seq$A_content <- str_count(Denrsh2_neg_translated_uORF_seq$sequence, "A")
Denrsh2_neg_translated_uORF_seq$T_content <- str_count(Denrsh2_neg_translated_uORF_seq$sequence, "T")
Denrsh2_neg_translated_uORF_seq$length <- nchar(Denrsh2_neg_translated_uORF_seq$sequence)
Denrsh2_neg_translated_uORF_seq$GC_content <- (Denrsh2_neg_translated_uORF_seq[,6]+Denrsh2_neg_translated_uORF_seq[,7])/Denrsh2_neg_translated_uORF_seq$length
hist(Denrsh2_neg_translated_uORF_seq$GC_content,breaks = 1500)
mean(Denrsh2_neg_translated_uORF_seq$GC_content)
median(Denrsh2_neg_translated_uORF_seq$GC_content)

Eif2dsh_neg_translated_uORF <- read.table('Eif2dsh_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Eif2dsh_neg_translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads')
Eif2dsh_neg_translated_uORF$Gene_ID <- substring(Eif2dsh_neg_translated_uORF[,1], 1, 18)
Eif2dsh_neg_translated_uORF$uORF_start <- Eif2dsh_neg_translated_uORF$uORF_start + 16
Eif2dsh_neg_translated_uORF$uORF_stop <- Eif2dsh_neg_translated_uORF$uORF_stop + 16
Eif2dsh_neg_translated_uORF_seq <- merge(seq,Eif2dsh_neg_translated_uORF, by = c('Gene_ID','uORF_start','uORF_stop'))
Eif2dsh_neg_translated_uORF_seq$G_content <- str_count(Eif2dsh_neg_translated_uORF_seq$sequence, "G")
Eif2dsh_neg_translated_uORF_seq$C_content <- str_count(Eif2dsh_neg_translated_uORF_seq$sequence, "C")
Eif2dsh_neg_translated_uORF_seq$A_content <- str_count(Eif2dsh_neg_translated_uORF_seq$sequence, "A")
Eif2dsh_neg_translated_uORF_seq$T_content <- str_count(Eif2dsh_neg_translated_uORF_seq$sequence, "T")
Eif2dsh_neg_translated_uORF_seq$length <- nchar(Eif2dsh_neg_translated_uORF_seq$sequence)
Eif2dsh_neg_translated_uORF_seq$GC_content <- (Eif2dsh_neg_translated_uORF_seq[,6]+Eif2dsh_neg_translated_uORF_seq[,7])/Eif2dsh_neg_translated_uORF_seq$length
hist(Eif2dsh_neg_translated_uORF_seq$GC_content,breaks = 1500)
mean(Eif2dsh_neg_translated_uORF_seq$GC_content)
median(Eif2dsh_neg_translated_uORF_seq$GC_content)

all_translated_uORF <- read.table('all_translated_uORFs.3.txt')
colnames(all_translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads')
all_translated_uORF$Gene_ID <- substring(all_translated_uORF[,1], 1, 18)
all_translated_uORF$uORF_start <- all_translated_uORF$uORF_start + 16
all_translated_uORF$uORF_stop <- all_translated_uORF$uORF_stop + 16
all_translated_uORF_seq <- merge(seq,all_translated_uORF, by = c('Gene_ID','uORF_start','uORF_stop'))
all_translated_uORF_seq$G_content <- str_count(all_translated_uORF_seq$sequence, "G")
all_translated_uORF_seq$C_content <- str_count(all_translated_uORF_seq$sequence, "C")
all_translated_uORF_seq$A_content <- str_count(all_translated_uORF_seq$sequence, "A")
all_translated_uORF_seq$T_content <- str_count(all_translated_uORF_seq$sequence, "T")
all_translated_uORF_seq$length <- nchar(all_translated_uORF_seq$sequence)
all_translated_uORF_seq$GC_content <- (all_translated_uORF_seq[,6]+all_translated_uORF_seq[,7])/all_translated_uORF_seq$length
hist(all_translated_uORF_seq$GC_content,breaks = 1500)
mean(all_translated_uORF_seq$GC_content)
median(all_translated_uORF_seq$GC_content)

GC_content_translated_uORFs <- as.matrix(all_translated_uORF_seq$GC_content,Denrsh2_neg_translated_uORF_seq$GC_content,Eif2dsh_neg_translated_uORF_seq$GC_content)
GC_content_translated_uORFs <- as.data.frame(GC_content_translated_uORFs)

#Kolmogorov-Smirnov test
ks.test(Eif2dsh_neg_translated_uORF_seq$GC_content, all_translated_uORF_seq$GC_content)
#box plot & violin plot
pdf('GC_content_all_negative_FDR_0.1_boxplot.pdf')
boxplot(all_translated_uORF_seq$GC_content,Denrsh2_neg_translated_uORF_seq$GC_content, Eif2dsh_neg_translated_uORF_seq$GC_content,col=c('lightgray','lightpink2','lightskyblue3'), names = c("all transcripts","Denr shRNA","Eif2d shRNA"))
dev.off()
pdf('GC_content_all_negative_FDR_0.1_violinplot.pdf')
vioplot(all_translated_uORF_seq$GC_content,Denrsh2_neg_translated_uORF_seq$GC_content, Eif2dsh_neg_translated_uORF_seq$GC_content,col=c('lightgray','lightpink2','lightskyblue3'), names = c("all transcripts","Denr shRNA","Eif2d shRNA"))
dev.off()

##5UTR GC content 

seq <- read.delim('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/annotations/Mus_musculus.GRCm38.100.5UTR.longest.seq.PC.txt', col.names = c('Gene_ID','sequence'))
seq$Gene_ID <- substring(seq[,1], 1, 18)
all_expressed <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/uORF_annotations/GeneID_transcriptID_selected_isoform.csv')
all_expressed  <- merge(seq,all_expressed, by ='Gene_ID')
all_expressed$G_content <- str_count(all_expressed$sequence, "G")
all_expressed$C_content <- str_count(all_expressed$sequence, "C")
all_expressed$A_content <- str_count(all_expressed$sequence, "A")
all_expressed$T_content <- str_count(all_expressed$sequence, "T")
all_expressed$length <- nchar(all_expressed$sequence)
all_expressed$GC_content <- (all_expressed[,4]+all_expressed[,5])/all_expressed$length
hist(all_expressed$GC_content,breaks = 1500)
mean(all_expressed$GC_content)
median(all_expressed$GC_content)
all_expressed$A_content <- (all_expressed[,6])/all_expressed$length
mean(all_expressed$A_content)
median(all_expressed$A_content)

Denrsh2_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Denrsh2_all_negative_TE_FC_FDR_0.1.csv')
Denrsh2_neg  <- merge(seq,Denrsh2_neg, by ='Gene_ID')
Denrsh2_neg$G_content <- str_count(Denrsh2_neg$sequence, "G")
Denrsh2_neg$C_content <- str_count(Denrsh2_neg$sequence, "C")
Denrsh2_neg$A_content <- str_count(Denrsh2_neg$sequence, "A")
Denrsh2_neg$T_content <- str_count(Denrsh2_neg$sequence, "T")
Denrsh2_neg$length <- nchar(Denrsh2_neg$sequence)
Denrsh2_neg$GC_content <- (Denrsh2_neg[,10]+Denrsh2_neg[,11])/Denrsh2_neg$length
mean(Denrsh2_neg$GC_content)
median(Denrsh2_neg$GC_content)
Denrsh2_neg$A_content <- (Denrsh2_neg[,12])/Denrsh2_neg$length
mean(Denrsh2_neg$A_content)
median(Denrsh2_neg$A_content)

Eif2dsh_neg <- read.csv('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/deltaTE/Eif2dsh_all_negative_TE_FC_FDR_0.1.csv')
Eif2dsh_neg  <- merge(seq,Eif2dsh_neg, by ='Gene_ID')
Eif2dsh_neg$G_content <- str_count(Eif2dsh_neg$sequence, "G")
Eif2dsh_neg$C_content <- str_count(Eif2dsh_neg$sequence, "C")
Eif2dsh_neg$A_content <- str_count(Eif2dsh_neg$sequence, "A")
Eif2dsh_neg$T_content <- str_count(Eif2dsh_neg$sequence, "T")
Eif2dsh_neg$length <- nchar(Eif2dsh_neg$sequence)
Eif2dsh_neg$GC_content <- (Eif2dsh_neg[,10]+Eif2dsh_neg[,11])/Eif2dsh_neg$length
mean(Eif2dsh_neg$GC_content)
median(Eif2dsh_neg$GC_content)
Eif2dsh_neg$A_content <- (Eif2dsh_neg[,12])/Eif2dsh_neg$length
mean(Eif2dsh_neg$A_content)
median(Eif2dsh_neg$A_content)

#Kolmogorov-Smirnov test
ks.test(Denrsh2_neg$GC_content, all_expressed$GC_content)
ks.test(Eif2dsh_neg$GC_content, all_expressed$GC_content)
ks.test(Denrsh2_neg$A_content, all_expressed$A_content)
ks.test(Eif2dsh_neg$A_content, all_expressed$A_content)

pdf('GC_content_5UTR_all_negative_FDR_0.1_boxplot.pdf')
boxplot(all_expressed$GC_content,Denrsh2_neg$GC_content, Eif2dsh_neg$GC_content,col=c('lightgray','lightpink2','lightskyblue3'), names = c("all transcripts","Denr shRNA","Eif2d shRNA"))
dev.off()

pdf('GC_content_5UTR_all_negative_FDR_0.1_violinplot.pdf')
vioplot(all_expressed$GC_content,Denrsh2_neg$GC_content, Eif2dsh_neg$GC_content,col=c('lightgray','lightpink2','lightskyblue3'), names = c("all transcripts","Denr shRNA","Eif2d shRNA"))
dev.off()

pdf('A_content_5UTR_all_negative_FDR_0.1_violinplot.pdf')
vioplot(all_expressed$A_content,Denrsh2_neg$A_content, Eif2dsh_neg$A_content,col=c('lightgray','lightpink2','lightskyblue3'), names = c("all transcripts","Denr shRNA","Eif2d shRNA"))
dev.off()

##Penultimate codons

seq <- read.table('GeneID_transcriptID_5UTR_single_isoform_sameCDS_longest_uORFs_seq.txt')
colnames(seq) <- c('Gene_ID','uORF_start','uORF_stop','sequence')
seq$uORF_length <- seq$uORF_stop-seq$uORF_start
seq <- seq[seq$uORF_length > 3,]

all_translated_uORF <- read.table('all_translated_uORFs.3.txt')
colnames(all_translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads')
all_translated_uORF$Gene_ID <- substring(all_translated_uORF[,1], 1, 18)
all_translated_uORF$uORF_start <- all_translated_uORF$uORF_start + 16
all_translated_uORF$uORF_stop <- all_translated_uORF$uORF_stop + 16
all_translated_uORF_seq <- merge(seq,all_translated_uORF, by = c('Gene_ID','uORF_start','uORF_stop'))
all_translated_uORF_seq$penultimate_codon <- substring(all_translated_uORF_seq$sequence, nchar(all_translated_uORF_seq$sequence)-5, nchar(all_translated_uORF_seq$sequence)-3)
all_penultimate_codon <- table(all_translated_uORF_seq$penultimate_codon)
write.table(all_penultimate_codon,'all_penultimate_codons.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
all_penultimate_codon <- read.table('all_penultimate_codons.txt')
total_all <- sum(all_penultimate_codon$V2)
all_penultimate_codon$ratio_all <- all_penultimate_codon$V2 / total_all

Denrsh2_translated_uORF <- read.table('Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Denrsh2_translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads')
Denrsh2_translated_uORF$Gene_ID <- substring(Denrsh2_translated_uORF[,1], 1, 18)
Denrsh2_translated_uORF$uORF_start <- Denrsh2_translated_uORF$uORF_start + 16
Denrsh2_translated_uORF$uORF_stop <- Denrsh2_translated_uORF$uORF_stop + 16
Denrsh2_translated_uORF_seq <- merge(seq,Denrsh2_translated_uORF, by = c('Gene_ID','uORF_start','uORF_stop'))
Denrsh2_translated_uORF_seq$penultimate_codon <- substring(Denrsh2_translated_uORF_seq$sequence, nchar(Denrsh2_translated_uORF_seq$sequence)-5, nchar(Denrsh2_translated_uORF_seq$sequence)-3)
Denrsh2_penultimate_codon <- table(Denrsh2_translated_uORF_seq$penultimate_codon)
write.table(Denrsh2_penultimate_codon,'Denrsh2_all_neg_FDR_0.1_translated_uORFs_penultimate_codons.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
Denrsh2_penultimate_codon <- read.table('Denrsh2_all_neg_FDR_0.1_translated_uORFs_penultimate_codons.txt')
total_Denrsh2 <- sum(Denrsh2_penultimate_codon$V2)
Denrsh2_penultimate_codon$ratio_Denrsh2 <- Denrsh2_penultimate_codon$V2 / total_Denrsh2
penultimate_codons <- merge(Denrsh2_penultimate_codon, all_penultimate_codon, by = 'V1')
penultimate_codons$ratio_Denrsh2_all <- penultimate_codons$ratio_Denrsh2 / penultimate_codons$ratio_all
penultimate_codons$log2_ratio_Denrsh2_all <- log2(penultimate_codons$ratio_Denrsh2 / penultimate_codons$ratio_all)
bt <- function(a, b, c) {binom.test(a, b, c, alternative= c("two.sided"), conf.level = 0.95)$p.value}
penultimate_codons$binom_test_Denrsh2 <- mapply(bt, penultimate_codons$V2.x, total_Denrsh2,(penultimate_codons$V2.y/total_all))
write.csv(penultimate_codons, 'Denr_sh2_penultimate_codons_log2.csv', row.names = FALSE)
pdf('penultimate_codon_Denrsh2_translated_uORFs_all_negative_FDR_0.1.pdf')
plot(log2(penultimate_codons$ratio_Denrsh2_all), -log10(penultimate_codons$binom_test_Denrsh2), pch=20,cex=2, col = ifelse(penultimate_codons$binom_test_Denrsh2 >= 0.05,'dark grey','orange'), xlab='fold enrichment', ylab='-log10(pvalue)',main='Denr shRNA2', ylim = c(0,3), xlim = c(-3.5,3.5))
text(penultimate_codons$V1, x=penultimate_codons$log2_ratio_Denrsh2_all, y = -log10(penultimate_codons$binom_test_Denrsh2), cex= 1, pos=3)
dev.off()

Eif2dsh_translated_uORF <- read.table('Eif2dsh_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Eif2dsh_translated_uORF) <- c('Gene_ID','uORF_start','uORF_stop','number_reads')
Eif2dsh_translated_uORF$Gene_ID <- substring(Eif2dsh_translated_uORF[,1], 1, 18)
Eif2dsh_translated_uORF$uORF_start <- Eif2dsh_translated_uORF$uORF_start + 16
Eif2dsh_translated_uORF$uORF_stop <- Eif2dsh_translated_uORF$uORF_stop + 16
Eif2dsh_translated_uORF_seq <- merge(seq,Eif2dsh_translated_uORF, by = c('Gene_ID','uORF_start','uORF_stop'))
Eif2dsh_translated_uORF_seq$penultimate_codon <- substring(Eif2dsh_translated_uORF_seq$sequence, nchar(Eif2dsh_translated_uORF_seq$sequence)-5, nchar(Eif2dsh_translated_uORF_seq$sequence)-3)
Eif2dsh_penultimate_codon <- table(Eif2dsh_translated_uORF_seq$penultimate_codon)
write.table(Eif2dsh_penultimate_codon,'Eif2dsh_all_neg_FDR_0.1_translated_uORFs_penultimate_codons.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
Eif2dsh_penultimate_codon <- read.table('Eif2dsh_all_neg_FDR_0.1_translated_uORFs_penultimate_codons.txt')
total_Eif2dsh <- sum(Eif2dsh_penultimate_codon$V2)
Eif2dsh_penultimate_codon$ratio_Eif2dsh <- Eif2dsh_penultimate_codon$V2 / total_Eif2dsh
penultimate_codons <- merge(Eif2dsh_penultimate_codon, all_penultimate_codon, by = 'V1')
penultimate_codons$ratio_Eif2dsh_all <- penultimate_codons$ratio_Eif2dsh / penultimate_codons$ratio_all
penultimate_codons$log2_ratio_Eif2dsh_all <- log2(penultimate_codons$ratio_Eif2dsh / penultimate_codons$ratio_all)
bt <- function(a, b, c) {binom.test(a, b, c, alternative= c("two.sided"), conf.level = 0.95)$p.value}
penultimate_codons$binom_test_Eif2dsh <- mapply(bt, penultimate_codons$V2.x, total_Eif2dsh, penultimate_codons$V2.y/total_all)
write.csv(penultimate_codons, 'Eif2d_sh_penultimate_codons_log2.csv')
pdf('penultimate_codon_Eif2dsh_translated_uORFs_all_negative_FDR_0.1.pdf')
plot(log2(penultimate_codons$ratio_Eif2dsh_all), -log10(penultimate_codons$binom_test_Eif2dsh), pch=20,cex=2, col = ifelse(penultimate_codons$binom_test_Eif2dsh > 0.05,'dark grey','orange'), xlab='fold enrichment', ylab='-log10(pvalue)',main='Eif2d shRNA',ylim = c(0,3), xlim = c(-3.5,3.5))
text(penultimate_codons$V1, x=penultimate_codons$log2_ratio_Eif2dsh_all, y = -log10(penultimate_codons$binom_test_Eif2dsh), cex= 1, pos=3)
dev.off()

##Stop codon plot translated uORFs

translated_uORFs <- read.table('all_translated_uORFs.3.txt')
colnames(translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
translated_uORFs$Gene_ID <- substring(translated_uORFs[,1], 1, 18)
translated_uORFs$Transcript_ID <- substring(translated_uORFs[,1], 20, 37)
translated_uORFs$uORF_stop <- translated_uORFs$uORF_stop + 16
translated_uORFs$uORF_start <- translated_uORFs$uORF_start + 16
translated_uORFs$stopCodon <- 'stopCodon'
translated_uORFs$strand <- '+'
translated_uORFs$number <- '1'
translated_uORFs <- translated_uORFs[(translated_uORFs$uORF_stop - translated_uORFs$uORF_start) >=6,]
translated_uORFs <- translated_uORFs[,c(6,7,3,8,9)]
write.table(translated_uORFs, 'all_translated_uORFs.stopCodon.min6.3.sga', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Denrsh2_neg_translated_uORFs <- read.table('Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Denrsh2_neg_translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
Denrsh2_neg_translated_uORFs$Gene_ID <- substring(Denrsh2_neg_translated_uORFs[,1], 1, 18)
Denrsh2_neg_translated_uORFs$Transcript_ID <- substring(Denrsh2_neg_translated_uORFs[,1], 20, 37)
Denrsh2_neg_translated_uORFs$uORF_stop <- Denrsh2_neg_translated_uORFs$uORF_stop + 16
Denrsh2_neg_translated_uORFs$uORF_start <- Denrsh2_neg_translated_uORFs$uORF_start + 16
Denrsh2_neg_translated_uORFs$stopCodon <- 'stopCodon'
Denrsh2_neg_translated_uORFs$strand <- '+'
Denrsh2_neg_translated_uORFs$number <- '1'
Denrsh2_neg_translated_uORFs <- Denrsh2_neg_translated_uORFs[(Denrsh2_neg_translated_uORFs$uORF_stop - Denrsh2_neg_translated_uORFs$uORF_start) >=6,]
Denrsh2_neg_translated_uORFs <- Denrsh2_neg_translated_uORFs[,c(6,7,3,8,9)]
write.table(Denrsh2_neg_translated_uORFs, 'Denrsh2_neg_FDR_0.1_translated_uORFs.stopCodon.min6.3.sga', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

setwd("/home/rmeurs/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/uORF_annotations/stop_codon_plots/")

ctrl <- read.table('translated_uORFs_stopCodon_vs_ctrl.Eif2dsh4_neg.min6.dat')
colnames(ctrl) <- c("position_relative_to_stop_codon", 'read_counts_ctrl')
Denrsh2 <- read.table('translated_uORFs_stopCodon_vs_Denrsh2.Eif2dsh4_neg.min6.dat')
colnames(Denrsh2) <- c("position_relative_to_stop_codon", 'read_counts_Denrsh2')
Eif2dsh3 <- read.table('translated_uORFs_stopCodon_vs_Eif2dsh3.Eif2dsh4_neg.min6.dat')
colnames(Eif2dsh3) <- c("position_relative_to_stop_codon", 'read_counts_Eif2dsh3')
Eif2dsh4 <- read.table('translated_uORFs_stopCodon_vs_Eif2dsh4.Eif2dsh4_neg.min6.dat')
colnames(Eif2dsh4) <- c("position_relative_to_stop_codon", 'read_counts_Eif2dsh4')
stop_codon <- merge(ctrl, Denrsh2, by = 'position_relative_to_stop_codon', all=TRUE)
stop_codon <- merge(stop_codon, Eif2dsh3, by = 'position_relative_to_stop_codon', all=TRUE)
stop_codon <- merge(stop_codon, Eif2dsh4, by = 'position_relative_to_stop_codon', all=TRUE)
stop_codon$read_counts_ctrl <- lapply(stop_codon[,2], function(x) (x*1000000)/102458233)
stop_codon$read_counts_Denrsh2 <- lapply(stop_codon[,3], function(x) (x*1000000)/66700613)
stop_codon$read_counts_Eif2dsh3 <- lapply(stop_codon[,4], function(x) (x*1000000)/78844780)
stop_codon$read_counts_Eif2dsh4 <- lapply(stop_codon[,5], function(x) (x*1000000)/55527239)
stop_codon$position_relative_to_stop_codon <- lapply(stop_codon[,1], function(x) (x+15))
stop_codon <- as.matrix(stop_codon)
write.csv(stop_codon, 'translated_uORFs_stopCodon.Eif2dsh4_neg.min6.csv', row.names = FALSE)

stop_codon <- read.csv('translated_uORFs_stopCodon.csv')
pdf('translated_uORFs_stopCodon_plot.pdf', width=80, height=40)
plot(stop_codon$position_relative_to_stop_codon, stop_codon$read_counts_ctrl, type = "l", xlim = c(-30,15), ylim = c(0,50), col='gray35', lwd=15)  
lines(stop_codon$position_relative_to_stop_codon, stop_codon$read_counts_Denrsh2, col = 'lightpink2', lwd=15)
lines(stop_codon$position_relative_to_stop_codon, stop_codon$read_counts_Eif2dsh3, col = 'lightseagreen', lwd=15)
lines(stop_codon$position_relative_to_stop_codon, stop_codon$read_counts_Eif2dsh4, col = 'lightskyblue3', lwd=15)
dev.off()

##Distance last uORF stop - CDS start

setwd("/home/rmeurs/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/uORF_annotations/")

CDS <- read.table('~/Documents/PhD/eIF2D_DENR_project/Ribosome_profiling/Analysis/annotations/Mus_musculus.GRCm38.100.startStop.IDS.bed')
colnames(CDS) <- c('IDs','CDS_start','CDS_stop','protein_coding')
translated_uORFs <- read.table('all_translated_uORFs.3.txt')
colnames(translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
uORF_stop_CDS_start <- merge(translated_uORFs, CDS, by = 'IDs')
uORF_stop_CDS_start <- uORF_stop_CDS_start[,c(1,3,5)]
uORF_stop_CDS_start$distance <- uORF_stop_CDS_start$CDS_start - uORF_stop_CDS_start$uORF_stop
uORF_stop_CDS_start_min <- uORF_stop_CDS_start %>% group_by(IDs) %>% top_n(-1, distance)
median(uORF_stop_CDS_start_min$distance)
hist(uORF_stop_CDS_start_min$distance,breaks = 1500, xlim = c(0,1500))

Dsh2_translated_uORFs <- read.table('Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Dsh2_translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
Dsh2_uORF_stop_CDS_start <- merge(Dsh2_translated_uORFs, CDS, by = 'IDs')
Dsh2_uORF_stop_CDS_start <- Dsh2_uORF_stop_CDS_start[,c(1,3,5)]
Dsh2_uORF_stop_CDS_start$distance <- Dsh2_uORF_stop_CDS_start$CDS_start - Dsh2_uORF_stop_CDS_start$uORF_stop
Dsh2_uORF_stop_CDS_start_min <- Dsh2_uORF_stop_CDS_start %>% group_by(IDs) %>% top_n(-1, distance)
median(Dsh2_uORF_stop_CDS_start_min$distance)
hist(Dsh2_uORF_stop_CDS_start_min$distance,breaks = 500, xlim = c(0,1500))
ks.test(Dsh2_uORF_stop_CDS_start_min$distance, uORF_stop_CDS_start_min$distance)

Esh_translated_uORFs <- read.table('Eif2dsh_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Esh_translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
Esh_uORF_stop_CDS_start <- merge(Esh_translated_uORFs, CDS, by = 'IDs')
Esh_uORF_stop_CDS_start <- Esh_uORF_stop_CDS_start[,c(1,3,5)]
Esh_uORF_stop_CDS_start$distance <- Esh_uORF_stop_CDS_start$CDS_start - Esh_uORF_stop_CDS_start$uORF_stop
Esh_uORF_stop_CDS_start_min <- Esh_uORF_stop_CDS_start %>% group_by(IDs) %>% top_n(-1, distance)
median(Esh_uORF_stop_CDS_start_min$distance)
hist(Esh_uORF_stop_CDS_start_min$distance,breaks = 500, xlim = c(0,1000))
ks.test(Esh_uORF_stop_CDS_start_min$distance, uORF_stop_CDS_start_min$distance)

pdf('distance_uORF_stop_CDS_start_all_negative_FDR_0.1_boxplot.pdf')
boxplot(uORF_stop_CDS_start_min$distance,Dsh2_uORF_stop_CDS_start_min$distance, Esh_uORF_stop_CDS_start_min$distance,ylim=c(0,1000),col=c('lightgray','lightpink2','lightskyblue3'), names = c("all transcripts","Denr shRNA","Eif2d shNRA"))
dev.off()

##Distance cap - first uORF

translated_uORFs <- read.table('all_translated_uORFs.3.txt')
colnames(translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
translated_uORFs_start_min <- translated_uORFs %>% group_by(IDs) %>% top_n(-1, uORF_start)
median(translated_uORFs_start_min$uORF_start)
hist(translated_uORFs_start_min$uORF_start,breaks = 1500, xlim = c(0,1500))

Dsh2_translated_uORFs <- read.table('Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Dsh2_translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
Dsh2_translated_uORFs_start_min <- Dsh2_translated_uORFs %>% group_by(IDs) %>% top_n(-1, uORF_start)
median(Dsh2_translated_uORFs_start_min$uORF_start)
hist(Dsh2_translated_uORFs_start_min$uORF_start,breaks = 500, xlim = c(0,1500))
ks.test(Dsh2_translated_uORFs_start_min$uORF_start, translated_uORFs_start_min$uORF_start)

Esh_translated_uORFs <- read.table('Eif2dsh_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Esh_translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
Esh_translated_uORFs_start_min <- Esh_translated_uORFs %>% group_by(IDs) %>% top_n(-1, uORF_start)
median(Esh_translated_uORFs_start_min$uORF_start)
hist(Esh_translated_uORFs_start_min$uORF_start,breaks = 500, xlim = c(0,1500))
ks.test(Esh_translated_uORFs_start_min$uORF_start, translated_uORFs_start_min$uORF_start)

pdf('distance_cap_start_uORF_all_negative_FDR_0.1_boxplot.pdf')
boxplot(translated_uORFs_start_min$uORF_start,Dsh2_translated_uORFs_start_min$uORF_start, Esh_translated_uORFs_start_min$uORF_start,ylim=c(0,1000),col=c('lightgray','lightpink2','lightskyblue3'), names = c("all transcripts","Denr shRNA","Eif2d shRNA"))
dev.off()

##Kozak context uORFs start (python merge_sequence.py)

Kozak <- read.table('all_translated_uORFs.3.Kozak_context.txt')
colnames(Kozak) <- c('IDs','uORF_start','uORF_stop','Kozak_context','number_reads')
Kozak$Kozak_score <- 'a'
for (n in 1:nrow(Kozak)) {
  Kozak_context <- Kozak[n,4]
  Kozak_score = 0
  one <- substr(Kozak_context,1,1)
  two <- substr(Kozak_context,2,2)
  three <- substr(Kozak_context,3,3)
  four <- substr(Kozak_context,4,4)
  five <- substr(Kozak_context,5,5)
  six <- substr(Kozak_context,6,6)
  ten <- substr(Kozak_context,10,10)
  if (one == 'G'){
    Kozak_score = Kozak_score + 3
  }
  if (two == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (three == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (four == 'A' || four == 'G') {
    Kozak_score = Kozak_score + 3
  }
  if (five == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (six == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (ten == 'G') {
    Kozak_score = Kozak_score + 3
  }
  Kozak[n,6] <- as.integer(Kozak_score)
}
median(as.integer(Kozak$Kozak_score))
mean(as.integer(Kozak$Kozak_score))
Kozak$Kozak_score <- as.numeric(Kozak$Kozak_score)

Kozak <- read.table('all_translated_uORFs.3.Kozak_context.txt')
colnames(Kozak) <- c('IDs','uORF_start','uORF_stop','Kozak_context','number_reads')
Dsh2 <- read.table('Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Dsh2) <- c('IDs','uORF_start','uORF_stop','number_reads')
Dsh2$uORF_start <- Dsh2$uORF_start + 16
Dsh2$uORF_stop <- Dsh2$uORF_stop + 16
Dsh2_Kozak <- merge(Dsh2, Kozak, by = c('IDs','uORF_start','uORF_stop','number_reads'))
Dsh2_Kozak$Kozak_score <- 'a'
for (n in 1:nrow(Dsh2_Kozak)) {
  Kozak_context <- Dsh2_Kozak[n,5]
  Kozak_score = 0
  one <- substr(Kozak_context,1,1)
  two <- substr(Kozak_context,2,2)
  three <- substr(Kozak_context,3,3)
  four <- substr(Kozak_context,4,4)
  five <- substr(Kozak_context,5,5)
  six <- substr(Kozak_context,6,6)
  ten <- substr(Kozak_context,10,10)
  if (one == 'G'){
    Kozak_score = Kozak_score + 3
  }
  if (two == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (three == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (four == 'A' || four == 'G') {
    Kozak_score = Kozak_score + 3
  }
  if (five == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (six == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (ten == 'G') {
    Kozak_score = Kozak_score + 3
  }
  Dsh2_Kozak[n,6] <- as.integer(Kozak_score)
}
median(as.integer(Dsh2_Kozak$Kozak_score))
mean(as.integer(Dsh2_Kozak$Kozak_score))
ks.test(as.integer(Dsh2_Kozak$Kozak_score), as.integer(Kozak$Kozak_score),alternative = c("two.sided"))
Dsh2_Kozak$Kozak_score <- as.numeric(Dsh2_Kozak$Kozak_score)

Esh <- read.table('Eif2dsh_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Esh) <- c('IDs','uORF_start','uORF_stop','number_reads')
Esh$uORF_start <- Esh$uORF_start + 16
Esh$uORF_stop <- Esh$uORF_stop + 16
Esh_Kozak <- merge(Esh, Kozak, by = c('IDs','uORF_start','uORF_stop','number_reads'))
Esh_Kozak$Kozak_score <- 'a'
for (n in 1:nrow(Esh_Kozak)) {
  Kozak_context <- Esh_Kozak[n,5]
  Kozak_score = 0
  one <- substr(Kozak_context,1,1)
  two <- substr(Kozak_context,2,2)
  three <- substr(Kozak_context,3,3)
  four <- substr(Kozak_context,4,4)
  five <- substr(Kozak_context,5,5)
  six <- substr(Kozak_context,6,6)
  ten <- substr(Kozak_context,10,10)
  if (one == 'G'){
    Kozak_score = Kozak_score + 3
  }
  if (two == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (three == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (four == 'A' || four == 'G') {
    Kozak_score = Kozak_score + 3
  }
  if (five == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (six == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (ten == 'G') {
    Kozak_score = Kozak_score + 3
  }
  Esh_Kozak[n,6] <- as.integer(Kozak_score)
}
median(as.integer(Esh_Kozak$Kozak_score))
mean(as.integer(Esh_Kozak$Kozak_score))
ks.test(as.integer(Esh_Kozak$Kozak_score), as.integer(Kozak$Kozak_score),alternative = c("two.sided"))
Esh_Kozak$Kozak_score <- as.numeric(Esh_Kozak$Kozak_score)

pdf('kozak_score_uORF_start_all_negative_FDR_0.1_boxplot.pdf')
boxplot(Kozak$Kozak_score,Dsh2_Kozak$Kozak_score, Esh_Kozak$Kozak_score, ylim=c(0,15),col=c('lightgray','lightgreen','skyblue2'), names = c("all transcripts","Denr shRNA2","Eif2d shRNA"))
dev.off()

pdf('kozak_score_uORF_start_all_negative_FDR_0.1_hist.pdf')
par(mfrow=c(2,2))
hist(as.integer(Kozak$Kozak_score), xlim = c(0,12), breaks = 12)
hist(as.integer(Dsh2_Kozak$Kozak_score), xlim = c(0,12), breaks = 12)
hist(as.integer(Esh_Kozak$Kozak_score), xlim = c(0,12), breaks = 12)
dev.off()



##Kozak context CDS start (python merge_sequence.py)

Kozak_CDS <- read.delim('Mus_musculus.GRCm38.100.startStop.Kozak_context.IDS.PC..bed', col.names = c("IDs","Kozak_context"))
Kozak_CDS$IDs <- substring(Kozak_CDS$IDs, 1, 37)
translated_uORFs <- read.table('all_translated_uORFs.3.txt')
colnames(translated_uORFs) <- c('IDs','uORF_start','uORF_stop','number_reads')
translated_uORFs_CDS_Kozak <- merge(translated_uORFs, Kozak_CDS, by = 'IDs')
translated_uORFs_CDS_Kozak <- translated_uORFs_CDS_Kozak[!duplicated(translated_uORFs_CDS_Kozak[,c(1,5)]), ]
translated_uORFs_CDS_Kozak$Kozak_score <- 'a'
for (n in 1:nrow(translated_uORFs_CDS_Kozak)) {
  Kozak_context <- translated_uORFs_CDS_Kozak[n,5]
  Kozak_score = 0
  one <- substr(Kozak_context,1,1)
  two <- substr(Kozak_context,2,2)
  three <- substr(Kozak_context,3,3)
  four <- substr(Kozak_context,4,4)
  five <- substr(Kozak_context,5,5)
  six <- substr(Kozak_context,6,6)
  ten <- substr(Kozak_context,10,10)
  if (one == 'G'){
    Kozak_score = Kozak_score + 3
  }
  if (two == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (three == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (four == 'A' || four == 'G') {
    Kozak_score = Kozak_score + 3
  }
  if (five == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (six == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (ten == 'G') {
    Kozak_score = Kozak_score + 3
  }
  translated_uORFs_CDS_Kozak[n,6] <- as.integer(Kozak_score)
}
median(as.integer(translated_uORFs_CDS_Kozak$Kozak_score))
mean(as.integer(translated_uORFs_CDS_Kozak$Kozak_score))
translated_uORFs_CDS_Kozak$Kozak_score <- as.numeric(translated_uORFs_CDS_Kozak$Kozak_score)

Dsh2 <- read.table('Denrsh2_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Dsh2) <- c('IDs','uORF_start','uORF_stop','number_reads')
Dsh2_Kozak_CDS <- merge(Dsh2, Kozak_CDS, by = 'IDs')
Dsh2_Kozak_CDS <- Dsh2_Kozak_CDS[!duplicated(Dsh2_Kozak_CDS[,c(1,5)]), ]
Dsh2_Kozak_CDS$Kozak_score <- 'a'
for (n in 1:nrow(Dsh2_Kozak_CDS)) {
  Kozak_context <- Dsh2_Kozak_CDS[n,5]
  Kozak_score = 0
  one <- substr(Kozak_context,1,1)
  two <- substr(Kozak_context,2,2)
  three <- substr(Kozak_context,3,3)
  four <- substr(Kozak_context,4,4)
  five <- substr(Kozak_context,5,5)
  six <- substr(Kozak_context,6,6)
  ten <- substr(Kozak_context,10,10)
  if (one == 'G'){
    Kozak_score = Kozak_score + 3
  }
  if (two == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (three == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (four == 'A' || four == 'G') {
    Kozak_score = Kozak_score + 3
  }
  if (five == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (six == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (ten == 'G') {
    Kozak_score = Kozak_score + 3
  }
  Dsh2_Kozak_CDS[n,6] <- as.integer(Kozak_score)
}
median(as.integer(Dsh2_Kozak_CDS$Kozak_score))
mean(as.integer(Dsh2_Kozak_CDS$Kozak_score))
ks.test(as.integer(Dsh2_Kozak_CDS$Kozak_score), as.integer(translated_uORFs_CDS_Kozak$Kozak_score),alternative = c("two.sided"))
Dsh2_Kozak_CDS$Kozak_score <- as.numeric( Dsh2_Kozak_CDS$Kozak_score)

Esh <- read.table('Eif2dsh_neg_TE_FC_FDR_0.1_translated_uORFs.3.txt')
colnames(Esh) <- c('IDs','uORF_start','uORF_stop','number_reads')
Esh_Kozak_CDS <- merge(Esh, Kozak_CDS, by = 'IDs')
Esh_Kozak_CDS <- Esh_Kozak_CDS[!duplicated(Esh_Kozak_CDS[,c(1,5)]), ]
Esh_Kozak_CDS$Kozak_score <- 'a'
for (n in 1:nrow(Esh_Kozak_CDS)) {
  Kozak_context <- Esh_Kozak_CDS[n,5]
  Kozak_score = 0
  one <- substr(Kozak_context,1,1)
  two <- substr(Kozak_context,2,2)
  three <- substr(Kozak_context,3,3)
  four <- substr(Kozak_context,4,4)
  five <- substr(Kozak_context,5,5)
  six <- substr(Kozak_context,6,6)
  ten <- substr(Kozak_context,10,10)
  if (one == 'G'){
    Kozak_score = Kozak_score + 3
  }
  if (two == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (three == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (four == 'A' || four == 'G') {
    Kozak_score = Kozak_score + 3
  }
  if (five == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (six == 'C') {
    Kozak_score = Kozak_score + 1
  }
  if (ten == 'G') {
    Kozak_score = Kozak_score + 3
  }
  Esh_Kozak_CDS[n,6] <- as.integer(Kozak_score)
}
median(as.integer(Esh_Kozak_CDS$Kozak_score))
mean(as.integer(Esh_Kozak_CDS$Kozak_score))
ks.test(as.integer(Esh_Kozak_CDS$Kozak_score), as.integer(translated_uORFs_CDS_Kozak$Kozak_score),alternative = c("two.sided"))
Esh_Kozak_CDS$Kozak_score <- as.numeric(Esh_Kozak_CDS$Kozak_score)

pdf('kozak_score_CDS_start_all_negative_FDR_0.1_boxplot.pdf')
boxplot(translated_uORFs_CDS_Kozak$Kozak_score,Dsh2_Kozak_CDS$Kozak_score, Esh_Kozak_CDS$Kozak_score, ylim=c(0,15),col=c('lightgray','lightgreen','skyblue2'), names = c("all transcripts","Denr shRNA2","Eif2d shRNA"))
dev.off()

pdf('kozak_score_CDS_start_all_negative_FDR_0.1_hist.pdf')
par(mfrow=c(2,2))
hist(as.integer(translated_uORFs_CDS_Kozak$Kozak_score), xlim = c(0,12), breaks = 12)
hist(as.integer(Dsh2_Kozak_CDS$Kozak_score), xlim = c(0,12), breaks = 12)
hist(as.integer(Esh_Kozak_CDS$Kozak_score), xlim = c(0,12), breaks = 12)
dev.off()