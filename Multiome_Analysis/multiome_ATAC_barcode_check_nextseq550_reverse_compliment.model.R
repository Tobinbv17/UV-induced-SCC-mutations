################################################################################
################################################################################
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
require(Biostrings)
#1.

QV1G12V_wt_atac1 = read.table('~/FASTQ/indexcheck/nextseq/0056/outs/fastq_path/HJM7NBGXJ/0056-QV1G12V-WT/QV1G12V-WT_S1_L001_R2.txt')
colnames(QV1G12V_wt_atac1) = c('Count','Barcode')
QV1G12V_wt_atac1$Lane = 'nextseq0056-Lane1'
write.csv(QV1G12V_wt_atac1, 'MultiomeSeq/multiome-atac/L1barcode.csv',row.names = F)
QV1G12V_wt_atac1 = read.csv('MultiomeSeq/multiome-atac/L1barcode.csv')

#make revese complement barcode for multiomeATAC-nextseq550
atac_wl = read.table('~/cellranger-v5_v6/Multiome/737K-arc-v1.txt')
colnames(atac_wl) = 'Barcode'
atac_wl$List = '10xATACbarcode'
atac_wl$revComp = sapply(atac_wl$Barcode, function(x) as.character(reverseComplement(DNAString(x))))
write.csv(atac_wl, 'MultiomeSeq/multiome-atac/737K-arc-ATAC-revComp-barcode.csv',row.names = F)
atac_wl = read.csv('MultiomeSeq/multiome-atac/737K-arc-ATAC-revComp-barcode.csv')

####### comparison between Barcode and revComp barcode
#barcode
df = left_join(QV1G12V_wt_atac1, atac_wl, by = 'Barcode')
df1 = df %>% group_by(Lane,List) %>% summarise(count = sum(Count))
#revComp barcode
QV1G12V_wt_atac1$revComp = QV1G12V_wt_atac1$Barcode
df_rc = left_join(QV1G12V_wt_atac1, atac_wl, by='revComp')
df2 = df_rc %>% group_by(Lane,List) %>% summarise(count = sum(Count))

#########Lane2,Lane3,Lane4
l2 = read.table('~/FASTQ/indexcheck/nextseq/0056/outs/fastq_path/HJM7NBGXJ/0056-QV1G12V-WT/QV1G12V-WT_S1_L002_R2.txt')
colnames(l2) = c('Count','Barcode')
l2$Lane = 'nextseq0056-Lane2'
write.csv(l2, 'MultiomeSeq/multiome-atac/L2barcode.csv',row.names = F)
l2 = read.csv('MultiomeSeq/multiome-atac/L2barcode.csv')

l2$revComp = l2$Barcode
df_rc = left_join(l2, atac_wl, by='revComp')
df3 = df_rc %>% group_by(Lane,List) %>% summarise(count = sum(Count))

#####
l3 = read.table('~/FASTQ/indexcheck/nextseq/0056/outs/fastq_path/HJM7NBGXJ/0056-QV1G12V-WT/QV1G12V-WT_S1_L003_R2.txt')
colnames(l3) = c('Count','Barcode')
l3$Lane = 'nextseq0056-Lane3'
write.csv(l3, 'MultiomeSeq/multiome-atac/L3barcode.csv',row.names = F)
l3 = read.csv('MultiomeSeq/multiome-atac/L3barcode.csv')

l3$revComp = l3$Barcode
df_rc = left_join(l3, atac_wl, by='revComp')
df4 = df_rc %>% group_by(Lane,List) %>% summarise(count = sum(Count))

#####
l4 = read.table('~/FASTQ/indexcheck/nextseq/0056/outs/fastq_path/HJM7NBGXJ/0056-QV1G12V-WT/QV1G12V-WT_S1_L004_R2.txt')
colnames(l4) = c('Count','Barcode')
l4$Lane = 'nextseq0056-Lane4'
write.csv(l4, 'MultiomeSeq/multiome-atac/L4barcode.csv',row.names = F)
l4 = read.csv('MultiomeSeq/multiome-atac/L4barcode.csv')

l4$revComp = l4$Barcode
df_rc = left_join(l4, atac_wl, by='revComp')
df5 = df_rc %>% group_by(Lane,List) %>% summarise(count = sum(Count))

#### combine l1,l2,l3,l4
df = rbind(df2,df3,df4,df5)
write.csv(df,'MultiomeSeq/multiome-atac/bc.comp.to.10Xmultiome.atacwl.csv')
grid.table(df)



