# import GTF files into dataframe
library(rtracklayer)
# http://bioconductor.org/packages/release/bioc/html/rtracklayer.html
gtf <- rtracklayer::import('/Users/yao/Desktop/czlabwork/Gannotation/gencode.vM22.annotation.gtf')
gtf_df=as.data.frame(gtf)
write.table(gtf_df, file='/Users/yao/Desktop/czlabwork/gtf_dfraw.csv',row.names = F, col.names = T)
# a = sort(table(gtf_df$type)) sort table according to gtf type features

# write exon level and transcript level table
library(dplyr) 
gtf_work = subset(gtf_df, select = c('seqnames', 'start', 'end', 'gene_id','strand', 'type', 'gene_type', 'gene_name'))
write.table(gtf_work, file='/Users/yao/Desktop/czlabwork/rwork/gtf_work.csv', row.names = F, col.names = T)

gtf_exon = subset(gtf_work, type == 'exon')
write.table(gtf_transcript, '/Users/yao/Desktop/czlabwork/rwork/gtf_exon.csv', row.names = F, col.names = T)

gtf_transcript = subset(gtf_work, type == 'transcript')
#data = read.csv('/Users/yao/Desktop/czlabwork/gtf_exon.csv', header = F, sep = '', stringAsFactors = F)
write.table(gtf_transcript, '/Users/yao/Desktop/czlabwork/rwork/gtf_transcript.csv', row.names = F, col.names = T)
# bed file contain factor "", can be edited by command line sed 's/"//g' < gtf2bed.bed > noquote.bed

# write bed file table by adding 1 in end coordinate
gtf_exon1 = gtf_exon
gtf_transcript1 = gtf_transcript
gtf_exon1$end = gtf_exon$end+1
gtf_transcript1$end = gtf_transcript$end+1

write.table(gtf_exon1, '/Users/yao/Desktop/czlabwork/rwork/gtf_exon.bed', sep = '\t', row.names = F, col.names = F)
write.table(gtf_transcript1, '/Users/yao/Desktop/czlabwork/rwork/gtf_transcript.bed', sep = '\t', row.names = F, col.names = F)

