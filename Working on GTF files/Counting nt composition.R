# import GTF files into dataframe
library(rtracklayer)
# http://bioconductor.org/packages/release/bioc/html/rtracklayer.html
gtf <- rtracklayer::import('/Users/yao/Desktop/czlabwork/Gannotation/gencode.vM20.annotation.gtf')
gtf_df=as.data.frame(gtf)
saveRDS(gtf_df, '/Users/yao/Desktop/czlabwork/GTFM20version/gtf.rds')
gtfread = readRDS('/Users/yao/Desktop/czlabwork/GTFM20version/gtf.rds')

library(dplyr)
# write subset of gtf data for downstream work
gtfwork = subset(gtfread, select = c ("seqnames", "start", "end", "width", "strand", "type", "gene_id", "gene_name", "transcript_id", "exon_id"))
saveRDS(gtfwork, '/Users/yao/Desktop/czlabwork/GTFM20version/gtfwork.rds')

# write exon level and transcript level table
# dplyr functions: select(), filter(), mutate(),  group_by(), and summarize(). To select columns of a data frame, use select(),To choose rows, use filter()

gtfexon = filter(gtfwork, type == 'exon')
saveRDS(gtfexon, './gtfexon.rds')

gtftrs = filter(gtfwork, type == 'transcript')
saveRDS(gtftrs, './gtftranscript.rds')

# calculate nt frequency by using 
# 1. oligonucleotideFrequency
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

GTFfile = "/Users/yao/Desktop/czlabwork/Gannotation/gencode.vM20.annotation.gtf"
FASTAfile = "/Users/yao/Desktop/czlabwork/Gannotation/GRCm38.p6.genome.fa"

# just look at exon level
GTF <- import.gff(GTFfile, format="gtf", feature.type="exon")
#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)
#Note: Dinucleotide counts for exons on the (-) strand are calculated from the reverse complement strand
# single nt & dint counts in exons
NUCcomp <- oligonucleotideFrequency(getSeq(FASTA, GTF), width = 1) # matrix
DINUCcomp = oligonucleotideFrequency(getSeq(FASTA, GTF), width = 2) # matrix
# data frames:columns (variables) can be different types (numeric/character/logical etc.). Matrices:same type.
ntcomp = bind_cols(gtfexon, as.data.frame(NUCcomp), as.data.frame(DINUCcomp))
saveRDS(ntcomp, './ntcomp_exon.rds')

# nt count group by transcript level
ntcompBytrpt = ntcomp %>% group_by(transcript_id) %>% summarise_each(funs(sum), "A", "C", "G", "T", "AA","AC","AG","AT","CA","CC", "CG","CT","GA","GC","GG","GT", "TA","TC","TG","TT")
# joint GTFtranscript with nt count by transcript 
ntBytrpt = gtftrs %>% left_join(ntcompBytrpt, by = "transcript_id")
ntBytrpt = ntBytrpt[!duplicated(ntBytrpt$transcript_id),]
# optional:
ntBytrpt = ntBytrpt %>% select (-c("exon_id"))
saveRDS(ntBytrpt, './ntBytranscript.rds')
#write.table(ntBytrpt, './ntBytrpt.dat', sep = "\t", quote = F, row.names = F)

#normalization by transcript length
ntBytrpt_norm = ntBytrpt %>% transmute(A = A/width, C = C/width, G = G/width, T = T/width, AA = AA/width, AC = AC/width, AG = AG/width, AT = AT/width, CA = CA/width, CC = CC/width, CG = CG/width, CT = CT/width, GA = GA/width, GC = GC/width, GG = GG/width, GT = GT/width, TA = TA/width, TC = TC/width, TG = TG/width, TT = TT/width)
ntBytrpt_norm = cbind(ntBytrpt[,c(1,5,9)], ntBytrpt_norm)
saveRDS(ntBytrpt_norm, './ntBytrpt_norm.rds')

