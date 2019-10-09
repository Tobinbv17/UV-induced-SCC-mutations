library(dplyr)
library(rtracklayer)


GTF = as.data.frame(rtracklayer::import('~/czlabwork/annotation/gencode.vM20.annotation.gtf', format = "gtf")) 
gtf1 = GTF %>% filter(type == 'transcript') %>% select (seqnames, start, end, strand, transcript_id, gene_id, gene_name)

#vcf_out_gtf <- NULL
#for (i in 1:nrow(vcf_out)) {
#  for (j in 1:nrow(gtf1)) {
#    if (vcf_out[i,]["chrom"] == gtf1[j,]["seqnames"]) { 
#      if (  vcf_out[i,]["pos"] > gtf1[j,]["start"] & vcf_out[i,]["pos"] < gtf1[j,]["end"]) { 
#        vcf_out_gtf = rbind(vcf_out_gtf, cbind ( vcf_out[i,], gtf1[j,]))
#        vcf_out_gtf <<- vcf_out_gtf
#      }
#    }
#  }
#}
#write.table(vcf_out_gtf, '~/czlabwork/vcfczold/20191003_R/vcf_out_gtf_r.dat', row.names = FALSE, sep = '\t', quote = FALSE)



vcf_gtf_cg_Nname = read.csv("~/czlabwork/vcfczold/20191007_py/vcf_gtf.dat", header = F, sep = "\t")
colnames(vcf_gtf_cg_Nname) = c("chrom", "pos", "ref", "alt", "strand", "sample", "GT", "AD", "gene_id", "transcript_id", "str_CBD")
vcf_gtf_cg_Nname$diff = c(diff(vcf_gtf_cg_Nname$pos),0)

colData = readRDS("~/czlabwork/vcfczold/20191007_py/colgroup.rds")
colData = data.frame(as.factor(colData$V1), as.factor(colData$V2))
colnames(colData) = c("group", "sample")


vcf_gtf_cg_Nname_diff1 = vcf_gtf_cg_Nname[vcf_gtf_cg_Nname$diff==1,] ### filter CC-TT/GG-AA dint variants only 
vcf_gtf_cg_Nname_diff1$group = left_join(as.data.frame(vcf_gtf_cg_Nname_diff1), colData, by = "sample")

a = vcf_gtf_cg_Nname_diff1$group
a_N_diCT = table(a$str_CBD, a$group)
a_N_diCT1 = t(table(a$str_CBD, a$group)) # transpose the dataframe
####         non-transcribed transcribed
#nonUV                70         117
#IgG                 807         252
#antiPD1             647         188
fisher.test(table(a$group, a$str_CBD))

saveRDS(vcf_gtf_cg_Nname, "~/czlabwork/vcfczold/20191007_py/vcf_gtf_cg_diff_all.rds")

################
# data frame, factorize each column
vcf_gtf_cg_Nname = as.data.frame(lapply(snp_ccgg_data_gname, as.factor))
# select duplicated rows
data_dup_row = vcf_gtf_cg_Nname[duplicated(vcf_gtf_cg_Nname$gene_id) | duplicated(vcf_gtf_cg_Nname$gene_id, fromLast = T),]
# add gene_name column
vcf_gtf_cg_Nname$gene_id = left_join(vcf_gtf_cg_Nname, gtf1$gene_name, by = "gene_name")






