library(dplyr)
library(rtracklayer)


GTF = as.data.frame(rtracklayer::import('~/czlabwork/annotation/gencode.vM20.annotation.gtf', format = "gtf")) 
gtf1 = GTF %>% filter(type == 'transcript') %>% select (seqnames, start, end, strand, transcript_id, gene_id, gene_name)

################# some tricks to edit df
# select duplicated rows
duplicated_row_all = snp[duplicated(snp$gene_id) | duplicated(snp$gene_id, fromLast = T),]
# data frame, factorize each column
snp_factorized = as.data.frame(lapply(snp, as.factor))
##################################################################################################


################################ Section 1.  read and save dataframe NO FACTOR!

# write.table(snp, "~/czlabwork/vcfczold/20191007_py/snp.dat", sep = "\t")
# snp1 = read.csv("~/czlabwork/vcfczold/20191007_py/snp.dat", header = T, sep = "\t", stringsAsFactors = F)

vcf_gtf = read.csv("~/czlabwork/vcfczold/20191007_py/vcf_gtf.dat", header = F, sep = "\t", stringsAsFactors = F)
gtf1 = read.csv("~/czlabwork/vcfczold/20191007_py/gtf1.dat", header = T, sep = "\t", stringsAsFactors = F)
colData = readRDS("~/czlabwork/vcfczold/20191003_R/colgroup.rds")

saveRDS(vcf_gtf, "~/czlabwork/vcfczold/20191007_py/vcf_gtf.rds")
saveRDS(gtf1, "~/czlabwork/vcfczold/20191007_py/gtf1.rds")
saveRDS(colData, "~/czlabwork/vcfczold/20191007_py/colData.rds")

########################### section 2.  snp_data_all data.frame

snp = readRDS("~/czlabwork/vcfczold/20191007_py/snp.rds")
colnames(snp) = c("chrom", "pos", "ref", "alt", "strand", "sample", "GT", "AD", "gene_id", "transcript_id", "str_CBD")

colData = readRDS("~/czlabwork/vcfczold/20191003_R/colgroup.rds")
colnames(colData) = c("group", "sample")
snp = left_join(snp, colData, by = "sample")
saveRDS(snp, "~/czlabwork/vcfczold/20191007_py/snp_diff_all.rds")

############################# section 3.  SNP_CCGG_UV_MUTATIONS ONLY
snp$chrom_pos = paste(snp$chrom, snp$pos, sep = "_")

#####1. clean and filter df at 3 levels

# filtering out snp mapped to overlapping genes in the opposite strand
snp = snp %>% group_by(chrom_pos) %>% filter(n_distinct(strand) == 1) 
# filter out snp mapped to multiple transcripts
snp = snp[!duplicated(snp$chrom_pos),] 
#!!! reorder df according to same sample + position ascending order
snp <- snp[order(snp$sample, snp$pos),] 


### 2. seperating C-T & G-A group respectively, incase C-T followed by G-A are counted as CC-TT/GG-AA count
snp_c = snp[snp$ref == "C", ]
snp_g = snp[snp$ref == "G", ]

snp_c$diff = c(diff(snp_c$pos),0)
snp_g$diff = c(diff(snp_g$pos),0)

# only extract C-T/G-A followed by C-T/G-A -------> UV induced CC-TT/GG-AA signatures
snp_g_diff1 = snp_g[snp_g$diff==1,]
snp_c_diff1 = snp_c[snp_c$diff==1,]

### 3. combind CC and GG togather and add gene_name
snp_UV_dint = rbind(snp_c_diff1, snp_g_diff1)
gtf1 = readRDS("~/czlabwork/vcfczold/20191007_py/gtf1.rds")
gtf1 = gtf1 %>% select(gene_id, gene_name)
gtf1 = gtf1[!duplicated(gtf1),]
snp_UV_dint = left_join(snp_UV_dint, gtf1, by= "gene_id")


#### 4. seperate AD into 2 columns sep by ","
library(tidyverse)
snp_UV_dint = separate(snp_UV_dint, AD, c("AD_ref", "AD_alt"), sep = ",")

saveRDS(snp_UV_dint, "~/czlabwork/vcfczold/20191007_py/snp_UV_dint.rds")
