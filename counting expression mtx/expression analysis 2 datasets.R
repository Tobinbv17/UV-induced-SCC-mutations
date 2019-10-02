load("~/Desktop/czlabwork/Expression/Expression.Rproj")
H7Gcount = readRDS('./H7GJWDRXX/counts.rds')
H7Gqc = readRDS('./H7GJWDRXX/QC.rds')
H7GcolData = readRDS('./H7GJWDRXX/coldata.rds')
H7Gsce = readRDS('./H7GJWDRXX/sce.rds')
Jstcount = readRDS('./GSE67602_Joost_et_al/counts.rds')
JstcolData = readRDS('./GSE67602_Joost_et_al/coldata.rds')

#####----- 1. combine 2 lanes get the mean
name = data.frame(matrix(nrow=length(H7Gqc$sample_id), ncol=1))
for (i in 1:length(H7Gqc$sample_id)) {name[i,] = strsplit(H7Gqc$sample_id, ".", fixed = T)[[i]][1]}
library(dplyr)
H7Gqc = bind_cols(name = name, H7Gqc)
H7Gqc = H7Gqc %>% group_by(H7Gqc[,1]) %>% summarise_each(mean, select = -c(1,2)) %>% mutate_if(is.numeric, round)

#####----- 2. QC for count matrix
#1.genes map to chrM
library(rtracklayer)
library(rtracklayer)
GTF = as.data.frame(rtracklayer::import('/Users/yao/Desktop/czlabwork/Gannotation/gencode.vM20.annotation.gtf'))
GTF = filter(GTF, seqnames == 'chrM' & type =='gene') %>% select(gene_id, gene_name) # filter out chrM genes

H7GgeneM = H7Gcount[GTF$gene_id, ]
H7GMcount = data.frame(row.names=colnames(H7Gcount), colSums(H7Gcount[GTF$gene_id,]))
H7GMcount = H7GMcount[match(H7Gqc$`H7Gqc[, 1]`, rownames(H7GMcount)),]
H7Gqc = bind_cols(H7Gqc, N_chrM = round(H7GMcount/H7Gqc$totalReads*100))

JstgeneM = Jstcount[GTF$gene_name, ]
Jstqc = data.frame(colnames(Jstcount), round(JstMcount/colSums(Jstcount)*100))
Jstercc = Jstcount[grepl("ERCC-", rownames(Jstcount)),]
Jstqc = bind_cols(Jstqc, N_ERCC= round(colSums(Jstercc)/Jstqc$colSums.Jstcount.*100))
Jstqc = bind_cols(Jstqc, expgene = round(colSums(Jstcount>0)/nrow(Jstcount)*100), th1 = colSums(Jstcount>1), th5 = colSums(Jstcount>5), th10 = colSums(Jstcount>10), totalReads = colSums(Jstcount))
colnames(Jstqc) = c("sample_id", "N_chrM", "N_ERCC", "N_gene_exp", "th1", "th5","th10","totalReads")






