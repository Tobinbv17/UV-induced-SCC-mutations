################################################################################
###1. Alignment workflow Using sctransform in Seurat --20210623
###In this vignette, we demonstrate how using sctransform based normalization enables recovering sharper biological distinction compared to log-normalization.
################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
#1.
wt_nuc = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0052-RNA-G12V-wt-stimulate/outs/filtered_feature_bc_matrix/')
pep_nuc = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0053-RNA-G12V-pep-stimulate/outs/filtered_feature_bc_matrix/')
wt_cyt = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0057-RNA-QV1G12V-WT-cyt/outs/filtered_feature_bc_matrix/')
pep_cyt = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0058-RNA-QV1G12V-pep-cyt/outs/filtered_feature_bc_matrix/')
wt_fb = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0059-RNA-QV1G12V-wt-cyt-batch2/outs/filtered_feature_bc_matrix/')
pep_fb = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0060-RNA-QV1G12V-pep-cyt-batch2/outs/filtered_feature_bc_matrix/')

# Set up control object
data = wt_nuc
data = pep_nuc
data = wt_cyt
data = pep_cyt
data = wt_fb
data = pep_fb

data = CreateSeuratObject(counts = data, project = '')
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')

ribo = grep('^RPS|^RPL',rownames(data))
data_ribo=data[ribo,]
ribo = grep('^RPS6[A-Z]', rownames(data_ribo))
data_ribo=data_ribo[-ribo]
ribo = grep('[-]',rownames(data_ribo))
data_ribo=data_ribo[-ribo]
ribo_num = Matrix::colSums(data_ribo@assays$RNA@counts>0)
data$ribo_num=ribo_num

wt_nuc = data
pep_nuc = data
wt_cyt = data
pep_cyt = data
wt_fb = data
pep_fb = data


#################################################################################################
wt_nuc$group='wt_nuc'
pep_nuc$group = 'pep_nuc'
wt_cyt$group = 'wt_cyt'
pep_cyt$group = 'pep_cyt'
wt_fb$group = 'wt_fb'
pep_fb$group = 'pep_fb'

metadata = rbind(wt_nuc@meta.data, pep_nuc@meta.data, wt_cyt@meta.data, pep_cyt@meta.data)

metadata = rbind(wt_fb@meta.data, pep_fb@meta.data)

metadata = subset(metadata, metadata$nFeature_RNA >= 3000)
metadata = subset(metadata, metadata$nFeature_RNA >= 4000)

p1 <- ggplot(metadata, aes(group, log(nCount_RNA), fill=group)) + 
  geom_jitter(size=.0001) + geom_violin(trim = F) + NoLegend() +
  geom_boxplot(fill='white', width=.1) + labs(x='', y='logUMI') + theme(axis.text.x = element_text(size = 7)) +
  scale_y_continuous(breaks=seq(6,12,1))

p2 <- ggplot(metadata, aes(group, percent.mt, fill=group)) + 
  geom_jitter(size=.0001) + geom_violin(trim = F) + NoLegend() +
  geom_boxplot(fill='white', width=.1) + labs(x='', y='percent.mt') +
  scale_y_continuous(breaks=seq(0,100,5))+ theme(axis.text.x = element_text(size = 7))

p3 <- ggplot(metadata, aes(group, nFeature_RNA, fill=group)) +
  geom_jitter(size=.0001) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Gene') + theme(axis.text.x = element_text(size = 7))+
  scale_y_continuous(breaks=seq(0,10000,2000))

p4 <- ggplot(metadata, aes(group, ribo_num, fill=group)) +
  geom_jitter(size=.0001) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='') + theme(axis.text.x = element_text(size = 7)) +
  scale_y_continuous(breaks=seq(0,100,5))
p1+p2+p3+p4

#######################################################
metadata = rbind(wt_nuc@meta.data, pep_nuc@meta.data, wt_cyt@meta.data, pep_cyt@meta.data)
metadata = rbind(wt_fb@meta.data, pep_fb@meta.data)

metadata$noFilter = 'noFilter'
metadata$Gene2.5K = ifelse(metadata$nFeature_RNA >=2500, '>=2500','<2500')
metadata$Gene3K = ifelse(metadata$nFeature_RNA >=3000, '>=3000','<3000')
metadata$Gene4K = ifelse(metadata$nFeature_RNA >=4000, '>=4000','<4000')

df = rbind(data.frame(table(metadata$noFilter,metadata$group)), 
           data.frame(table(metadata$Gene2.5K, metadata$group)),
           data.frame(table(metadata$Gene3K, metadata$group)))

df = rbind(data.frame(table(metadata$noFilter,metadata$group)), 
           data.frame(table(metadata$Gene3K, metadata$group)),
           data.frame(table(metadata$Gene4K, metadata$group)))

df1 = df %>% group_by(Var2, Var1) %>% summarise(count = Freq)
colnames(df1) = c('sample','filter','count')

library(reshape)
df2 = cast(df1, sample~filter)
write.table(df2,'MultiomeSeq/numberbyfilter_2featurebarcode.csv',quote = F, sep = ',',row.names = F)




