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
library(reshape2)

#1.
wt_nuc = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0052-RNA-G12V-wt-stimulate/outs/filtered_feature_bc_matrix/')
pep_nuc = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0053-RNA-G12V-pep-stimulate/outs/filtered_feature_bc_matrix/')
wt_cyt = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0057-RNA-QV1G12V-WT-cyt/outs/filtered_feature_bc_matrix/')
pep_cyt = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0058-RNA-QV1G12V-pep-cyt/outs/filtered_feature_bc_matrix/')
wt_fb = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0059-RNA-QV1G12V-wt-cyt-batch2/outs/filtered_feature_bc_matrix/')
pep_fb = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0060-RNA-QV1G12V-pep-cyt-batch2/outs/filtered_feature_bc_matrix/')

# Set up control object
####1.
data = wt_nuc
data = CreateSeuratObject(counts = data, project = '')
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
data <- data %>% NormalizeData(.,verbose = FALSE) %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(., features = rownames(.)) %>% RunPCA(., features = VariableFeatures(object = .)) %>%
  FindNeighbors(., dims = 1:10) %>% FindClusters(., resolution = .5) %>% 
  RunTSNE(., dims = 1:10) %>% RunUMAP(., dims = 1:10) 
wt_nuc = data

data = pep_nuc
data = CreateSeuratObject(counts = data, project = '')
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
data <- data %>% NormalizeData(.,verbose = FALSE) %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(., features = rownames(.)) %>% RunPCA(., features = VariableFeatures(object = .)) %>%
  FindNeighbors(., dims = 1:10) %>% FindClusters(., resolution = .5) %>% 
  RunTSNE(., dims = 1:10) %>% RunUMAP(., dims = 1:10) 
pep_nuc = data

data = wt_cyt
data = CreateSeuratObject(counts = data, project = '')
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
data <- data %>% NormalizeData(.,verbose = FALSE) %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(., features = rownames(.)) %>% RunPCA(., features = VariableFeatures(object = .)) %>%
  FindNeighbors(., dims = 1:10) %>% FindClusters(., resolution = .5) %>% 
  RunTSNE(., dims = 1:10) %>% RunUMAP(., dims = 1:10) 
wt_cyt = data

data = pep_cyt
data = CreateSeuratObject(counts = data, project = '')
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
data <- data %>% NormalizeData(.,verbose = FALSE) %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(., features = rownames(.)) %>% RunPCA(., features = VariableFeatures(object = .)) %>%
  FindNeighbors(., dims = 1:10) %>% FindClusters(., resolution = .5) %>% 
  RunTSNE(., dims = 1:10) %>% RunUMAP(., dims = 1:10) 
pep_cyt = data

data = wt_fb
data = CreateSeuratObject(counts = data, project = '')
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
data <- data %>% NormalizeData(.,verbose = FALSE) %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(., features = rownames(.)) %>% RunPCA(., features = VariableFeatures(object = .)) %>%
  FindNeighbors(., dims = 1:10) %>% FindClusters(., resolution = .5) %>% 
  RunTSNE(., dims = 1:10) %>% RunUMAP(., dims = 1:10) 
wt_fb = data

data = pep_fb
data = CreateSeuratObject(counts = data, project = '')
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
data <- data %>% NormalizeData(.,verbose = FALSE) %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(., features = rownames(.)) %>% RunPCA(., features = VariableFeatures(object = .)) %>%
  FindNeighbors(., dims = 1:10) %>% FindClusters(., resolution = .5) %>% 
  RunTSNE(., dims = 1:10) %>% RunUMAP(., dims = 1:10) 
pep_fb = data

#####5.
wt_nuc$group='wt_nuc'
pep_nuc$group = 'pep_nuc'
wt_cyt$group = 'wt_cyt'
pep_cyt$group = 'pep_cyt'
wt_fb$group = 'wt_fb'
pep_fb$group = 'pep_fb'

#######################################20210628 UPDATE TNF,IFNG,GZMB,LAG3 UMI in cutoff cells
data = wt_nuc %>% subset(nFeature_RNA>=3000)
data = pep_nuc %>% subset(nFeature_RNA>=3000)
data = wt_cyt %>% subset(nFeature_RNA>=3000)
data = pep_cyt %>% subset(nFeature_RNA>=3000)
data = wt_fb %>% subset(nFeature_RNA>=4000)
data = pep_fb %>% subset(nFeature_RNA>=4000)

feature=c('TNF','IFNG','GZMB','LAG3')
metadata = cbind(data@meta.data,
                 as.matrix(GetAssayData(data,assay = 'RNA',slot = 'counts'))[feature,] 
                 %>% t(),FetchData(data,feature),Embeddings(data,'umap'))
colnames(metadata)[12:15]=paste(colnames(metadata)[12:15],'logNorm',sep = '_')
##############################################################################
p1 = ggplot(metadata, aes(group, (TNF),fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Expression of TNF') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p2 = ggplot(metadata, aes(group, log2(TNF), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='log2(Expression of TNF)') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p1/p2

#1. by calculating TNF IFNG GZMB expression levels
quantile(metadata$TNF)
metadata <- metadata %>% mutate(tbl=cut(TNF,breaks = c(seq(0,45,by=5)), include.lowest = F))
df = rbind(data.frame(Var1='=0',Freq=sum(metadata$TNF==0)), data.frame(table(metadata$tbl)))%>% 
  dplyr::rename(TNF_UMIs='Var1',Cell_Nbr='Freq')
grid.table(df)
###########################################################
p1 = ggplot(metadata, aes(group, (IFNG), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Expression of IFNG') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p2 = ggplot(metadata, aes(group, log2(IFNG), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='log2(Expression of IFNG)') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p1/p2

#1. by calculating TNF IFNG GZMB expression levels
quantile(metadata$IFNG)
metadata <- metadata %>% mutate(tbl=cut(IFNG,breaks = c(seq(0,50,by=5),seq(100,1000,by=200)), include.lowest = F))
df = rbind(data.frame(Var1='=0',Freq=sum(metadata$IFNG==0)), data.frame(table(metadata$tbl)))%>% 
  dplyr::rename(IFNG_UMIs='Var1',Cell_Nbr='Freq')
grid.table(df)
######################################################
p1 = ggplot(metadata, aes(group, (GZMB), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Expression of GZMB') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p2 = ggplot(metadata, aes(group, log2(GZMB), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='log2(Expression of GZMB)') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p1/p2

#1. by calculating TNF IFNG GZMB expression levels
quantile(metadata$GZMB)
metadata <- metadata %>% mutate(tbl=cut(GZMB,breaks = c(seq(0,300,by=30),seq(500,2000,by=500)), include.lowest = F))
df = rbind(data.frame(Var1='=0',Freq=sum(metadata$GZMB==0)), data.frame(table(metadata$tbl)))%>% 
  dplyr::rename(GZMB_UMIs='Var1',Cell_Nbr='Freq')
grid.table(df)
#############################################################
p1 = ggplot(metadata, aes(group, (LAG3), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Expression of LAG3') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p2 = ggplot(metadata, aes(group, log2(LAG3), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='log2(Expression of LAG3)') + 
  scale_fill_manual(values = c('#56B4E9','gold'))
p1/p2
#1. by calculating TNF IFNG GZMB expression levels
quantile(metadata$LAG3)
metadata <- metadata %>% mutate(tbl=cut(LAG3,breaks = c(seq(0,30,by=3),seq(50,120,by=10)), include.lowest = F))
df = rbind(data.frame(Var1='=0',Freq=sum(metadata$LAG3==0)), data.frame(table(metadata$tbl)))%>% 
  dplyr::rename(LAG3_UMIs='Var1',Cell_Nbr='Freq')
grid.table(df)

################################20210628 CD4 and CD8 UMI in cutoff cells in 6 samples
data = wt_nuc %>% subset(nFeature_RNA>=3000)
data = pep_nuc %>% subset(nFeature_RNA>=3000)
data = wt_cyt %>% subset(nFeature_RNA>=3000)
data = pep_cyt %>% subset(nFeature_RNA>=3000)
data = wt_fb %>% subset(nFeature_RNA>=4000)
data = pep_fb %>% subset(nFeature_RNA>=4000)

feature=c('CD3D','CD3E','CD4','CD8A','CD8B')
metadata = cbind(data@meta.data,
                 as.matrix(GetAssayData(data,assay = 'RNA',slot = 'counts'))[feature,]%>% t(),
                 Embeddings(data,'umap')) %>% select(7:12)
#wide to long table
library(reshape2)
# Specify id.vars: the variables to keep but not split apart on
df = melt(metadata, id.vars=c("group")) %>% mutate(sample_gene = paste(group,variable,sep = '_'))

#############################################################
p1 = ggplot(df, aes(sample_gene, (value), fill=variable)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Expression of CD3/4/8') 
p2 = ggplot(df, aes(sample_gene, log2(value), fill=variable)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='log2(Expression of CD3/4/8)') 
p1/p2
#1. by calculating CD3D,CD3E,CD4,CD8A,CD8B UMIA expression levels
metadata$CD4.positive = ifelse(metadata$CD4>0,'CD4>0','CD4=0')
metadata$CD8A.positive = ifelse(metadata$CD8A>0,'CD8A>0','CD8A=0')
metadata$CD8B.positive = ifelse(metadata$CD8B>0,'CD8B>0','CD8B=0')
metadata$CD8AB.positive = ifelse(metadata$CD8A>0 & metadata$CD8B>0, 'CD8A+CD8B+','non_double+')
df = rbind(data.frame(table(metadata$CD4.positive,metadata$group)),
           data.frame(table(metadata$CD8A.positive,metadata$group)),
           data.frame(table(metadata$CD8B.positive,metadata$group)),
           data.frame(table(metadata$CD8AB.positive,metadata$group)))
colnames(df) = c('group','sample','Cell_Nbr')
df$pct = round(df$Cell_Nbr/nrow(metadata)*100,2)
grid.table(df)

################################20210630 update CD4+,CD8+,CD4+CD8+,CD4-CD8- UMI in cutoff cells in 6 samples
data = wt_nuc %>% subset(nFeature_RNA>=3000)
data = pep_nuc %>% subset(nFeature_RNA>=3000)
data = wt_cyt %>% subset(nFeature_RNA>=3000)
data = pep_cyt %>% subset(nFeature_RNA>=3000)
data = wt_fb %>% subset(nFeature_RNA>=4000)
data = pep_fb %>% subset(nFeature_RNA>=4000)

feature=c('CD3D','CD3E','CD8A','CD8B','CD4')
metadata = cbind(data@meta.data,
                 as.matrix(GetAssayData(data,assay = 'RNA',slot = 'counts'))[feature,]%>% t(),
                 Embeddings(data,'umap')) %>% select(7:14)
####################
#1.featureplot for each gene expression level
p1=FeaturePlot(data,features = feature)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p2=DimPlot(data,pt.size = .0001,split.by = 'group')&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7),
        legend.key.size = unit(0.1,'cm'))
p1+p2
#2. calculating single and double positive , double negative CD4,CD8 cells
metadata$CD4.positive = ifelse(metadata$CD4>0,'CD4+','CD4-')
metadata$CD8.positive = ifelse(metadata$CD8A>0 | metadata$CD8B>0,'CD8+','CD8-')
metadata$CD4.CD8.DP = ifelse(metadata$CD4.positive=='CD4+' & metadata$CD8.positive=='CD8+', 
                             'CD4+CD8+','not_double+')
metadata$CD4.CD8.DN = ifelse(metadata$CD4==0 & metadata$CD8A==0 & metadata$CD8B==0, 
                             'CD4-CD8-','not_double-')
df = rbind(data.frame(table(metadata$CD4.positive,metadata$group)),
           data.frame(table(metadata$CD8.positive,metadata$group)),
           data.frame(table(metadata$CD4.CD8.DP,metadata$group)),
           data.frame(table(metadata$CD4.CD8.DN,metadata$group)))
colnames(df) = c('group','sample','Cell_Nbr')
df$pct = round(df$Cell_Nbr/nrow(metadata)*100,2)
df1 = dcast(df,sample~group,value.var = 'Cell_Nbr')
grid.table(df1)
#3.show single double positive,negative cells on umap
data = AddMetaData(data,metadata = metadata)
data@meta.data$CD4.CD8.DP = factor(metadata$CD4.CD8.DP,levels = c('not_double+','CD4+CD8+'))
data@meta.data$CD4.CD8.DN = factor(metadata$CD4.CD8.DN,levels = c('not_double-','CD4-CD8-'))

plot.list = list()
for (i in unique(x = colnames(data@meta.data[15:18]))) {
  plot.list[[i]] <- DimPlot(data,pt.size = .0001,group.by= i)+
    scale_color_manual(values = c('grey','blue')) +
    theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
          plot.title = element_text(size=10),legend.text = element_text(size=7),
          legend.key.size = unit(0.1,'cm'))}
patchwork::wrap_plots(plots = plot.list, ncol = 2)

################################20210701 update CD8+ cells with TNF,IFNG,GZMB,LAG3 UMI in cutoff cells in 6 samples
data = wt_nuc %>% subset(nFeature_RNA>=3000)
data = pep_nuc %>% subset(nFeature_RNA>=3000)
data = wt_cyt %>% subset(nFeature_RNA>=3000)
data = pep_cyt %>% subset(nFeature_RNA>=3000)
data = wt_fb %>% subset(nFeature_RNA>=4000)
data = pep_fb %>% subset(nFeature_RNA>=4000)

feature=c('CD8A','CD8B','TNF','IFNG','GZMB','LAG3')
metadata = data.frame(as.matrix(GetAssayData(data,assay = 'RNA',slot = 'counts'))[feature,]%>% t(),
                      data$group)
####################
#2. calculating single and double positive , double negative CD4,CD8 cells
metadata$CD8.positive = ifelse(metadata$CD8A>0 | metadata$CD8B>0,'CD8+','CD8-')
metadata = metadata %>% mutate(group = paste(data.group,CD8.positive,sep = '.')) %>% 
  select(feature[3:6],group)
df = melt(metadata,id.vars = 'group')
#3. draw ggplots
p1 = ggplot(df, aes(group, value, fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Gene Expression (UMI count)') + 
  scale_fill_manual(values = c('#56B4E9','gold')) + facet_wrap(~variable,nrow = 1)
p2 = ggplot(df, aes(group, log(value), fill=group)) +
  geom_jitter(size=.05) + geom_violin(trim=F) + NoLegend() +
  geom_boxplot(fill='white',width=.1) + labs(x='',y='Gene Expression log(UMI count)') + 
  scale_fill_manual(values = c('#56B4E9','gold')) + facet_wrap(~variable,nrow = 1)
p1/p2
#4. calculate UMI by ranges
quantile(df$value)
df1 <- df %>% mutate(tbl=cut(value,breaks = c(seq(0,30,by=5),seq(100,1900,by=400)), include.lowest = F))
df2 = rbind(data.frame(table(df1$value==0,df1$variable,df1$group)),
            data.frame(table(df1$tbl,df1$variable,df1$group)))%>%filter(Var1 != FALSE)
df2$Var1 = gsub('TRUE','=0',df2$Var1)
df3 = dcast(df2, Var1~Var2+Var3, value.var = 'Freq') %>% rename('UMIs'='Var1')
df3$UMIs = factor(df3$UMIs, levels = unique(df2$Var1))
df3 = df3[order(df3$UMIs,levels=unique(df2$Var1)),]
grid.table(df3[,1:5])
grid.table(df3[,c(1,6:9)])










