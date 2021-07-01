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
wt_fb = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0059-RNA-QV1G12V-wt-cyt-batch2/outs/filtered_feature_bc_matrix/')
pep_fb = Read10X('/home/yao.zhan/analysis-CRv4/nextseq0060-RNA-QV1G12V-pep-cyt-batch2/outs/filtered_feature_bc_matrix/')

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

wt_fb$group = 'wt_fb'
pep_fb$group = 'pep_fb'

#filter gene < 4K cells from cutoff UMI analysis performed today 20210625
wt_fb = subset(wt_fb, nFeature_RNA >=4000)
pep_fb = subset(pep_fb, nFeature_RNA >=4000)
#################################################################################################
#### 1. seurat normalized data integration
#################################################################################################
ifnb.list = list(wt_fb,pep_fb)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- immune.combined %>% ScaleData(., verbose = FALSE) %>% RunPCA(., npcs = 30, verbose = FALSE) %>% 
  RunUMAP(., reduction = "pca", dims = 1:30) %>% FindNeighbors(., reduction = "pca", dims = 1:30) %>% 
  FindClusters(., resolution = 0.5)
# Visualization
p1 <-DimPlot(immune.combined, reduction = "umap", group.by = "group")
p2 <-DimPlot(immune.combined, reduction = "umap", label = TRUE)
p3 <-DimPlot(immune.combined, reduction = 'umap', split.by = 'group')
(p1+p2)/p3
saveRDS(immune.combined,'MultiomeSeq/fb.lognorm.integrate.rds')
#Identify differential expressed genes across conditions
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- "group"
DEG <- FindMarkers(immune.combined, ident.1 = "wt_fb", ident.2 = "pep_fb", verbose = FALSE) #,logfc.threshold = log(2),min.pct = 0.5)
write.csv(DEG,'MultiomeSeq/cytoplasmic_RNA_feature_barcode/fb_DEG_wt_pep.logNormIntegrate.csv')

#################### DE gene analysis and visualization
group_marker = DEG

group_marker$annotation = ifelse(group_marker$avg_log2FC>1 & group_marker$p_val_adj < .05 | 
                                   group_marker$avg_log2FC < -1 & group_marker$p_val_adj < .05, '2 Fold change','')
#1)
ggplot(group_marker, aes(avg_log2FC, -log10(p_val_adj), color=annotation)) + geom_point(size=.2) + 
  ggtitle('fb_DEG in QV1G12Vwt Vs pep ')
#2)
df = data.frame('fb_DEG_2X_up_PepVSWt' = sum(group_marker$avg_log2FC< -1),
                'fb_DEG_2X_down_PepVSWt' = sum(group_marker$avg_log2FC> 1))
grid.table(t(df))
#3)
DefaultAssay(immune.combined) <- 'RNA'
immune.combined <- immune.combined %>% NormalizeData(.) %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% ScaleData(., verbose = FALSE) %>% 
  RunPCA(., npcs = 30, verbose = FALSE) %>% RunUMAP(., reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(., reduction = "pca", dims = 1:30) %>% FindClusters(., resolution = 0.5)

DEG_wt_nuc = FetchData(object = subset(immune.combined,group %in% c('wt_nuc','pep_nuc')), 
                          slot = 'data', c(rownames(DEG),'group'))
ggplot()



#################################################################################################
#### 2. SCTransform
#################################################################################################
#1. Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
#   1)During normalization, also remove confounding sources of variation, for example, mitochondrial mapping percentage
#2. As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
#3. Run the PrepSCTIntegration() function prior to identifying anchors
#4. When running FindIntegrationAnchors(), and IntegrateData(), set the normalization.method parameter to the value SCT.
#5. When running sctransform-based workflows, including integration, do not run the ScaleData() function
wt_fb <- PercentageFeatureSet(wt_fb, pattern = "^MT-", col.name = "percent.mt")
pep_fb <- PercentageFeatureSet(pep_fb, pattern = '^MT-',col.name = 'percent.mt')
sct_wt_fb <- SCTransform(wt_fb, vars.to.regress = "percent.mt", verbose = FALSE)
sct_pep_fb <- SCTransform(pep_fb, vars.to.regress = "percent.mt", verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = list(sct_wt_fb,sct_pep_fb), 
                                      nfeatures = 6000)
data <- PrepSCTIntegration(object.list = list(sct_wt_fb,sct_pep_fb), anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = data, normalization.method = "SCT", anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- immune.combined.sct %>% RunPCA(., verbose = FALSE) %>% 
  RunUMAP(., reduction = "pca", dims = 1:30) %>% FindNeighbors(., reduction = "pca", dims = 1:30) %>% 
  FindClusters(., resolution = 0.5)

#step 7 visualization
p1 = DimPlot(immune.combined.sct, reduction = "umap", group.by = "group",pt.size = .0001)+
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.3,'cm')) + 
  guides(color = guide_legend(override.aes = list(size=2)))
p2 = DimPlot(immune.combined.sct, reduction = 'umap', label = T,pt.size = .0001)+
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.3,'cm')) + 
  guides(color = guide_legend(override.aes = list(size=2)))
p3 = DimPlot(immune.combined.sct, reduction = "umap", split.by = "group",pt.size = .0001)+
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.3,'cm')) + 
  guides(color = guide_legend(override.aes = list(size=2)))
(p1+p2)/p3
saveRDS(immune.combined.sct,'MultiomeSeq/cytoplasmic_RNA_feature_barcode/fb_wt_pep.sct.integrated4.rds')
#DEgene analysis
DefaultAssay(immune.combined.sct) <- "SCT"
Idents(immune.combined.sct) <-'group'
SCT_deg1 = FindMarkers(immune.combined.sct, ident.1 = 'wt_fb', ident.2 = 'pep_fb', assay = 'SCT',slot = 'data')
write.csv(SCT_deg1,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_nuc_wt_pep_SCTintegrate.csv')
###########################################










