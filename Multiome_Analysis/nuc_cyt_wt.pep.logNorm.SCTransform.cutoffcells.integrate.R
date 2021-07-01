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

# Set up control object
wt_nuc <- CreateSeuratObject(counts = wt_nuc, min.cells = 5)
wt_nuc@meta.data$group <- "wt_nuc"
pep_nuc <- CreateSeuratObject(counts = pep_nuc, min.cells = 5)
pep_nuc@meta.data$group <- "pep_nuc"
wt_cyt <- CreateSeuratObject(counts = wt_cyt, min.cells = 5)
wt_cyt@meta.data$group <- "wt_cyt"
pep_cyt <- CreateSeuratObject(counts = pep_cyt, min.cells = 5)
pep_cyt@meta.data$group <- "pep_cyt"

#filter gene < 3K cells from cutoff UMI analysis performed today 20210624
wt_nuc = subset(wt_nuc, nFeature_RNA >=3000)
pep_nuc = subset(pep_nuc, nFeature_RNA >=3000)
wt_cyt = subset(wt_cyt, nFeature_RNA >=3000)
pep_cyt = subset(pep_cyt, nFeature_RNA >=3000)
#################################################################################################
#### 1. seurat normalized data integration
#################################################################################################
ifnb.list = list(wt_nuc,pep_nuc,wt_cyt,pep_cyt)
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
saveRDS(immune.combined,'MultiomeSeq//4sample.lognorm.integrate.rds')
#Identify differential expressed genes across conditions
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- "group"
DEG <- FindMarkers(immune.combined, ident.1 = "wt_nuc", ident.2 = "pep_nuc", verbose = FALSE) #,logfc.threshold = log(2),min.pct = 0.5)
DEG2 <- FindMarkers(immune.combined, ident.1 = 'wt_cyt', ident.2 = 'pep_cyt')
write.csv(DEG,'MultiomeSeq/multiome-nuclear_RNA/DEG_nuc_wt_pep.logNormIntegrate.csv')
write.csv(DEG2,'MultiomeSeq/multiome-nuclear_RNA/DEG_cyt_wt_pep.logNormIntegrate.csv')

DEG3 <- FindMarkers(immune.combined, ident.1 = "wt_nuc", ident.2 = "wt_cyt", verbose = FALSE) #,logfc.threshold = log(2),min.pct = 0.5)
DEG4 <- FindMarkers(immune.combined, ident.1 = 'pep_nuc', ident.2 = 'pep_cyt')
write.csv(DEG3,'MultiomeSeq/multiome-nuclear_RNA/DEG_nuc_cyt_wt.logNormIntegrate.csv')
write.csv(DEG4,'MultiomeSeq/multiome-nuclear_RNA/DEG_nuc_cyt_pep.logNormIntegrate.csv')

#################### DE gene analysis and visualization
group_marker = DEG
group_marker = DEG2
group_marker = DEG3
group_marker = DEG4

group_marker$annotation = ifelse(group_marker$avg_log2FC>1 & group_marker$p_val_adj < .05 | 
                                   group_marker$avg_log2FC < -1 & group_marker$p_val_adj < .05, '2 Fold change','')
#1)
ggplot(group_marker, aes(avg_log2FC, -log10(p_val_adj), color=annotation)) + geom_point(size=.2) + 
  ggtitle('Cyt_DEG in QV1G12Vwt Vs pep ')
#2)
df = data.frame('Cyt_DEG_2X_up_PepVSWt' = sum(group_marker$avg_log2FC< -1),
                'Cyt_DEG_2X_down_PepVSWt' = sum(group_marker$avg_log2FC> 1))
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
wt_nuc <- PercentageFeatureSet(wt_nuc, pattern = "^MT-", col.name = "percent.mt")
pep_nuc <- PercentageFeatureSet(pep_nuc, pattern = '^MT-',col.name = 'percent.mt')
wt_cyt <- PercentageFeatureSet(wt_cyt, pattern = "^MT-", col.name = "percent.mt")
pep_cyt <- PercentageFeatureSet(pep_cyt, pattern = '^MT-',col.name = 'percent.mt')

sct_wt_nuc <- SCTransform(wt_nuc, vars.to.regress = "percent.mt", verbose = FALSE)
sct_pep_nuc <- SCTransform(pep_nuc, vars.to.regress = "percent.mt", verbose = FALSE)
sct_wt_cyt <- SCTransform(wt_cyt, vars.to.regress = "percent.mt", verbose = FALSE)
sct_pep_cyt <- SCTransform(pep_cyt, vars.to.regress = "percent.mt", verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = list(sct_wt_nuc,sct_pep_nuc,sct_wt_cyt,sct_pep_cyt), 
                                      nfeatures = 6000)
data <- PrepSCTIntegration(object.list = list(sct_wt_nuc,sct_pep_nuc,sct_wt_cyt,sct_pep_cyt), anchor.features = features)
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
saveRDS(immune.combined.sct,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/nuc_cyt_wt_pep.sct.integrated4.rds')
#DEgene analysis
DefaultAssay(immune.combined.sct) <- "SCT"
Idents(immune.combined.sct) <-'group'
SCT_deg1 = FindMarkers(immune.combined.sct, ident.1 = 'wt_nuc', ident.2 = 'pep_nuc', assay = 'SCT',slot = 'data')
SCT_deg2 = FindMarkers(immune.combined.sct, ident.1 = 'wt_cyt', ident.2 = 'pep_cyt', assay = 'SCT',slot = 'data')
SCT_deg3 = FindMarkers(immune.combined.sct, ident.1 = 'wt_nuc', ident.2 = 'wt_cyt', assay = 'SCT',slot = 'data')
SCT_deg4 = FindMarkers(immune.combined.sct, ident.1 = 'pep_nuc', ident.2 = 'pep_cyt', assay = 'SCT',slot = 'data')
write.csv(SCT_deg1,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_nuc_wt_pep_SCTintegrate.csv')
write.csv(SCT_deg2,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_cyt_wt_pep_SCTintegrate.csv')
write.csv(SCT_deg3,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_nuc_cyt_wt_SCTintegrate.csv')
write.csv(SCT_deg4,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_nuc_cyt_pep_SCTintegrate.csv')

###########################################










