################################################################################
###1. Using Seurat with multimodal data --20210628
################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)

#1.# Load in the cellranger count for ADT assay matrix (input both RNA + ADT fastqs in libraries.csv)
wt_fb_citeseq = Read10X('/home/yao.zhan/cellranger-v5_v6/feature_barcode/nextseq0061-fb-wt-0628/outs/filtered_feature_bc_matrix/')
rownames(x = wt_fb_citeseq[['Antibody Capture']]) <- gsub(pattern = ".1", replacement = "",x = rownames(x = wt_fb_citeseq[["Antibody Capture"]]))

pep_fb_citeseq = Read10X('/home/yao.zhan/cellranger-v5_v6/feature_barcode/nextseq0061-fb-pep-0628/outs/filtered_feature_bc_matrix/')
rownames(x = pep_fb_citeseq[['Antibody Capture']]) <- gsub(pattern = ".1", replacement = "",x = rownames(x = pep_fb_citeseq[["Antibody Capture"]]))

#create seurat object for RNA UMI assay and append protein ADT assay in the same object
#sample1
wt_fb = CreateSeuratObject(counts = wt_fb_citeseq$`Gene Expression`, project = '')
Assays(wt_fb)
wt_fb[['Protein']] <- CreateAssayObject(counts = wt_fb_citeseq$`Antibody Capture`)
Assays(wt_fb)
rownames(wt_fb[['Protein']])
wt_fb$group = 'wt_fb'

#sample2
pep_fb = CreateSeuratObject(counts = pep_fb_citeseq$`Gene Expression`, project = '')
Assays(pep_fb)
pep_fb[['Protein']] <- CreateAssayObject(counts = pep_fb_citeseq$`Antibody Capture`)
Assays(pep_fb)
rownames(pep_fb[['Protein']])
pep_fb$group = 'pep_fb'

all.equal(colnames(wt_fb@assays$RNA), colnames(wt_fb@assays$Protein))
all.equal(colnames(pep_fb@assays$RNA), colnames(pep_fb@assays$Protein))

#2.Cluster cells on the basis of their scRNA-seq profiles
DefaultAssay(wt_fb) <- "RNA"
DefaultAssay(wt_fb)
wt_fb <- wt_fb %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA(.) %>% 
  FindNeighbors(., dims = 1:30) %>% FindClusters(., resolution = 0.8) %>% RunUMAP(., dims = 1:30)
DimPlot(wt_fb, label = TRUE)+ggtitle('wt_fb')

DefaultAssay(pep_fb) <- "RNA"
DefaultAssay(pep_fb)
pep_fb <- pep_fb %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA(.) %>% 
  FindNeighbors(., dims = 1:30) %>% FindClusters(., resolution = 0.8) %>% RunUMAP(., dims = 1:30)
DimPlot(pep_fb, label = TRUE)+ggtitle('pep_fb')

#filter gene < 4K cells from cutoff UMI analysis performed today 20210625
wt_fb = subset(wt_fb, nFeature_RNA >=4000)
pep_fb = subset(pep_fb, nFeature_RNA >=4000)

#3. Visualize multiple modalities side-by-side
# Normalize ADT data,
wt_fb <- wt_fb %>% NormalizeData(., normalization.method = "CLR", margin = 2, assay = "Protein") %>%
  FindVariableFeatures(.,assay = 'Protein') %>% ScaleData(.,assay = 'Protein')
pep_fb <- pep_fb %>% NormalizeData(., normalization.method = "CLR", margin = 2, assay = "Protein") %>%
  FindVariableFeatures(.,assay = 'Protein') %>% ScaleData(.,assay = 'Protein')

# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
Key(wt_fb[["RNA"]])
Key(wt_fb[["Protein"]])
gene = c('CD3D','CD4','CD8A','CD28','CD274','CTLA4')
atd=paste('protein_',gene,'.1',sep='')
rna=paste('rna_',gene,sep='')
p1 <- FeaturePlot(wt_fb, atd[1:3], cols = c("snow2", "darkgreen"),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p2 <- FeaturePlot(wt_fb, rna[1:3],ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p3 <- FeaturePlot(wt_fb, atd[4:6], cols = c("snow2", "darkgreen"),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p4 <- FeaturePlot(wt_fb, rna[4:6],ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p1|p2|p3|p4

p1 <- FeaturePlot(pep_fb, atd[1:3], cols = c("snow2", "darkgreen"),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p2 <- FeaturePlot(pep_fb, rna[1:3],ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p3 <- FeaturePlot(pep_fb, atd[4:6], cols = c("snow2", "darkgreen"),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p4 <- FeaturePlot(pep_fb, rna[4:6],ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p1|p2|p3|p4

p1 <- FeaturePlot(wt_fb, c('protein_CD19.1','protein_NCAM1.1'), cols = c("snow2", "darkgreen"),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p2 <- FeaturePlot(wt_fb, c('rna_CD19','rna_NCAM1'),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p3 <- FeaturePlot(pep_fb, c('protein_CD19.1','protein_NCAM1.1'), cols = c("snow2", "darkgreen"),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p4 <- FeaturePlot(pep_fb, c('rna_CD19','rna_NCAM1'),ncol = 1)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
(p1|p2)/(p3|p4)

#4.Identify cell surface markers for scRNA-seq clusters
#We can leverage our paired CITE-seq measurements to help annotate clusters derived from scRNA-seq, and to identify both protein and RNA markers.
#4.1 using known B cell and T cell markers CD19,CD4,CD8,CD3D
p1 = DimPlot(wt_fb,pt.size = .0001,label = T)+ggtitle('wt_fb')&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),axis.title.x = element_blank(),
        plot.title = element_text(size=10),legend.text = element_text(size=7),
        legend.key.size = unit(0.3,'cm'))
p3 = VlnPlot(wt_fb, c(atd[1:3],'protein_CD19.1','protein_NCAM1.1'),ncol=2,pt.size = .0001)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),axis.title.x = element_blank(),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
(p3+p1)

p2 = DimPlot(pep_fb,pt.size = .0001,label = T)+ggtitle('pep_fb')&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),axis.title.x = element_blank(),
        plot.title = element_text(size=10),legend.text = element_text(size=7),
        legend.key.size = unit(0.3,'cm'))
p4 = VlnPlot(pep_fb, c(atd[1:3],'protein_CD19.1','protein_NCAM1.1'),ncol=2,pt.size = .0001)&
  theme(axis.text = element_text(size=7), axis.title = element_text(size=8),axis.title.x = element_blank(),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
(p4+p2)
#4.2# we can also identify alternative protein and RNA markers for this cluster through
# differential expression
adt_markers_wt <- FindAllMarkers(wt_fb, min.pct = 0.25, assay = "Protein")
rna_markers_wt <- FindAllMarkers(wt_fb, min.pct = .25, assay = "RNA")
write.csv(adt_markers_wt,'MultiomeSeq/cytoplasmic_RNA_feature_barcode/adt_markers_wt.csv')
write.csv(rna_markers_wt,'MultiomeSeq/cytoplasmic_RNA_feature_barcode/rna_markers_wt.csv')

#protein marker in each cluster
adt_markers_wt %>% filter(p_val_adj<.05, avg_log2FC>1) %>% group_by(cluster)-> top10
DoHeatmap(wt_fb, features = top10$gene,assay = 'Protein',size=3,angle=90,draw.lines = T) + NoLegend() +
  scale_fill_gradientn(colors = c("blue", "black", "yellow"),na.value = 'white') + 
  theme(text = element_text(size=8))
#gene marker in each cluster
rna_markers_wt %>% filter(p_val_adj<.05, avg_log2FC>1) %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(wt_fb, features = top10$gene,assay = 'RNA',size=3,angle=90,draw.lines = T) + NoLegend() +
  scale_fill_gradientn(colors = c("blue", "white", "red"),na.value = 'white') + 
  theme(text = element_text(size=8))
##################################
adt_markers_pep <- FindAllMarkers(pep_fb, min.pct = 0.25, assay = "Protein")
rna_markers_pep <- FindAllMarkers(pep_fb, min.pct = .25, assay = "RNA")
write.csv(adt_markers_pep,'MultiomeSeq/cytoplasmic_RNA_feature_barcode/adt_markers_pep.csv')
write.csv(rna_markers_pep,'MultiomeSeq/cytoplasmic_RNA_feature_barcode/rna_markers_pep.csv')

#protein marker in each cluster
adt_markers_pep %>% filter(p_val_adj<.05, avg_log2FC>.6) %>% group_by(cluster) -> top10
DoHeatmap(pep_fb, features = top10$gene,assay = 'Protein',size=3,angle=90,draw.lines = T) + NoLegend() +
  scale_fill_gradientn(colors = c("blue", "black", "yellow"),na.value = 'white') + 
  theme(text = element_text(size=8))
#gene marker in each cluster
rna_markers_pep %>% filter(p_val_adj<.05, avg_log2FC>1) %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pep_fb, features = top10$gene,assay = 'RNA',size=3,angle=90,draw.lines = T) + NoLegend() +
  scale_fill_gradientn(colors = c("blue", "white", "red"),na.value = 'white') + 
  theme(text = element_text(size=7.5))



