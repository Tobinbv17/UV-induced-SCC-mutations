################################################################################
###1. Alignment workflow Using sctransform in Seurat --20210428
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
QV1G12V_wt = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0052-RNA-G12V-wt-stimulate/outs/filtered_feature_bc_matrix/')
QV1G12V_pep = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0053-RNA-G12V-pep-stimulate/outs/filtered_feature_bc_matrix/')

QV1G12V_wt = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0052-RNA-G12V-wt-stimulate-force//outs/filtered_feature_bc_matrix/')
QV1G12V_pep = Read10X('/home/yao.zhan/cellranger-v5_v6/Multiome/nextseq0053-RNA-G12V-pep-stimulate-force//outs/filtered_feature_bc_matrix/')

# Set up control object
QV1G12V_wt <- CreateSeuratObject(counts = QV1G12V_wt, min.cells = 5)
QV1G12V_wt@meta.data$group <- "QV1G12V_wt"
# Set up treatment1 object
QV1G12V_pep <- CreateSeuratObject(counts = QV1G12V_pep, min.cells = 5)
QV1G12V_pep@meta.data$group <- "QV1G12V_pep"

#################################################################################################
#### 1. seurat normalized data integration
#################################################################################################
bm280k.list <- lapply(X = list(QV1G12V_wt,QV1G12V_pep), FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)})

features <- SelectIntegrationFeatures(object.list = bm280k.list, nfeatures = 5000)
bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = bm280k.list, reduction = "rpca",dims = 1:50, anchor.features = features)
bm280k.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

bm280k.integrated <- ScaleData(bm280k.integrated, verbose = FALSE)
bm280k.integrated <- RunPCA(bm280k.integrated, verbose = FALSE)
bm280k.integrated <- RunUMAP(bm280k.integrated, dims = 1:50)

DimPlot(bm280k.integrated, group.by = "group")


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

# cluster cell count in each dataset
p1/p2
df = as.data.frame.matrix(table(immune.combined.sct$seurat_clusters, immune.combined.sct$group))
b = round(sweep(df,2,colSums(df),'/')*100,2)
df = cbind(df,b)
df = cbind(df, cluster=as.data.frame(table(immune.combined.sct$seurat_clusters))[,1])
colnames(df)[1:4]=c(colnames(df)[1:2],paste(colnames(df)[1:2],'%',sep = ''))
grid.table(df)

p1/DimPlot(subset(immune.combined.sct, clonotype %in% seq), reduction = 'umap',group.by ='clonotype', pt.size = .01)
data=subset(immune.combined.sct,clonotype%in%seq)@meta.data
df = as.data.frame.matrix(table(data$clonotype, data$group))
b = round(sweep(df,2,table(immune.combined.sct$group),'/')*100,2)
df = cbind(df,b)
df = cbind(df, clonotype=as.data.frame(table(data$clonotype))[,1])
colnames(df)[1:6]=c(colnames(df)[1:3],paste(colnames(df)[1:3],'%',sep = ''))
grid.table(df)

###############################################################################
###2. secondary analysis -- differential expression analysis in SCTransformed object
###############################################################################
#1. Exploring known cell type markers -- https://hbctraining.github.io/scRNA-seq/lessons/08_SC_clustering_quality_control.html
#With the cells clustered, we can explore the cell type identities by looking for known markers. The UMAP plot with clusters marked is shown, followed by the different cell types expected.
#For example if we were interested in exploring known immune cell markers, such as:
#CD14+ monocytes	CD14, LYZ
#FCGR3A+ monocytes	FCGR3A, MS4A7
#Conventional dendritic cells	FCER1A, CST3
#Plasmacytoid dendritic cells	IL3RA, GZMB, SERPINF1, ITM2C
#B cells	CD79A, MS4A1
#T cells	CD3D
#CD4+ T cells	CD3D, IL7R, CCR7
#CD8+ T cells	CD3D, CD8A
#NK cells	GNLY, NKG7
#Megakaryocytes	PPBP
#Erythrocytes	HBB, HBA2
#########################################################################################
# The results of sctransfrom are stored in the “SCT” assay
#pbmc[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA. 
#The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts. we also convert Pearson residuals back to ‘corrected’ UMI counts. You can interpret these as the UMI counts we would expect to observe if all cells were sequenced to the same depth.
#We store log-normalized versions of these corrected counts in pbmc[["SCT"]]@data, which are very helpful for visualization.
#You can use the corrected log-normalized counts for differential expression and integration. However, in principle, it would be most optimal to perform these calculations directly on the residuals (stored in the scale.data slot) themselves. This is not currently supported in Seurat v3, but will be soon.
##########################################################################################
#20210621 update
DefaultAssay(bm280k.integrated) <- "RNA"
# same cluster in 2 groups need to reach more than 50 cells = 0,1,4,7,8,10
Idents(bm280k.integrated) <- bm280k.integrated$group
group_marker = FindMarkers(bm280k.integrated, ident.1 = 'QV1G12V_wt', ident.2 = 'QV1G12V_pep', 
                           min.cells.group = 100, 
                           logfc.threshold = log(1.5),min.pct = 0.5)
group_marker = group_marker %>% filter(p_val_adj < .05) %>% arrange(avg_log2FC) %>% mutate(gene=rownames(.))

ggplot(group_marker, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(size=.2) + 
  ggtitle('DEG in QV1G12V-WT Vs pep')
write.csv(group_marker, 'MultiomeSeq/rpca.wt.pep.deg.forcedcell.csv')

#marker_use = group_marker %>% filter(pct.1>.1 & pct.2>.1)
p1 = FeaturePlot(immune.combined.sct, features = head(group_marker$gene,4), split.by = "group",
            cols = c("grey", "red"))
p2 = FeaturePlot(immune.combined.sct, features = tail(group_marker$gene,4), split.by = "group",
                 cols = c("grey", "red"))
p1|p2

group_marker$annotation = ifelse(group_marker$avg_log2FC>1 & group_marker$p_val_adj < .05 | 
                                   group_marker$avg_log2FC < -1 & group_marker$p_val_adj < .05, '2 Fold change','')
write.csv(group_marker, 'MultiomeSeq/QV1G12V-wt_vs_pep.csv')

ggplot(group_marker, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(size=.2) + 
  ggtitle('DEG in QV1G12V-WT Vs pep')
sum(group_marker$avg_log2FC< -1)
sum(group_marker$avg_log2FC> 1)
df = data.frame('DEG_2X_up_PepVSWt' = sum(group_marker$avg_log2FC< -1),
                'DEG_2X_down_PepVSWt' = sum(group_marker$avg_log2FC> 1))
grid.table(t(df))
#### 20210618 update 2x up and 2x down gene
gene.list = group_marker[group_marker$avg_log2FC< -1,]$gene
DoHeatmap(immune.combined.sct, features =gene.list,size = 4,angle=0,draw.lines = F)+ guides(color=F)+
  theme(axis.text.y = element_text(size = 5)) +
  scale_fill_gradientn(colors = c('blue','white','red'), limits=c(-3,3))

df1 = as.matrix(GetAssayData(immune.combined.sct, slot = 'scale.data'))
gene = intersect(gene.list,rownames(df1))
df1 = df1[gene,]

pdf('MultiomeSeq/2X_up_DEG_pepVSwt.pdf')
hmap <- Heatmap((df1),name="Expression",col <- colorRamp2(c(-2,0,2),c("blue","white", "red")),
                heatmap_legend_param=list(color_bar="continuous", legend_direction="vertical", 
                                          legend_width=unit(8,"mm"), title_position="topleft", 
                                          title_gp=gpar(fontsize=10, fontface="bold")),
                top_annotation = HeatmapAnnotation(sample = rep(c('QV1G12V-Wt','QV1G12V-Pep'),each=5000),
                                                   col = list(sample = c("QV1G12V-Wt" = "yellow", "QV1G12V-Pep" = "green"))),
                
                cluster_rows=F,show_row_dend=F,row_title_side="left",row_title_gp=gpar(fontsize=8),
                show_row_names=TRUE,row_names_side="left",row_title_rot=0,row_gap = unit(1,'mm'),row_names_gp = gpar(fontsize=6),
                
                cluster_columns = F,
                column_title = NULL,column_title_side="top",column_title_rot=0,show_column_names=F,
                column_names_gp = gpar(fontsize=6),column_names_side = 'top',
                
                clustering_distance_columns=function(x) as.dist(1-cor(t(x))),clustering_method_columns="ward.D2",
                clustering_distance_rows="euclidean",clustering_method_rows="ward.D2",row_dend_width=unit(5,"mm"),
                column_dend_height=unit(5,"mm"))
draw(hmap, heatmap_legend_side="left")
dev.off()

#### 2x down gene top 100
gene.list = group_marker %>% filter(avg_log2FC > 1) %>% arrange(-avg_log2FC) %>% rownames(.)

DoHeatmap(immune.combined.sct, features =gene.list[1:100],size = 4,angle=0,draw.lines = F)+ guides(color=F)+
  theme(axis.text.y = element_text(size = 4)) +
  scale_fill_gradientn(colors = c('blue','white','red'), limits=c(-3,3))

df1 = as.matrix(GetAssayData(immune.combined.sct, slot = 'scale.data'))
gene = intersect(gene.list[1:100],rownames(df1))
df1 = df1[gene,]
pdf('MultiomeSeq/2X_down_DEG_pepVSwt.pdf')
hmap <- Heatmap((df1),name="Expression",col <- colorRamp2(c(-2,0,2),c("blue","white", "red")),
                heatmap_legend_param=list(color_bar="continuous", legend_direction="vertical", 
                                          legend_width=unit(8,"mm"), title_position="topleft", 
                                          title_gp=gpar(fontsize=10, fontface="bold")),
                
                top_annotation = HeatmapAnnotation(sample = rep(c('QV1G12V-Wt','QV1G12V-Pep'),each=5000),
                                                   col = list(sample = c("QV1G12V-Wt" = "yellow", "QV1G12V-Pep" = "green"))),
                 
                cluster_rows=F,show_row_dend=F,row_title_side="left",row_title_gp=gpar(fontsize=8),
                show_row_names=TRUE,row_names_side="left",row_title_rot=0,row_gap = unit(1,'mm'),row_names_gp = gpar(fontsize=6),
                
                cluster_columns = F,
                column_title="",column_title_side="top",column_title_rot=0,show_column_names=F,
                column_names_gp = gpar(fontsize=6),column_names_side = 'top',
                
                clustering_distance_columns=function(x) as.dist(1-cor(t(x))),clustering_method_columns="ward.D2",
                clustering_distance_rows="euclidean",clustering_method_rows="ward.D2",row_dend_width=unit(5,"mm"),
                column_dend_height=unit(5,"mm")) 

draw(hmap, heatmap_legend_side="left")
dev.off()

####### 







# 2. Identify differential expressed genes across conditions
#20210510 update --https://hbctraining.github.io/scRNA-seq/lessons/09_merged_SC_marker_identification.html
#2.0 Since we have samples representing different conditions in our dataset, our best 
#option is to find conserved markers This function internally separates out cells by sample group/condition, 
#and then performs differential gene expression testing for a single specified cluster against all other clusters 
#(or a second cluster, if specified). Gene-level p-values are computed for each condition and then combined across 
#groups using meta-analysis methods from the MetaDE R package.

#2.1 Before we start our marker identification we will explicitly set our default assay, we want to use the original counts and not the integrated data.
DefaultAssay(immune.combined.sct) <- "RNA"
Idents(immune.combined.sct) <- immune.combined.sct$seurat_clusters
#ident.1: this function only evaluates one cluster at a time; here you would specify the cluster of interest.
#grouping.var: the variable (column header) in your metadata which specifies the separation of cells into groups
cluster1_conserved_markers <- FindConservedMarkers(immune.combined.sct,
                                                   ident.1 = 1,
                                                   grouping.var = "group",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

cluster3_conserved_markers <- FindConservedMarkers(immune.combined.sct,
                                                   ident.1 = 3,
                                                   grouping.var = "group",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)



