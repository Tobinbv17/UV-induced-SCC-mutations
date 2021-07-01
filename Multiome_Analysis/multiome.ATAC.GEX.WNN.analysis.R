################################################################################
##https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
##WNN analysis of 10x Multiome, RNA + ATAC
################################################################################
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
require(Biostrings)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
#1.
# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("~/cellranger-v5_v6/Multiome/nextseq0062-arc-QV1G12V-WT/outs/filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
wt <- CreateSeuratObject(counts = rna_counts,project = 'wt_ATAC-GEX')
wt[["percent.mt"]] <- PercentageFeatureSet(wt, pattern = "^MT-")
# Now add in the ATAC-seq data, we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "~/cellranger-v5_v6/Multiome/nextseq0062-arc-QV1G12V-WT/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations)
wt[["ATAC"]] <- chrom_assay

#2.We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.
p1=VlnPlot(wt, features = c("nCount_ATAC", "nCount_RNA","percent.mt",'nFeature_ATAC','nFeature_RNA'), 
        ncol = 5,log = TRUE,pt.size = .000001) + NoLegend()&
  theme(axis.text = element_text(size=7),axis.text.x = element_text(angle = 0,hjust =.5), 
        axis.title = element_text(size=8),axis.title.x = element_blank(),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
wt_fil <- subset(x = wt,subset = nFeature_RNA>=3000)
p2 = VlnPlot(wt_fil, features = c("nCount_ATAC", "nCount_RNA","percent.mt",'nFeature_ATAC','nFeature_RNA'), 
             ncol = 5,log = TRUE,pt.size = .000001) + NoLegend()&
  theme(axis.text = element_text(size=7),axis.text.x = element_text(angle = 0,hjust =.5), 
        axis.title = element_text(size=8),axis.title.x = element_blank(),
        plot.title = element_text(size=10),legend.text = element_text(size=7))
p1/p2

#3.We next perform pre-processing and dimensional reduction on both assays independently, using standard approaches for RNA and ATAC-seq data.
# RNA analysis
DefaultAssay(wt_fil) <- "RNA"
wt_fil <- SCTransform(wt_fil, verbose = FALSE) %>% RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# ATAC analysis ,We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(wt_fil) <- "ATAC"
wt_fil <- RunTFIDF(wt_fil) %>% FindTopFeatures(.,min.cutoff='q0') %>% RunSVD(.) %>%
  RunUMAP(.,reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#4.We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
wt_fil <- FindMultiModalNeighbors(wt_fil, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
wt_fil <- RunUMAP(wt_fil, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
wt_fil <- FindClusters(wt_fil, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
DimPlot(wt_fil)
# perform sub-clustering on cluster 6 to find additional structure
wt_fil <- FindSubCluster(wt_fil, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(wt_fil) <- "sub.cluster"
Idents(wt_fil) <- 'seurat_clusters'
#We can visualize clustering based on gene expression, ATAC-seq, or WNN analysis. The differences are more subtle than in the previous analysis (you can explore the weights, which are more evenly split than in our CITE-seq example), but we find that WNN provides the clearest separation of cell states.
p1 <- DimPlot(wt_fil, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(wt_fil, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(wt_fil, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

#5.Visualize cell marker by ATACseq assistance. 
#For example, the ATAC-seq data assists in the separation of CD4 and CD8 T cell states. This is due to the presence of multiple loci that exhibit differential accessibility between different T cell subtypes. For example, we can visualize ‘pseudobulk’ tracks of the CD8A locus alongside violin plots of gene expression levels, using tools in the Signac visualization vignette.
## to make the visualization easier, subset T cell clusters
celltype.names <- Idents(wt_fil)
tcell.names <- celltype.names#grep("CD4|CD8|Treg", celltype.names,value = TRUE)
tcells <- subset(wt_fil, idents = tcell.names)
CoveragePlot(tcells, region = 'CD8A', features = 'CD8A', assay = 'ATAC', expression.assay = 'SCT', peaks = T)

#6.Examine the accessible regions of each cell for enriched motifs
#Next, we will examine the accessible regions of each cell to determine enriched motifs. As described in the Signac motifs vignette, there are a few ways to do this, but we will use the chromVAR package from the Greenleaf lab. This calculates a per-cell accessibility score for known motifs, and adds these scores as a third assay (chromvar) in the Seurat object.
#chromVAR for the analysis of motif accessibility in scATAC-seq; presto for fast differential expression analyses.;TFBSTools for TFBS analysis
#JASPAR2020 for JASPAR motif models; motifmatchr for motif matching; BSgenome.Hsapiens.UCSC.hg38 for chromVAR
remotes::install_github("immunogenomics/presto")
BiocManager::install(c("chromVAR", "TFBSTools", "JASPAR2020", "motifmatchr", "BSgenome.Hsapiens.UCSC.hg38"))  
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(pbmc) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(pbmc), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
pbmc <- SetAssayData(pbmc, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)





