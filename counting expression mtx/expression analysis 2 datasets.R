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
#.genes map to chrM, th1, 5,10
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

#####------ 3. Single-Cell Analysis Toolkit for Gene Expression Data in R
## creating singlecellexperiment object
BiocManager::install("scater")
library(scater)
#assay colnames() must be NULL or identical to colData rownames()
H7Gcol= data.frame(row.names = H7Gcol$sampleName, group = H7Gcol[,2])
H7Gsce = SingleCellExperiment(assays = list(counts = as.matrix(H7Gcount)), colData = as.matrix(H7GcolData))
JstcolData = data.frame(row.names = JstcolData$title, JstcolData)
Jstsce = SingleCellExperiment(assays = list(counts = as.matrix(Jstcount)), colData = as.matrix(JstcolData))
H7Gsce = H7Gsce[rowSums(counts(H7Gsce)>0)>0, ]
Jstsce = Jstsce[rowSums(counts(Jstsce)>0)>0, ]

# some usages
identical(exprs(H7Gsce), logcounts(H7Gsce))
cpm(H7Gsce) = calculateCPM(H7Gsce)
assay(H7Gsce, "is_expr") = counts(H7Gsce) > 0
assayNames(H7Gsce)  
#The calcAverage function will compute the average count for each gene after scaling each cell’s counts by its size factor
head(calcAverage(H7Gsce))

## expression plot
Jstsce = normalize(Jstsce)
H7Gsce = normalize(H7Gsce)

# x = colData group, plot gene express against column metadata
colnames(colData(H7Gsce)) = c("group", "whee")
plotExpression(H7Gsce, rownames(H7Gsce)[1:6], x = "group", exprs_values = "logcounts")
colnames(colData(Jstsce))
plotExpression(Jstsce, rownames(Jstsce)[1:6], x = "type", exprs_values = "logcounts")
# x = gene, plot gene express against a certain gene
plotExpression(H7Gsce, rownames(H7Gsce)[1:6], x = "ENSMUSG00000026096.14", exprs_values = "logcounts")
plotExpression(Jstsce, rownames(Jstsce)[1:6], x = "Krt14", exprs_values = "logcounts")
# colour_by, size_by, shape_by
plotExpression(H7Gsce, rownames(H7Gsce)[1:6], colour_by = "group", shape_by = "group", size_by = "ENSMUSG00000026096.14")
plotExpression(Jstsce, rownames(Jstsce)[1:6], colour_by = "Krt14", size_by = "Krt14", shape_by = "type")
#draw median
plotExpression(H7Gsce, rownames(H7Gsce)[88:93], x = "group", exprs_value = "counts", colour = "group", show_median = T, xlab = "group", log = T)
plotExpression(Jstsce, rownames(Jstsce)[88:95], x = "Krt5", exprs_value = "counts", colour = "Krt5", show_median = T, xlab = "Krt5", log = T)
#draw per gene 
plotExpression(H7Gsce, rownames(H7Gsce)[88:93])
plotExpression(Jstsce, rownames(Jstsce)[1:3])

## dimentionality reduction
H7Gsce = runPCA(H7Gsce)
Jstsce = runPCA(Jstsce)
plotReducedDim(H7Gsce, use_dimred = "PCA", colour_by = "group", shape_by = "group")
plotReducedDim(Jstsce, use_dimred = "PCA", colour_by = "cell.type.level.1.ch1")
plotReducedDim(H7Gsce, use_dimred = "PCA", colour_by = "group", size_by = "ENSMUSG00000045545.8")
plotReducedDim(Jstsce, use_dimred = "PCA", colour_by = "cell.type.level.1.ch1", size_by = "Krt14")


# if pre-PCA is calculated by reducedim : plotPCA(example_sce) 
H7Gsce2 <- runPCA(H7Gsce, feature_set = rowData(H7Gsce)$is_feature_control)
plotPCA(H7Gsce2)
# calculate multiple components
H7Gsce = runPCA(H7Gsce, ncomponents = 20)
plotPCA(H7Gsce, ncomponents = 4, colour_by = "group")

## t-SNE plots
set.seed(1000)
H7Gsce = runTSNE(H7Gsce, perplexity = 50)
plotTSNE(H7Gsce, colour_by = "group", size_by = "ENSMUSG00000045545.8")

set.seed(1000)
Jstsce = runTSNE(Jstsce, perplexity = 50)
plotTSNE(Jstsce, colour_by = "cell.type.level.1.ch1", size_by = "Krt14")

## Diffusion plot
BiocManager::install("destiny")
library(destiny)
H7Gsce = runDiffusionMap(H7Gsce)
plotDiffusionMap(H7Gsce, colour_by = 'group', size_by = "ENSMUSG00000045545.8")
plotDiffusionMap(Jstsce, colour_by = 'cell.type.level.1.ch1', size_by = "Krt14")

##########-----4. QC metrics
H7Gsceqc = SingleCellExperiment(assays = list(counts = as.matrix(H7Gcount)), colData = as.matrix(H7GcolData))
H7Gsceqc = calculateQCMetrics(H7Gsceqc)
colnames(colData(H7Gsceqc))
colnames(rowData(H7Gsceqc))

Jstsceqc = SingleCellExperiment(assays = list(counts = as.matrix(Jstcount)), colData = as.matrix(JstcolData) )
Jstsceqc = calculateQCMetrics(Jstsceqc)
colnames(colData(Jstsceqc))
colnames(rowData(Jstsceqc))

## QC for each control sets
H7Gsceqc = calculateQCMetrics(H7Gsceqc, feature_controls = list(ERCC = 1:20, mito = 500:1000), cell_controls = list(empty = 1:5, damaged = 31:40))
all_col_qcH7G = colnames(colData(H7Gsceqc))
all_col_qcH7G_ERCC = all_col_qcH7G[grep("ERCC", all_col_qcH7G)]

Jstsceqc = calculateQCMetrics(Jstsceqc, feature_controls = list(ERCC = 1:20, mito = 500:1000), cell_controls = list(empty = 1:5, damaged = 31:40))

## Examing the most expressed features
plotHighestExprs(H7Gsceqc, exprs_values = "counts", n = 50)
plotHighestExprs(Jstsceqc, exprs_values = "counts", n = 50)

## Frequency of expression as a function of the mean
plotExprsFreqVsMean(H7Gsceqc)
plotExprsFreqVsMean(Jstsceqc)

## Percentage of counts assigned to feature controls
plotColData(H7Gsceqc, x = "total_features_by_counts", y = "pct_counts_feature_control", colour = "H7GcolData...2.") + theme(legend.position = "top") + stat_smooth(method = "lm", size = 2)
plotColData(Jstsceqc, x = "total_features_by_counts", y = "pct_counts_feature_control", colour = "type") + theme(legend.position = "top") + stat_smooth(method = "lm", size = 2)

## Cumulative expression plot
plotScater(H7Gsceqc, block1 = "H7GcolData...2.", colour_by = "H7GcolData...2.", nfeatures = 300, exprs_values = "counts")
plotScater(Jstsceqc, block1 = "cell.type.level.1.ch1", colour_by = "cell.type.level.1.ch1", nfeatures = 300, exprs_values = "counts")

## other quality control plots

p1 = plotHighestExprs(H7Gsceqc[, H7Gsceqc$H7GcolData...2. == "nonUV"])
p2 = plotHighestExprs(H7Gsceqc[, H7Gsceqc$H7GcolData...2. == "UVIgG"])
p3 = plotHighestExprs(H7Gsceqc[, H7Gsceqc$H7GcolData...2. == "UVantiPD1"])
multiplot(p1,p2,p3, cols = 3)

#########------ 6. differential expressed gene
library(scran)
var.fit.nospike = trendVar(H7Gsceqc, parametric=TRUE, use.spikes=FALSE, loess.args=list(span=0.2))
var.out.nospike = decomposeVar(H7Gsceqc, var.fit.nospike)
plot(var.out.nospike$mean, var.out.nospike$total, pch=16,cex=0.6, xlab="Mean log-exprs", ylab='Variance of log-exprs')
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=T)

is.spike.Jst = grepl("^ERCC", rownames(Jstsceqc))
isSpike(Jstsceqc, "ERCC") <- is.spike.Jst
summary(is.spike.Jst)
var.fit.Jst = trendVar(Jstsceqc, parametric=TRUE, loess.args=list(span=0.3))
var.out.Jst = decomposeVar(Jstsceqc, var.fit.Jst)
plot(var.out.Jst$mean, var.out.Jst$total, pch=16, cex = 0.6, xlab = "Mean log-expression", ylab = "Variance of log-expression")
curve(var.fit.Jst$trend(x), col="dodgerblue", lwd=2, add=T)
cur.spike = isSpike(Jstsceqc)
points(var.out.Jst$mean[cur.spike], var.out.Jst$total[cur.spike], col="red", pch=16)
hvg.out.Jst = var.out.Jst[which(var.out.Jst$FDR <=0.05), ]
hvg.out.Jst = as.data.frame(hvg.out.Jst)
hvg.out.Jst = hvg.out.Jst[order(hvg.out.Jst$FDR, decreasing = F),]





