library(DESeq2)

countdata = readRDS("./20201222_Scc_htseqcounts.rds")
coldata = readRDS("./20201222_scc_QC.rds")
treatment = c(rep(c("cal27_1mM","cal27_con","scc_9_1mM","scc_9_con"),each=3))
order = c("cal27_con","cal27_1mM","scc_9_con","scc_9_1mM")
coldata$treatment = factor(treatment, levels = order)
#coldata$treatment = relevel(factor(treatment), "con_2h")

ddsMat = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
#1). transformation of the counts for visually exploration sample relationship
#2). go back to original count data for calculating statisticaal tests

#1. pre-filtering of zeros
ddsMat = ddsMat[rowSums(counts(ddsMat)) > 50, ]
nrow(ddsMat)
#2.rlog transformation ( regularized-logarithm transformation or rlog  stabilize the variance across the mean count)
#(rlog on high counts similar to log2, on low counts are shrunken towards genes' averages across all samples)
rld = rlog(ddsMat, blind = F)
head(assay(rld), 5)
plot(assay(rld)[,1:2], pch=16, cex=.3)

#3.sample distance (R function dist to calculate the Euclidean distance between samples, expects the different samples to be rows of its argument, and different dimensions (here,
#genes) to be columns
sampleDists = dist(t(assay(rld))) #(Euclidean distance between samples)
sampleDists
library(pheatmap)
library(RColorBrewer)
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld$barcode)
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)
#or Poisson Distance measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration
#BiocManager::install("PoiClaClu")
library(PoiClaClu)
poisd = PoissonDistance(t(counts(ddsMat)))
samplePoisDistMatrix = as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) = paste(rld$barcode)
colnames(samplePoisDistMatrix) = NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
#4.PCA plot
library(ggplot2)
#The best way to customize the plot is to use plotPCA to return a small data.frame and then use ggplot2 to customize the graph.
df <- plotPCA(rld, intgroup =c('barcode',"treatment"), return=T) 
ggplot(df, aes(PC1,PC2,shape=treatment,color=barcode)) + geom_point(size=3,alpha=.8) +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.3,'cm')) + 
  guides(color = guide_legend(override.aes = list(size=2)))
#plotPCA(rld, intgroup = 'treatment')

#5.Differential expression analysis
#DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there is zero effect of the treatment on the gene and that the observed difference between treatment and control was merely caused by experimental variability
#the set of genes with adjusted p value less than 0.1 should contain no more than 10% false positives.
dds = DESeq(ddsMat)
(res = results(dds))
mcols(res, use.names = T)
summary(res)

res.05 = results(dds, alpha = .05)
table(res.05$padj < .05)
resLFC1 = results(dds, lfcThreshold = 1) #(gene counts more than doubling or less than halving, because 2^1 = 2)
table(resLFC1$padj < .1)
#p values in NA, is reporting all counts for this gene were zero, and hence no test was applied. or,it contained an extreme count outlie
#library(rtracklayer)
#gtf = rtracklayer::import('Homo_sapiens.GRCh38.102.gtf')
#gtf = as.data.frame(gtf)
#gtf <- gtf[,c(1:3,5,7,10,12,14)]
#gtf = gtf %>% filter(type=='gene')
saveRDS(gtf, "Homo_sapiens.GRCh38.102.gtf.rds")
gtf = readRDS("Homo_sapiens.GRCh38.102.gtf.rds")
gtf = gtf %>% select(5:8)

library(dplyr)
library(tidyr)
#gtf = gtf %>% filter(V3 == "gene")
#gtf$V9 = gsub("gene_id | gene_version | gene_name | gene_source | gene_biotype", "", gtf$V9)
#gtf = gtf %>% separate (V9, c('1','2','3','4','5'), sep=";")
#gtf = gtf %>% select('1','3','5') %>% rename('gene_id' = '1', 'gene_name' = '3', 'gene_type' = "5")
#saveRDS(gtf, "./human_gene.gtf.rds")

#6.other comparison (results for a comparison of any two levels of a variable can be extracted using the contrast )
#specify three values: variable, the level of numerator, and the level for the denominator.
res_cal27_phen = results(dds, contrast = c("treatment", "cal27_1mM", "cal27_con"))
sum(res_cal27_phen$padj < .05, na.rm = T)
#res_phen10 = results(dds, contrast = c("treatment", "phen_10h", "con_10h"), alpha = .05)
DEG_cal27_phen = data.frame(subset(res_cal27_phen, padj < .05))
#res_atf3_oe = subset(res_atf3_oe, log2FoldChange > .6 | log2FoldChange < -.6)
DEG_cal27_phen = data.frame(DEG_cal27_phen[order(DEG_cal27_phen$log2FoldChange, decreasing = T),])
DEG_cal27_phen$gene_id = paste(rownames(DEG_cal27_phen))
DEG_cal27_phen = left_join(DEG_cal27_phen, gtf, by = "gene_id")
saveRDS(DEG_cal27_phen, "./DEG_cal27_phen1mMVSctl.rds")

res_scc9_phen = results(dds, contrast = c("treatment", "scc_9_1mM", "scc_9_con"))
sum(res_scc9_phen$padj < .05, na.rm = T)
#res_phen10 = results(dds, contrast = c("treatment", "phen_10h", "con_10h"), alpha = .05)
DEG_scc9_phen = data.frame(subset(res_scc9_phen, padj < .05))
#res_atf3_oe = subset(res_atf3_oe, log2FoldChange > .6 | log2FoldChange < -.6)
DEG_scc9_phen = data.frame(DEG_scc9_phen[order(DEG_scc9_phen$log2FoldChange, decreasing = T),])
DEG_scc9_phen$gene_id = paste(rownames(DEG_scc9_phen))
DEG_scc9_phen = left_join(DEG_scc9_phen, gtf, by = "gene_id")
saveRDS(DEG_scc9_phen, "./DEG_scc9_phen1mMVSctl.rds")

summary(res_cal27_phen)

#7.Plotting results
plotMA(res_cal27_phen, ylim=c(-5,5),main = 'Cal27 phenformin 1mM VS Cal27 control')
plotMA(res_scc9_phen,ylim=c(-5,5),main='Scc9 phenformin 1mM VS Scc9 control')

#7.1 volcano plot
library(ggplot2)
library(ggrepel)

data = data.frame(res_cal27_phen)
#data = data.frame(res_scc9_phen)

data$gene_id = paste(rownames(data))
data = data %>% left_join(gtf, by = 'gene_id') %>% arrange(-log2FoldChange) %>% 
  filter(gene_biotype == 'protein_coding' & baseMean > 10) 
# absolute fold change greater than 1.5 fold 2^0.6 = 1.5
data$threshold = ifelse(data$log2FoldChange > .6 & data$padj <.05 | data$log2FoldChange < -.6 & data$padj < .05, '>1.5fold','')
write.table(data, './DEGtotal_cal27phenVSctl.csv', quote = F, sep = ',',row.names = F)
#write.table(data, './DEGtotal_sccphenVSctl.csv', quote = F, sep = ',',row.names = F)

mini_data = subset(data, threshold == ">1.5fold" & baseMean > 100)
mini_data = rbind(head(mini_data,10), tail(mini_data, 10))
ggplot(data, aes(log2FoldChange, -log10(padj),  color = threshold)) + 
  geom_point(size = 1) +
  geom_text_repel(data = mini_data, aes(label = gene_name), size = 4, color = 'black') + 
  scale_color_manual(values = c("grey","red")) +
  labs(title = "scc9 1mM phenformin vs control", color = "") + theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = .5,size = 16))
#geom_hline(yintercept = -log2(0.05), color = 'black', linetype = "dashed") + 
#geom_vline(xintercept = .6, color = "black", linetype = "dashed") + 
#geom_vline(xintercept = -.6, linetype = "dashed", color = "black") + 
#geom_text(aes(label = ifelse(hgnc_symbol %in% gene , as.character(hgnc_symbol), '')), hjust=0, vjust=0, size = 2.5) + 

a = subset(data, data$log2FoldChange>0.6 & data$padj < 0.05)
a = subset(data, data$log2FoldChange< -0.6 & data$padj < 0.05)

topGene = DEG_cal27_phen$gene_id[which.min(DEG_cal27_phen$padj)]
plotCounts(dds, gene = topGene, intgroup = "treatment")

#hist(res_cal27_phen$pvalue[res_cal27_phen$baseMean >10000], breaks = 0:20/20, col='grey50', border = 'white', 
#     xlab = 'p_Value[basemean >10000]', main = 'p value of genes in DE analysis' )
########################
#7.2 venn diagram
library(ggVennDiagram)

data_cal27 = read.csv('DEGtotal_cal27phenVSctl.csv')
data_scc9 = read.csv('DEGtotal_sccphenVSctl.csv')

a = subset(data_cal27, data_cal27$log2FoldChange < -.6 & data_cal27$padj<.05)$gene_id
b = subset(data_scc9, data_scc9$log2FoldChange < -.6 & data_scc9$padj < .05)$gene_id
x = list(cal27_down=a, scc9_down=b)
ggVennDiagram(x[1:2],label='count') + scale_fill_gradient(low = 'white',high = 'blue')

a = subset(data_cal27, data_cal27$log2FoldChange < -.6 & data_cal27$padj<.05)$gene_id
b = subset(data_scc9, data_scc9$log2FoldChange < -.6 & data_scc9$padj < .05)$gene_id
down_share = intersect(a,b)
down_share_scc9 = data_scc9[data_scc9$gene_id %in% down_share,]
down_share_cal27 = data_cal27[data_cal27$gene_id %in% down_share,]

a = subset(data_cal27, data_cal27$log2FoldChange > .6 & data_cal27$padj<.05)$gene_id
b = subset(data_scc9, data_scc9$log2FoldChange > .6 & data_scc9$padj < .05)$gene_id
up_share = intersect(a,b)
up_share_scc9 = data_scc9[data_scc9$gene_id %in% up_share,]
up_share_cal27 = data_cal27[data_cal27$gene_id %in% up_share,]

write.csv(down_share_cal27,'down_shared_gene_cal27.csv',row.names = F)
write.csv(down_share_scc9,'down_shared_gene_scc9.csv',row.names = F)
write.csv(up_share_cal27,'up_shared_gene_cal27.csv',row.names = F)
write.csv(up_share_scc9,'up_shared_gene_scc9.csv',row.names = F)

#########################GO classification
library(clusterProfiler)
library(org.Hs.eg.db)

geneList = read.csv('DEGtotal_cal27phenVSctl.csv') %>% filter(threshold == '>1.5fold' )
geneList = geneList %>% select(gene_name, log2FoldChange)
gene.df <- bitr(geneList$gene_id, fromType = "ENSEMBL",
                toType = c("SYMBOL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

ranks = geneList$log2FoldChange
#ranks = gseaDat$logFC
names(ranks) = geneList$gene_name

ego3 <- gseGO(geneList     = ranks,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              nPerm = 10000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = T, 
              keyType = "SYMBOL")

ego3 = ego3@result
#not working, results not meaningful
###########################################
#try fgsea
###### 1)fgsea
library(fgsea)
library(data.table)
library(ggplot2)
data_cal27
data_scc9

geneList = read.csv('DEGtotal_cal27phenVSctl.csv') %>% filter(threshold == '>1.5fold' )
ranks = geneList$log2FoldChange
#ranks = gseaDat$logFC
names(ranks) = geneList$gene_name
barplot(sort(ranks, decreasing = T))
#KEGG pathway gmt data
fgsea_kegg = fgsea(pathways = gmtPathways("c2.cp.kegg.v7.2.symbols.gmt"), 
                   ranks, nperm=1000) %>% arrange(padj,-size)
for (i in 1:nrow(fgsea_kegg)) {
  a = length(fgsea_kegg$leadingEdge[[i]])
  fgsea_kegg$count[i]=a
}
fgsea_kegg$GeneRatio = fgsea_kegg$count/fgsea_kegg$size
ggplot(head(fgsea_kegg,15), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = pval < .05)) +
  labs(x='Pathway', y = 'Normalized Enrichment Score', title = 'KEGG pathways NES from GSEA') +
  theme_minimal() + coord_flip()
fwrite(fgsea_kegg, file = 'cal27_DEG_kegg_nofilter.csv')
fwrite(subset(fgsea_kegg,size>9),file = 'cal27_DEG_kegg_pwyGeneMoreThanTen.csv')

fgsea_kegg = subset(fgsea_kegg,size>9)
fgsea_kegg$order = rownames(fgsea_kegg)
fgsea_kegg = fgsea_kegg %>% arrange(-GeneRatio)
pdf("Cal27_KEGG.pdf", width = 8, height = 5)
ggplot(fgsea_kegg, aes(GeneRatio, factor(order,levels = fgsea_kegg$order[23:1]), size = count, color = -log10(pval))) +
  geom_point() + scale_size(range = c(3,5), name = 'Count') + 
  theme(axis.text.y = element_text(size = 8)) +scale_color_gradient(low='green', high='red') + 
  labs(y='', title = 'Cal27 phenformin vs control KEGG pwy')+
  scale_y_discrete(labels=fgsea_kegg$pathway[23:1]) 
dev.off()


#GO annotation gmt data
fgsea_go = fgsea(pathways = gmtPathways("c5.go.v7.2.symbols.gmt"), 
                 ranks, minSize=15, maxSize=500,nperm=1000) %>% arrange(pval,-size)
fgsea_go = subset(fgsea_go, pval <.05)
for (i in 1:nrow(fgsea_go)) {
  a = length(fgsea_go$leadingEdge[[i]])
  fgsea_go$count[i]=a
}
fgsea_go$GeneRatio = fgsea_go$count/fgsea_go$size

ggplot(head(fgsea_go,15), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < .05)) +
  labs(x='Pathway', y = 'Normalized Enrichment Score', title = 'GO term NES from GSEA') +
  theme_minimal() + coord_flip()
pathways = gmtPathways("~/Documents/xunwei_RNAsesq/c5.go.all.v7.1.symbols.gmt")
plotEnrichment(pathways[['GO_CELL_CHEMOTAXIS']], ranks) + labs(title = 'cheotaxis phen_1_10h')
plotEnrichment(pathways[['GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY']], ranks) + labs(title = 'cytokine phen_1_10h')

fwrite(fgsea_go, file = 'cal27_DEG_GO_pwyGeneMoreThan15.csv')
fgsea_go$order = rownames(fgsea_go)
fgsea_go = fgsea_go %>% arrange(-GeneRatio)
pdf("Cal27_GO.pdf", width = 8, height = 5)
ggplot(fgsea_go[1:20,], aes(GeneRatio, factor(order,levels = fgsea_go$order[20:1]), size = count, color = -log10(pval))) +
  geom_point() + scale_size(range = c(3,5), name = 'Count') + 
  theme(axis.text.y = element_text(size = 8)) +scale_color_gradient(low='green', high='red') + 
  labs(y='', title = 'Cal27 phenformin vs control GO pwy')+
  scale_y_discrete(labels=fgsea_go$pathway[20:1]) 
dev.off()

##########
geneList = read.csv('DEGtotal_sccphenVSctl.csv') %>% filter(threshold == '>1.5fold' )
ranks = geneList$log2FoldChange
#ranks = gseaDat$logFC
names(ranks) = geneList$gene_name
barplot(sort(ranks, decreasing = T))
#KEGG pathway gmt data
fgsea_kegg = fgsea(pathways = gmtPathways("c2.cp.kegg.v7.2.symbols.gmt"), 
                   ranks, nperm=1000) %>% arrange(padj,-size)
for (i in 1:nrow(fgsea_kegg)) {
  a = length(fgsea_kegg$leadingEdge[[i]])
  fgsea_kegg$count[i]=a
}
fgsea_kegg$GeneRatio = fgsea_kegg$count/fgsea_kegg$size
ggplot(head(fgsea_kegg,15), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = pval < .05)) +
  labs(x='Pathway', y = 'Normalized Enrichment Score', title = 'KEGG pathways NES from GSEA') +
  theme_minimal() + coord_flip()
fgsea_kegg = fgsea_kegg %>% filter(pval<.05, size >9)
fwrite(fgsea_kegg, file = 'scc9_DEG_kegg_nofilter.csv')

fgsea_kegg$order = rownames(fgsea_kegg)
fgsea_kegg = fgsea_kegg %>% arrange(-GeneRatio)
pdf("Scc9_KEGG.pdf", width = 8, height = 5)
ggplot(fgsea_kegg, aes(GeneRatio, factor(order,levels = fgsea_kegg$order[9:1]), size = count, color = -log10(pval))) +
  geom_point() + scale_size(range = c(3,5), name = 'Count') + 
  theme(axis.text.y = element_text(size = 8)) +scale_color_gradient(low='green', high='red') + 
  labs(y='', title = 'Scc9 phenformin vs control KEGG pwy')+
  scale_y_discrete(labels=fgsea_kegg$pathway[9:1]) 
dev.off()


#GO annotation gmt data
fgsea_go = fgsea(pathways = gmtPathways("c5.go.v7.2.symbols.gmt"), 
                 ranks, minSize=15, maxSize=500,nperm=1000) %>% arrange(pval,-size)
fgsea_go = subset(fgsea_go, pval <.05)
for (i in 1:nrow(fgsea_go)) {
  a = length(fgsea_go$leadingEdge[[i]])
  fgsea_go$count[i]=a
}
fgsea_go$GeneRatio = fgsea_go$count/fgsea_go$size

ggplot(head(fgsea_go,15), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < .05)) +
  labs(x='Pathway', y = 'Normalized Enrichment Score', title = 'GO term NES from GSEA') +
  theme_minimal() + coord_flip()
pathways = gmtPathways("~/Documents/xunwei_RNAsesq/c5.go.all.v7.1.symbols.gmt")
plotEnrichment(pathways[['GO_CELL_CHEMOTAXIS']], ranks) + labs(title = 'cheotaxis phen_1_10h')
plotEnrichment(pathways[['GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY']], ranks) + labs(title = 'cytokine phen_1_10h')

fwrite(fgsea_go, file = 'scc9_DEG_GO_pwyGeneMoreThan15.csv')
fgsea_go$order = rownames(fgsea_go)
fgsea_go = fgsea_go %>% arrange(-GeneRatio)
pdf("scc9_GO.pdf", width = 8, height = 5)
ggplot(fgsea_go[1:20,], aes(GeneRatio, factor(order,levels = fgsea_go$order[20:1]), size = count, color = -log10(pval))) +
  geom_point() + scale_size(range = c(3,5), name = 'Count') + 
  theme(axis.text.y = element_text(size = 8)) +scale_color_gradient(low='green', high='red') + 
  labs(y='', title = 'Scc9 phenformin vs control GO pwy')+
  scale_y_discrete(labels=fgsea_go$pathway[20:1]) 
dev.off()

#8.clustering (only on highly variable genes)
library(genefilter)
#ATF3_SG1 in TSCC
topVarGenes = head(order(rowVars(assay(rld[,-c(3,4)])), decreasing=T), 100)
mat = assay(rld[,-c(3,4)])[topVarGenes,]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld[,-c(3,4)])[,'barcode'])
rownames(df) = colnames(mat)
colnames(df) = 'sample'
dfr = data.frame(gene_id = row.names(mat))
dfr = left_join(dfr,gtf, by = "gene_id") %>% select(gene_name)
row.names(mat) = dfr$gene_name
pheatmap(mat, annotation_col = df, fontsize_row = 4, fontsize_col = 4, angle_col = 0)

#9 independent filtering
qs = c(0, quantile(res2_phen10$baseMean[res2_phen10$baseMean>0], 0:6/6))
bins = cut(res2_phen10$baseMean, qs)
levels(bins) = paste0('~', round(signif(.5*qs[-1] +.5*qs[-length(qs)],2)))
ratios = tapply(res2_phen10$pvalue, bins, function(p) mean(p < .05, na.rm=T))
barplot(ratios, xlab = 'mean normalized count', ylab='ratio of samll p values')

#10 exporting data
DEG_atf3_oe = data.frame(res_atf3_oe) %>% mutate(gene_id = row.names(res_atf3_oe)) %>% 
  left_join(gtf, by = "gene_id") %>% filter(padj < .05)

#BiocManager::install('ReportingTools')
library(ReportingTools)
htmlRep = HTMLReport(shortName = 'report', title = 'myreport', reportDirectory = './report')
publish(DEG_atf3_oe, htmlRep)
url = finish(htmlRep)
browseURL(url)

#11 Time course experiments .....
coldata2 = coldata2 %>% mutate(time = c("10h",'2h','10h','2h','10h','2h')) %>% mutate(tx = c(rep("AICAR",2),rep('con',2),rep('phen',2)))
dds2Mat2 = DESeqDataSetFromMatrix(countData = count2, colData = coldata2, design = ~ tx + time + tx:time)
dds2tc = DESeq(dds2Mat2)
resultsNames(dds2tc)

#12. Functional analysis of DEGs
library(goseq)
library(GenomicFeatures)
library(GO.db)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# get gene length if not in goseq db ########################
txdb = makeTxDbFromGFF("./Homo_sapiens.GRCh38.99.gtf", format = "gtf")
saveRDS(txdb, "./txdb.hg38.gtf.rds")
exon.per.gene = exonsBy(txdb, by = 'gene') #collect exons per gene id
exonic.gene.size = sum(width(reduce(exon.per.gene)))
lengthfilt = exonic.gene.size[match(row.names(gene_res_phen10),names(exonic.gene.size))]
#############################################################

###### 1)fgsea
library(fgsea)
library(data.table)
library(ggplot2)
data_cal27
data_scc9

geneList = read.csv('DEGtotal_cal27phenVSctl.csv') %>% filter(threshold == '>1.5fold' )

ranks = geneList$log2FoldChange
#ranks = gseaDat$logFC
names(ranks) = geneList$gene_name
barplot(sort(ranks, decreasing = T))

#load pathway from downloaded MsigDB gmt file, hall marker pathway
pathways = gmtPathways("~/Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_RNAsesq/h.all.v7.1.symbols.gmt")
#pathways %>% head() %>% lapply(head) 
fgseaRes = fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500,nperm = 1000)
fgseaRes = fgseaRes[order(padj),]
plotEnrichment(pathways[['HALLMARK_G2M_CHECKPOINT']], ranks) + labs(title = 'HALLMARK_G2M_CHECKPOINT ATF3 OE')
plotEnrichment(pathways[['HALLMARK_E2F_TARGETS']], ranks) + labs(title = 'HALLMARK_E2F_TARGETS ATF3 OE')
plotEnrichment(pathways[['HALLMARK_MYC_TARGETS_V1']], ranks) + labs(title = 'HALLMARK_MYC_TARGETS_V1 ATF3 OE')
plotEnrichment(pathways[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']], ranks) + labs(title = 'HALLMARK_OXIDATIVE_PHOSPHORYLATION ATF3 OE')

ggplot(head(fgseaRes,15), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < .05)) +
  labs(x='Pathway', y = 'Normalized Enrichment Score', title = 'Hallmark pathways NES from GSEA ') +
  theme_minimal() + coord_flip()
library(data.table)
fwrite(fgseaRes[1:30,], file = './Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_ATF3/atf3_SGb_fgsea.csv', row.names = T, col.names = T)




###################### compare replicate_1 and replicate_2 
#piechart on DEgene #
data = data.frame(group = c("up", "down"), value = c(sum(res_phen10$log2FoldChange >.6), sum(res_phen10$log2FoldChange< -.6)))
data = data.frame(group = c("up", "down"), value = c(sum(res2_phen10$log2FoldChange>.6), sum(res2_phen10$log2FoldChange< -.6)))
data = data.frame(group = c("up", "down"), value = c(sum(res3_phen10$log2FoldChange>.6), sum(res3_phen10$log2FoldChange< -.6)))

data = data %>% mutate(midpoint = cumsum(value)-value/2)

ggplot(data, aes(x='', y=value, fill=group)) + 
  geom_bar(stat = 'identity', color = "white",size=2, fill = c("yellow",'purple')) + coord_polar('y',start = 0) +
  geom_text(aes( y = midpoint, label = paste(group,'gene','\n',value)), size=4, color='black') +
  theme_void() + labs(title='phen_2_10h VS con_2', caption = '(padj < .05 |FC| > 1.5)') + 
  theme(plot.title = element_text(hjust = .5, vjust = -5,size = 12), 
        plot.caption = element_text(hjust = .5, vjust = 10))

#vignet chart on shared gene in categories
library(VennDiagram)
#) up genes
set1 = row.names(subset(res_phen10, log2FoldChange < -.6))
set2 = row.names(subset(res2_phen10, log2FoldChange< -.6))
set3 = row.names(subset(res3_phen10, log2FoldChange < -.6))

vd = venn.diagram(x=list(set1,set2,set3),
             category.names = c("phen10h_ttl","phen10_rep1","phen10_rep2"),
             filename = 'up gene venn diagram.png', output = T,
             
             imagetype = 'png',height = 480, width = 480,
             resolution = 300, compression = 'lzw',
             
             lwd=2, col=c("#440154ff", '#21908dff', 'orange'), 
             fill=c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('orange',0.3)),
             cex=.6, cat.cex=.6,
             fontfamily = "sans",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 175),
             cat.dist = c(0.085, 0.085, 0.055),
             cat.fontfamily = "sans",
             cat.col = c("#440154ff", '#21908dff', 'orange'),
             rotation = 1,
             main='down genes', main.cex =.8)

##########################################
########3# mouse to human gene conversion
library(biomaRt)

Mouse2Human <- function(MouseGenes){
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesMousetoHuman = getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), 
                             filters = "mgi_symbol", 
                             values = MouseGenes , 
                             mart = mouse, 
                             attributesL = c("ensembl_gene_id", "hgnc_symbol"), 
                             martL = human, 
                             uniqueRows = TRUE)
  
  colnames(genesMousetoHuman) <- c("Mouse.Gene_ID", "MGI", "Human.Gene_ID", "HGNC")
  
  return(genesMousetoHuman) 
  
}

## Get mouse genes
mmusculus_genes <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),  
                         mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                         useCache = FALSE)

## create the conversion table
Mouse2HumanTable <- Mouse2Human(MouseGenes = mmusculus_genes$mgi_symbol)







