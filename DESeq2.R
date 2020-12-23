library(DESeq2)

countdata = readRDS("./Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_PGENETICS/1st_cal27_htseqcounts.rds")
coldata = readRDS("./Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_PGENETICS/1st_cal27_QC.rds")
treatment = c("CAL27_SG1","CAL27_SG2","CAL27_CTL")
order = c('CAL27_CTL','CAL27_SG1','CAL27_SG2')
coldata$treatment = factor(treatment, levels = order)
#coldata$treatment = relevel(factor(treatment), "con_2h")

ddsMat = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
#1). transformation of the counts for visually exploration sample relationship
#2). go back to original count data for calculating statisticaal tests

#1. pre-filtering of zeros
ddsMat = ddsMat[rowSums(counts(ddsMat)) > 1, ]
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
plotPCA(rld, intgroup ='barcode')
#plotPCA(rld, intgroup = 'treatment')

#5.Differential expression analysis
#DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there is zero effect of the treatment on the gene and that the observed difference between treatment and control was merely caused by experimental variability
#he set of genes with adjusted p value less than 0.1 should contain no more than 10% false positives.
dds = DESeq(ddsMat)
(res = results(dds))
mcols(res, use.names = T)
summary(res)

res.05 = results(dds, alpha = .05)
table(res.05$padj < .05)
resLFC1 = results(dds, lfcThreshold = 1) #(gene counts more than doubling or less than halving, because 2^1 = 2)
table(resLFC1$padj < .1)
#p values in NA, is reporting all counts for this gene were zero, and hence no test was applied. or,it contained an extreme count outlie

gtf = readRDS("~/Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei202006_phenformin/human_gene.gtf.rds")
library(dplyr)
library(tidyr)
#gtf = gtf %>% filter(V3 == "gene")
#gtf$V9 = gsub("gene_id | gene_version | gene_name | gene_source | gene_biotype", "", gtf$V9)
#gtf = gtf %>% separate (V9, c('1','2','3','4','5'), sep=";")
#gtf = gtf %>% select('1','3','5') %>% rename('gene_id' = '1', 'gene_name' = '3', 'gene_type' = "5")
#saveRDS(gtf, "./human_gene.gtf.rds")

#6.other comparison (results for a comparison of any two levels of a variable can be extracted using the contrast )
#specify three values: variable, the level of numerator, and the level for the denominator.
res_atf3_oe = results(dds, contrast = c("treatment", "CAL27_SG1", "CAL27_CTL"))
sum(res_atf3_oe$padj < .05, na.rm = T)
#res_phen10 = results(dds, contrast = c("treatment", "phen_10h", "con_10h"), alpha = .05)
DEG_atf3_oe = data.frame(subset(res_atf3_oe, padj < .05))
#res_atf3_oe = subset(res_atf3_oe, log2FoldChange > .6 | log2FoldChange < -.6)
DEG_atf3_oe = data.frame(DEG_atf3_oe[order(DEG_atf3_oe$log2FoldChange, decreasing = T),])
DEG_atf3_oe$gene_id = row.names(DEG_atf3_oe)
DEG_atf3_oe = left_join(DEG_atf3_oe, gtf, by = "gene_id")
saveRDS(res_atf3_oe, "./res_atf3_oe.rds")

#7.Plotting results
topGene = rownames(res_atf3_oe)[which.min(res_atf3_oe$padj)]
plotCounts(dds, gene = topGene, intgroup = "treatment")
plotMA(res_atf3_oe, ylim=c(-5,5))
hist(res_atf3_oe$pvalue[res_atf3_oe$baseMean >10000], breaks = 0:20/20, col='grey50', border = 'white', 
     xlab = 'p_Value[basemean >10000]', main = 'p value of genes in DE analysis' )

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
gene_res_phen10 = as.integer(row.names(dds) %in% row.names(res_phen10))
names(gene_res_phen10) = row.names(dds)
#fitting proability weighting function depending on gene length
pwf = nullp(gene_res_phen10, 'hg38', 'ensGene', bias.data = lengthfilt)
#calculate the over and under expressed GO among DE genes
go.wall = goseq(pwf,'hg38','ensGene')
enriched.GO = go.wall$category[p.adjust(go.wall$over_represented_pvalue, method = 'BH')<.05]
for (go in enriched.GO[1:10]) {
  print(GOTERM[[go]]) 
  #cat("--------\n")
}
enriched.GO.1 = as.data.frame(enriched.GO)
go.wall.filt = go.wall[match(enriched.GO.1$enriched.GO, go.wall$category),]

#KEGG pathway analysis
#library(org.Hs.eg.db)
#ens2enz = as.list(org.Hs.egENSEMBL2EG)
#enz2kegg = as.list(org.Hs.egPATH)
#grepKEGG = function(id,mapkeys){unique(unlist(mapkeys[id], use.names = F))}
#kegg=lapply(ens2enz,grepKEGG,enz2kegg)
#KEGG=goseq(pwf,gene2cat=kegg)
KEGG = goseq(pwf,'hg38','ensGene',test.cats = 'KEGG')

###### 1)fgsea
library(fgsea)
library(data.table)
library(ggplot2)
# converting ensembl to entrezgene
#library(biomaRt)
#mart = useDataset('hsapiens_gene_ensembl',useMart('ensembl'))
#genes = getBM(filters = 'ensembl_gene_id', 
#              attributes = c('ensembl_gene_id','entrezgene_id'), 
#               values=row.names(dds),mart = mart)

#res_atf3_oe = read.table("./Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_ATF3/DEG.atf3_SGa.csv", 
#                         sep = ",", header = T, row.names = 1)
res_atf3_oe_de = data.frame(res_atf3_oe) %>% mutate(gene_id = row.names(res_atf3_oe)) %>%
  left_join(gtf, by = "gene_id") 

#res_atf3_oe_de = data.frame(res_atf3_oe) %>% mutate(gene_name = row.names(res_atf3_oe)) %>%
#    left_join(gtf, by = "gene_name") 
  
gseaDat = filter(res_atf3_oe_de, !is.na(gene_name))
ranks = gseaDat$log2FoldChange
#ranks = gseaDat$logFC
names(ranks) = gseaDat$gene_name
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


#KEGG pathway gmt data
fgsea_kegg = fgsea(pathways = gmtPathways("~/Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_RNAsesq/c2.cp.kegg.v7.1.symbols.gmt"), 
                   ranks, nperm=1000) %>% arrange(padj)
ggplot(head(fgsea_kegg,15), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < .05)) +
  labs(x='Pathway', y = 'Normalized Enrichment Score', title = 'KEGG pathways NES from GSEA') +
  theme_minimal() + coord_flip()
pathways = gmtPathways("~/Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_RNAsesq/c2.cp.kegg.v7.1.symbols.gmt")
plotEnrichment(pathways[['KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION']], ranks) + labs(title = 'cytokine_signaling phen_1_10h')
plotEnrichment(pathways[['KEGG_CHEMOKINE_SIGNALING_PATHWAY']], ranks) + labs(title = 'chemokine_signaling phen_1_10h')

fwrite(fgsea_kegg[1:30,], file = './Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_ATF3/atf3_SGb_kegg_pathway.csv', row.names = T, col.names = T)


#GO annotation gmt data
fgsea_go = fgsea(pathways = gmtPathways("~/Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_RNAsesq/c5.go.all.v7.1.symbols.gmt"), 
                   ranks, minSize=15, maxSize=500,nperm=1000) %>% arrange(padj)
fgsea_go = subset(fgsea_go, padj <.05)
ggplot(head(fgsea_go,15), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < .05)) +
  labs(x='Pathway', y = 'Normalized Enrichment Score', title = 'GO term NES from GSEA') +
  theme_minimal() + coord_flip()
pathways = gmtPathways("~/Documents/xunwei_RNAsesq/c5.go.all.v7.1.symbols.gmt")
plotEnrichment(pathways[['GO_CELL_CHEMOTAXIS']], ranks) + labs(title = 'cheotaxis phen_1_10h')
plotEnrichment(pathways[['GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY']], ranks) + labs(title = 'cytokine phen_1_10h')

fwrite(fgsea_go[1:30,], file = './Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_ATF3/atf3_SGb_goterm_pathway.csv', row.names = T, col.names = T)


library(ggplot2)
library(ggrepel)

data = data.frame(results(dds, contrast = c("treatment", "ATF3_OE", "NEO")))
#######
data$gene_id = row.names(data)
data = left_join(data, gtf, by = "gene_id") %>% filter(gene_type %in% c(' protein_coding' , " lncRNA"), !is.na(padj))
# absolute fold change greater than 1.5 fold 2^0.6 = 1.5
data$threshold = ifelse(data$log2FoldChange > .6 & data$padj <.05 | data$log2FoldChange < -.6 & data$padj < .05, '>1.5fold','')
data = data[order(data$log2FoldChange, decreasing = T),][-1,]
write.table(data, './DEG.atf3.OE.csv', quote = F, sep = ',')

mini_data = subset(data, threshold == ">1.5fold")
mini_data = rbind(head(mini_data,10), tail(mini_data, 10))

pdf("./Volcano plots.atf3_oe VS neo.pdf")
#gene = c("CXCL2", "CXCL1", "CCL2", "CXCL6", "IL24", "CXCL8", "MYC",'IL33','IL27RA')
ggplot(data, aes(log2FoldChange, -log(padj),  color = threshold)) + 
  geom_point(size = 0.5) +
  geom_text_repel(data = mini_data, aes(label = gene_name), size = 2, color = 'black') + 
                    scale_color_manual(values = c("grey","red")) +
                    labs(title = "ATF3_OE vs NEO \n", color = "") +
                    theme(legend.position = "none")
dev.off()
  #geom_hline(yintercept = -log2(0.05), color = 'black', linetype = "dashed") + 
  #geom_vline(xintercept = .6, color = "black", linetype = "dashed") + 
  #geom_vline(xintercept = -.6, linetype = "dashed", color = "black") + 
  #geom_text(aes(label = ifelse(hgnc_symbol %in% gene , as.character(hgnc_symbol), '')), hjust=0, vjust=0, size = 2.5) + 

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







