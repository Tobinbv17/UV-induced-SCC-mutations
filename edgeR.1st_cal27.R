library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

countdata = readRDS("./Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_PGENETICS/1st_cal27_htseqcounts.rds")
coldata = readRDS("./Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_PGENETICS/1st_cal27_QC.rds")
treatment = c("CAL27_SG1","CAL27_SG2","CAL27_CTL")
order = c('CAL27_CTL','CAL27_SG1','CAL27_SG2')
coldata$treatment = factor(treatment, levels = order)

#1. The DGEList data class
y = DGEList(counts=countdata, samples = coldata)

#2. pre-filtering of zeros
keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=F]

#3. Normalization
#1) sequencing depth This is part ofthe basic modeling procedure and flows automatically into fold-change or p-value calculations.It is always present, and doesn’t require any user intervention.
#2).If a small proportion of highly expressed genes consume a substantial proportion of the total library size for a particular sample, this will cause the remaining genes to be under-sampled for that sample. Unless this effect is adjusted for, the remaining genes may falsely appear to be down-regulated in that sample
#The calcNormFactors function normalizes the library sizes by finding a set of scaling factors for the library sizes that minimizes the log-fold changes between the samples for most genes. The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples

y = calcNormFactors(y)
y$samples

#4. Negative binomial models The classic edgeR pipeline: pairwise comparisons between two or more groups
# qCML method is only applicable on datasets with a single factor design, replicated data
#5. 1) fitting the model more complex experiments glm Generalized linear models functionality 2) estimate dispersion
#e.g. group <- factor(c(1,1,2,2,3,3))
#e.g. design <- model.matrix(~group)
y <- estimateDisp(y, design)

#6. if you have no replicates 
#Typical values for the common BCV (square-rootdispersion) for datasets arising from well-controlled experiments are 0.4 for human data,0.1 for data on genetically identical model organisms or 0.01 for technical replicates
bcv = 0.4
y = DGEList(counts = countdata, group = coldata$treatment, genes = rownames(countdata))
gtf = readRDS( "./Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei202006_phenformin/human_gene.gtf.rds")

keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=F]
m = match(y$genes$genes, gtf$gene_id)
y$genes$gene_name = gtf$gene_name[m]
d = duplicated(y$genes$gene_name)
y = y[!d,]
rownames(y$counts) = rownames(y$genes) <- y$genes$gene_name

et = exactTest(y, dispersion = bcv^2)
#et_6_inh = exactTest(y, dispersion = bcv^2, pair = c('ctl _ 6h', 'mir_21_inh _ 6h'))
#et_6_mim = exactTest(y, dispersion = bcv^2, pair = c('ctl _ 6h', 'mir_21_mimic _ 6h'))
et_atf3_sg = exactTest(y, dispersion = bcv^2, pair = c('CAL27_CTL', 'CAL27_SG1'))
et_atf3_sg1 = data.frame(et_atf3_sg)

et_atf3_sg2 = exactTest(y, dispersion = bcv^2, pair = c('CAL27_CTL', 'CAL27_SG2'))
et_atf3_sg2 = data.frame(et_atf3_sg2)
#et_6_inh_sig = data.frame(subset(et_6_inh$table, PValue < .05))
#et_6_mim_sig = data.frame(subset(et_6_mim$table, PValue < .05))
et_atf3_sg_sig = data.frame(subset(et_atf3_sg$table, PValue < .05))

data = et_atf3_sg$table
mini_data = et_atf3_sg_sig

ggplot(data, aes(logFC, -log(PValue),  color = PValue<.05)) + 
  geom_point(size = 0.5) +
  geom_text_repel(data = mini_data, aes(label = row.names(mini_data)), size = 2, color = 'black') + 
  scale_color_manual(values = c("grey","red")) +
  labs(title = "CAL27_SG1 VS CTL", color = "") +
  theme(legend.position = "none")

#6.1.2.clustering (only on highly variable genes)
library(genefilter)
library(pheatmap)
library(dplyr)
#ATF3_SG1 VS ctl in cal27
y1 = y[,c(1,3)]
topVarGenes =head(order(rowVars(log2(y1$counts+1)), decreasing=T),350)
mat = log2(y1$counts[topVarGenes,])
mat = mat - rowMeans(mat)
mat = mat[-1,] # remove inf and NaN the first row
df = as.data.frame(y1$samples[,'group'])
rownames(df) = colnames(mat)
colnames(df) = 'sample'

dfr = data.frame(gene_name = row.names(mat)) %>% left_join(gtf, by = 'gene_name') %>% filter(gene_type == ' protein_coding')

#dfr = left_join(dfr,gtf, by = "gene_id") %>% select(gene_name)
mat1 = mat[rownames(mat) %in% dfr$gene_name,]
row.names(mat1) == dfr$gene_name
pheatmap(mat1, annotation_col = df, fontsize_row = 4, fontsize_col = 4, angle_col = 0)

#ATF3_SG2 VS ctl in cal27
y1 = y[,c(2,3)]
topVarGenes =head(order(rowVars(log2(y1$counts+1)), decreasing=T),350)
mat = log2(y1$counts[topVarGenes,])
mat = mat - rowMeans(mat)
mat = mat[-c(1,2,11),] # remove inf and NaN the rows
df = as.data.frame(y1$samples[,'group'])
rownames(df) = colnames(mat)
colnames(df) = 'sample'

dfr = data.frame(gene_name = row.names(mat)) %>% left_join(gtf, by = 'gene_name') %>% filter(gene_type == ' protein_coding')

#dfr = left_join(dfr,gtf, by = "gene_id") %>% select(gene_name)
mat1 = mat[rownames(mat) %in% dfr$gene_name,]
row.names(mat1) == dfr$gene_name
pheatmap(mat1, annotation_col = df, fontsize_row = 4, fontsize_col = 4, angle_col = 0)


#6.2 treat similar datasets as 'biological replicates'
y1 = y[,4:9]
treatment = factor(c('s_12_ctl','s_12_Mkd','s_12_ctl', 's_24_ctl','s_24_Mkd','s_24_ctl'))
design = model.matrix(~treatment)
rownames(design) = colnames(y1)
y1 = estimateGLMCommonDisp(y1, design, method = 'deviance', robust = T, subset = NULL)
y1 = estimateGLMTrendedDisp(y1, design)
y1 = estimateGLMTagwiseDisp(y1, design)

#testing for DE genes likelihood ratio test
# 12h mir_inh VS ctl_12h
fit = glmFit(y1, design)
lrt.12_mkdVS12_ctl = glmLRT(fit, coef = 2)
summary(decideTests(lrt.12_mkdVS12_ctl))
plotMD(lrt.12_mkdVS12_ctl)
abline(h=c(-1,1),col='blue')

data = data.frame(topTags(lrt.12_mkdVS12_ctl, n=nrow(lrt.12_mkdVS12_ctl), adjust.method = 'BH', sort.by = 'PValue', p.value = 1))
mini_data = subset(data, FDR < .05 & logFC < -0.6 | FDR < .05 & logFC > .6)
ggplot(data, aes(logFC, -log(PValue),  color = FDR < .05 & logFC < -0.6 | FDR < .05 & logFC > .6)) + 
  geom_point(size = 0.5) +
  geom_text_repel(data = mini_data, aes(label = row.names(mini_data)), size = 2, color = 'black') + 
  scale_color_manual(values = c("grey","red")) +
  labs(title = "mir_21_12h_Mkd VS ctl_12h", color = "") +
  theme(legend.position = "none")
pheatmap(cpm(y1,log = T)[mini_data$gene_name,1:3], scale = 'row')
write.table(data, "./mir21_12h_Mkd VS ctl_12h.csv", sep= ",", row.names = F)

# 24h mir_inh VS ctl_24h
fit = glmFit(y1, design)
lrt.24_mkdVS24_ctl = glmLRT(fit, contrast = c(0,0,-1,1))
summary(decideTests(lrt.24_mkdVS24_ctl))
plotMD(lrt.24_mkdVS24_ctl)
abline(h=c(-1,1),col='blue')

data = data.frame(topTags(lrt.24_mkdVS24_ctl, n=nrow(lrt.24_mkdVS24_ctl), adjust.method = 'BH', sort.by = 'PValue', p.value = 1))
mini_data = subset(data, FDR < .05 & logFC < -0.6 | FDR < .05 & logFC > .6) %>% arrange(logFC)
mini_data = rbind(head(mini_data,10), tail(mini_data,10))
ggplot(data, aes(logFC, -log(PValue),  color = FDR < .05 & logFC < -0.6 | FDR < .05 & logFC > .6)) + 
  geom_point(size = 0.5) +
  geom_text_repel(data = mini_data, aes(label = row.names(mini_data)), size = 2, color = 'black') + 
  scale_color_manual(values = c("grey","red")) +
  labs(title = "mir_21_24h_Mkd VS ctl_24h", color = "") +
  theme(legend.position = "none")
pheatmap(cpm(y1,log = T)[mini_data$gene_name,4:6], scale = 'row')
write.table(data, "./mir21_24h_Mkd VS ctl_24h.csv", sep= ",", row.names = F)


library(fgsea)
#ATF3_SG1 VS ctl in cal27
ranks = et_atf3_sg1$logFC
names(ranks) = et_atf3_sg1$gene_name
barplot(sort(ranks, decreasing = T))

#ATF3_SG2 VS ctl in cal27
ranks = et_atf3_sg2$logFC
names(ranks) = et_atf3_sg2$gene_name
barplot(sort(ranks, decreasing = T))

#load pathway from downloaded MsigDB gmt file, hall marker pathway
pathways = gmtPathways('~/Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_RNAsesq/c2.cp.kegg.v7.1.symbols.gmt')
#pathways %>% head() %>% lapply(head) 
fgseaRes = fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500,nperm = 1000)
fgseaRes = fgseaRes[order(padj),]
library(data.table)
fwrite(fgseaRes[1:50,], file = './Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_PGENETICS/atf3_SG2_ctl_kegg.csv', row.names = T, col.names = T)

#GO annotation gmt data
fgsea_go = fgsea(pathways = gmtPathways("~/Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_RNAsesq/c5.go.all.v7.1.symbols.gmt"), 
                 ranks, minSize=15, maxSize=500,nperm=1000) %>% arrange(padj)
fwrite(fgsea_go[1:50,], file = './Documents/Bioinformatics projects/RNAseq/xunwei_202006_single end/xunwei_202008_PGENETICS/atf3_SG2_ctl_goterm.csv', row.names = T, col.names = T)

#bubblechart
a = sapply(fgseaRes$leadingEdge, function(x) {length(x)})
fgseaRes$Count = a
fgseaRes$GeneRatio = fgseaRes$Count/fgseaRes$size
fgseaRes$order = rownames(fgseaRes)

library(ggplot2)
ggplot(fgseaRes[1:20,], aes(GeneRatio, pathway, size = Count, color = -log10(pval))) +
  geom_point() + scale_size(range = c(2,5), name = 'Count') + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_color_gradient(low='green', high='red') + labs(y='', title = 'Atf3_Sg1 vs. Control')

#bubblechart
############20201006 update
fgsea_go = read.csv("./atf3_SG1_ctl_goterm.csv")
fgsea_go = data.frame(fgsea_go)

a = c(24,9,7,3,33,40,28,69,87,7,5,19,97,7,54)
fgsea_go$Count = a
fgsea_go$GeneRatio = fgsea_go$Count/fgsea_go$size
fgsea_go$order = rownames(fgsea_go)

library(ggplot2)
pdf("Atf3_Sg1_vs_ctl_20201006_top15.pdf", width = 8, height = 5)
ggplot(fgsea_go[1:15,], aes(GeneRatio, factor(order,levels = fgsea_go$X[15:1]), size = Count, color = -log10(pval))) +
  geom_point() + scale_size(range = c(3,5), name = 'Count') + 
  theme(axis.text.y = element_text(size = 8)) +scale_color_gradient(low='green', high='red') + 
  labs(y='', title = 'Atf3_Sg1 vs. Control')+
  scale_y_discrete(labels=fgsea_go$pathway[15:1]) 
dev.off()
write.table(fgsea_go,"fgsea_go_20201006_top15.csv", quote = F, sep = ",")
######################################
#10-14 update
library(dplyr)
fgsea_go = read.csv("fgsea_go_20201006_top15.csv")
fgsea_go = fgsea_go %>% arrange(-GeneRatio, -Count)

library(ggplot2)
pdf("Atf3_Sg1_vs_ctl_20201014_top15.pdf", width = 8, height = 5)
ggplot(fgsea_go, aes(GeneRatio, factor(order,levels = fgsea_go$X[15:1]), size = Count, color = -log10(pval))) +
  geom_point() + scale_size(range = c(3,5), name = 'Count') + 
  theme(axis.text.y = element_text(size = 8)) +scale_color_gradient(low='green', high='red') + 
  labs(y='', title = 'Atf3_Sg1 vs. Control')+
  scale_y_discrete(labels=fgsea_go$pathway[15:1]) 
dev.off()





#####################################
#enrichplot
library(enrichplot)
library(clusterProfiler)

#data(geneList, package="DOSE")
#de <- names(geneList)[abs(geneList) > 2]
#ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
#goplot(ego)

#7. case study 
y = DGEList(counts = countdata, genes = rownames(countdata))
#7.1 annotation
gtf = readRDS( "../xunwei202006_phenformin/human_gene.gtf.rds")
m = match(y$genes$genes, gtf$gene_id)
y$genes$gene_name = gtf$gene_name[m]
head(y$genes)
#7.2 filtering and normalization
o = order(rowSums(y$counts), decreasing = T)
y = y[o,]
d = duplicated(y$genes$gene_name)
y = y[!d,]
nrow(y)
y = y[rowSums(y$counts<50)==0,]
y$samples$lib.size = colSums(y$counts)
rownames(y$counts) = rownames(y$genes) <- y$genes$gene_name
y = calcNormFactors(y)
plotMDS(y)

treatment = rep(c('ctl', 'mir_21_inh', 'mir_21_mimic'),3)
time = rep(c('6h','12h','24h'),each =3)
data.frame(sample = colnames(y), treatment, time)

design = model.matrix(~colnames(y))
rownames(design) = colnames(y)
y = estimateDisp(y, design, robust = T)
y$common.dispersion
plotBCV(y)

fit = glmFit(y,design)
lrt = glmLRT(fit)
topTags(lrt)

o = order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1,1),col='blue')











