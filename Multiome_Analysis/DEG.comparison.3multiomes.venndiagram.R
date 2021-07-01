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
library(VennDiagram)
library(ggVennDiagram)
#1.

DEG_nuc_wp = read.csv('MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_nuc_wt_pep.logNormIntegrate.csv',row.names = 1)
DEG_cyt_wp = read.csv('MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_cyt_wt_pep.logNormIntegrate.csv',row.names = 1)
DEG_fb_wp = read.csv('MultiomeSeq/cytoplasmic_RNA_feature_barcode/fb_DEG_wt_pep.logNormIntegrate.csv',row.names = 1)

x = list('DEG_nuclear_2Xup' = rownames(DEG_nuc_wp[DEG_nuc_wp$avg_log2FC< -1 & DEG_nuc_wp$p_val_adj<.05,]),
         'DEG_cyto\n_2Xup' = rownames(DEG_cyt_wp[DEG_cyt_wp$avg_log2FC< -1 & DEG_cyt_wp$p_val_adj<.05,]),
         'DEG_fb\n_2Xup' = rownames(DEG_fb_wp[DEG_fb_wp$avg_log2FC< -1 & DEG_fb_wp$p_val_adj<.05,]))
ggVennDiagram(x) + scale_fill_gradient(low="white",high = "red")

y = list('DEG_nuclear_2Xdown' = rownames(DEG_nuc_wp[DEG_nuc_wp$avg_log2FC>1 & DEG_nuc_wp$p_val_adj<.05,]),
         'DEG_cyto\n_2Xdown' = rownames(DEG_cyt_wp[DEG_cyt_wp$avg_log2FC>1 & DEG_cyt_wp$p_val_adj<.05,]),
         'DEG_fb\n_2Xdown' = rownames(DEG_fb_wp[DEG_fb_wp$avg_log2FC>1 & DEG_fb_wp$p_val_adj<.05,]))
ggVennDiagram(y) + scale_fill_gradient(low="white",high = "blue")

DEG_nuc_wp$annotation = ifelse(DEG_nuc_wp$avg_log2FC>1 & DEG_nuc_wp$p_val_adj<.05 | DEG_nuc_wp$avg_log2FC < -1 & DEG_nuc_wp$p_val_adj<.05,
                               '2 fold change','')
DEG_cyt_wp$annotation = ifelse(DEG_cyt_wp$avg_log2FC>1 & DEG_cyt_wp$p_val_adj<.05 | DEG_cyt_wp$avg_log2FC < -1 & DEG_cyt_wp$p_val_adj<.05,
                               '2 fold change','')
DEG_fb_wp$annotation = ifelse(DEG_fb_wp$avg_log2FC>1 & DEG_fb_wp$p_val_adj<.05 | DEG_fb_wp$avg_log2FC < -1 & DEG_fb_wp$p_val_adj<.05,
                               '2 fold change','')
write.csv(DEG_nuc_wp,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_nuc_wt_pep.logNormIntegrate.csv')
write.csv(DEG_cyt_wp,'MultiomeSeq/nuclearRNA_cytoplasmicRNA/DEG_cyt_wt_pep.logNormIntegrate.csv')
write.csv(DEG_fb_wp,'MultiomeSeq/cytoplasmic_RNA_feature_barcode/fb_DEG_wt_pep.logNormIntegrate.csv')

#################20210628 update, extract the overlapping genes in both 2xup and 2xdown in all the 6 samples
x = list('DEG_nuclear_2Xup' = rownames(DEG_nuc_wp[DEG_nuc_wp$avg_log2FC< -1 & DEG_nuc_wp$p_val_adj<.05,]),
         'DEG_cyto\n_2Xup' = rownames(DEG_cyt_wp[DEG_cyt_wp$avg_log2FC< -1 & DEG_cyt_wp$p_val_adj<.05,]),
         'DEG_fb\n_2Xup' = rownames(DEG_fb_wp[DEG_fb_wp$avg_log2FC< -1 & DEG_fb_wp$p_val_adj<.05,]))
up = Reduce(intersect, x)
df_up = data.frame(gene = up,annotation='2xUP_Pep_vs_Wt',
                   nuc_wp_avg.log2.FC = DEG_nuc_wp[up,2],nuc_wp_p.val.adj=DEG_nuc_wp[up,5],
                   cyt_wp_avg.log2.FC = DEG_cyt_wp[up,2],cyt_wp_p.val.adj=DEG_cyt_wp[up,5],
                   fb_wp_avg.log2.FC = DEG_fb_wp[up,2],fb_wp_p.val.adj=DEG_fb_wp[up,5])

y = list('DEG_nuclear_2Xdown' = rownames(DEG_nuc_wp[DEG_nuc_wp$avg_log2FC>1 & DEG_nuc_wp$p_val_adj<.05,]),
         'DEG_cyto\n_2Xdown' = rownames(DEG_cyt_wp[DEG_cyt_wp$avg_log2FC>1 & DEG_cyt_wp$p_val_adj<.05,]),
         'DEG_fb\n_2Xdown' = rownames(DEG_fb_wp[DEG_fb_wp$avg_log2FC>1 & DEG_fb_wp$p_val_adj<.05,]))
down = Reduce(intersect, y)
df_down = data.frame(gene = down,annotation='2xDown_Pep_vs_Wt',
                   nuc_wp_avg.log2.FC = DEG_nuc_wp[down,2],nuc_wp_p.val.adj=DEG_nuc_wp[down,5],
                   cyt_wp_avg.log2.FC = DEG_cyt_wp[down,2],cyt_wp_p.val.adj=DEG_cyt_wp[down,5],
                   fb_wp_avg.log2.FC = DEG_fb_wp[down,2],fb_wp_p.val.adj=DEG_fb_wp[down,5])
df = rbind(df_up,df_down)
write.csv(df,'MultiomeSeq/DEG_comparison_6samples/DEG_2xUp.2xDown.3comparison.csv')

################################################### making feature_ref.csv
ablist =read.csv('~/cellranger-v5_v6/feature_barcode/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.csv')

fb = data.frame(id = ablist$Gene.name, name = ablist$Gene.name, read = 'R2',
                pattern = '5PNNNNNNNNNN(BC)', sequence = ablist$Barcode, feature_type='Antibody Capture')
write.csv(ablist, '~/cellranger-v5_v6/feature_barcode/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.csv',row.names = F)
write.csv(fb, '~/cellranger-v5_v6/feature_barcode/feature_ref_totalseqC.csv',row.names = F)




