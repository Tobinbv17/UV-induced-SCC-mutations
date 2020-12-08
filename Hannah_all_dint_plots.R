Skip to content
Search or jump to…

Pull requests
Issues
Marketplace
Explore
 
@Tobinbv17 
Learn Git and GitHub without any code!
Using the Hello World guide, you’ll start a branch, write comments, and open a pull request.


etaijacob
/
MOUSE_MELANOMA
Private
2
10
Code
Issues
Pull requests
Actions
Projects
Security
Insights
MOUSE_MELANOMA/latest_scripts/dint_analysis_plots_new.R
@hannah-almubarak
hannah-almubarak all dint plots
…
Latest commit d1e52a0 on Dec 13, 2019
 History
 1 contributor
313 lines (225 sloc)  12 KB
  
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)

#read in the dint_table including sample, strand, and CC GG counts per transcript
all_dint_in <- read.table("/Users/hannah/Desktop/DFCI/UV_proj/vcf/high_depth/all_dint_high_depth.dat", stringsAsFactors = FALSE, header = T)

#read in gene TPM values from rsem output
gene_tpm <- read.table("/Users/hannah/Desktop/rsem/rsem_gene_TPM.dat", stringsAsFactors = F, header = T, check.names = F)
gene_tpm <- gene_tpm %>%  separate(gene_id, c("gene_id", "gene_name"), sep = "_")

#read in cell type info from the classification analysis (including sample name and cell type)
cell_type <- read.table("/Users/hannah/Desktop/DFCI/UV_proj/UVh7g_cellAnnotation.dat", stringsAsFactors = F, header = T, check.names = F)
cell_type <- cell_type %>% select(1, 4)
colnames(cell_type) <- c("sample", "cell_type")

#groups
groups <- read.table("/Users/hannah/Desktop/DFCI/UV_proj/vcf/groups.txt", stringsAsFactors = FALSE, header = T)


#get TPM values per gene per sample
all_dint_tpm <- NULL
for (i in unique(all_dint_in$sample)){
  all_dint2 <- all_dint_in
  all_dint2 <- all_dint2[all_dint2$sample == i,]
  
  gene_tpm_sub <- gene_tpm
  gene_tpm_sub <- as.data.frame(gene_tpm_sub[,c("gene_id", i)])
  colnames(gene_tpm_sub) <- c("gene_id", "TPM")
  
  all_dint2 <- all_dint2 %>% left_join(gene_tpm_sub, by = "gene_id")
  all_dint_tpm <- all_dint_tpm %>% rbind(all_dint2)
}


all_dint_tpm2 <- all_dint_tpm %>% left_join(cell_type, by="sample")

#exclude immune cells
all_dint_tpm2 <- all_dint_tpm2[! is.na(all_dint_tpm2$cell_type),]
all_dint_tpm2 <- all_dint_tpm2[! (all_dint_tpm2$cell_type == "TC" | all_dint_tpm2$cell_type == "LH"),]

#exclude any shared mutations
all_dint_tpm2 <- all_dint_tpm2 %>% distinct(chrom_pos, .keep_all = T)

#apply general filters
all_dint_tpm2 <- all_dint_tpm2 %>% mutate(DP = ref_AD + alt_AD)
all_dint_tpm2 <- all_dint_tpm2[(all_dint_tpm2$TPM >= 10 & all_dint_tpm2$TPM <= 1000),]
all_dint_tpm2 <- all_dint_tpm2[(all_dint_tpm2$DP >= 10 & all_dint_tpm2$DP <= 100),]


###PLOTS###

##PLOT-1##
#plot the difference between number of mutations on transcribed and non-transcribed strand per gene

#save all plots for different thresholds in a single pdf file
pdf("/Users/hannah/Desktop/DFCI/UV_proj/difference_plots.pdf")
#apply specific filters based on reads supporting variant allele
for(threshold in seq(1,50,2)){
#filter based on alt_AD
all_dint_tpm3 <- all_dint_tpm2[all_dint_tpm2$alt_AD >= threshold, ]

#take number of available substrates per transcript into account
all_dint_tpm3 <- all_dint_tpm3 %>% mutate(n_substrate = case_when(ref== 'C' ~ CC_count, ref == 'G' ~ GG_count))
all_dint_tpm4 <- all_dint_tpm3 %>% group_by(ref, sign, strand, group, gene_id, gene_name) %>% summarise(mt_count = n(), mean_TPM = mean(TPM)) 

all_dint_ccgg <- all_dint_tpm3 %>% select(gene_id, ref, n_substrate) %>% mutate(tmp = paste(gene_id, ref, sep = '_'))
all_dint_ccgg <- all_dint_ccgg[! duplicated(all_dint_ccgg$tmp),]
all_dint_ccgg <- all_dint_ccgg[, -4]

all_dint_tpm4 <- all_dint_tpm4 %>% left_join(all_dint_ccgg, by = c("gene_id", "ref"))
all_dint_tpm4 <- all_dint_tpm4 %>% mutate(fraction = mt_count/n_substrate)


#to get the difference between the transcribed and non-transcribed mt fraction
temp <- all_dint_tpm4 %>%
  group_by(strand, group, gene_id) %>%
  summarize(fraction = fraction) %>%
  ungroup() %>%
  complete(strand, group, gene_id,
           fill = list(fraction = 0))

save_mean <- all_dint_tpm4[, c(4,5,8)]

temp1 <- temp %>% group_by(group, gene_id) %>% summarise(diff = diff(fraction)) 
temp1 <- temp1 %>% left_join(save_mean, by = c("group", "gene_id"))

#ratio of genes with more non-transcribed mt/genes with more transcribed mt
#antiPD1:
antipd1_r <- length(which(temp1[temp1$group=="antiPD1",]$diff <0))/length(which(temp1[temp1$group=="antiPD1",]$diff >0))
#IgG:
Igg_r <- length(which(temp1[temp1$group=="IgG",]$diff <0))/length(which(temp1[temp1$group=="IgG",]$diff >0))
#nonUV:
nonuv_r <- length(which(temp1[temp1$group=="nonUV",]$diff <0))/length(which(temp1[temp1$group=="nonUV",]$diff >0))

#create dataframe to store ratios
dat_text <- data.frame(
  label = c(paste("ratio= ", round(antipd1_r,3)), paste("ratio= ", round(Igg_r,3)), paste("ratio= ", round(nonuv_r,3))),
  group   = c("antiPD1", "IgG", "nonUV"))


p <- ggplot(temp1, aes(x=log(mean_TPM), y=diff, color=group)) +
  facet_grid(.~ factor(group)) +
  labs( title =paste("Dinucleotide CC>TT Mutations (alt_AD>=", threshold, ")"), y = "mutation fraction diff (transcribed - non-transcribed)") +
  geom_point(size=0.7) + 
  geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.1,
    vjust   = -1,
    col = "black"
  ) 

plot(p)

}
dev.off() 


#additional plots


#BARPLOTS#

##PLOT-2##
#Total Number of CC>TT Mutations

pdf("/Users/hannah/Desktop/DFCI/UV_proj/overall_barplots.pdf")
#apply specific filters based on reads supporting variant allele
for(threshold in seq(1,50,2)){
  
all_dint_tpm3 <- all_dint_tpm2[all_dint_tpm2$alt_AD >= threshold, ]
  
plot <- as.data.frame(table(all_dint_tpm3$group, all_dint_tpm3$strand))
colnames(plot) <- c("group", "strand", "count")

p <- ggplot(plot , aes(group, count, fill=strand, group=strand)) +
  geom_bar(stat='identity', position='dodge') + labs( title =paste("Total Number of CC>TT Mutations (alt_AD>=", threshold, ")")) +
  theme(plot.title = element_text(hjust = 0.1)) + ylim(c(0,800))

plot(p)

}
dev.off() 


##PLOT-3##
#Total Number of CC>TT Mutations (sign specific)
t1 <- as.data.frame(table(all_dint$group, all_dint$strand, all_dint$sign))
colnames(t1) <- c("group", "strand", "sign","count")

ggplot(t1 , aes(group, count, fill=strand, group=strand)) +
  geom_bar(stat='identity', position='dodge') + labs( title ="Total Number of CC>TT Mutations") +
  facet_grid(.~factor(sign)) +
  theme(plot.title = element_text(hjust = 0.1))


##PLOT-4##
# number of cells with a given number of mutations (by group)
all_dint2 <- all_dint_tpm2 %>% group_by(sample, group) %>% summarise(CC_count=n())

all_dint3 <- all_dint2 %>% group_by(CC_count, group) %>% summarise(sample_count = n())

ggplot(all_dint3, aes(x=factor(CC_count), y = sample_count)) + 
  facet_grid(factor(group)~.) +
  geom_bar(stat="identity", width =0.7)+ theme_bw() + labs(x="Number of CC>TT Mutations", y="Number of Cells") + 
  theme(axis.text.x = element_text(angle=60, size = 3.5, hjust=1)) + theme(legend.position="top")


##PLOT-5##
#get mutations per cell with strand info (per treatment group)
all_dint4 <-  all_dint_tpm2 %>% group_by(sample, group, strand) %>% summarise(CC_count=n())

test <- all_dint4  %>%
  group_by(sample,strand) %>%
  summarize(CC_count = sum(CC_count)) %>%
  ungroup() %>%
  complete(sample, strand,
           fill = list(CC_count = 0))

#re-join with group names
test <- test %>% left_join(groups, by = "sample")

#specify group name 
ggplot(subset(test,  group == "antiPD1"), aes(x=reorder(sample,-CC_count), y=CC_count, fill=strand)) + geom_bar(stat="identity", width =0.7)+ 
  theme_bw() + labs(x="antiPD1 Cells", y="Number of CC>TT Mutations") + theme(axis.text.x = element_text(angle=60, size = 3.5, hjust=1)) +
  theme(legend.position="top") + labs(fill="") +ylim(0,100)

ggplot(subset(test,  group == "IgG"), aes(x=reorder(sample,-CC_count), y=CC_count, fill=strand)) + geom_bar(stat="identity", width =0.7)+ 
  theme_bw() + labs(x="IgG Cells", y="Number of CC>TT Mutations") + theme(axis.text.x = element_text(angle=60, size = 3.5, hjust=1)) +
  theme(legend.position="top") + labs(fill="") +ylim(0,100)

ggplot(subset(test,  group == "nonUV"), aes(x=reorder(sample,-CC_count), y=CC_count, fill=strand)) + geom_bar(stat="identity", width =0.7)+ 
  theme_bw() + labs(x="nonUV Cells", y="Number of CC>TT Mutations") + theme(axis.text.x = element_text(angle=60, size = 3.5, hjust=1)) +
  theme(legend.position="top") + labs(fill="") +ylim(0,65)


#PLOT-6#
#density of transcribed mutations
mutations_per_cell <-  all_dint_tpm2 %>% group_by(sample, group, strand) %>% summarise(count=n())
transcribed <- mutations_per_cell[ mutations_per_cell$strand == "transcribed",]
tbl_trans <- transcribed[,c(2,4)]

ggplot(tbl_trans, aes(x=count, fill=group)) + geom_density(alpha=0.4) +
  labs( title ="Transcribed Strand", x = "Number of CC>TT mutations per Cell") +
  theme(plot.title = element_text(hjust = 0.5)) +ylim(c(0,0.3))

#density of non_transcribed mutations
non_transcribed <- mutations_per_cell[ mutations_per_cell$strand == "non-transcribed",]
tbl_non_trans <- non_transcribed[,c(2,4)]

ggplot(tbl_non_trans, aes(x=count, fill=group)) + geom_density(alpha=0.4) +
  labs( title ="non-transcribed Strand", x = "Number of CC>TT mutations per Cell") +
  theme(plot.title = element_text(hjust = 0.5)) +ylim(c(0,0.3))

#input either tbl_trans or tbl_non_trans
ggplot(tbl_trans, aes(x=log(count), fill=group)) + geom_density(alpha=0.4) +
  labs( title ="Transcribed Strand", x = "log(Number of CC>TT mutations per Cell)") +
  theme(plot.title = element_text(hjust = 0.5)) +ylim(c(0,0.7))

ggplot(tbl_non_trans, aes(x=log(count), fill=group)) + geom_density(alpha=0.4) +
  labs( title ="non-transcribed Strand", x = "log(Number of CC>TT mutations per Cell)") +
  theme(plot.title = element_text(hjust = 0.5)) +ylim(c(0,0.7))

##PLOT-7##

#Ratio
ratio <- mutations_per_cell  %>%
  group_by(sample,strand) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  complete(sample, strand,
           fill = list(count = 0))

ratio_wide <- spread(ratio, strand, count)
ratio_wide <- ratio_wide %>% left_join(groups, by = "sample")
ratio_wide <- ratio_wide[,c(1,4,2,3)]

#calculate ratio and add 1 to both mutation counts
ratio_wide$ratio <- (ratio_wide$transcribed+1)/(ratio_wide$`non-transcribed`+1)

tbl_ratio <- ratio_wide[,c(2,5)]

ggplot(tbl_ratio, aes(x=log(ratio), fill=group)) + geom_density(alpha=0.4) +
  labs( title ="Ratio", x = "log(ratio) transcribed/non-transcribed") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,0.6)


ggplot(tbl_ratio, aes(x=ratio, fill=group)) + geom_density(alpha=0.4) +
  labs( title ="Ratio", x = "ratio transcribed/non-transcribed") +
  theme(plot.title = element_text(hjust = 0.5)) 

##PLOT-8##

##FISHER##

#Fisher.test plot
all_dint2 <- all_dint_tpm2 %>% group_by(sample, group) %>% summarise(CC_count=n())

fisher_out <- NULL
for (i in 2:max(all_dint2$CC_count)){
  all_dint2 <- all_dint2 %>% mutate(count = case_when(CC_count >= i ~ "high",CC_count < i ~ "low"))
  table(all_dint2$group, all_dint2$count)[c(1,2),] -> test1
  test1 <- test1[c(2,1), c(2,1)]
  res <- fisher.test(test1)
  fisher_out <- rbind(fisher_out, c(i, res$p.value))
  fisher_out <- as.data.frame(fisher_out)
  fisher_out <<- fisher_out
} 

colnames(fisher_out) <- c("count", "p_val")
fisher_out <- fisher_out[order(fisher_out$count),]

ggplot(fisher_out, aes(x=count, y = -log(p_val))) + geom_line() + geom_point() +
  labs( title ="Fisher Test (IgG vs antiPD1)", x = "mutation count threshold") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_x_continuous(breaks = round(seq(min(fisher_out$count), max(fisher_out$count), by = 1),1))+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 6)) 


##PLOT-9##-- additional

#number of cells with a specific gene CC>TT mutation
all_genes1 <-  all_dint_tpm2 %>% group_by(gene_name) %>% summarise(count=n())
all_genes1 <- all_genes1[all_genes1$count >2,]

all_genes <-  all_dint_tpm2 %>% group_by(gene_name, sample) %>% summarise(count=n())
all_genes <- all_genes %>% left_join(groups, by = "sample")

genes <- all_genes  %>%
  group_by(gene_name,group) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  complete(gene_name, group,
           fill = list(count = 0))

genes1 <- genes[genes$gene_name %in% all_genes1$gene_name,]

genes2 <- genes1  %>%
  group_by(gene_name,group) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  complete(gene_name, group,
           fill = list(count = 0))

ggplot(genes2, aes(x=reorder(gene_name,-count), y=count, fill=group)) + geom_bar(stat="identity", width =0.7)+ 
  theme_bw() + labs(x="gene name", y="Number of cells") + theme(axis.text.x = element_text(angle=60, size = 4, hjust=1)) +
  theme(legend.position="top") + labs(fill="") +ylim(c(0,23))




© 2020 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
