### author: yao zhan
### date: 02/04/2020

#install.packages("gdata")
library(gdata)
mela = read.xls("../Melanocyte_shXLOC006189.xlsx", sheet = 1, stringsAsFactors = F, header = T)
UACC257 = read.xls("..//UACC257_shXLOC006189.xlsx", sheet = 1, stringsAsFactors = F, header = T)

mela$threshold = ifelse(mela$logFC > 1.5 & mela$FDR<0.05 | mela$logFC < (-1.5) & mela$FDR<0.05 , "true", "false" )
UACC257$threshold = ifelse(UACC257$logFC > 1.5 & UACC257$FDR<0.05 | UACC257$logFC < (-1.5) & UACC257$FDR<0.05 , "true", "false" )
saveRDS(mela, "../Melanocyte_shXLOC006189.rds")
saveRDS(UACC257, "../UACC257_shXLOC006189.rds")
mela = readRDS("../Melanocyte_shXLOC006189.rds")
UACC257 = readRDS("../UACC257_shXLOC006189.rds")

library(dplyr)
mela_DEG = mela %>% filter(threshold == "true")
UACC257_DEG = UACC257 %>% filter(threshold == "true")
write.table(mela_DEG, "../melanocyte.DEG.cvs", quote = F, sep = "\t")
write.table(UACC257_DEG, "../UACC257_DEG.cvs", quote = F, sep = "\t")
# Valcano plot
library(ggplot2)
library(ggrepel)
pdf("../Volcano plots.melanotytic markers.pdf")
gene = c("TYR", "TRPM1", "DCT", "TYRP1", "CDK2", "BCL2", "BCL2A", "PMEL")
ggplot(mela, aes(logFC, -log(FDR), label=hgnc_symbol, color = threshold)) + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = -log(0.05), color = 'red', linetype = "dashed") + 
  geom_vline(xintercept = 1.5, color = "green", linetype = "dashed") + 
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "green") + 
  #geom_text(aes(label = ifelse(hgnc_symbol %in% gene , as.character(hgnc_symbol), '')), hjust=0, vjust=0, size = 2.5) + 
  geom_text_repel(aes(label = ifelse(hgnc_symbol %in% gene , as.character(hgnc_symbol), '')), hjust=0, vjust=1.1, size =2.5) + 
  scale_color_manual(values = c("grey","black")) +
  labs(title = "Melanocyte_shXLOC006189\n", color = "") +
  theme(legend.position = "none")
  
#install.packages("ggrepel")
ggplot(UACC257, aes(logFC, -log(FDR), label=hgnc_symbol, color = threshold)) + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = -log(0.05), color = 'red', linetype = "dashed") + 
  geom_vline(xintercept = 1.5, color = "green", linetype = "dashed") + 
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "green") + 
  scale_color_manual(values = c("grey","black")) +
  labs(title = "UACC257_shXLOC006189\n", color = "") +
  theme(legend.position = "none") +
  geom_text_repel(aes(label = ifelse(hgnc_symbol %in% gene, as.character(hgnc_symbol), '')), hjust=0, vjust=1.1, size=2.5)
  #geom_text_repel(aes(label = ifelse(logFC>1.5 & -log(FDR) > 30 , as.character(hgnc_symbol), '')), hjust=1, vjust=1, size = 2.5) + 
  #geom_text_repel(aes(label = ifelse(logFC < (-1.5) & -log(FDR) > 30 , as.character(hgnc_symbol), '')), hjust=0, vjust=1.1, size = 2.5) 
dev.off()
