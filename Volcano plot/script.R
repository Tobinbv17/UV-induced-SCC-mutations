### author: yao zhan
### date: 02/04/2020



#install.packages("gdata")
library(gdata)
mela = read.xls("Documents/2020 Phillip/Melanocyte_shXLOC006189.xlsx", sheet = 1, stringsAsFactors = F, header = T)
UACC257 = read.xls("Documents/2020 Phillip/UACC257_shXLOC006189.xlsx", sheet = 1, stringsAsFactors = F, header = T)

mela$threshold = ifelse(mela$logFC > 1.5 & -log(mela$FDR) >2 | mela$logFC < (-1.5) & -log(mela$FDR) >2 , "true", "false" )
UACC257$threshold = ifelse(UACC257$logFC > 1.5 & -log(UACC257$FDR) >2 | UACC257$logFC < (-1.5) & -log(UACC257$FDR) >2 , "true", "false" )
saveRDS(mela, "Documents/2020 Phillip/Melanocyte_shXLOC006189.rds")
saveRDS(UACC257, "Documents/2020 Phillip/UACC257_shXLOC006189.rds")

mela_DEG = mela %>% filter(threshold == "true")
UACC257_DEG = UACC257 %>% filter(threshold == "true")
write.table(mela_DEG, "Documents/2020 Phillip/melanocyte.DEG.cvs", quote = F, sep = "\t")
write.table(UACC257_DEG, "Documents/2020 Phillip/UACC257_DEG.cvs", quote = F, sep = "\t")
# Valcano plot
pdf("Documents/2020 Phillip/Volcano plots.pdf")
ggplot(mela, aes(logFC, -log(FDR), label=hgnc_symbol, color = threshold)) + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = 2, color = 'red', linetype = "dashed") + 
  geom_vline(xintercept = 1.5, color = "green", linetype = "dashed") + 
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "green") + 
  geom_text(aes(label = ifelse(logFC>2.5 & -log(FDR) > 50 , as.character(hgnc_symbol), '')), hjust=0, vjust=0, size = 2.5) + 
  geom_text_repel(aes(label = ifelse(logFC < (-3) & -log(FDR) > 50 , as.character(hgnc_symbol), '')), hjust=0, vjust=1.1, size = 2.5) + 
  scale_color_manual(values = c("grey","black")) +
  labs(title = "Melanocyte_shXLOC006189\n", color = "") +
  theme(legend.position = "none")
  
#install.packages("ggrepel")
#library(ggrepel)
ggplot(UACC257, aes(logFC, -log(FDR), label=hgnc_symbol, color = threshold)) + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = 2, color = 'red', linetype = "dashed") + 
  geom_vline(xintercept = 1.5, color = "green", linetype = "dashed") + 
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "green") + 
  scale_color_manual(values = c("grey","black")) +
  labs(title = "UACC257_shXLOC006189\n", color = "") +
  theme(legend.position = "none") + 
  geom_text_repel(aes(label = ifelse(logFC>1.5 & -log(FDR) > 30 , as.character(hgnc_symbol), '')), 
                  hjust=1, vjust=1, size = 2.5) + 
  geom_text_repel(aes(label = ifelse(logFC < (-1.5) & -log(FDR) > 30 , as.character(hgnc_symbol), '')), 
                  hjust=0, vjust=1.1, size = 2.5) 
dev.off()

