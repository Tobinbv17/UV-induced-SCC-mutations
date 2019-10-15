require(data.table)
require(ggplot2)

####### binomial modelling    mod2 <- glm(diabetic ~ age + frame + chol, data=diabetes2, family="binomial")
####### linear modelling    mod1 <- lm(glyhb ~ age + frame + chol, data=diabetes2)

#################################### Numbers of CC-TT mutations in all groups
ggplot(snp_UV_dint, aes(group, fill= str_CBD)) + 
  geom_bar(position = "dodge") + 
  labs(x = "Numbers of UV mutations", y = "Counts", fill = "Strand")
##################################################################################################
##### Fisher test on # of cells with 1:30 UV mutations 

df = data.frame(threshold = numeric(), fisher_test = numeric())
for (i in 1:29) {
  Occurrance = snp_UV_dint %>% group_by(sample, group) %>% summarize(count = n())
  test = Occurrance %>% mutate(count = case_when(count > i ~ "high", count <= i ~ " low"))
  test1 = as.data.frame.matrix(table(test$group, test$count))
  
  df[i,] = data.frame(threshold = i, p_value = fisher.test(test1[2:3,])$p.value)
}

################### Number of cells pass threshold
df_number1 = data.frame(var1 = character(), var2 = character(), var3 = numeric(), var4 = numeric())
for (i in 1:30) {
  Occurrance = snp_UV_dint %>% group_by(sample, group) %>% summarize(count = n())
  test = Occurrance %>% mutate(count = case_when(count > i ~ "high", count <= i ~ " low"))
  test1 = as.data.frame(table(test$group, test$count))
  test1$var4 = rep(i, nrow(test1))
  
  df_number1 = rbind(df_number1, test1)
}

colnames(df_number1) = c("group", "categories", "count", "threshold")

####################### log2 transform in ggplot
ggplot(df, aes(threshold, fisher_test)) + 
  geom_bar(stat = "identity") + 
  scale_x_continuous(breaks = c(1:30)) + 
  scale_y_continuous(trans = "log2")

#################### -log2 transformation in dataframe first
df$log2_f.test = -log2(df$fisher_test)

ggplot(df, aes(threshold, log2_f.test)) + 
  geom_bar(stat = "identity") + 
  scale_x_continuous(breaks = c(1:30)) + 
  labs(x = "Threshold (Number of UV mutations per Cell)" , y = "-log2(p_Value)(Fisher.test)") + 
  geom_hline(yintercept = 4.3, linetype = "dashed", size = 1, color = "red") + 
  ylim(0,7)

########################### log 10 transformation
ggplot(df, aes(threshold, -log10(fisher_test))) + 
  geom_bar(stat = "identity", fill = "orange") +
  scale_x_continuous(breaks = c(1:28)) + 
  labs(x = "Number of UV mutations per Cell" , y = "-log10(p_Value)(Fisher.test)") + 
  geom_hline(yintercept = 1.3, linetype = "dashed", size = 1, color = "black") + 
  ylim(0,2.5)


Occurrance10 = snp_UV_dint %>% group_by(sample) %>% summarize(count = n() >= 10) %>% filter(count == TRUE)
df10 = snp_UV_dint[snp_UV_dint$sample %in% Occurrance10$sample,]
table(df10$group, df10$str_CBD)


######## N of CC-TT in nonUV
ggplot(a1[(a1$group=="XPCeenonUV"),], aes(sample, fill = str_CBD)) + geom_bar() + theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust = 0, size = 8))
ggplot(a1[(a1$group=="XPCeeUVIgG"),], aes(sample, fill = str_CBD)) + geom_bar() + theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust = 0, size = 5))
ggplot(a1[(a1$group=="XPCeeUVantiPD1"),], aes(sample, fill = str_CBD)) + geom_bar() + theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust = 0, size = 5))


