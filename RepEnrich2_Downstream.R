

rm(list=ls())

# Input data ####

files = list.files(pattern = "*__fraction_counts.txt")
myfiles = do.call(cbind, 
                  lapply(files, 
                         function(x) read.delim(x, 
                                                header = FALSE,
                                                col.names = rep(substr(x,1,nchar(x)-30),4),
                                                stringsAsFactors = FALSE)))#1267 36
info <- myfiles[,c(1:3)]
rownames(info) <- info[,1]
info <- info[,-1]
colnames(info) <- c('name1','name2')

rownames(myfiles) <- myfiles[,1]
myfiles <- myfiles[,c(seq(from=4,to=24,by=4))]

#####

# DEG analysis - DESeq ####

library('DESeq2')
input_data <- myfiles[,c(1,3,5,4,6)]
condition <- factor(c(rep("DMSO",3),rep("MS023_5d",2)))
coldata <- data.frame(row.names = colnames(input_data),condition)
#DESeq input matrix
dds <- DESeqDataSetFromMatrix(countData=input_data,colData=coldata, 
                              design=~condition)

dds <- DESeq(dds)
res <- results(dds,contrast = c("condition","MS023_5d","DMSO"))
summary(res)
#order：根据padj列排序
res <- res[order(res$padj),]
res <- as.data.frame(res)

out <- merge(res,input_data,by = 'row.names')
rownames(out) <- out$Row.names
out <- merge(info,out,by = 'row.names')
out <- out[,-4]
write.csv(out,"Repeats_MS023_5d_DMSO_DESeq2.csv",row.names = F)

#####

# Heatmap ####
out <- read.csv('Repeats_MS023_5d_DMSO_DESeq2.csv',row.names = 1,check.names = F)#1267
out <- out[order(out$name1),]
out$name1 <- gsub('\\?','',out$name1)
table(out$name1)

hm <- out[which(out$name1 %in% c('LINE','SINE','LTR','Satellite')),]#837
colnames(hm)[1] <- 'Repeat Element'
hm$DMSO <- 0
hm$MS023 <- 0
for (i in 1:nrow(hm)) {
  hm$DMSO[i] = (mean(as.numeric(hm[i,9:11]))+1)
  hm$MS023[i] = (mean(as.numeric(hm[i,12:13]))+1)
}
library(pheatmap)
pdf(file = 'Heatm_Repeat Element expr of MDAMB468 MS023 5d.pdf',width = 6,height = 8)
pheatmap(hm[,14:15],
         show_colnames = T,show_rownames = F,
         main = 'Repeat Element expr of MDAMB468 MS023 5d',
         annotation_row = hm[1],gaps_col = 1,
         cluster_rows = T,cluster_cols = F)
dev.off()
#####

# Volcanol plot ####
out <- read.csv('Repeats_MS023_5d_DMSO_DESeq2.csv',row.names = 1,check.names = F)#1267
out <- out[order(out$name1),]
out$name1 <- gsub('\\?','',out$name1)
out$significant <- as.factor(out$padj<0.05 & abs(out$log2FoldChange) > 1)
out$symbol <- rownames(out)

library(ggplot2)
library(ggrepel)
pdf(file = 'VolcP_Repeats Elements of MDAMB468 MS023.pdf',width = 6,height = 6)
ggplot(data=out, aes(x=log2FoldChange, y =-log10(padj),color =significant)) +
  geom_point() +
  scale_color_manual(values =c("black","red"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme_classic()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  )+
  labs(title="Repeats Elements of MDAMB468 MS023", x="log2 (fold change)",y="-log10 (padj)")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#####

# Author : Ji Yishuai .
