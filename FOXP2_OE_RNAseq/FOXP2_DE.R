########################################

## FOXP2 overexpression in U2OS (Sacroma)

## William Reay (2023)

########################################

library(edgeR)
library(data.table)
library(dplyr)
library(tidyverse)

setwd("~/Desktop/Retinol_GWAS/Gene_priortisation/FOXP2_differential_expression/")

DF <- fread("GSE138938_annotated_combined.counts.txt", header = T)
DF_filt <- DF %>% select(id, Control_1, Control_2,
                         Control_3, Control_4, Control_5,
                         FOXP2_WT_1, FOXP2_WT_2, FOXP2_WT_3,
                         FOXP2_WT_4, FOXP2_WT_5)

x <- DF_filt %>% remove_rownames %>% column_to_rownames(var="id")

x <- as.matrix(x)

## Define groups for subsequent analyses

group<-c(rep("C",5),rep("T",5))
y<-DGEList(counts=x, group=group)

## MDS plot

cn.col="blue"
tr.col="green"
cols=c(rep(cn.col,5),rep(tr.col,5))
plotMDS(y,pch = 16,col=cols,las=1)
plotMDS(y,pch = 16,col=cols,las=1, labels = colnames(y$counts))

## Remove genes with < 10 counts in the lowest sample 

## Lowest library size was 15755842, therefore divide 10/15.76 (in millions) = 

keep <- rowSums(cpm(y)>0.4108463) >= 5
y2 <- y[keep,,keep.lib.sizes=FALSE]
y2$samples

plotMDS(y2,pch = 16,col=cols,las=1, labels = colnames(y$counts))
 

## Scale the raw library sizes and estimate neg binomial dispersion (empirical Bayes)
y3 <-calcNormFactors(y2)
y3 <- estimateDisp(y3)

plotBCV(y3)

## Perform differential expression (exact test for data ~ neg binomial)

et.1 <- exactTest(y3, c("C", "T"))
topTags(et.1)
results <- topTags(et.1, n=dim(y3)[1]+1, adjust.method = "BH",sort.by = "PValue")

write.table(results, file="FOXP2_DE_all_samples.txt", sep="\t", row.names = T, quote = F)

DE <- fread("FOXP2_DE_all_samples.txt")

## Volcano plot

DE$threshold <- ifelse((DE$FDR < 0.01 & abs(DE$logFC) > 0.6), 1, 0)

DE$threshold <- as.factor(DE$threshold)

ggplot(data = DE,
       aes(x=logFC, y = -log10(FDR), colour = threshold)) +
  geom_point(size = 0.9) +
  labs(x = expression("log FC"), y = expression(paste("-log"[10], "FDR"))) +
  theme_bw() +
  theme(legend.position ="none") +
  theme(legend.position ="none") +
  theme(axis.title = element_text(face="bold", size=12)) +
  geom_vline(xintercept = -0.6 , linetype="longdash") +
  geom_vline(xintercept = 0.6 , linetype="longdash") +
  geom_hline(yintercept=2, linetype="longdash") +
  scale_color_manual(name="threshold", values=c('#999999', 'mediumpurple1')) + 
  ggtitle("FOXP2 OE vs CTRL")

## Stricter VP

DE$threshold_2 <- ifelse((DE$PValue < 2.854044e-06 & abs(DE$logFC) > 2), 1, 0)

DE$threshold_2 <- as.factor(DE$threshold_2)

ggplot(data = DE,
       aes(x=logFC, y = -log10(PValue), colour = threshold_2)) +
  geom_point(size = 0.9) +
  labs(x = expression("log FC"), y = expression(paste("-log"[10], "P-value"))) +
  theme_bw() +
  theme(legend.position ="none") +
  theme(legend.position ="none") +
  theme(axis.title = element_text(face="bold", size=12)) +
  geom_vline(xintercept = -2 , linetype="longdash") +
  geom_vline(xintercept = 2 , linetype="longdash") +
  geom_hline(yintercept=5.54, linetype="longdash") +
  scale_color_manual(name="threshold_2", values=c('#999999', 'dodgerblue1')) + 
  ggtitle("FOXP2 OE vs CTRL")

## Above to submit to gprofiler

DE$threshold_3 <- ifelse((DE$PValue < 2.854044e-06 & abs(DE$logFC) > 1.5), 1, 0)

DE$threshold_3 <- as.factor(DE$threshold_3)

DE_thresh_2 <- DE %>% filter(threshold_3 == 1) %>%
  select(V1)

write.table(DE_thresh_4, file="Bonf_log2FC_1.5_DE_FOXP2_all_Samples.txt",
            row.names = F, quote = F)

DE$threshold_4 <- ifelse((DE$PValue < 2.854044e-06 & abs(DE$logFC) > 5), 1, 0)

DE$threshold_4 <- as.factor(DE$threshold_4)

DE_thresh_3 <- DE %>% filter(threshold_4 == 1) %>%
  select(V1)

write.table(DE_thresh_3, file="Bonf_log2FC_5_DE_FOXP2_all_Samples.txt",
            row.names = F, quote = F)

## Remove FOXP2 replicates 1 and 3 as less well clustered

x_filt <- x[,c(1:5,7,9:10)]

group_2<-c(rep("C",5),rep("T",3))
y_2<-DGEList(counts=x_filt, group=group_2)

## MDS plot

cn.col="blue"
tr.col="green"
cols=c(rep(cn.col,5),rep(tr.col,3))
plotMDS(y_2,pch = 16,col=cols,las=1)
plotMDS(y_2,pch = 16,col=cols,las=1, labels = colnames(y_2$counts))

## Remove genes with < 10 counts in the lowest sample 

## Lowest library size was 15755842, therefore divide 10/15.76 (in millions) = 

keep_2 <- rowSums(cpm(y)>0.4108463) >= 5
y2_2 <- y_2[keep_2,,keep.lib.sizes=FALSE]
y2_2$samples

plotMDS(y2_2,pch = 16,col=cols,las=1, labels = colnames(y_2$counts))


## Scale the raw library sizes and estimate neg binomial dispersion (empirical Bayes)
y3_2 <-calcNormFactors(y2_2)
y3_2 <- estimateDisp(y3_2)

plotBCV(y3_2)

et.1_2 <- exactTest(y3_2, c("C", "T"))
topTags(et.1_2)
results_2 <- topTags(et.1_2, n=dim(y3_2)[1]+1, adjust.method = "BH",sort.by = "PValue")

