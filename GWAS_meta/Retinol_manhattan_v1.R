###################################

## Manhattan plot - retinol (METSIM+INTERVAL) & All three meta

## William Reay (2022)

##################################

library(CMplot)
library(data.table)
library(dplyr)

setwd("~/Desktop/Retinol_GWAS/")

MET_INT_input <- fread("Meta_analyses/Common_var/Stouffers_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz", header = T)

MET_INT_input <- MET_INT_input %>% select(rsids, chrom, pos, `P-value`) 

MET_INT_input <-  rename(MET_INT_input, "Chromosome"="chrom", "SNP"="rsids",
                      "Position"="pos", "P"="P-value")

SNPs <- as.vector(c("rs1260326", "rs34898035", "rs11762406", "rs6601299",
                    "rs10882283", "rs12149203", "rs1667226", "rs6029188"))


Closest_genes <- as.vector(c("GCKR", "TFCP2L1", "FOXP2", "PP1R3B", "RBP4",
                             "MAF", "TTR", "MAFB"))

CMplot(MET_INT_input, plot.type="m",LOG10=TRUE,col=c("grey30","steelblue1"),highlight=SNPs, highlight.type="p",
       highlight.col=c("firebrick1","firebrick1", "firebrick1", "firebrick1", "firebrick1",
                       "firebrick1", "firebrick1", "firebrick1"),
       ,highlight.text=Closest_genes,      
       highlight.text.col=c("firebrick1","firebrick1", "firebrick1", "firebrick1", "firebrick1",
                       "firebrick1", "firebrick1", "firebrick1"),threshold=5e-08,threshold.lty=2,   
       amplify=FALSE,file="tiff",memo="",dpi=1000,file.output=TRUE,
       verbose=TRUE,width=16,height=6)

## All three 

All_three <- fread("Meta_analyses/Common_var/Stouffers_Full_METSIM_INTERVAL_PLCO_ATBC1.tbl")

All_three <-  rename(All_three, "rsids"="MarkerName", "P"="P-value")

MET_INT_input <- fread("Meta_analyses/Common_var/Stouffers_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz", header = T)

MET_INT_input <- MET_INT_input %>% select(rsids, chrom, pos) 

Merged <- merge(All_three, MET_INT_input, by = "rsids")
Merged$SNP <- Merged$rsids

Merged <- Merged %>% select(SNP, chrom, pos, P) 
Merged <- rename(Merged, "CHR"="chrom", "BP"="pos")

## Make Manhattan

SNPs <- as.vector(c("rs1260326", "rs6601299", "rs11187547", "rs11865979",
                    "rs4799581", "rs6029188"))


Closest_genes <- as.vector(c("GCKR", "PPP1R3B", "RBP4",
                             "MAF", "TTR", "MAFB"))

CMplot(Merged, plot.type="m",LOG10=TRUE,col=c("grey30","steelblue1"),highlight=SNPs, highlight.type="p",
       highlight.col=c("firebrick1","firebrick1", "firebrick1", "firebrick1"),
       highlight.text=Closest_genes,      
       highlight.text.col=c("firebrick1","firebrick1", "firebrick1", "firebrick1"),threshold=5e-08,threshold.lty=2,   
       amplify=FALSE,file="tiff",memo="",dpi=1000,file.output=TRUE,
       verbose=TRUE,width=16,height=6)

