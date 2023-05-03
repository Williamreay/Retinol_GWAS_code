######################

## Intersecting genome-wide sig SNPs with GTEX eQTL finemapping results (DAP-G, PIP > 0.95)

## William Reay (2022)

#######################

setwd("~/Desktop/Retinol_GWAS/")

library(data.table)
library(dplyr)

## METSIM+INTERVAL - DAP-G

INTERVAL_METSIM_GW_sig <- fread("Meta_analyses/METSIM_INTERVAL_stouffer_GTEx_IDS_GW_sig.txt", header = F)
INTERVAL_METSIM_GW_sig$V3 <- INTERVAL_METSIM_GW_sig$V1

Finemapped_eQTL <- fread("Meta_analyses/FUMA/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.CS95.txt")

Merged_INT_MET <- merge(Finemapped_eQTL, INTERVAL_METSIM_GW_sig, by = "V3")

write.table(Merged_INT_MET, file="Meta_analyses/FUMA/GTEx_v8_finemapping_DAPG/INTERVAL_METSIM_CS95.txt",
            sep = "\t", row.names = F, quote = F)


