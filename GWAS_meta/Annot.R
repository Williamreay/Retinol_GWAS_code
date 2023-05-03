##########################################

## Retinol GWAS meta - annotate METSIM+INTERVAL with dbSNP

## William Reay (2022)

##########################################

library(dplyr)
library(data.table)

setwd("~/data/users/william/Meta_retinol_GWAS")

MET_IVW <- fread("~/data/users/william/Meta_retinol_GWAS/Meta_analysis_output/IVW_METSIM_INTERVAL1.tbl")
MET_IVW <- MET_IVW %>% select(-dbSNP)

INT <- fread("~/data/users/william/Meta_retinol_GWAS/Common_var/INTERVAL_METSIM_meta/METSIM_harmonised_INTERVAL_METSIM_COMMON.txt")
INT <- INT %>% select(CHR_BP_REF_ALT, rsids, chrom, pos)
INT <- rename(INT, "MarkerName"="CHR_BP_REF_ALT")

Merge_IVW <- merge(MET_IVW, INT, by = "MarkerName")

write.table(Merge_IVW, file="Meta_analysis_output/Annotated_sumstats/IVW_METSIM_INTERVAL_sumstats_retinol.txt",
            sep = "\t", row.names = F, quote = F)

MET_ST <- fread("~/data/users/william/Meta_retinol_GWAS/Meta_analysis_output/Stouffers_METSIM_INTERVAL1.tbl")
MET_ST <- MET_ST %>% select(-dbSNP)

Merge_ST <- merge(MET_ST, INT, by = "MarkerName")

write.table(Merge_ST, file="Meta_analysis_output/Annotated_sumstats/Stouffers_METSIM_INTERVAL_sumstats_retinol.txt",
            sep = "\t", row.names = F, quote = F)
