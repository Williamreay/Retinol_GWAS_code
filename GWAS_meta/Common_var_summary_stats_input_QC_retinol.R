##########################################

## Retinol GWAS meta - summary statistics QC and harmonisation

## COMMON VARIANTS

## William Reay (2022)

##########################################

setwd("~/data/users/william/Meta_retinol_GWAS/")

library(dplyr)
library(data.table)

## INTERVAL raw - 11132 is max number of participants

INTERVAL_raw <- fread("Raw_sumstats/INT-WGS-gwas-metabolon-metabolon_1806.tsv.gz")

## Retain only common variants in INTERVAL

INTERVAL_common <- INTERVAL_raw %>% filter(AF > 0.01 & AF < 0.99)
## Make a locus_ref_alt column for merging
INTERVAL_common$CHR_BP_REF_ALT <- paste(INTERVAL_common$locus,":", 
                                      INTERVAL_common$REF, ":", INTERVAL_common$ALT, sep="")

## Check for any missing test statistics

INTERVAL_common %>% filter(is.na(beta))

## Read in METSIM data and filter for common variants

METSIM_raw <- fread("Raw_sumstats/METSIM_retinol_GWAS_C498.txt.gz")

## Retain only common variants in METSIM

METSIM_common <- METSIM_raw %>% filter(maf > 0.01 & maf < 0.99)
## Make a locus_ref_alt column for merging
METSIM_common$CHR_BP_REF_ALT <- paste(METSIM_common$chrom,":", METSIM_common$pos, ":", 
                                      METSIM_common$ref, ":", METSIM_common$alt, sep="")
METSIM_common$CHR_BP_REF_ALT <- paste("chr", METSIM_common$CHR_BP_REF_ALT, sep="")

## Get variant ID list between two studies with matching IDs and alleles - remove SNPs without an ID

INTERVAL_SNP_info <- INTERVAL_common %>% select(CHR_BP_REF_ALT, rsid, REF, ALT)

METSIM_SNP_info <- METSIM_common %>% select(CHR_BP_REF_ALT, rsids, ref, alt)

## Merge variant lists together

METSIM_INTERVAL_merge <- merge(INTERVAL_SNP_info, METSIM_SNP_info, by = "CHR_BP_REF_ALT")

## Check if any discordant rsids

METSIM_INTERVAL_merge %>% filter(rsid != rsids)
METSIM_INTERVAL_merge %>% filter(is.na(rsid))

## Select METSIM rsids as precedent rsIDs

METSIM_INTERVAL_merge <- METSIM_INTERVAL_merge %>% select(CHR_BP_REF_ALT, rsids)
METSIM_INTERVAL_merge %>% filter(!duplicated(CHR_BP_REF_ALT))

## Extract these ids from INTERVAL

INTERVAL_common_harmonised <- merge(METSIM_INTERVAL_merge, INTERVAL_common, by = "CHR_BP_REF_ALT")
INTERVAL_common_harmonised <- INTERVAL_common_harmonised %>% select(-rsid)

## Output
write.table(INTERVAL_common_harmonised, file="Common_var/INTERVAL_METSIM_meta/INTERVAL_harmonised_INTERVAL_METSIM_COMMON.txt",
            sep = "\t", row.names = F, quote = F)


## Extract these ids from METSIM

METSIM_common_harmonised <- merge(METSIM_INTERVAL_merge, METSIM_common, by = "CHR_BP_REF_ALT")
METSIM_common_harmonised <- METSIM_common_harmonised %>% select(-rsids.y)
METSIM_common_harmonised <- rename(METSIM_common_harmonised, "rsids"="rsids.x")

## Output
write.table(METSIM_common_harmonised, file="Common_var/INTERVAL_METSIM_meta/METSIM_harmonised_INTERVAL_METSIM_COMMON.txt",
            sep = "\t", row.names = F, quote = F)

## Next meta analysis - all three sumstats ATBC/PLCO included ##

PLCO_ATBC_raw <- fread("Raw_sumstats/Over_0.8_FIZI_new_1000G_march_PLCO_ATBC.txt")
PLCO_ATBC_raw <- rename(PLCO_ATBC_raw, "rsids"="SNP")

PLCO_ATBC_SNPs <- PLCO_ATBC_raw %>% select(rsids, A_EFFECT, A_ALT)

## Check for overlap with METSIM_INTERVAL harmonised set

Merged_all_SNPs <- merge(PLCO_ATBC_SNPs, METSIM_INTERVAL_merge, by ="rsids")

## Extract this harmonised SNP set from all three

Merged_all_SNPs <- Merged_all_SNPs %>% select(rsids)
Merged_all_SNPs <- Merged_all_SNPs %>% filter(!duplicated(rsids))

## Extract from PLCO_ATBC

PLCO_ATBC_all_harmonised <- merge(Merged_all_SNPs, PLCO_ATBC_raw, by = "rsids")

write.table(PLCO_ATBC_all_harmonised, file="Common_var/Full_meta/Full_meta_common_PLCO_ATBC.txt",
            sep = "\t", row.names = F, quote = F)

## Extract from METSIM

METSIM_all_harmonised <- merge(Merged_all_SNPs, METSIM_common_harmonised, by = "rsids")

write.table(METSIM_all_harmonised, file="Common_var/Full_meta/Full_meta_common_METSIM.txt",
            sep = "\t", row.names = F, quote = F)

## Extract from INTERVAL

INTERVAL_all_harmonised <- merge(Merged_all_SNPs, INTERVAL_common_harmonised, by = "rsids")

write.table(INTERVAL_all_harmonised, file="Common_var/Full_meta/Full_meta_common_INTERVAL.txt",
            sep = "\t", row.names = F, quote = F)
