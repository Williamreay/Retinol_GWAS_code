##########################################

## Retinol GWAS meta - summary statistics QC and harmonisation

## RARE VARIANTS

## William Reay (2022)

##########################################

setwd("~/data/users/william/Meta_retinol_GWAS/")

library(dplyr)
library(data.table)

## INTERVAL raw - 11132 is max number of participants

INTERVAL_raw <- fread("Raw_sumstats/INT-WGS-gwas-metabolon-metabolon_1806.tsv.gz")

## Retain only rare variants in INTERVAL

INTERVAL_rare <- INTERVAL_raw %>% filter(AF < 0.01 | AF > 0.99)
## Make a locus_ref_alt column for merging
INTERVAL_rare$CHR_BP_REF_ALT <- paste(INTERVAL_rare$locus,":", 
                                        INTERVAL_rare$REF, ":", INTERVAL_rare$ALT, sep="")
INTERVAL_rare <- rename(INTERVAL_rare, "rsids"="rsid")

## Check for any missing test statistics

INTERVAL_rare %>% filter(is.na(beta))

## Read in METSIM data and filter for common variants

METSIM_raw <- fread("Raw_sumstats/METSIM_retinol_GWAS_C498.txt.gz")

## Retain only rare variants in METSIM
METSIM_rare <- METSIM_raw %>% filter(maf < 0.01 | maf > 0.99)
## Make a locus_ref_alt column for merging
METSIM_rare$CHR_BP_REF_ALT <- paste(METSIM_rare$chrom,":", METSIM_rare$pos, ":", 
                                      METSIM_rare$ref, ":", METSIM_rare$alt, sep="")
METSIM_rare$CHR_BP_REF_ALT <- paste("chr", METSIM_rare$CHR_BP_REF_ALT, sep="")

## Get variant ID list between two studies with matching IDs and alleles - remove SNPs without an ID

INTERVAL_SNP_info <- INTERVAL_rare %>% select(CHR_BP_REF_ALT)

METSIM_SNP_info <- METSIM_rare %>% select(CHR_BP_REF_ALT, rsids)

## Merge variant lists together

METSIM_INTERVAL_merge <- merge(INTERVAL_SNP_info, METSIM_SNP_info, by = "CHR_BP_REF_ALT")

## Extract from METSIM

METSIM_rare <- METSIM_rare %>% select(-rsids)

METSIM_rare_harmonised <- merge(METSIM_INTERVAL_merge, METSIM_rare, by="CHR_BP_REF_ALT")
METSIM_rare_harmonised <- METSIM_rare_harmonised %>% filter(!duplicated(CHR_BP_REF_ALT))

write.table(METSIM_rare_harmonised, file="Rare_var/Rare_var_METSIM_harmonised.txt",
            sep = "\t", row.names = F, quote = F)

## Extract from INTERVAL

INTERVAL_rare_harmonised <- merge(METSIM_INTERVAL_merge, INTERVAL_rare, by="CHR_BP_REF_ALT")
INTERVAL_rare_harmonised <- INTERVAL_rare_harmonised %>% filter(!duplicated(CHR_BP_REF_ALT))

write.table(INTERVAL_rare_harmonised, file="Rare_var/Rare_var_INTERVAL_harmonised.txt",
            sep = "\t", row.names = F, quote = F)

