####################################

## Reverse MR - tier #2 and tier #3 traits

## William Reay (2023)

####################################

library(dplyr)
library(data.table)
library(TwoSampleMR)

setwd("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/")

Tier_3 <- fread("Processed_output/Tier_3/Passing_traits_tier_3.txt")

Tier_3_IDs <- as.data.frame(unique(Tier_3$id.outcome))

Exp_IVs <- extract_instruments(outcomes = Tier_3_IDs$`unique(Tier_3$id.outcome)`)

table(Exp_IVs$id.exposure)

## Read in retinol GWAS

Ret_GWAS <- fread("../Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

Ret_GWAS$Allele1 <- toupper(Ret_GWAS$Allele1)
Ret_GWAS$Allele2 <- toupper(Ret_GWAS$Allele2)

Formatted_ret <- format_data(Ret_GWAS,
                            type = "outcome",
                            snps = Exp_IVs$SNP,
                            beta_col = "Effect",
                            se_col = "StdErr",
                            effect_allele_col = "Allele1", 
                            other_allele_col = "Allele2",
                            pval_col = "P-Value",
                            snp_col = "rsids")

## Harmonise and MR

Harm_df <- harmonise_data(Exp_IVs, Formatted_ret, action =3)

MR_all <- mr(Harm_df, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_wald_ratio"))

mr_forest_plot(MR_all)

write.table(MR_all, file = "~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/Reverse_causality_for_cis_acting_SNPs.txt",
            sep = "\t", row.names = F)

MR_all %>% filter(pval < 0.05)

SS <- mr_leaveoneout(Harm_df)
