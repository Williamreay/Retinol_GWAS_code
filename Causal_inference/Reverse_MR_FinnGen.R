####################################

## Reverse MR - FinnGen

## William Reay (2023)

####################################

library(dplyr)
library(data.table)
library(TwoSampleMR)

setwd("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/FinnGen_r8_MR/")

Ret_GWAS <- fread("../../Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

Ret_GWAS$Allele1 <- toupper(Ret_GWAS$Allele1)
Ret_GWAS$Allele2 <- toupper(Ret_GWAS$Allele2)


## Function to format

IV_func <- function(df, Name) {
  df_dat <- df
  
  ## Approximate SE
  
  df_dat$Z <- sign(df_dat$beta) * sqrt(qchisq(df_dat$pval, 1, lower =F))
  df_dat$SE <- abs(df_dat$beta/df_dat$Z)
  
  ## Add name
  
  df_dat$Name <- Name
  
  ## Format as exposure
  
  Exp <- format_data(df_dat, type = "exposure", phenotype_col = "Name",
                      beta_col = "beta", snp_col = "rsids",
                     se_col = "SE", pval_col = "pval",
                     effect_allele_col = "alt", other_allele_col = "ref")
  
  Clump <- clump_data(Exp, clump_p1 = 5e-08)
  
  return(Clump)
}

## Read in IVs for binary exposures

DM <- fread("Reverse_causality/finngen_R8_E4_DM2COMA_lead.tsv")

DM_IV <- IV_func(DM, "Diabetes_with_coma")

LIV <- fread("Reverse_causality/finngen_R8_K11_OTHINFLIV_lead.tsv")

LIV_IV <- IV_func(LIV, "Inflammatory_liver")

Erys <- fread("Reverse_causality/finngen_R8_AB1_ERYSIPELAS_lead.tsv")

Erys_IV <- IV_func(Erys, "Erysipelas")

Congen <- fread("Reverse_causality/finngen_R8_CONGEN_HEART_ARTER_lead.tsv")

Congen_IV <- IV_func(Congen, "Congential_artery_malformations")

Merge <- rbind(DM_IV, LIV_IV, Erys_IV, Congen_IV)

## Format outcome

Formatted_ret <- format_data(Ret_GWAS,
                             type = "outcome",
                             snps = Merge$SNP,
                             beta_col = "Effect",
                             se_col = "StdErr",
                             effect_allele_col = "Allele1", 
                             other_allele_col = "Allele2",
                             pval_col = "P-Value",
                             snp_col = "rsids")

Harm_df <- harmonise_data(Merge, Formatted_ret, action = 2)

## Run MR IVW-FE and single snp

MR_all_Fin <- mr(Harm_df, method_list = c("mr_ivw_mre"))

write.table(MR_all_Fin, file = "Reverse_causality/GW_sig_reverse_causality.txt",
            sep = "\t", row.names = F, quote = F)
