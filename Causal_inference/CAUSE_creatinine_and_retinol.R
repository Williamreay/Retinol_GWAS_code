###########################

## CAUSE creatinine and retinol

## Retinol IVW

## Dr William Reay (2023)

###########################

library(data.table)
library(dplyr)
library(cause)
library(genetics.binaRies)
library(tidyr)
library(ggplot2)
library(ieugwasr)

setwd("~/Desktop/Retinol_GWAS/")

Retinol_df <- fread("Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz",
                    sep="\t")
Retinol_df$Allele1 <- toupper(Retinol_df$Allele1)
Retinol_df$Allele2 <- toupper(Retinol_df$Allele2)

## Creatnine

Creatinine_df <- fread("~/Downloads/met-d-Creatinine.vcf.gz", skip="CHROM\tPOS",
                    stringsAsFactors=FALSE, data.table=FALSE)

NM <- colnames(Creatinine_df)

COLNAME <- NM[10]

## Split VCF fields

Creatinine_df_formatted <- separate(data = Creatinine_df, col = COLNAME, into = c("Beta", "SE", "LogP", "AF", "SS", "rsID"), sep = ":")

## Convert Beta and SE to numeric

Creatinine_df_formatted$Beta <- as.numeric(Creatinine_df_formatted$Beta)
Creatinine_df_formatted$SE <- as.numeric(Creatinine_df_formatted$SE) 

## Merge and harmonise using inbuilt CAUSE function

Merged <- gwas_merge(Creatinine_df_formatted, Retinol_df,
                     snp_name_cols = c("ID", "rsids"),
                     beta_hat_cols = c("Beta", "Effect"),
                     se_cols = c("SE", "StdErr"),
                     A1_cols = c("ALT", "Allele1"),
                     A2_cols = c("REF", "Allele2"))

set.seed(1500)

varlist <-  with(Merged, sample(snp, size=1000000, replace=FALSE))

nusiance_param <- est_cause_params(Merged, varlist)

## LD clumping

Merged_pval <- Merged %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail = F))

Clump_input <- rename(Merged_pval, "pval"="pval1", "rsid"="snp")

Clumped <- ld_clump(dplyr::tibble(rsid=Clump_input$rsid, pval=Clump_input$pval), 
                    clump_r2 =  0.01, clump_p = 1e-03, plink_bin = genetics.binaRies::get_plink_binary(),
                    bfile = "/Users/williamreay/cloudstor/Cross_disorder_PES_2019/GEUVADIS_gene_expression/g1000_eur")

## Fit CAUSE model

Mod <- cause(X=Merged, variants = Clumped$rsid, param_est = nusiance_param)

loo::pareto_k_table(Mod$loos[[2]])

Mod$elpd

summary(Mod, ci_size = 0.95, digits = 3)

PLT <- plot(Mod, type="data")

