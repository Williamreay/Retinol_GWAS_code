#################################

## RBP4 retinol IVs - MR-PheWAS: IEUGWAS: v6.9.4 - 2023-09-01

## William Reay (2022)

#################################

library(dplyr)
library(data.table)
library(ggplot2)
library(psych)
library(ggpubr)
library(ieugwasr)
library(coloc)
library(viridis)

setwd("~/Desktop/Retinol_GWAS/RBP4_MR_protein_retinol/")

## Effect size of lead SNP on RBP4 protein abundance and retinol

## Retinol: rs10882283, EA = A , NEA = C, Beta = 0.1086 SD

## RBP4 protein: rs10882283, EA = A, NEA = C, Beta = 0.18099 SD


## Use adapted previous function developed for this purpose - https://github.com/Williamreay/Pneumonia_meta_GWAS/blob/master/MR/MR_pheWAS_IEUGWAS.R
## Do not include multiple-esting corr step to check for any spurious estimates

IEU_mr_pheWAS <- function(SNP, EA, NEA, IV_beta) {
  
  ## Extract IV-outcome effect sizes from the IEUGWAS db
  
  SNP_pheWAS <- phewas(SNP, pval = 1)
  
  ## Check for mismatching effect and non-effect alleles
  
  mismatch_SNP <- which(SNP_pheWAS$ea != EA & SNP_pheWAS$nea != NEA,arr.ind=TRUE)
  SNP_pheWAS[mismatch_SNP,]$beta <- SNP_pheWAS[mismatch_SNP,]$beta*-1
  SNP_pheWAS[mismatch_SNP,]$ea <- EA
  SNP_pheWAS[mismatch_SNP,]$nea <- NEA
  
  ## Perform MR
  
  SNP_pheWAS$MR_beta <- SNP_pheWAS$beta/IV_beta
  SNP_pheWAS$MR_SE <- abs(SNP_pheWAS$se/IV_beta)
  SNP_pheWAS$MR_pval <- pnorm(abs(SNP_pheWAS$MR_beta)/SNP_pheWAS$MR_SE, lower.tail=FALSE) * 2

  return(SNP_pheWAS)
}

## Retinol pheWAS

Retinol_MR <- IEU_mr_pheWAS("rs10882283", "A", "C", 0.1086)

## Exclude four IEUGWAS entries with incorrect formatted beta/SE estimates

Retinol_MR <- Retinol_MR %>% filter(id != "ebi-a-GCST007236" & id != "ebi-a-GCST007800" & id != "ebi-a-GCST007799" &
                                      id != "ieu-b-5070")

Retinol_MR$MR_FWER <- p.adjust(Retinol_MR$MR_pval, method="bonferroni")
Retinol_MR$MR_FDR <- p.adjust(Retinol_MR$MR_pval, method="fdr")
Retinol_MR$Exposure <- "Retinol"

write.table(Retinol_MR, file="Retinol_RBP_IV_IEUGWAS_UKBB_MR_PheWAS.txt", quote = F, row.names = F, sep = "\t")

## RBP4 pheWAS

RBP4_MR <- IEU_mr_pheWAS("rs10882283", "A", "C", 0.18099)

## Exclude four IEUGWAS entries with incorrect formatted beta/SE estimates

RBP4_MR <- RBP4_MR %>% filter(id != "ebi-a-GCST007236" & id != "ebi-a-GCST007800" & id != "ebi-a-GCST007799" &
                                      id != "ieu-b-5070")

RBP4_MR$MR_FWER <- p.adjust(RBP4_MR$MR_pval, method="bonferroni")
RBP4_MR$MR_FDR <- p.adjust(RBP4_MR$MR_pval, method="fdr")
RBP4_MR$Exposure <- "RBP4_protein_expression"

write.table(RBP4_MR, file="RBP4_pQTL_RBP_IV_IEUGWAS_UKBB_MR_PheWAS.txt", quote = F, row.names = F, sep = "\t")

## Colocalisation for select number of phenotypes that are at least MR-phenome-wide significant (P < 1e-05)

## Define function to do the following:
## 1. Retrieve SNP effect sizes from IEUGWASdb
## 2. Merge with retinol GWAS
## 3. Perform colocalisation

Ret_GWAS <- fread("../Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

## Hg38 locus for RBP4 signal

Ret_GWAS <- Ret_GWAS %>% filter(chrom == 10) %>% filter(pos > 93536118 & pos < 93601208)
Ret_GWAS <- rename(Ret_GWAS, "rsid"="rsids")
Ret_GWAS$SNPVar <- (Ret_GWAS$StdErr)^2

Coloc_ieugwas_func <- function(ieugwas_id) {
  
  ## Pull region in hg19 from IEUGWAS
  
  Trait2 <- associations(variants="10:95295876-95360965", id=c(ieugwas_id), proxies=0)
  Trait2$SNPVar <- (Trait2$se)^2
  
  ## Merge by SNP ID with retinol GWAS
  Merged_sumstats <- merge(Ret_GWAS, Trait2, by = "rsid")
  
  ## Perform remaning colocalisation steps
  Retinol_df <- Merged_sumstats %>% select(rsid, Effect, SNPVar.x)
  Trait2_df <- Merged_sumstats %>% select(rsid, beta, SNPVar.y)
  
  ## Define coloc list input
  Retinol=list(beta=Retinol_df$Effect, varbeta=Retinol_df$SNPVar.x, 
               sdY=1, type="quant")
  Trait2_df=list(beta=Trait2_df$beta, varbeta=Trait2_df$SNPVar.y, 
               sdY=1, type="quant")
  
  # Perform colocalisation under default priors
  coloc_df <- coloc.abf(Retinol, Trait2_df, MAF=NULL)
  
  return(coloc_df)
  
}

## White Blood Cell Count (UKBB) - ukb-d-30000_irnt, units = SD

WBC <- Coloc_ieugwas_func("ukb-d-30000_irnt")
sensitivity(WBC, rule="H4 > 0.8", plot.manhattans = T)
write.table(WBC$summary, file = "Colocalisation_RBP4_region/WBC_coloc_results.txt",
            sep = "\t", row.names = T, quote = F)

## Triglycerides (UKBB) - ieu-b-111, units = SD

TG <- Coloc_ieugwas_func("ieu-b-111")
sensitivity(TG, rule="H3 > 0.8", plot.manhattans = T)
write.table(TG$summary, file = "Colocalisation_RBP4_region/Triglycerides_coloc_results.txt",
            sep = "\t", row.names = T, quote = F)

## Optic disc area - ebi-a-GCST004076, units = SD

OD <- Coloc_ieugwas_func("ebi-a-GCST004076")
sensitivity(OD, rule="H4 > 0.8", plot.manhattans = T)
write.table(OD$summary, file = "Colocalisation_RBP4_region/Optic_disc_area_coloc_results.txt",
            sep = "\t", row.names = T, quote = F)

## NET100 0364 (UKBB) - esting-state fMRI network edge, ubm-a-1520, units = SD

fMRI <- Coloc_ieugwas_func("ubm-a-1520")
sensitivity(fMRI, rule="H4 > 0.8", plot.manhattans = T)
write.table(fMRI$summary, file = "Colocalisation_RBP4_region/NET100_0364_fMRI_network_edge_coloc_results.txt",
            sep = "\t", row.names = T, quote = F)

## Make forest plot for the above and indicate whether colocalisation occured

FP_input <- Retinol_MR %>% filter(id == "ubm-a-1520" | id == "ebi-a-GCST004076" |
                                    id == "ieu-b-111" | id == "ukb-d-30000_irnt")

FP_input$Phenotype <- c("Triglycerides", "Leukocyte Count", "Optic Disc Area", "Resting state fMRI network edge 0364")
FP_input$Colocalised <- c("No", "Yes", "Yes", "Yes")
FP_input$L_CI <- FP_input$MR_beta - (1.96*FP_input$MR_SE)
FP_input$U_CI <- FP_input$MR_beta + (1.96*FP_input$MR_SE)

FP <- ggplot(data = FP_input, aes(x=Phenotype, y=MR_beta,
                            ymax=U_CI, ymin=L_CI, 
                            colour=Colocalised)) +
  coord_flip() +
  theme_bw() +
  geom_pointrange() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Effect on outcome (SD) per SD increase in serum retinol") +
  xlab("Outcome")

FP2 <- FP + scale_color_manual(values = c("grey30", "orange"))  


