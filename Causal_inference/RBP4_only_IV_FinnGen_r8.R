##############################

## RBP4 only IV

## FinnGen r8

## Dr William Reay (2023)

##############################

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/FinnGen_r8_MR/")

## Harmonise data to generate heterogeneity statistics

Clumped_IVs <- fread("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/Clumped_retinol_IVs.txt")

## Read each phenome-wide statistics for each IV from Finngen, concatenate, then approximate SE

N1 <- fread("10_93601207_A_C_phenotype_associations.tsv")
N1$SNP <- "rs10882283"
N1$ref <- "A"
N1$alt <- "C"


Concat <- N1

## Approximate SE

Concat$Z <- sign(Concat$beta) * sqrt(qchisq(Concat$pval, 1, lower =F))
Concat$SE <- abs(Concat$beta/Concat$Z)

## Filter with cases >1000

Concat <- Concat %>% filter(n_case > 1000)
## Format as outcome

Out_df <- format_data(Concat,
                      type = "outcome",
                      phenotype_col = "phenocode",
                      snps = Clumped_IVs$SNP,
                      snp_col = "SNP",
                      beta_col = "beta", se_col = "SE",
                      pval_col = "pval",
                      effect_allele_col = "alt",
                      other_allele_col = "ref",
                      ncase_col = "n_case", ncontrol_col = "n_control")

Harm_df <- harmonise_data(Clumped_IVs, Out_df, action =3)

MR_RBP4_finnGen <- mr(Harm_df)

MR_RBP4_finnGen$FDR <- p.adjust(MR_RBP4_finnGen$pval,method="fdr")
Out <- generate_odds_ratios(MR_RBP4_finnGen)

Out <- rename(Out, "phenocode"="outcome")

Names <- fread("10_93601207_A_C_phenotype_associations.tsv")
Names <- Names %>% select(phenocode, phenostring, category)

Merge <- merge(Out, Names, by="phenocode")

write.table(Merge, "RBP4_IV_With_names_FinnGen_MR_results.txt", sep="\t", row.names = F, quote = F)

Merge$Z <- Merge$b/Merge$se

## Correlate with the all results

IVW <- fread("FinnGen_r8_retinol_MR_Ncase_over_1000.txt")
IVW <-IVW %>% filter(method == "Inverse variance weighted (multiplicative random effects)")
IVW$Z <- IVW$b/IVW$se
IVW$phenocode <- IVW$outcome

Merge_corr <- merge(Merge, IVW, by="phenocode")

cor.test(Merge_corr$Z.x, Merge_corr$Z.y)

##  r = 0.4821001 
library(ggpubr)


CR <- ggplot(Merge_corr, aes(x=Z.x, y=Z.y)) +
  geom_point(alpha = 0.6) +
  stat_smooth(method = "lm") +
  theme_bw() +
  ylab("MR Z-score - all IVs (IVW-MRE)") +
  xlab("MR Z-score - RBP4 IV (Wald Ratio)")

CR + stat_cor(method="pearson")
