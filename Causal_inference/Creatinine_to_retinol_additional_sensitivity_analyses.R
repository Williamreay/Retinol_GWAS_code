################################

## Creatinine to retinol

## Additional sensitivity analyses

## Dr William Reay (2023)

################################

library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
library(data.table)
library(ggplot2)
library(MendelianRandomization)
library(MVMR)

## Penalised and robust estimates ##

## Read in creatinine exposures

Creat_raw <- extract_instruments("met-d-Creatinine", clump = T)

## Read in retinol

Out_df <- fread("~/Desktop/Retinol_GWAS/Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

Out_df$Pheno <- "Retinol"
Out_df$N <- 17268
Out_df$Allele1 <- toupper(Out_df$Allele1)
Out_df$Allele2 <- toupper(Out_df$Allele2)

Ret_df <- format_data(Out_df,
                      type = "outcome",
                      snps = Creat_raw$SNP,
                      beta_col = "Effect",
                      se_col = "StdErr",
                      snp_col = "rsids",
                      effect_allele_col = "Allele1",
                      other_allele_col = "Allele2",
                      phenotype_col = "Pheno",
                      samplesize_col = "N")

Harmonised <- harmonise_data(Creat_raw, Ret_df, action = 3)

## Robust and penalised tests

Creat_to_ret_input <- mr_input(bx = Harmonised$beta.exposure,
                            bxse = Harmonised$se.exposure,
                            by = Harmonised$beta.outcome,
                            byse = Harmonised$se.outcome,
                            snps = Harmonised$SNP,
                            effect_allele = Harmonised$effect_allele.exposure,
                            other_allele = Harmonised$other_allele.exposure,
                            exposure = "Creatinine (SD)",
                            outcome = "Retinol (SD)")

Creat_all <- mr_allmethods(Creat_to_ret_input, method = "all")

mr_plot(Creat_all)

## Multivariable MR - Creatinine with lipids (HDL, LDL, Triglycerides)

## Define multivariable IDs

mv_ids <- c("met-d-Creatinine", "ieu-a-299", "ieu-a-300", "ieu-a-302")

## Mutlvariable IVs - strict threshold for all

MV_IVs_STRICT <- mv_extract_exposures(mv_ids, clump_r2 = 0.001, pval_threshold = 5e-08)

MV_Ret_df <- format_data(Out_df,
                      type = "outcome",
                      snps = MV_IVs_STRICT$SNP,
                      beta_col = "Effect",
                      se_col = "StdErr",
                      snp_col = "rsids",
                      effect_allele_col = "Allele1",
                      other_allele_col = "Allele2",
                      phenotype_col = "Pheno",
                      samplesize_col = "N")

MV_harmonised_STRICT <- mv_harmonise_data(MV_IVs_STRICT, MV_Ret_df, harmonise_strictness = 3)

## Convert to MendelianRandomisation package format

MRMV_STRICT <- mr_mvinput(bx = MV_harmonised_STRICT$exposure_beta,
                                  bxse = MV_harmonised_STRICT$exposure_se,
                                  by = MV_harmonised_STRICT$outcome_beta,
                                  byse = MV_harmonised_STRICT$outcome_se,
                                  exposure = c("HDL","LDL", "Triglycerides", "Creatinine"),
                                  outcome="Retinol")


## Assess instrument strength, set pairwise covariance to zero, even though lipids from same sample so may be underestimate https://doi.org/10.1093/ije/dyy262

STRICT_MV_IV_strength <- strength_mvmr(MRMV_STRICT, gencov = 0)

## Conditional F-statistics for instrument strength

#exposure1 exposure2 exposure3 exposure4
#F-statistic  44.98258  39.42955  19.09354  9.733294

## Perform mutlivariable methods - IVW, median, Egger, and LASSO

mvIVW <- mr_mvivw(MRMV_STRICT)

mvLASSO <- mr_mvlasso(MRMV_STRICT)

mvEgger <- mr_mvegger(MRMV_STRICT)

mvMedian <- mr_mvmedian(MRMV_STRICT)

## Format output and plot forest plot

Output_IVW <- as.data.frame(cbind(mvIVW@Exposure, mvIVW@Estimate, mvIVW@StdError))
Output_IVW$V4 <- "mvIVW"

Output_LASSO <- as.data.frame(cbind(mvLASSO@Exposure, mvLASSO@Estimate, mvLASSO@StdError))
Output_LASSO$V4 <- "mvLASSO"

Output_Egger <- as.data.frame(cbind(mvEgger@Exposure, mvEgger@Estimate, mvEgger@StdError.Est))
Output_Egger$V4 <- "mvEgger"

Output_Median <- as.data.frame(cbind(mvMedian@Exposure, mvMedian@Estimate, mvMedian@StdError))
Output_Median$V4 <- "mvMedian"

Combined_df <- rbind(Output_IVW, Output_LASSO, Output_Median, Output_Egger)
Combined_df <- rename(Combined_df, "Exposure"="V1", "Beta"="V2", "SE"="V3", "MVMR_method"="V4")

write.table(Combined_df, file = "~/Desktop/Retinol_GWAS/MR_retinol_as_outcome/Creatinine_HDL_LDL_TG_MVMR.txt",
            sep = "\t", row.names = F, quote = F)

## FP
Combined_df$Beta <- as.numeric(Combined_df$Beta)
Combined_df$SE <- as.numeric(Combined_df$SE)
Combined_df$UCI <- Combined_df$Beta + (1.96*Combined_df$SE)
Combined_df$LCI <- Combined_df$Beta - (1.96*Combined_df$SE)

ggplot(data = Combined_df, aes(x=Exposure, y=Beta, ymin=LCI, ymax=UCI, colour=Exposure)) +
  geom_pointrange() +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~MVMR_method,strip.position="top",nrow=5,scales = "free_y") +
  coord_flip() +
  ylab("Effect on retinol [95% CI] per SD in exposure increase") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank()) +
  scale_colour_brewer(palette = "Dark2") +
  ylim(-0.3,0.6)
 