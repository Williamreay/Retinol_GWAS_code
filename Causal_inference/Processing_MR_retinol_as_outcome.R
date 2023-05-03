#############################

## Processing outcome of MR automation script 

## Retinol as outcome

## Dr William Reay (2023)

#############################

setwd("~/Desktop/Retinol_GWAS/MR_retinol_as_outcome/")

library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(TwoSampleMR)

## Read in concatenated MR data

Raw_df <- fread("Processed_output/MR_retinol_as_outcome_MR_tests.txt")

## MRE based (5 or more IVs) ##

## Extract MRE estimates and ensure no duplication

MRE_select <- Raw_df %>% filter(method == "Inverse variance weighted (multiplicative random effects)") %>%
  filter(!duplicated(id.exposure))

## Filter traits with 5 or greater IVs

MRE_filt <- MRE_select %>% filter(nsnp > 4)

## Apply FDR and FWER

MRE_filt$FWER <- p.adjust(MRE_filt$pval, method = "bonferroni")
MRE_filt$FDR <- p.adjust(MRE_filt$pval, method = "fdr")

write.table(MRE_filt, file="Multiple_testing_correction_IVW_MRE_5_or_more_IVs.txt",
            sep = "\t", row.names = F, quote = F)

## Select FDR < 0.01 traits

MRE_filt_FDR <- MRE_filt %>% filter(FDR < 0.01)

PLT_FDR <- generate_odds_ratios(MRE_filt_FDR)
PLT_FDR$Trait <- c("Triglycerides", "Creatinine (UKBB)", "Sphingomyelins", "Creatinine (Kettunen)",
                   "Free cholesterol:lipids large HDL", "Cholesteryl esters:total lipids medium VLDL",
                   "Choelsterol:total lipids medium VLDL", "Free cholesterol:total medium lipids",
                   "Phospholipids:total lipids medium VLDL", "Triglycerides in medium VLDL", "Triglycerides:total lipids medium VLDL",
                   "Phospholipids:total lipids in small VLDL", "Triglycerides:total lipids in small VLDL",
                   "Average diameter for VLDL particles", "Cholesteryl esters:total lipids very large VLDL",
                   "Cholesterol:total lipids very large VLDL", "Free cholesterol:total lipids very large VLDL" ,"Triglycerides in very large VLDL",
                   "Triglycerides:total lipids very large VLDL", "T2* swMRI left putamen",
                   "Frequency of solarium/sun lamp use")
PLT_FDR$Category <- c("Lipids", "Metabolite", "Lipids", "Metabolite",
                   "Lipids", "Lipids", "Lipids", 
                   "Lipids", "Lipids", 
                   "Lipids", "Lipids", "Lipids", 
                   "Lipids", "Lipids", 
                   "Lipids", "Lipids", 
                   "Lipids", "Lipids", 
                   "Lipids", "Neuroimaging",
                   "Behavioural")

ggplot(data = PLT_FDR, aes(x = Trait, y = b, ymin = lo_ci, ymax = up_ci, colour = Category)) +
  geom_pointrange() +
  theme_bw() +
  theme(legend.position = "right") +
  coord_flip() +
  xlab(" ") +
  ylab("Effect on retinol (95% CI) per SD increase in exposure") +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_color_brewer(palette = "Paired") +
  ggtitle("IVW-MRE, FDR < 0.01")

## Plot FDR < 0.01 traits

## Tier system - as in exposure ##

List_of_traits <- MRE_filt_FDR %>% select(id.outcome)

Selected_traits <- merge(Raw_df, List_of_traits, by = "id.outcome")

Nom_selected <- Selected_traits %>% filter(pval < 0.05) %>% select(id.outcome)

## Tier 1
Tier_1 <- Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 4)

Tier_1_traits <- merge(Tier_1, Selected_traits, by = "id.outcome")

## Tier 2
Tier_2 <- Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 3)

Tier_2_traits <- merge(Tier_2, Selected_traits, by = "id.outcome")

## Tier 3
Tier_3 <- Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 2)

Tier_3_traits <- merge(Tier_3, Selected_traits, by = "id.outcome")

## Check the following -
## i) No heterogeneity between IVs, ii) Egger intercept non-significant, iii) no leave-one-out estimate non-significant

## Heterogeneity - import and retain only traits without significant heterogeneity
Het_raw <- fread("Processed_output/MR_retinol_as_outcome_HET_tests.txt", sep = "\t")

Het_filt <- Het_raw %>% filter(Q_pval > 0.05) %>% select(id.outcome)

## Egger intercept - import and retain only traits without significant heterogeneity
Egger_raw <- fread("Processed_output/MR_retinol_as_outcome_EI_tests.txt", sep= "\t")

Egger_filt <- Egger_raw %>% filter(pval > 0.05) %>% select(id.outcome)

## LOO - import and retain only traits without any evidence of a LOO estimate being > 0.05
LOO_raw <- fread("Processed_output/MR_retinol_as_outcome_LOO_tests.txt", sep= "\t")

LOO_filt <- LOO_raw %>% filter(p < 0.05) %>% filter(SNP != "All") %>% select(id.outcome)

LOO_count <- LOO_filt %>% group_by_all() %>% summarise(LOO_COUNT = n()) %>% arrange(desc(LOO_COUNT))  %>% filter(LOO_COUNT > 5)

## See how many Tier 1 traits pass the above
List_cond <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "id.outcome"),
                      list(Het_filt, Egger_filt, LOO_count))

List_cond <- List_cond %>% filter(!duplicated(id.outcome))

Tier_1_Merge_pass_het_and_Egger_test <- merge(Tier_1_traits, List_cond, by = "id.outcome")

write.table(Tier_1_Merge_pass_het_and_Egger_test, file="Passing_traits_tier_1.txt",
            sep = "\t", row.names = F, quote = F)

## Tier 2

Tier_2_Merge_pass_het_and_Egger_test <- merge(Tier_2_traits, List_cond, by = "id.outcome")

## Just creatinine again

## Tier 3

Tier_3_Merge_pass_het_and_Egger_test <- merge(Tier_3_traits, List_cond, by = "id.outcome")

## Just creatinine again

library(TwoSampleMR)
PLT <- generate_odds_ratios(Tier_1_Merge_pass_het_and_Egger_test)
PLT$METHOD <- c("IVW-FE", "Weighted Median", "MR Egger", "Weighted Mode", "IVW-MRE")

CR_FP <- ggplot(data = PLT, aes(x = METHOD, y = b, ymin = lo_ci, ymax = up_ci, colour = METHOD)) +
  geom_pointrange() +
  theme_bw() +
  theme(legend.position = "none") +
  coord_flip() +
  xlab(" ") +
  ylab("Effect on retinol (95% CI) per SD increase in serum creatinine") +
  ylim(-0.2, 1) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_color_brewer(palette = "Paired") +
  ggtitle("Creatinine to retinol")

## Leave one out density plot

Creat_LOO <- LOO_raw %>% filter(id.exposure == "met-d-Creatinine")

Creat_LOO <- Creat_LOO %>% filter(SNP != "All")

Creat_LOO$Z <- Creat_LOO$b/Creat_LOO$se

Creat_LOO$Num <- c(1:61)

CP <- ggplot(Creat_LOO, aes(x=Z)) +
  geom_density(fill = "orange", alpha = 0.4) +
  xlim(1, 7) +
  theme_bw() +
  geom_vline(xintercept = 1.96, lty = "dashed") +
  xlab("MR IVW Z score - leave-one-out analysis") +
  ylab("Density") +
  ggtitle("Distribution of leave-one-out results")

ggarrange(CR_FP, CP)

## Proceed to traits with 4 or less IVs

FE_filt <- Raw_df %>% filter(method == "Inverse variance weighted (fixed effects)") %>%
  filter(!duplicated(id.exposure)) %>% filter(nsnp < 5)

FE_filt$FDR <- p.adjust(FE_filt$pval, method = "fdr")
FE_filt$FWER <- p.adjust(FE_filt$pval, method = "bonferroni")

FDR_FE <- FE_filt %>% filter(FDR < 0.01)

write.table(FDR_FE, file="FDR_sig_Fixed_effects_less_than_5_IVs_Ret_as_outcome.txt",
            sep = "\t", row.names = F)

## Check single SNP

Alanine_SSNP <- fread("Raw_output/met-a-469_IVW_METSIM_INTERVAL_retinol_SSNP.txt")

## Driven by effect of GCKR

Alanine_DT <- fread("Raw_output/met-a-469_IVW_METSIM_INTERVAL_retinol_DT.txt")

## Look for bidirectional effects