#############################

## Bidirectional retinol MR effects

## Dr William Reay (2023)

#############################

library(data.table)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Retinol_GWAS/")

## Read in ret as exposure (IVW-MRE)

Ret_exp <- fread("MR_all_cis_acting_SNPs/All_MRE_IVW_results_ret_as_exp.txt")
Ret_exp$ID_merge <- Ret_exp$ID
## Read in ret as outcome

Ret_outcome <- fread("MR_retinol_as_outcome/Multiple_testing_correction_IVW_MRE_5_or_more_IVs.txt")
Ret_outcome$ID_merge <- Ret_outcome$id.exposure

## Identify anything shared (P < 0.05) for both

Merge_ret <- merge(Ret_exp, Ret_outcome, by ="ID_merge")


Filt <- Merge_ret %>% filter(pval.x < 0.05 & pval.y < 0.05)
Filt$Z.1 <- Filt$b.x/Filt$se.x
Filt$Z.2 <- Filt$b.y/Filt$se.y

Filt$Name <- c("Creatinine (UKBB)", "Age hypertension diagnosed",
               "Age of menarche", "Body fat %",
               "Haematocrit %")

EXP <- ggplot(Filt, aes(x=Name, y=Z.1)) +
  geom_bar(stat = "identity", position=position_dodge(), colour = "black", fill = "steelblue1") +
  coord_flip() +
  theme_bw() +
  ylab("MR Z score - IVW-MRE") +
  ggtitle("Retinol as exposure") +
  geom_hline(yintercept = 0, lty="dashed") +
  ylim(-7, 7)

OUT <- ggplot(Filt, aes(x=Name, y=Z.2)) +
  geom_bar(stat = "identity", position=position_dodge(), colour = "black", fill = "grey72") +
  coord_flip() +
  theme_bw() +
  ylab("MR Z score - IVW-MRE") +
  ggtitle("Retinol as outcome") +
  geom_hline(yintercept = 0, lty="dashed") +
  ylim(-7, 7)

ggarrange(EXP, OUT, nrow = 2)
