#############################

## Processing outcome of MR automation script 

## Retinol as exposure

## Dr William Reay (2023)

#############################

setwd("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/")

library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(psych)
library(ggpubr)

## Read in raw concatenated data

Raw_df <- fread("Concatenated_MR_results_RAW_retinol_as_exposure.txt", sep = "\t")

## Select outcomes that only have at least 6 IVs of the 7

Step_1 <- Raw_df %>% filter(nsnp > 5)

## Extract multiplicative random effects estimates for multiple testing correction

MRE_select <- Step_1 %>% filter(method == "Inverse variance weighted (multiplicative random effects)")

## Ensure no duplicates

MRE_select <- MRE_select %>% filter(!duplicated(id.outcome))

## Apply FDR and FWER correction

MRE_select$FWER <- p.adjust(MRE_select$pval, method = "bonferroni")
MRE_select$FDR <- p.adjust(MRE_select$pval, method = "fdr")

write.table(MRE_select, file="All_MRE_IVW_results_ret_as_exp.txt",
            sep = "\t", row.names = F, quote = F)
## Select FDR < 0.01 sig traits

FDR_sig_MRE <- MRE_select %>% filter(FDR < 0.01)

write.table(FDR_sig_MRE, file = "Processed_output/FDR_sig_0.01_IVW_MRE_retinol_as_exp_6_IVs_or_more.txt",
            sep = "\t", row.names = F, quote = F)

## Look at consistency across different MR methods for FDR sig traits

List_of_traits <- FDR_sig_MRE %>% select(id.outcome)

Selected_traits <- merge(Step_1, List_of_traits, by = "id.outcome")

## Tier #1 - all methods at least nominally significant - none satisfy this

## Tier #2 - four out of five methods at least nominally significant

Nom_selected <- Selected_traits %>% filter(pval < 0.05) %>% select(id.outcome)

Tier_2 <- Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 3)

Tier_2_traits <- merge(Tier_2, Selected_traits, by = "id.outcome")

## Tier 3 - three out of five methods at least nominally significant

Tier_3 <- Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 2)

Tier_3_traits <- merge(Tier_3, Selected_traits, by = "id.outcome")

## Check the following -
## i) No heterogeneity between IVs, ii) Egger intercept non-significant, iii) no leave-one-out estimate non-significant

## Heterogeneity - import and retain only traits without significant heterogeneity
Het_raw <- fread("Heterogeneity_concatenated_ret_as_exp.txt", sep = "\t")

Het_filt <- Het_raw %>% filter(Q_pval > 0.05) %>% select(id.outcome)

## Egger intercept - import and retain only traits without significant heterogeneity
Egger_raw <- fread("Egger_intercept_concatenated_ret_as_exp.txt", sep= "\t")

Egger_filt <- Egger_raw %>% filter(pval > 0.05) %>% select(id.outcome)

## LOO - import and retain only traits without any evidence of a LOO estimate being > 0.05
LOO_raw <- fread("LOO_concatenated_ret_as_exp.txt", sep= "\t")

LOO_filt <- LOO_raw %>% filter(p < 0.05) %>% filter(SNP != "All") %>% select(id.outcome)

LOO_count <- LOO_filt %>% group_by_all() %>% summarise(LOO_COUNT = n()) %>% arrange(desc(LOO_COUNT))  %>% filter(LOO_COUNT > 5)

## See how many Tier 2 traits pass the above
List_tier_2 <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "id.outcome"),
                                               list(Het_filt, Egger_filt, LOO_count))

List_tier_2 <- List_tier_2 %>% filter(!duplicated(id.outcome))

Tier_2_Merge_pass_het_and_Egger_test <- merge(Tier_2_traits, List_tier_2, by = "id.outcome")

## Ensure number of significant LOO matches number of IVs
Tier_2_Merge_pass_het_and_Egger_test <- Tier_2_Merge_pass_het_and_Egger_test %>% filter(LOO_COUNT == nsnp)

write.table(Tier_2_Merge_pass_het_and_Egger_test, file="Processed_output/Tier_2/Passing_traits_tier_2.txt",
            sep = "\t", row.names = F, quote = F)

## Plot the Tier 2 results

Plot_input_tier_2 <- Tier_2_Merge_pass_het_and_Egger_test %>% filter(method == "Inverse variance weighted (multiplicative random effects)") 

Plot_input_tier_2$Plot_name <- c("IgG seropositivity (Herpesvirus 6 IE1A antigen)", "Ruminococcaceae prevalence (OTU97_105)", 
                                 "Ruminococcaceae prevalence (OTU97_120)", "FAM177A1 protein expression",
                                 "Right rostral anterior cingulate thickness", "Trunk fat percentage",
                                 "Age hay fever, rhinitis or eczema diagnosed")

Plot_input_tier_2$Category <- c("Immune", "Microbiome", "Microbiome", "Protein", "Neuroimaging", "Anthropometric", "Immune")
Plot_input_tier_2$L_CI <- Plot_input_tier_2$b - (1.96*Plot_input_tier_2$se)
Plot_input_tier_2$U_CI <- Plot_input_tier_2$b + (1.96*Plot_input_tier_2$se)
Plot_input_tier_2$Z <- Plot_input_tier_2$b/Plot_input_tier_2$se

BP1 <- ggplot(Plot_input_tier_2, aes(x=Z, y=Plot_name, fill=Category)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  theme_bw() +
  xlab("MR Z score - standardised effect of retinol (IVW-MRE)") +
  ylab(" ") +
  geom_vline(xintercept = 0, lty="dashed") +
  scale_fill_manual(values = c("#A6CEE3", "#E31A1C", "#FF7F00", 
                          "#CAB2D6", "#FFFF99")) +
  theme(legend.position = "none")

## Correlate with individual RBP4 Z scores

RBP4_effects_single_IV <- fread("../RBP4_MR_protein_retinol/Retinol_RBP_IV_IEUGWAS_UKBB_MR_PheWAS.txt")
RBP4_effects_single_IV$Z <- (RBP4_effects_single_IV$MR_beta/RBP4_effects_single_IV$MR_SE)
RBP4_effects_single_IV <- RBP4_effects_single_IV %>% select(id, Z)
RBP4_effects_single_IV$id.outcome <- RBP4_effects_single_IV$id

MRE_select$Z <- (MRE_select$b/MRE_select$se)
Z_MRE <- MRE_select %>% select(id.outcome, Z)

Merge_RBP4 <- merge(RBP4_effects_single_IV, Z_MRE,  by = "id.outcome")
Merge_RBP4 <- Merge_RBP4 %>% filter(!duplicated(id.outcome))

cor.test(Merge_RBP4$Z.x, Merge_RBP4$Z.y)

CP <- ggplot(Merge_RBP4, aes(x=Z.x, y=Z.y)) +
  geom_point(alpha = 0.6) +
  stat_smooth(method = "lm") +
  theme_bw() +
  ylab("MR Z-score - all IVs (IVW-MRE)") +
  xlab("MR Z-score - RBP4 IV (Wald Ratio)")

CP_out <- CP + stat_cor(method = "pearson")

Merge_plot_comp <- merge(RBP4_effects_single_IV, Plot_input_tier_2, by = "id.outcome")
Merge_plot_comp$MR_Z <- Merge_plot_comp$b/Merge_plot_comp$se

Tier_2_com <- ggplot(Merge_plot_comp, aes(x=Z.y, y=Z.x, colour=Category)) +
  geom_point() +
  theme_bw() +
  xlab("MR Z-score - all IVs (IVW-MRE)") +
  ylab("MR Z-score - RBP4 IV (Wald Ratio)") +
  xlim(-7,7) +
  ylim(-7,7) +
  geom_hline(yintercept = 1.96, lty = "dashed") +
  geom_hline(yintercept = -1.96, lty = "dashed") +
  geom_vline(xintercept = 1.96, lty = "dashed") +
  geom_vline(xintercept = -1.96, lty = "dashed") +
  scale_colour_manual(values = c("#A6CEE3", "#E31A1C", "#FF7F00", 
                               "#CAB2D6", "#FFFF99"))

ggarrange(BP1, Tier_2_com)

## Repeat for Tier 3 traits

Tier_3_Merge_pass_het_and_Egger_test <- merge(Tier_3_traits, List_tier_2, by = "id.outcome")

## Ensure number of significant LOO matches number of IVs
Tier_3_Merge_pass_het_and_Egger_test <- Tier_3_Merge_pass_het_and_Egger_test %>% filter(LOO_COUNT == nsnp)

write.table(Tier_3_Merge_pass_het_and_Egger_test, file="Processed_output/Tier_3/Passing_traits_tier_3.txt",
            sep = "\t", row.names = F, quote = F)

## Plot the top 10 most significant tier 3 results

Tier_3_MRE <- Tier_3_Merge_pass_het_and_Egger_test %>% filter(method == "Inverse variance weighted (multiplicative random effects)")

Tier_3_MRE$Plot_name <- c("IgG seropositivity (Herpesvirus 6 IE1A antigen)", "G. Oscillibacter Abundance",
                          "Bacteroides abundance (OTU99_197)", "Ruminococcaceae prevalence (OTU97_105)", 
                          "Ruminococcaceae prevalence (OTU97_120)", "Ruminococcaceae prevalence (OTU99_121)",
                          "Faecalibacterium prevalence (OTU99_45)", "Sequalae of tuberculosis",
                          "Rash/non-specific skin erruption", "Sitting height ratio", "FAM177A1 protein expression",
                          "PEAR-1 protein expression", "ADPGK protein expression", "rfMRI connectivity ICA25 edge 0161",
                          "rfMRI connectivity ICA100 edge 13", "Left pallidum volume", "rfMRI connectivity ICA100 edge 0306",
                          "Interm prim-Jensen syrus surface area", "Left cuneus thickness", "Right caudal anterior cingulate thickness",
                          "Right rostral anterior cingulate thickness", "Right percular inferior frontal gyrus thickness",
                          "Right pericallosal thicnkess", "Body fat percentage", "Left leg fat percentage", "Trunk fat percentage", "Dental problems (none of above)",
                          "Keratometry:3mm weak meridian (left)", "Keratometry:3mm weak meridian (right)", "Keratometry:6mm weak meridian (left)",
                          "Keratometry:6mm weak meridian (left)", "Other savoury snack intake", "Other bread intake", "Home area pop. density - hamlet/isolated dwelling",
                          "Medication related adverse effects", "Coxarthrosis [arthrosis of hip](FG)", "Medication related adverse effects (Asthma/COPD)",
                          "Cereal Type", "Age hay fever/rhinitis/eczema diagnosed")

Tier_3_MRE$Category <-  c("Immune", "Microbiome", "Microbiome", "Microbiome", "Microbiome", "Microbiome", "Microbiome",
                          "Immune", "Dermatological", "Anthropometric", "Protein", "Protein", "Protein", "Neuroimaging",
                          "Neuroimaging","Neuroimaging","Neuroimaging","Neuroimaging","Neuroimaging","Neuroimaging","Neuroimaging",
                          "Neuroimaging","Neuroimaging", "Anthropometric", "Anthropometric","Anthropometric", "Dental",
                          "Opthamological","Opthamological", "Opthamological", "Opthamological", "Food intake", "Food intake",
                          "Demographic", "Medication", "Rheumatological", "Medication", "Food intake", "Immune")

## Plot

Tier_3_MRE$L_CI <- Tier_3_MRE$b - (1.96*Tier_3_MRE$se)
Tier_3_MRE$U_CI <- Tier_3_MRE$b + (1.96*Tier_3_MRE$se)
Tier_3_MRE$Z <- Tier_3_MRE$b/Tier_3_MRE$se

BP_tier3 <- ggplot(Tier_3_MRE, aes(x=Z, y=Plot_name, fill=Category)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  theme_bw() +
  xlab("MR Z score - standardised effect of retinol (IVW-MRE)") +
  ylab(" ") +
  geom_vline(xintercept = 0, lty="dashed") +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "none")

## Comparison plot

Merge_plot_com2p <- merge(RBP4_effects_single_IV, Tier_3_MRE, by = "id.outcome")

Com <- ggplot(Merge_plot_com2p, aes(x=Z.y, y=Z.x, colour=Category)) +
  geom_point() +
  theme_bw() +
  xlab("MR Z-score - all IVs (IVW-MRE)") +
  ylab("MR Z-score - RBP4 IV (Wald Ratio)") +
  xlim(-7,7) +
  ylim(-7,7) +
  geom_hline(yintercept = 1.96, lty = "dashed") +
  geom_hline(yintercept = -1.96, lty = "dashed") +
  geom_vline(xintercept = 1.96, lty = "dashed") +
  geom_vline(xintercept = -1.96, lty = "dashed") +
  scale_color_brewer(palette = "Paired")

ggarrange(BP_tier3, Com)

