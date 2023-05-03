#############################
## Automated MR of all retinol IVs

## FinnGen r8 data - 1141 binary phenotypes with N cases > 1000

## William Reay (2023)
#############################

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/FinnGen_r8_MR/")

Clumped_IVs <- fread("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/Clumped_retinol_IVs.txt")

## Read each phenome-wide statistics for each IV from Finngen, concatenate, then approximate SE

N1 <- fread("10_93601207_A_C_phenotype_associations.tsv")
N1$SNP <- "rs10882283"
N1$ref <- "A"
N1$alt <- "C"
N2 <- fread("16_79684284_C_G_phenotype_associations.tsv")
N2$SNP <- "rs12149203"
N2$ref <- "C"
N2$alt <- 'G'
N3 <- fread("18_31593832_A_G_phenotype_associations.tsv")
N3$SNP <- "rs1080094"
N3$ref <- "A"
N3$alt <- "G"
N4 <- fread("2_121326709_G_A_phenotype_associations.tsv")
N4$SNP <- "rs34898035"
N4$ref <- "G"
N4$alt <- "A"
N5 <- fread("2_27508073_T_C_phenotype_associations.tsv")
N5$SNP <- "rs1260326"
N5$ref <- "T"
N5$alt <- "C"
N6 <- fread("20_40572274_A_G_phenotype_associations.tsv")
N6$SNP <- "rs6029188"
N6$ref <- "A"
N6$alt <- "G"
N7 <- fread("7_114540237_C_A_phenotype_associations.tsv")
N7$SNP <- "rs11762406"
N7$ref <- "C"
N7$alt <- "A"
N8 <- fread("8_9327181_T_C_phenotype_associations.tsv")
N8$SNP <- "rs6601299"
N8$ref <- "T"
N8$alt <- "C"

Concat <- rbind(N1, N2, N3, N4, N5, N6, N7, N8)

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

## Harmonise

Harmonise_df <- harmonise_data(Clumped_IVs, Out_df, action = 3)

## MR

MR <- mr(Harmonise_df, method_list=c("mr_ivw_mre", "mr_ivw_fe", "mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

write.table(MR, file="FinnGen_r8_retinol_MR_Ncase_over_1000.txt",
            sep = "\t", row.names = F, quote = F)

## IVW-MRE FDR < 0.01

FDR_sig <- MR %>% filter(method == "Inverse variance weighted (multiplicative random effects)")
FDR_sig$FDR <- p.adjust(FDR_sig$pval, method="fdr")

N1$outcome <- N1$phenocode

FDR_sig_name <- merge(N1, FDR_sig, by = "outcome")

FDR_sig <- FDR_sig %>% filter(FDR < 0.01)

## Look at consistency across different MR methods for FDR sig traits

List_of_traits <- FDR_sig %>% select(outcome)

Selected_traits <- merge(MR, List_of_traits, by = "outcome")

## Tier system based on significance across multiple methods

Nom_selected <- Selected_traits %>% filter(pval < 0.05) %>% select(outcome)

## Tier 1 - no tier 1
Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 4)

## Tier 2  - no tier 2

Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 3)

## Tier 3 - 5 tier 3

Tier_3 <- Nom_selected %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))  %>% filter(COUNT > 2)

Tier_3_traits <- merge(Tier_3, Selected_traits, by = "outcome")

## Heterogeneity

Het <- mr_heterogeneity(Harmonise_df)

Het_filt <- Het %>% filter(Q_pval > 0.05) %>% select(outcome)

## Egger

Egger_raw <- mr_pleiotropy_test(Harmonise_df)

Egger_filt <- Egger_raw %>% filter(pval > 0.05) %>% select(outcome)

## Leave-one-out

LOO_raw <- mr_leaveoneout(Harmonise_df)

LOO_filt <- LOO_raw %>% filter(p < 0.05) %>% filter(SNP != "All") %>% select(outcome)

LOO_count <- LOO_filt %>% group_by_all() %>% summarise(LOO_COUNT = n()) %>% arrange(desc(LOO_COUNT))  %>% filter(LOO_COUNT > 6)

## Check for tier 3 traits that pass the above

List_tier_3 <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "outcome"),
                      list(Het_filt, Egger_filt, LOO_count))

List_tier_3 <- List_tier_3 %>% filter(!duplicated(outcome))

Tier_3_Merge_pass_het_and_Egger_test <- merge(Tier_3_traits, List_tier_3, by = "outcome")


Tier_3_Merge_pass_het_and_Egger_test <- Tier_3_Merge_pass_het_and_Egger_test %>% filter(LOO_COUNT == nsnp)

write.table(Tier_3_Merge_pass_het_and_Egger_test, file="FinnGen_Passing_traits_tier_3.txt",
           sep = "\t", row.names = F, quote = F)

Tier_3_plot <- merge(Tier_3_Merge_pass_het_and_Egger_test, N1, by = "outcome")

## Make forest plot

Tier_3_plot <- generate_odds_ratios(Tier_3_plot)
Tier_3_plot$METHOD <- c("IVW-FE", "MR-Egger",
                        "Weighted median", "Weighted mode", "IVW-MRE",
                        "IVW-MRE", "Weighted mode", "MR-Egger", "Weighted median",
                        "IVW-FE", "IVW-MRE", "Weighted median", "Weighted mode",
                        "IVW-FE", "MR-Egger")

IVW_MRE <- Tier_3_plot %>% filter(METHOD == "IVW-MRE")

Other <- Tier_3_plot %>% filter(METHOD != "IVW-MRE")

FP_1 <- ggplot(data = IVW_MRE, aes(x=phenostring, y=or,
                                       ymax=or_uci95, ymin=or_lci95, 
                                       colour=phenostring)) +
  facet_wrap(~METHOD, nrow = 3, ncol = 2) + 
  coord_flip() +
  theme_bw() +
  geom_pointrange() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 1, lty = "dashed") +
  ylab("OR [95% CI] per SD increase in circulating retinol") +
  xlab(" ") +
  scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")

FP_2 <- ggplot(data = Other, aes(x=phenostring, y=or,
                           ymax=or_uci95, ymin=or_lci95, 
                           colour=phenostring)) +
  facet_wrap(~METHOD, nrow = 3, ncol = 2) + 
  coord_flip() +
  theme_bw() +
  geom_pointrange() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 1, lty = "dashed") +
  ylab("OR [95% CI] per SD increase in circulating retinol") +
  xlab(" ") +
  scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")

library(ggpubr)

ggarrange(FP_1, FP_2, nrow = 2)
## Exploratory IVW-MRE in matched category

table(Tier_3_plot$category)

FDR_sig_name <- generate_odds_ratios(FDR_sig_name)

FDR_sig_name$Z <- FDR_sig_name$b/FDR_sig_name$se

OTHER <- FDR_sig_name %>% filter(category == "Cardiometabolic endpoints" | category == "Certain infectious and parasitic diseases (AB1_)" | category == "Diseases of the digestive system (K11_)" |
category == "Endocrine, nutritional and metabolic diseases (E4_)")

PSYCH <- FDR_sig_name %>% filter(category == "Psychiatric endpoints from Katri Räikkönen" |
                                   category == "Mental and behavioural disorders (F5_)")

CANCER <- FDR_sig_name %>% filter(category == "Neoplasms, from cancer register (ICD-O-3)")

ggplot(data = OTHER, aes(x=Z, fill = category)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~category) +
  theme_bw() +
  theme(legend.position = "right") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("Density") +
  xlab("MR Z score - effect of retinol") +
  scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")

## Psych endpoints

ggplot(data = CANCER, aes(x=phenostring, y=or,
                               ymax=or_uci95, ymin=or_lci95, 
                               colour=phenocode)) +
  coord_flip() +
  theme_bw() +
  geom_pointrange() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 1, lty = "dashed") +
  ylab("OR [95% CI] per SD increase in circulating retinol") +
  xlab(" ") +
  theme(legend.position = "none") +
  ggtitle("Neoplasms - cancer registry")

LOO_raw %>% filter(outcome == "C3_NSCLC_SQUAM_EXALLC")

## Lung cancer leave-one-out, driven by RBP4
