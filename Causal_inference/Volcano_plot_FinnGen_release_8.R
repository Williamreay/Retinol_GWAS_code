###############################

## FinnGen forest plot (release 8)

## Dr William Reay (2023)

###############################

library(ggplot2)
library(data.table)
library(dplyr)
library(TwoSampleMR)

setwd("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/FinnGen_r8_MR/")

IVW <- fread("FinnGen_r8_retinol_MR_Ncase_over_1000.txt")

N3 <- fread("18_31593832_A_G_phenotype_associations.tsv")
N3$N <- paste("Number of cases = ", N3$n_case, sep = "")
N3$Plot_name <- paste(N3$phenostring, " (",N3$N,")", sep="")

N3 <- N3 %>% select(phenocode, phenostring, category, Plot_name)
N3$outcome <- N3$phenocode
Merge <- merge(IVW, N3, by = "outcome")

FDR_sig <- Merge %>% filter(method == "Inverse variance weighted (multiplicative random effects)")
FDR_sig$FDR <- p.adjust(FDR_sig$pval, method="fdr")

FDR_sig <- generate_odds_ratios(FDR_sig)

FDR_plt <- FDR_sig %>% filter(FDR < 0.01)

ggplot(data = FDR_plt, aes(x=Plot_name, y=or,
                               ymax=or_uci95, ymin=or_lci95, 
                               colour=category)) +
  coord_flip() +
  theme_bw() +
  geom_pointrange() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 1, lty = "dashed") +
  ylab("OR [95% CI] per SD increase in circulating retinol") +
  xlab(" ") +
  scale_color_brewer(palette = "Paired") +
  theme(legend.position = "right") +
  ggtitle("FDR < 0.01 - IVW-MRE")

## Volcano plot - first paste in formatted input to rename the categories

FDR_sig$category <- as.factor(FDR_sig$category)

levels(FDR_sig$category)[levels(FDR_sig$category)=="Certain infectious and parasitic diseases (AB1_)"] <- "Infectious disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the digestive system (K11_)"] <- "Diseases of digestive system"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Neurological endpoints"] <-"Neurological and Psychiatric"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Mental and behavioural disorders (F5_)"] <-"Neurological and Psychiatric"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Alcohol related diseases"] <-"Neurological and Psychiatric"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Asthma and related endpoints"] <-"Respiratory disorders"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Comorbidities of Asthma"] <-"Respiratory disorders"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Comorbidities of Neurological endpoints"] <-"Neurological and Psychiatric"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases marked as autimmune origin"] <-"Inflammatory/Automimmune related"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Endocrine, nutritional and metabolic diseases (E4_)"] <-"Metabolic/Endocrine disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Miscellaneous, not yet classified endpoints"] <-"Miscellaneous"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Neoplasms, from cancer register (ICD-O-3)"] <-"Neoplasms"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Comorbidities of COPD"] <-"Respiratory disorders"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Neoplasms from hospital discharges (CD2_)"] <-"Neoplasms"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Gastrointestinal endpoints"] <- "Diseases of digestive system"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Cardiometabolic endpoints"] <-"Metabolic/Endocrine disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="COPD and related endpoints"] <- "Respiratory disorders"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Interstitial lung disease endpoints"] <-"Respiratory disorders"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Comorbidities of Diabetes"] <- "Metabolic/Endocrine disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)"] <-"Inflammatory/Automimmune related"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Common endpoint"] <-"Miscellaneous"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Dental endpoints"] <- "Dental"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diabetes endpoints"] <- "Metabolic/Endocrine disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the eye and adnexa (H7_)"] <- "Opthalamological"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the nervous system (G6_)"] <- "Neurological and Psychiatric"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the circulatory system (I9_)"] <-"Cardiovascular Disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Rheuma endpoints"] <- "Rheumatological endpoints"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the ear and mastoid process (H8_)"] <-"Ear disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the respiratory system (J10_)"] <-"Respiratory disorders"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the musculoskeletal system and connective tissue (M13_)"] <- "Musculoskeletal/Connective Tissue disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Psychiatric endpoints from Katri Räikkönen"] <- "Neurological and Psychiatric"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the skin and subcutaneous tissue (L12_)"] <-"Skin/subcutaneous disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Other, not yet classified endpoints (same as #MISC)"] <-"Miscellaneous"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Diseases of the genitourinary system (N14_)"] <-"Genitourinary disease"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Comorbidities of Gastrointestinal endpoints"] <-"Diseases of digestive system"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Pregnancy, childbirth and the puerperium (O15_)"] <-"Pregnancy, childbirth and the puerperium"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Congenital malformations, deformations and chromosomal abnormalities (Q17)"] <- "Congential Malformations"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R18_)"] <- "Miscellaneous"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Drug purchase endpoints"] <-"Drug purchase endpoints"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Comorbidities of Interstitial lung disease endpoints"] <- "Respiratory disorders"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Injury, poisoning and certain other consequences of external causes (ST19_)"] <-"Injury/poisoning"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Codes for special purposes (U22_)"] <-"Miscellaneous"
levels(FDR_sig$category)[levels(FDR_sig$category)=="Factors influencing health status and contact with health services (Z21_)"] <-"Miscellaneous"

FDR_sig <- rename(FDR_sig, "Category"="category")

V<- ggplot(data = FDR_sig,
       aes(x=b, y = -log10(pval),
           colour = Category)) +
  geom_point(alpha = 0.8, size = 1.75) +
  geom_hline(yintercept = -log10(8.191719e-05), lty = "dashed") +
  labs(x = expression("Log odds per SD in circulating retinol"), y = expression(paste("-log"[10], "P-value"))) +
  theme_bw() +
  theme(legend.position ="bottom") +
  guides(colour=guide_legend(ncol=3)) 
  theme(axis.title = element_text(face="bold", size=12))
  ggtitle("FVC:Class B/2 secretin family receptors") +
  xlim(c(-4,4)) +
  ylim(c(0, 5))

V + scale_color_manual(values = c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1",
    "blue1", "steelblue4",
    "darkturquoise", "green1"))


