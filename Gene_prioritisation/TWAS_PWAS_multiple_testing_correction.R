################################################

## TWAS/PWAS multiple-testing correction and filtering

## Dr William Reay (2023)

################################################

library(dplyr)
library(data.table)

setwd("~/Desktop/Retinol_GWAS/TWAS_PWAS/")

Raw_df <- fread("NEW_full_coloc_METSIM_INTERVAL_TWAS_blood_adipose_pancreas_breast_intestine_liver.txt")

## Remove NA TWAS.Z

Filt_df <- Raw_df %>% filter(!is.na(TWAS.Z))

## FDR and FWER

Filt_df$FWER <- p.adjust(Filt_df$TWAS.P, method="bonferroni")
Filt_df$FDR <- p.adjust(Filt_df$TWAS.P, method = "fdr")

## Annotate with names

Gene_names <- rtracklayer::import("OLD/Multiple_testing_corr/gencode.v26.GRCh38.genes.gtf")
Gene_names <- as.data.frame(Gene_names)
Gene_names <- Gene_names %>% filter(!duplicated(transcript_id))

Gene_names <- rename(Gene_names, "ID"="gene_id")

Merged_filt <- merge(Filt_df, Gene_names, by= "ID")

write.table(Merged_filt, file="TWAS_METSIM_INTERVAL_annot_multiple_testing_corr.txt",
            sep = "\t", row.names = F, quote = F)

FDR_coloc <- Merged_filt %>% filter(FDR < 0.05 & COLOC.PP4 > 0.8)

write.table(FDR_coloc, file="FDR_sig_coloc_H4_over_0.8_METSIM_INTERVAL_TWAS.txt",
            sep = "\t", row.names = F)

## Colocalised transcriptome-wide (PPh4 > 0.8)

Coloc_transcriptome_wide <- Merged_filt %>% filter(COLOC.PP4 > 0.8)

Coloc_transcriptome_wide <- Coloc_transcriptome_wide[order(Coloc_transcriptome_wide$ID, -abs(Coloc_transcriptome_wide$TWAS.Z) ), ]

Coloc_transcriptome_wide <- Coloc_transcriptome_wide[ !duplicated(Coloc_transcriptome_wide$ID), ]

## Suggestive colocalisation transcriptome-wide (PPh4 > 0.4)

Sug_coloc_transcriptome_wide <- Merged_filt %>% filter(COLOC.PP4 > 0.4 & COLOC.PP1 < 0.4 & COLOC.PP2 < 0.4)

Sug_coloc_transcriptome_wide <- Sug_coloc_transcriptome_wide[order(Sug_coloc_transcriptome_wide$ID, -abs(Sug_coloc_transcriptome_wide$TWAS.Z) ), ]

Sug_coloc_transcriptome_wide <- Sug_coloc_transcriptome_wide[ !duplicated(Sug_coloc_transcriptome_wide$ID), ]

## Repeat above for PWAS and then merge the output for cMAP

PWAS_Raw_df <- fread("NEW_blood_PWAS_METSIM_INTERVAL.txt")

## Remove NA TWAS.Z

PWAS_Filt_df <- PWAS_Raw_df %>% filter(!is.na(PWAS.Z))

## FDR and FWER

PWAS_Filt_df$FWER <- p.adjust(PWAS_Filt_df$PWAS.P, method="bonferroni")
PWAS_Filt_df$FDR <- p.adjust(PWAS_Filt_df$PWAS.P, method = "fdr")

write.table(PWAS_Filt_df, file="PWAS_METSIM_INTERVAL_annot_multiple_testing_corr.txt",
            sep = "\t", row.names = F, quote = F)

PWAS_FDR_coloc <- PWAS_Filt_df %>% filter(FDR < 0.05 & COLOC.PP4 > 0.8)

write.table(PWAS_FDR_coloc, file="FDR_sig_coloc_H4_over_0.8_METSIM_INTERVAL_PWAS.txt",
            sep = "\t", row.names = F)

## Colocalised proteome wide (PPh4 > 0.8)

Coloc_proteome_wide <- PWAS_Filt_df %>% filter(COLOC.PP4 > 0.8)

Coloc_proteome_wide <- Coloc_proteome_wide[order(Coloc_proteome_wide$ID, -abs(Coloc_proteome_wide$PWAS.Z) ), ]

Coloc_proteome_wide <- Coloc_proteome_wide[ !duplicated(Coloc_proteome_wide$ID), ]

## Suggestive colocalisation proteome-wide (PPh4 > 0.4) && PP.H1 (only associated with expression) or H2 (only GWAS) < 0.4

Sug_coloc_proteome_wide <- PWAS_Filt_df %>% filter(COLOC.PP4 > 0.4 & COLOC.PP1 < 0.4 & COLOC.PP2 < 0.4)

Sug_coloc_proteome_wide <- Sug_coloc_proteome_wide[order(Sug_coloc_proteome_wide$ID, -abs(Sug_coloc_proteome_wide$PWAS.Z) ), ]

Sug_coloc_proteome_wide <- Sug_coloc_proteome_wide[ !duplicated(Sug_coloc_proteome_wide$ID), ]


## Merge outcome - colocalised transcriptome/proteome wide and suggestive colocalisation transcriptome or proteome-wide

TWAS_coloc <- Coloc_transcriptome_wide %>% select(gene_name, TWAS.Z, TWAS.P, COLOC.PP4)
TWAS_coloc <- rename(TWAS_coloc, "ID"="gene_name", "Z"="TWAS.Z", "P"="TWAS.P")
PWAS_coloc <- Coloc_proteome_wide %>% select(ID, PWAS.Z, PWAS.P, COLOC.PP4)
PWAS_coloc <- rename(PWAS_coloc, "Z"="PWAS.Z", "P"="PWAS.P")

Colocalised_PWAS_TWAS <- rbind(TWAS_coloc, PWAS_coloc)

## Add RBP4 pQTL MR analyses

RBP4 <- c("RBP4", "9.7837", "1.23847737511899e-19", "9.999925e-01")

Colocalised_PWAS_TWAS_RBP4 <- rbind(Colocalised_PWAS_TWAS, RBP4)

write.table(Colocalised_PWAS_TWAS_RBP4, file="Colocalised_H4_0.8_TWAS_PWAS_pQTL_MR_RBP4.txt",
            sep = "\t", row.names = F, quote = F)

