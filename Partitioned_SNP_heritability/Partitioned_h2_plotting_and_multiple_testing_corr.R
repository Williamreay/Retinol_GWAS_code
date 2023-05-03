#################################

## Partitioned h2 - GTEx/Franke tissues and cell types

## William Reay (2022)

#################################

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

setwd("~/Desktop/Retinol_GWAS/Partitioned_h2/")

## Read in all three Stouffer meta partitioned h2

All_three_raw <- fread("Stouffer_all_three_partitioned_h2.cell_type_results.txt")

All_three_raw$FWER <- p.adjust(All_three_raw$Coefficient_P_value, method = "bonferroni")
All_three_raw$FDR <- p.adjust(All_three_raw$Coefficient_P_value, method = "fdr")

write.csv(All_three_raw, file="Multiple_testing_corr_all_three_Stouffer.csv",
          row.names = F, quote = F)

## Make plot

Input_1 <- All_three_raw %>% filter(Coefficient_P_value < 0.01)

P1 <- ggplot(data = Input_1, aes(x=Name, y=-log10(Coefficient_P_value), fill=Name)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  coord_flip() + theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size=10)) +
  geom_hline(yintercept=1.30103, lty=2) +
  ylab("-log10 P-value") + xlab(" ") +
  ggtitle("All three - Heritability Enrichment (P < 0.01)") +
  ylim(0,9)

## Read in METSIM+INTERVAL

MET_INT <- fread("Stouffer_INTERVAL_METSIM_partitioned_h2.cell_type_results.txt")

MET_INT$FWER <- p.adjust(MET_INT$Coefficient_P_value, method = "bonferroni")
MET_INT$FDR <- p.adjust(MET_INT$Coefficient_P_value, method = "fdr")

write.csv(MET_INT, file="Multiple_testing_corr_METISM_INTERVAL_Stouffer.csv",
          row.names = F, quote = F)


Input_2 <- MET_INT %>% filter(Coefficient_P_value < 0.01)

P2 <- ggplot(data = Input_2, aes(x=Name, y=-log10(Coefficient_P_value), fill=Name)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  coord_flip() + theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size=10)) +
  geom_hline(yintercept=1.30103, lty=2) +
  ylab("-log10 P-value") + xlab(" ") +
  ggtitle("METSIM+INTERVAL - Heritability Enrichment (P < 0.01)") +
  ylim(0,9)

ggarrange(P1, P2)
