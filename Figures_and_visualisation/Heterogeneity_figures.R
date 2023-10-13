##################################

## Heterogeneity figures - retinol lead SNPs

## Dr William Reay (2023)

###################################

setwd("~/Desktop/Retinol_GWAS/Manuscript/Nat_comm_resubmission/New_figures_or_analyses/")

library(dplyr)
library(data.table)
library(ggplot2)

## METSIM + INTERVAL

MET_1 <- fread("METSIM_lead_SNPs_meta_1.txt")
MET_1$Study <- "METSIM"
MET_1 <- MET_1 %>% select(rsids, ref, alt, beta, sebeta, Study)
MET_1$Z <- MET_1$beta/MET_1$sebeta

INT_1 <- fread("INTERVAL_lead_SNPs_meta_1.txt")
INT_1$Study <- "INTERVAL"
INT_1 <- INT_1 %>% select(rsid, REF, ALT, beta, standard_error, Study)
INT_1$Z <- INT_1$beta/INT_1$standard_error
INT_1 <- rename(INT_1, "rsids"="rsid", "ref"="REF", "alt"="ALT", "sebeta"="standard_error")

MET_INT_1 <- rbind(MET_1, INT_1)
MET_INT_1$LCI <- MET_INT_1$beta - (1.96*MET_INT_1$sebeta)
MET_INT_1$UCI <- MET_INT_1$beta + (1.96*MET_INT_1$sebeta)

P1 <- ggplot(MET_INT_1, aes(x=rsids, y=beta, ymin=LCI, ymax=UCI, colour=Study)) +
  geom_pointrange(position = position_dodge(width = 0.25)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  ylab("Effect on retinol (SD) per allele") +
  ggtitle("Per study effects: METSIM + INTERVAL (beta)") +
  geom_hline(yintercept = 0, lty ="dashed") +
  xlab("SNP")

P2 <- ggplot(MET_INT_1, aes(x=rsids, y=Z, fill=Study)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2") +
  ylab("Retinol Z score") +
  ggtitle("Per study effects: METSIM + INTERVAL (Z-score)") +
  geom_hline(yintercept = 0) +
  xlab("SNP")

library(ggpubr)

ggarrange(P1, P2, common.legend = T, legend = c("bottom"))

## All three

MET_2 <- fread("METSIM_lead_SNPs_meta_2.txt")
MET_2$Study <- "METSIM"
MET_2$Z <- MET_2$beta/MET_2$sebeta
MET_2 <- MET_2 %>% select(rsids, ref, alt, Z, Study)

INT_2 <- fread("INTERVAL_lead_SNPs_meta_2.txt")
INT_2$Study <- "INTERVAL"
INT_2$Z <- INT_2$beta/INT_2$standard_error
INT_2 <- INT_2 %>% select(rsid, REF, ALT, Z, Study)
INT_2 <- rename(INT_2, "rsids"="rsid", "ref"="REF", "alt"="ALT")

A_P <- fread("ATBC_PLCO_lead_SNPs_meta_2.txt")
A_P <- A_P %>% select(SNP, A_EFFECT, A_ALT, Z)
A_P <- rename(A_P, "rsids"="SNP", "alt"="A_EFFECT", "ref"="A_ALT")
A_P$Study <- "ATBC+PLCO"

All_three <- rbind(MET_2, INT_2, A_P)
All_three$Z <- All_three$Z*-1

## Align to effect allele in table1

All_three <- All_three %>% mutate(case_when(rsids == "rs1260326" ~ Z*-1, 
                                            rsids == "rs6601299" ~ Z*-1,
                                            
                                            TRUE ~ Z))

ggplot(All_three, aes(x=rsids, y=Z, fill = Study)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired") +
  ylab("Retinol Z score") +
  ggtitle("Per study effects: METSIM + INTERVAL + ATBC/PLCO (Z-score)") +
  geom_hline(yintercept = 0) +
  xlab("SNP")

## Heterogeneity across all HapMap SNPs

MET_INT <- fread("../../../Meta_analyses/Common_var/Stouffers_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

HapMap <- fread("w_hm3.noMHC.snplist")
HapMap <- rename(HapMap, "rsids"="SNP")

Merge_MET_INT <- merge(MET_INT, HapMap, by="rsids")

H1 <- ggplot(Merge_MET_INT, aes(y=-log10(HetPVal), x=-log10(`P-value`))) +
  geom_point() +
  geom_hline(yintercept = 5, lty ="dashed", colour = "red") +
  geom_vline(xintercept = 5, lty ="dashed", colour = "red") +
  theme_bw() +
  xlab("-log10 GWAS P-value") +
  ylab("-log10 Heterogeneity P-value") +
  ggtitle("METSIM+INTERVAL")

All_three <- fread("../../../Meta_analyses/Common_var/Stouffers_Full_METSIM_INTERVAL_PLCO_ATBC1.tbl.gz")
All_three <- rename(All_three, "rsids"="MarkerName")

Merge_all <- merge(All_three, HapMap, by="rsids")

ggplot(Merge_all, aes(y=-log10(HetPVal), x=-log10(`P-value`))) +
  geom_point() +
  geom_hline(yintercept = 5, lty ="dashed", colour = "red") +
  geom_vline(xintercept = 5, lty ="dashed", colour = "red") +
  theme_bw() +
  xlab("-log10 GWAS P-value") +
  ylab("-log10 Heterogeneity P-value") +
  ggtitle("METSIM+INTERVAL+ATBC+PLCO")
