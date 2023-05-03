####################################

## AshR per bin of LD score

## METSIM+INTERVAL meta

## William Reay (2022)

####################################

suppressMessages(library(ashr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

## Load GWAS

METSIM_INTERVAL <- fread("~/data/users/william/Meta_retinol_GWAS/Meta_analysis_output/Annotated_sumstats/IVW_METSIM_INTERVAL_sumstats_retinol.txt")
METSIM_INTERVAL$SNP <- METSIM_INTERVAL$rsids

## Load LD scores (HapMap3)

HapMap3_LD_score <- fread("unannotated_LDscore.l2.ldsc", header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))

Merged_df <- merge(HapMap3_LD_score, METSIM_INTERVAL, by="SNP")

## Sort by LD score

Sorted_df <- Merged_df[order(Merged_df[,L2]),]

Sorted_df$L2 <- as.numeric(Sorted_df$L2)

## Load Hm3 alleles and merge with allele list

Hm3_alleles <- fread("../Cytokine_pneumonia_lung_function/ldsc/w_hm3.noMHC.snplist",
                     header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))

Sorted_df <- merge(Sorted_df, Hm3_alleles, by = "SNP")

## Check if any mismatches
Sorted_df <- Sorted_df %>% filter(A1 != Allele1 | A1 != Allele2)

## Remove any duplicates

Sorted_df <- Sorted_df %>% filter(duplicated(SNP)==FALSE)

Sorted_df$Z <- (Sorted_df$Effect/Sorted_df$StdErr)

## Method adapted from Boyle et al, omnigenic paper

## Genome-wide ashR

ash(Sorted_df$Effect, Sorted_df$StdErr)$fitted_g$pi

## 0.1867836 proportion of SNPs with non-null effects (both causal and tagging SNPs)

## 1000 LD score bins

EB_ASH_1000 <- Sorted_df %>% 
  mutate(ldBin = ntile(L2, 1000)) %>% group_by(ldBin) %>% 
  summarize(mean.ld = mean(L2), se.ld=sd(L2)/sqrt(n()), 
            mean.chisq = mean(Z**2, na.rm=T), 
            se.chisq=sd(Z**2, na.rm=T)/sqrt(sum(!is.na(Z))), 
            prop.null = ash(Effect, StdErr)$fitted_g$pi[1], n=n())

Global_mean_1000_bins <- mean((1-EB_ASH_1000$prop.null))
Global_median_1000_bins <- median((1-EB_ASH_1000$prop.null))

## Mean = 0.02344151 , ## Median = 0.003283514

## 5000 LD score bins

EB_ASH_5000 <- Sorted_df %>% 
  mutate(ldBin = ntile(L2, 5000)) %>% group_by(ldBin) %>% 
  summarize(mean.ld = mean(L2), se.ld=sd(L2)/sqrt(n()), 
            mean.chisq = mean(Z**2, na.rm=T), 
            se.chisq=sd(Z**2, na.rm=T)/sqrt(sum(!is.na(Z))), 
            prop.null = ash(Effect, StdErr)$fitted_g$pi[1], n=n())


Global_mean_5000_bins <- mean((1-EB_ASH_5000$prop.null))
Global_median_5000_bins <- median((1-EB_ASH_5000$prop.null))

LD_1000G_1 <- ggplot(EB_ASH_1000, aes(x=ldBin, y=(1-prop.null))) +
  geom_point(cex = 0.9) +
  geom_smooth() +
  #geom_text(x=380, y=0.7, label=paste("Mean fraction = ", Global_mean_1000_bins, sep = "")) +
  ggtitle("Retinol (1000G EUR LD)") +
  xlab("LD score bin (N=1000)") +
  ylab("Proportion of non-zero effect sizes") +
  ylim(0, 1) +
  theme_bw()

LD_1000G_2 <- ggplot(EB_ASH_5000, aes(x=ldBin, y=(1-prop.null))) +
  geom_point(cex=0.9) +
  geom_smooth() +
  #geom_text(x=2000, y = 0.7, label=paste("Mean fraction = ", Global_mean_5000_bins, sep = "")) +
  ggtitle("Retinol (1000G EUR LD)") +
  xlab("LD score bin (N=5000)") +
  ylab("Proportion of non-zero effect sizes") +
  ylim(0, 1) +
  theme_bw()


## UKBB LD (EUR)

LDscores_UKBB <- fread("../Meta_retinol_GWAS/Meta_analysis_output/UKBB.EUR.rsid.l2.ldscore.gz")

Merged_df_2 <- merge(LDscores_UKBB, METSIM_INTERVAL, by="SNP")

## Sort by LD score

Sorted_df_2 <- Merged_df_2[order(Merged_df_2[,L2]),]

Sorted_df_2$L2 <- as.numeric(Sorted_df_2$L2)

## Load Hm3 alleles and merge with allele list

Sorted_df_2 <- merge(Sorted_df_2, Hm3_alleles, by = "SNP")

## Check if any mismatches
Sorted_df_2 <- Sorted_df_2 %>% filter(A1 != Allele1 | A1 != Allele2)

## Remove any duplicates

Sorted_df_2 <- Sorted_df_2 %>% filter(duplicated(SNP)==FALSE)

Sorted_df_2$Z <- (Sorted_df_2$Effect/Sorted_df_2$StdErr)

## Repeat the above

ash(Sorted_df_2$Effect, Sorted_df_2$StdErr)$fitted_g$pi

## 0.1534897 non_null

UKBB_EB_ASH_1000 <- Sorted_df_2 %>% 
  mutate(ldBin = ntile(L2, 1000)) %>% group_by(ldBin) %>% 
  summarize(mean.ld = mean(L2), se.ld=sd(L2)/sqrt(n()), 
            mean.chisq = mean(Z**2, na.rm=T), 
            se.chisq=sd(Z**2, na.rm=T)/sqrt(sum(!is.na(Z))), 
            prop.null = ash(Effect, StdErr)$fitted_g$pi[1], n=n())

UKBB_Global_mean_1000_bins <- mean((1-UKBB_EB_ASH_1000$prop.null))
UKBB_Global_median_1000_bins <- median((1-UKBB_EB_ASH_1000$prop.null))

## 5000 LD score bins

UKBB_EB_ASH_5000 <- Sorted_df_2 %>% 
  mutate(ldBin = ntile(L2, 5000)) %>% group_by(ldBin) %>% 
  summarize(mean.ld = mean(L2), se.ld=sd(L2)/sqrt(n()), 
            mean.chisq = mean(Z**2, na.rm=T), 
            se.chisq=sd(Z**2, na.rm=T)/sqrt(sum(!is.na(Z))), 
            prop.null = ash(Effect, StdErr)$fitted_g$pi[1], n=n())


UKBB_Global_mean_5000_bins <- mean((1-UKBB_EB_ASH_5000$prop.null))
UKBB_Global_median_5000_bins <- median((1-UKBB_EB_ASH_5000$prop.null))

LD_UKBB_1 <- ggplot(UKBB_EB_ASH_1000, aes(x=ldBin, y=(1-prop.null))) +
  geom_point(cex = 0.9) +
  geom_smooth(col="purple") +
  #geom_text(x=380, y=0.7, label=paste("Mean fraction = ", UKBB_Global_mean_1000_bins, sep = "")) +
  ggtitle("Retinol (UKBB EUR LD)") +
  xlab("LD score bin (N=1000)") +
  ylab("Proportion of non-zero effect sizes") +
  ylim(0, 1) +
  theme_bw()

LD_UKBB_2 <- ggplot(UKBB_EB_ASH_5000, aes(x=ldBin, y=(1-prop.null))) +
  geom_point(cex = 0.9) +
  geom_smooth(col="purple") +
  #geom_text(x=2000, y=0.7, label=paste("Mean fraction = ", UKBB_Global_mean_5000_bins, sep = "")) +
  ggtitle("Retinol (UKBB EUR LD)") +
  xlab("LD score bin (N=5000)") +
  ylab("Proportion of non-zero effect sizes") +
  ylim(0, 1) +
  theme_bw()

ggarrange(LD_1000G_1, LD_1000G_2, LD_UKBB_1, LD_UKBB_2, nrow = 2, ncol = 2)

## Zoom in and plot as both GAM and LOESS
 ggplot(UKBB_EB_ASH_1000, aes(x=ldBin, y=(1-prop.null))) + geom_smooth(col="red") + geom_smooth(col="purple", method = "loess") + geom_smooth(col="blue", method = "lm") +
   theme_bw()
  


