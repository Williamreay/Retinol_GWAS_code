#################################################

## Retinol RapidoPGS style ABF calculation for PGS weighting

## William Reay (2023)

#################################################

library(dplyr)
library(data.table)
library(RapidoPGS)

## Load ABF function from GitHub - https://github.com/Williamreay/Pneumonia_meta_GWAS/blob/master/Finemapping/ABF_function.R
source("ABF_function.R")

setwd("/home/control/data/users/william/Meta_retinol_GWAS/")

## Import Sumstats 

Ret_raw <- fread("IVW_input_for_PGS_calc_MET_INT.txt", header = T)

## Approximate W using h2 - W - h2/v(1/Nsigma)

Ret_formatted <- rename(Ret_raw, "SNPID"="CHR_BP", "ALT"="Allele1", "REF"="Allele2",
                        "BETA"="Effect", "SE"="StdErr")

sd.prior.est(Ret_formatted, 0.0668, 17268)

## W = 0.01853208

Hg19_LDblocks <- fread("RapidoPGS_approach/fourier_ls-all.bed")

## Write function to deploy ABF function for all loci

ABF_automate <- function(chr, BP1, BP2, df, W) {
  Subset <- df %>% filter(CHR == chr & BP > BP1 & BP < BP2)
  ABF <- Finemapping_abf(Subset$CHR_BP, Subset$Effect, Subset$StdErr, W = W)
  return(ABF)
}

chr_list <- as.list(Hg19_LDblocks$chr)
BP1_list <- as.list(Hg19_LDblocks$start)
BP2_list <- as.list(Hg19_LDblocks$stop)

Apply_ABF <- mapply(ABF_automate, chr = chr_list, BP1 = BP1_list, BP2 = BP2_list,
                    MoreArgs = list(Ret_raw, 0.01853208),
                    SIMPLIFY = F)

Extract_ABF <- data.frame()

for (i in 1:length(BP1_list)) {
  Extract_ABF <- rbind(Extract_ABF, Apply_ABF[[i]])
}

## Merge and multiply by PP

Extract_ABF <- rename(Extract_ABF, "CHR_BP"="SNP")

Merged_inferred_prior_SD <- merge(Extract_ABF, Ret_raw, by = "CHR_BP")

## Repeat with 0.15 prior SD

Apply_ABF_2 <- mapply(ABF_automate, chr = chr_list, BP1 = BP1_list, BP2 = BP2_list,
                      MoreArgs = list(Ret_raw, 0.15),
                      SIMPLIFY = F)

Extract_ABF_2 <- data.frame()

for (i in 1:length(BP1_list)) {
  Extract_ABF_2 <- rbind(Extract_ABF_2, Apply_ABF_2[[i]])
}

Extract_ABF_2 <- rename(Extract_ABF_2, "CHR_BP"="SNP")

Merged_default_prior <- merge(Extract_ABF_2, Ret_raw, by = "CHR_BP")

Merged_default_prior_SD$PP_Beta <- Merged_default_prior_SD$PP*Merged_default_prior_SD$Effect

write.table(Merged_default_prior_SD, file="PP_weighted_PGS_default_SD_prior_IVW.txt", sep = "\t",
            row.names = F, quote = F)
