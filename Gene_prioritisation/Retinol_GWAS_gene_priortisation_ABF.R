#################################################

## Retinol gene priortisation - ABF finemapping

## William Reay (2022)

#################################################

library(dplyr)
library(data.table)


## Load ABF function from GitHub - https://github.com/Williamreay/Pneumonia_meta_GWAS/blob/master/Finemapping/ABF_function.R
source("~/Desktop/Retinol_GWAS/Scripts/ABF_function.R")

## Read in Sumstats in hg38 - IVW

Retinol_INTERVAL_METSIM <- fread("~/Desktop/Retinol_GWAS/Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

## Loci (hg19)

#2:27598097:27752871
#2:122078406:122084285
#7:114014488:114286611
#8:9167797:9224907
#10:95295876:95360964
#16:79696939:79756197
#18:29134171:29190174
#20:39142516:39234223

## Post Liftover loci bounaries to hg38

Hg38_boundaries <- fread("~/Desktop/Retinol_GWAS/Gene_priortisation/Hg38_loci.txt")

## Write function to deploy ABF function for all loci - W = 0.01853208 from RapidoPGS comparing h2 to N see that paper

ABF_automate <- function(chr, BP1, BP2, df) {
  Subset <- df %>% filter(chrom == chr & pos > BP1 & pos < BP2)
  ABF <- Finemapping_abf(Subset$rsids, Subset$Effect, Subset$StdErr, W = 0.15)
  ABF$Locus_hg38 <- paste(chr,":",BP1,":", BP2, sep="")
  return(ABF)
}

chr_list <- as.list(Hg38_boundaries$V1)
BP1_list <- as.list(Hg38_boundaries$V2)
BP2_list <- as.list(Hg38_boundaries$V3)

Apply_ABF <- mapply(ABF_automate, chr = chr_list, BP1 = BP1_list, BP2 = BP2_list,
                    MoreArgs = list(Retinol_INTERVAL_METSIM),
                    SIMPLIFY = F)

Extract_ABF <- data.frame()

for (i in 1:length(BP1_list)) {
  Extract_ABF <- rbind(Extract_ABF, Apply_ABF[[i]])
}

write.table(Extract_ABF, file="~/Desktop/Retinol_GWAS/Gene_priortisation/New_prior_0.15_ABF_finemapping_IVW_METSIM_INTERVAL.txt",
            sep = "\t", row.names = F, quote = F)
