##################################

## Effect of body fat % on retinol associated brain regions

## William Reay (2023)

##################################

library(data.table)
library(dplyr)
library(TwoSampleMR)

setwd("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/")

## Import IVs for body fat %

BF_IVs <- extract_instruments("ebi-a-GCST003435", clump = T)

## Read in brain measures to test

Brain_measures_raw <- fread("Processed_output/Tier_3/Passing_traits_tier_3.txt")

Brain_measures_filt <- Brain_measures_raw %>% select(id.outcome)
  
Brain_measures_filt <- unique(filter(Brain_measures_filt, grepl("ubm", id.outcome)))

List <- as.list(Brain_measures_filt$id.outcome)
## Extract IVs for those outcomes

Out_brain <- extract_outcome_data(snps = BF_IVs$SNP,
                                  outcomes = Brain_measures_filt$id.outcome)

## Harmonise

Harm_brain <- harmonise_data(BF_IVs, Out_brain, action = 3)

## MR

MR_brain <- mr(Harm_brain, method_list = c("mr_ivw_mre"))

## Output

write.table(MR_brain, "Body_fat_percentage_on_retinol_assoc_brain_Regions.txt",
            sep = "\t", row.names = F)
