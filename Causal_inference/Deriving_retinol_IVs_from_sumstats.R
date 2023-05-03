library(dplyr)
library(data.table)

SNPs <- fread("~/Desktop/Retinol_GWAS/Meta_analyses/FUMA/Stouffer_METSIM_INTERVAL/rs.txt")
SNPs <- rename(SNPs, "rsids"="rsID")

Ret <- fread("~/Desktop/Retinol_GWAS/Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

Merge <- merge(SNPs, Ret, by = "rsids")

Merge <- Merge %>% select(rsids, Allele1, Allele2, Effect, StdErr, `P-value`)

Merge$Allele1 <- toupper(Merge$Allele1)
Merge$Allele2 <- toupper(Merge$Allele2)

Merge <- rename(Merge, "beta"="Effect", "se"="StdErr", 
                "P"="P-value")

write.table(Merge, file="~/Desktop/Retinol_GWAS/Retinol_exp_all_IVs/All_IVs_retinol_METSIM_INTERVAL.txt",
            sep = "\t", row.names = F, quote = F)
