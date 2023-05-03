#############################

## Effect of lead SNPs on estimated retinol intake (UKBB - IEUGWAS)

## Dr William Reay (2023)

#########################

library(ieugwasr)
library(dplyr)
library(ggplot2)
library(viridis)

Ret_intake <- associations(c("rs1260326", "rs34898035", "rs11762406",
                             "rs6601299", "rs10882283", "rs12149203",
                             "rs1667226", "rs2207132", "rs6029188"), id = "ukb-b-17406")

write.table(Ret_intake, file="~/Desktop/Retinol_GWAS/UKBB_retinol_intake_lead_SNPs/Lead_snps_retinol_intake.txt",
            sep = "\t", row.names = F, quote = F)

Ret_intake$LCI <- (Ret_intake$beta - (1.96*Ret_intake$se))
Ret_intake$UCI <- (Ret_intake$beta + (1.96*Ret_intake$se))

ggplot(Ret_intake, aes(x=target_snp, y=beta, ymin = LCI, ymax = UCI, colour = target_snp)) +
  geom_pointrange() +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Effect on estimated retinol intake (SD) per effect allele") +
  xlab(" ") +
  theme(legend.position = "none") +
  scale_color_viridis(discrete = T)
