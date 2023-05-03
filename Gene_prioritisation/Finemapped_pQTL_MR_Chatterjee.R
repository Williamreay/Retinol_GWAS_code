####################################

## Retinol - Finemapped pQTL MR

## William Reay (2022)

#####################################

library(TwoSampleMR)
library(dplyr)
library(data.table)
library(hudson)
library(ggplot2)
library(viridis)
library(ggpubr)

setwd("~/Desktop/Retinol_GWAS/ARIC_MR_and_coloc/")

## Read in finemapped Chatterjee pQTL (PIP > 0.5, SuSiE)

Raw_IVs <- fread("MR_input_finemapped_PIP_over_0.5.txt")
Raw_IVs$A2 <- ifelse(Raw_IVs$A1 == Raw_IVs$ALT, Raw_IVs$REF, Raw_IVs$ALT)

Formatted_IVs <- format_data(dat = Raw_IVs, type="exposure",
                             effect_allele_col = "A1",
                             other_allele_col = "A2",
                             beta_col = "BETA", se_col = "SE",
                             phenotype_col = "Protein_ID",
                             snp_col = "ID")

## Check LD clumping

Clumped_IVs <- clump_data(Formatted_IVs)

## Save Clumped IVs as Rdata frame

saveRDS(Clumped_IVs, file="Finemapped_clumped_chatterjee_EA.rds")

## Read in retinol GWAS

Ret_GWAS <- fread("../Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")
Ret_GWAS$Allele1 <- toupper(Ret_GWAS$Allele1)
Ret_GWAS$Allele2 <- toupper(Ret_GWAS$Allele2)

Outcome_formatted <- format_data(dat = Ret_GWAS, type = "outcome",
                                 snps = Clumped_IVs$SNP, snp_col = "rsids",
                                 beta_col = "Effect", se_col = "StdErr",
                                 effect_allele_col = "Allele1",
                                 other_allele_col = "Allele2")

## Harmonise - remove palindromic SNPs

Harmonised_pQTL_retinol <- harmonise_data(Clumped_IVs, Outcome_formatted, action =  3)

## MR tests - Wald ratio and IVW-FE (when > 1 SNP)

pQTL_MR <- mr(Harmonised_pQTL_retinol, method_list= c("mr_ivw_fe", "mr_wald_ratio"))

## Output data

pQTL_MR$FWER <- p.adjust(pQTL_MR$pval, method = "bonferroni")
pQTL_MR$FDR <- p.adjust(pQTL_MR$pval, method = "fdr")

Gene_IDs <- fread("seqid.txt")
Gene_IDs <- rename(Gene_IDs, "exposure"="seqid_in_sample")

pQTL_MR_gene_annot <- merge(pQTL_MR, Gene_IDs, by = "exposure")

write.table(pQTL_MR_gene_annot, file="MR_finemapped_SNPs/pQTL_chatterjee_finemapped_SNPs_MR.txt",
            sep = "\t", row.names = F, quote = F)

## Make forest plot and Miami plot

Miami_plot_input <- pQTL_MR_gene_annot %>% select(entrezgenesymbol, chromosome_name,
                                                  transcription_start_site, pval, b)
Miami_plot_input <- rename(Miami_plot_input, "SNP"="entrezgenesymbol",
                           "CHR"="chromosome_name", "POS"="transcription_start_site",
                           "pvalue"="pval")
Up <- Miami_plot_input %>% filter(b > 0)
Down <- Miami_plot_input %>% filter(b < 0)

Miami <- gmirror(top=Up, bottom=Down, tline=5.336179e-05, bline=5.336179e-05,
        highlight_p = c(5.336179e-05, 5.336179e-05), highlighter = "blue",
        toptitle = "Finemapped plasma pQTLs as IVs")

FP_input <- pQTL_MR_gene_annot %>% filter(FWER < 0.05)
FP_input$LB <- FP_input$b - (1.96*FP_input$se)
FP_input$UB <- FP_input$b + (1.96*FP_input$se)

FP <- ggplot(data = FP_input, aes(x=entrezgenesymbol, y=b,
                            ymax=UB, ymin=LB, 
                            colour=entrezgenesymbol)) +
  coord_flip() +
  theme_bw() +
  geom_pointrange() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Effect on retinol (SD) per SD increase in protein abundance") +
  xlab("Protein") +
  scale_colour_viridis(discrete = T, option = "ggplot default")

ggarrange(Miami, FP,
          nrow = 1)
