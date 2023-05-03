#############################
## Automated MR of all retinol IVs

## Adapted from script by Dr Dylan Kiltschewskij (2022)

## William Reay (2023)
#############################

## Load dependencies
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(MRInstruments))
suppressMessages(library(R.utils))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ggplot2))

# Specify command line inputs
option_list=list(
  make_option("--out_id", action="store",default=NA,type="character",help="IEU-GWAS id for the outcome trait [required]."),
  make_option("--exp", action="store",default=NA,type="character",help="Exposure file for retinol IVs"),
  make_option("--exp_name", action="store", defaul=NA,type="character",help="Exposure name [required]."),
  make_option("--oc_name", action="store", defaul=NA,type="character",help="Outcome name [required]."),
  make_option("--path", action="store", defaul=NA,type="character",help="Path to working directory [required].")
)
opt=parse_args(OptionParser(option_list=option_list))
setwd(paste(opt$path))

## Import IVs for ret per and clump for LD
cat("####################")
cat("\n")
cat("Reading in IVs for retinol")
cat("\n")

Exp_df <- fread(file=paste("",opt$exp,"",sep=""))
Exp_df$Pheno <- "Retinol"
Exp_df$N <- 17268

ivs <- format_data(Exp_df,
                   type = "exposure",
                   beta_col = "beta",
                   se_col = "se",
                   snp_col = "rsids",
                   effect_allele_col = "Allele1",
                   other_allele_col = "Allele2",
                   phenotype_col = "Pheno",
                   samplesize_col = "N",
                   pval_col = "P")

ivs<-clump_data(ivs)

## Extract instrument variables for IEUGWASdb GWAS id-code supplied
cat("####################")
cat("\n")
cat("Extracting instruments from GWAS")
cat("\n")
outcome2 <- extract_outcome_data(
  snps = ivs$SNP,
  outcomes = opt$out_id,
  proxies = F)

## Harmonise data and exclude palindromic SNPs
cat("####################")
cat("\n")
cat("Harmonising data")
cat("\n")
harmonised <- harmonise_data(exposure_dat = ivs, 
                             outcome_dat = outcome2, 
                             action = 3)

## Perform MR (TwoSampleMR included tests - IVW-MRE/FE, Weighted Median, Weighted, Mode, MR-Egger, robust and penalised)
cat("####################")
cat("\n")
cat("Performing MR")
cat("\n")
mr <- mr(harmonised ,method_list=c("mr_ivw_mre", "mr_ivw_fe", "mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
mr$ID <- opt$out_id
mr$Name <- opt$oc_name
write.table(mr, file=(paste("",opt$exp_name,"_",opt$out_id,"_MR.txt",sep="")), sep="\t",quote=F,row.names = F)

## Heterogeneity test via Cochran's Q
cat("####################")
cat("\n")
cat("Performing heterogeneity test")
cat("\n")
het <- mr_heterogeneity(harmonised)
write.table(het, file=(paste("",opt$exp_name,"_",opt$out_id,"_HET.txt",sep="")), sep="\t",quote=F,row.names = F)


## Test if Egger intercept is significantly different from zero
cat("####################")
cat("\n")
cat("Performing pleiotropy test")
cat("\n")
Egger_intercept <- mr_pleiotropy_test(harmonised)
write.table(Egger_intercept, file=(paste("",opt$exp_name,"_",opt$out_id,"_EI.txt",sep="")), sep="\t",quote=F,row.names = F)

## Leave one out analysis to test for IVs which disproprtionately contribute to the causal estimate
cat("####################")
cat("\n")
cat("Performing leave-one-out test")
cat("\n")
LOO <- mr_leaveoneout(harmonised)
write.table(LOO, file=(paste("",opt$exp_name,"_",opt$out_id,"_LOO.txt",sep="")), sep="\t",quote=F,row.names = F)
LOO_plot <- mr_leaveoneout_plot(LOO)
#ggsave(LOO_plot, file=(paste("",opt$exp_name,"_",opt$oc_name,"_LOO.png",sep="")),dpi=300)
pdf(paste("",opt$exp_name,"_",opt$out_id,"_LOO.pdf",sep=""))
LOO_plot[1]
dev.off()


# single SNP analysis
cat("####################")
cat("\n")
cat("Performing single SNP analysis")
cat("\n")
single_SNP <- mr_singlesnp(harmonised)
write.table(single_SNP, file=(paste("",opt$exp_name,"_",opt$out_id,"_SSNP.txt",sep="")), sep="\t",quote=F,row.names = F)
SSFP<-mr_forest_plot(single_SNP)
#ggsave(SSFP[1], file=(paste("",opt$exp_name,"_",opt$oc_name,"_SSNP.png",sep="")),dpi=300)
pdf(paste("",opt$exp_name,"_",opt$out_id,"_SSNP.pdf",sep=""))
SSFP[1]
dev.off()

## Scatter plot 
cat("####################")
cat("\n")
cat("Saving scatter plot")
cat("\n")
sp <- mr_scatter_plot(mr, harmonised)
#ggsave(sp[1], file=(paste("",opt$exp_name,"_",opt$oc_name,"_SP.png",sep="")),dpi=300)
pdf(paste("",opt$exp_name,"_",opt$out_id,"_SP.pdf",sep=""))
sp[1]
dev.off()

# directionality test
cat("####################")
cat("\n")
cat("Performing directionality test")
cat("\n")
dt<-directionality_test(harmonised)
write.table(dt, file=(paste("",opt$exp_name,"_",opt$out_id,"_DT.txt",sep="")), sep="\t",quote=F,row.names = F)