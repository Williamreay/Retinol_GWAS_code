####################################

## Retinol ~ protein colocalisation

## William Reay (2022)

#####################################

setwd("~/Desktop/Retinol_GWAS/ARIC_MR_and_coloc/")

library(dplyr)
library(data.table)
library(readxl)
library(coloc)

## Read in retinol GWAS

Ret_GWAS <- fread("../Meta_analyses/Common_var/IVW_METSIM_INTERVAL_sumstats_retinol_hg38.txt.gz")

## Square se to estimate per variant variance

Ret_GWAS$SNPVAR <- (Ret_GWAS$StdErr^2)
Ret_GWAS$ID <- Ret_GWAS$rsids

## List of proteins to test
## RBP4 - SeqId_15633_6
## INHBC - SeqId_15686_49
## GCKR - SeqId_5223_59
## DSG2 - SeqId_9484_75

## Define function to test colocalisation

Coloc_func_protein <- function(Protein_df) {
  Protein_df$SNPVAR <- Protein_df$SE^2
  Merged_sumstats <- merge(Ret_GWAS, Protein_df, by = "ID")
  pQTL_df <- Merged_sumstats %>% select(ID, BETA, SNPVAR.y)
  Retinol_df <- Merged_sumstats %>% select(ID, Effect, SNPVAR.x)
  ## Define coloc list input
  pQTL=list(beta=pQTL_df$BETA, varbeta=pQTL_df$SNPVAR.y, 
            sdY=1, type="quant")
  Retinol=list(beta=Retinol_df$Effect, varbeta=Retinol_df$SNPVAR.x, 
            sdY=1, type="quant")
  
  # Perform colocalisation under default priors
  coloc_df <- coloc.abf(pQTL, Retinol, MAF=NULL)
  
  return(coloc_df)
}

##RBP4
RBP4 <- fread("EA/SeqId_15633_6.PHENO1.glm.linear")

Coloc_RBP4 <- Coloc_func_protein(RBP4)

sensitivity(Coloc_RBP4, rule = "H4 > 0.8")

##INHBC
INHBC <- fread("EA/SeqId_15686_49.PHENO1.glm.linear")

Coloc_INHBC <- Coloc_func_protein(INHBC)

sensitivity(Coloc_INHBC, rule = "H4 > 0.8")

##GCKR
GCKR <- fread("EA/SeqId_5223_59.PHENO1.glm.linear")

Coloc_GCKR <- Coloc_func_protein(GCKR)

sensitivity(Coloc_GCKR, rule = "H4 > 0.8")

##DSG2
DSG2 <- fread("EA/SeqId_9484_75.PHENO1.glm.linear")

Coloc_DSG2 <- Coloc_func_protein(DSG2)

sensitivity(Coloc_DSG2, rule = "H3 > 0.8")
