#############################
## Strength of retinol IVs (F-statistic and I^2 stat)

## William Reay (2023)
#############################

library(data.table)
library(TwoSampleMR)
library(dplyr)

ivs <- fread("~/Desktop/Retinol_GWAS/MR_all_cis_acting_SNPs/Clumped_retinol_IVs.txt")

## F statistic

F_statistic <- function(r,n,k){r*(n-1-k)/(1-r)/k}

Approximated_r2 <- sum(get_r_from_bsen(ivs$beta.exposure, ivs$se.exposure,17268)^2)

F_statistic(Approximated_r2, 17268, 8)

## F statistic = 47.07724 all IVS

## RBP4 on its own

RBP4 <- ivs %>% filter(SNP == "rs10882283")

RBP4_Approximated_r2 <- sum(get_r_from_bsen(RBP4$beta.exposure, RBP4$se.exposure,17268)^2)

F_statistic(RBP4_Approximated_r2, 17268, 1)

## 95.72243

## I^2

Isq = function(y,s){
  k = length(y)
  w = 1/s^2; sum.w = sum(w)
  mu.hat = sum(y*w)/sum.w
  Q = sum(w*(y-mu.hat)^2)
  Isq = (Q - (k-1))/Q
  Isq = max(0,Isq)
  return(Isq)
}

All_IVs_Isq <-Isq(ivs$beta.exposure, 
                         ivs$se.exposure)

Isq(RBP4$beta.exposure, RBP4$se.exposure)
