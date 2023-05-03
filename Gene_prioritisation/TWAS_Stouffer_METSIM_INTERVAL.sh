#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=36:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## TWAS - METSIM+INTERVAL TWAS (GTEx v8 weights)

while read TISSUE CHROMOSOME;
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats  ~/data/users/william/Meta_retinol_GWAS/Meta_munged/FULL_munged_Stouffer_METSIM_INTERVAL.sumstats.gz \
--weights ~/data/users/william/Meta_retinol_GWAS/TWAS_PWAS/GTEx_v8_weights/GTExv8.ALL.${TISSUE}.pos \
--weights_dir ~/data/users/william/Meta_retinol_GWAS/TWAS_PWAS/GTEx_v8_weights/ \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--coloc_P 0.05 --PANELN ~/data/users/william/Meta_retinol_GWAS/TWAS_PWAS/GTEx_v8_weights/GTExv8.ALL.${TISSUE}.pos \
--GWASN 17268 \
--chr $CHROMOSOME --out ~/data/users/william/Meta_retinol_GWAS/TWAS_PWAS/TWAS_out_Stouffer_METSIM_INTERVAL/METSIM_INTERVAL_Stouffer_retinol_${TISSUE}_${CHROMOSOME};
done < /home/control/data/users/william/Meta_retinol_GWAS/TWAS_PWAS/Retinol_METSIM_INTERVAL_input.txt
