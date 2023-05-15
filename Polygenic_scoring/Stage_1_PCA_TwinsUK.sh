#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## PCA processing - STAGE 1

for chr in $(seq 1 22);
  do
    /home/wrr239/c3202333/plink2 --bfile /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/0.8_rsq_PGS_calc/Chr${chr}_0.8_r2_HWE_MAF \
    --maf 0.05 --extract /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/Genotyped_twinsuk.txt \
    --exclude /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/long-range-ld.txt \
    --make-bed --keep-allele-order --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Chr${chr}_Stage_1_PCA;
  done
