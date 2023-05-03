#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## First tranche of INFO filtering 
## MAF, missing, and HWE

for chr in $(seq 1 23);
  do
    /home/wrr239/c3202333/plink2  \
    --pfile /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/Raw_pgen/Raw_chr${chr}_TwinsUK \
    --mach-r2-filter 0.8 \
    --maf 0.01 --hwe 0.0000000001 \
    --make-bed --keep-allele-order --allow-extra-chr \
    --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/0.8_rsq_PGS_calc/Chr${chr}_0.8_r2_HWE_MAF;
  done
