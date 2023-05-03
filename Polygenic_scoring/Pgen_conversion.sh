#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## Convert TwinsUK vcf into pgen

for chr in $(seq 1 23);
  do
    /home/wrr239/c3202333/plink2  \
    --vcf /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/1000G_P3V5_PID_chr${chr}.dose.vcf.gz dosage=DS \
    --make-pgen --keep-allele-order --allow-extra-chr \
    --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/Raw_pgen/Raw_chr${chr}_TwinsUK;
  done
