#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## PCA processing - STAGE 2

for chr in $(seq 1 22);
  do
    /home/wrr239/c3202333/plink2 --bfile /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Chr${chr}_Stage_1_PCA \
    --keep-allele-order --indep-pairwise 1000 50 0.05 --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Chr${chr}_LD_pruned_Stage_2_PCA;
  done
