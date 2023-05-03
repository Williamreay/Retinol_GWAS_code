#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## KING relatedness estimates in each subcohort

for i in $(seq 1 2);
  do
    /home/wrr239/c3202333/plink2 \
    --bfile /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Merged_stage_4/Subset_${i}_all_chr_PCA_input \
    --make-king-table --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Merged_stage_4/KING_relatedness_Subset_${i}_all_chr_PCA_input;
  done
