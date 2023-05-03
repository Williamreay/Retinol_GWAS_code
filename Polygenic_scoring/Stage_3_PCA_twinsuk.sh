#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## PCA processing - STAGE 3

for chr in $(seq 1 22);
  do
    /home/wrr239/c3202333/plink2 --bfile /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Chr${chr}_Stage_1_PCA\
    --extract /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Chr${chr}_LD_pruned_Stage_2_PCA.prune.in --make-bed \
    --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Stage_3_output/Chr${chr}_ready_for_merging_stage_3;
  done
