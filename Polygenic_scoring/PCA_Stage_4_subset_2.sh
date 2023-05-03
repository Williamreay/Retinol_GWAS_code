#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## PCA processing - STAGE 4: merging and extracting each half of the cohort - half one

/home/wrr239/c3202333/plink2 \
    --pmerge-list /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/Input_for_merge.txt bfile --merge-max-allele-ct 2  \
    --keep /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/Subset2_IDs.txt --make-bed \
    --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/PCA_calculations_in_either_TwinSubset/Merged_stage_4/Subset_2_all_chr_PCA_input
