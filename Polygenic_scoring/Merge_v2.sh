#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

/home/wrr239/c3202333/plink2 --pmerge-list /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/To_merge.txt --merge-max-allele-ct 2 \
--pmerge-list-dir /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/0.3_rsq_assoc \
--make-pgen --out /home/wrr239/TwinsUK/DTR_Genotypes_1000G_V3P5_Collab_PublicID/0.3_rsq_assoc/Merged_all_chr_0.3 
