#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=36:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## PWAS - METSIM+INTERVAL TWAS (Chatterjee - Blood)

for chr in $(seq 1 22);
do
	Rscript ~/data/users/william/PWAS_weights/Blood/scripts/PWAS.assoc_test.R \
	--sumstats  ~/data/users/william/Meta_retinol_GWAS/Meta_munged/FULL_munged_Stouffer_ALL_THREE.sumstats.gz \
	--weights ~/data/users/william/PWAS_weights/Blood/PWAS_EA/Plasma_Protein_EA_hg38.pos \
	--weights_dir ~/data/users/william/PWAS_weights/Blood/PWAS_EA/Plasma_Protein_weights_EA_enet \
	--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
	--coloc_P 0.05 --PANELN ~/data/users/william/PWAS_weights/Blood/PWAS_EA/Plasma_Protein_EA_hg38.pos \
	--GWASN 22274 \
	--chr $chr --out ~/data/users/william/Meta_retinol_GWAS/TWAS_PWAS/Blood_PWAS/Full_meta_all_three_blood_PWAS_chr${chr};
done
