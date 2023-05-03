#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --time=13:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

./ldsc.py --h2-cts ~/data/users/william/Meta_retinol_GWAS/Meta_munged/Hm3_munged_Stouffer_INTERVAL_METSIM.sumstats.gz \
--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
--out ~/data/users/william/Meta_retinol_GWAS/Part_h2/Stouffer_INTERVAL_METSIM_partitioned_h2 \
--ref-ld-chr-cts Multi_tissue_gene_expr.ldcts --w-ld-chr weights_hm3_no_hla/weights.
