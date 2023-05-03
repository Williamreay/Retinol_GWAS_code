#!/bin/bash
#SBATCH --cpus-per-task=3
#SBATCH --time=90:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## FIZI Retinol sumstats imputation

eval "$(conda shell.bash hook)"
conda activate fizi_env

cd /home/control/data/users/william/Retinol_GWAS_summ_stats_imputation/FIZI_imputation

for chr in $(seq 1 22);
  do
    fizi impute FIZI_munged_ATBC_PLCO_retinol.sumstats.gz \
    ../../1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --min-prop 0.01 --output New_1000G_panel/Chr_${chr}_FIZI_default_retinol_imputed_ATBC_PLCO;
  done
