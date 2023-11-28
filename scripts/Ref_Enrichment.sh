#!/bin/bash

#SBATCH --job-name=Ref_Enrich_21
#SBATCH --time=48:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -o ../slurm_out/%j.out
#SBATCH --mail-user=erkin.alacamli@ut.ee
#SBATCH --mail-type=ALL

module load snakemake
#module load gatk

snakemake -s Ref_Enrichment.smk --profile config_dir/ --use-envmodules 