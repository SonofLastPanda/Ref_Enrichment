#!/bin/bash

#SBATCH --job-name=Vcf_operations
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -o ../slurm_out/%j.out
#SBATCH --mail-user=erkin.alacamli@ut.ee
#SBATCH --mail-type=ALL

module load snakemake
#module load gatk

#snakemake -s Vcf_Operations.smk --profile config_dir/ --use-envmodules --forcerun IndStats
snakemake -s Vcf_Operations.smk  --profile config_dir/ --use-envmodules --forceall