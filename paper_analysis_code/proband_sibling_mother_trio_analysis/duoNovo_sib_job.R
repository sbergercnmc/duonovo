#!/bin/bash
#SBATCH --job-name=duoNovo_proband_sibling_mother
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-100

module load r

Rscript run_duoNovo_sib.R ${SLURM_ARRAY_TASK_ID}