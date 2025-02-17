#!/bin/bash
#SBATCH --job-name=process_duoNovo_outputs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-42

module load r

Rscript process_duoNovo_output.R ${SLURM_ARRAY_TASK_ID}
