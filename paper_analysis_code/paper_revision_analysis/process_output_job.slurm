#!/bin/bash
#SBATCH --job-name=process_duoNovo_outputs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-100

DIR_FILE="trio_directories.txt"  
INDEX=$SLURM_ARRAY_TASK_ID
module load r



RUNDIR=$( head -n $INDEX | tail -n 1)

Rscript process_duoNovo_output.R $RUNDIR
Rscript duonovo_performance_metrics.R $RUNDIR
Rscript de_novo_variant_granges.R $RUNDIR
