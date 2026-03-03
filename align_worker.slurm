#!/bin/bash
#SBATCH --job-name=medium_queue_job
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --nodelist=rothhpc409

# Load conda
module load conda
source activate alignment-env

# Your SLURM job commands go here
echo "Running alignment for the sample..."
# alignment_command
