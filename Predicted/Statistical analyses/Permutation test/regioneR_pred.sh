#!/bin/bash
#SBATCH -J regioneR_pred
#SBATCH --time=0-05:00
#SBATCH --chdir=/work/hartmaca/regioneR_pred
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=15G

# output files
#SBATCH -o /work/%u/job_logs/%x-%A-%a.out
#SBATCH -e /work/%u/job_logs/%x-%A-%a.err

#loading modules
module load GCC/13.3.0
module load R/4.5.1

#command
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_LIBS_USER=/data/laub-rna/hiwi/cajus/thesis/R_packages
Rscript --vanilla regioneR_cluster_pred.R