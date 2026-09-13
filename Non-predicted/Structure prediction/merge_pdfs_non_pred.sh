#!/bin/bash
#SBATCH -J merge_pdfs_non_pred
#SBATCH --time=0-01:00
#SBATCH --chdir=/work/hartmaca/RNAfold_non_pred/merged_pdfs
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=15G

# output files
#SBATCH -o /work/%u/job_logs/%x-%A-%a.out
#SBATCH -e /work/%u/job_logs/%x-%A-%a.err

#loading modules
module load GCC/13.3.0
module load R/4.5.1
module load ImageMagick/7.1.1-38
module load poppler/25.07.0

#command
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_LIBS_USER=/data/laub-rna/hiwi/cajus/thesis/R_packages
Rscript --vanilla merge_pdfs_non_pred.R