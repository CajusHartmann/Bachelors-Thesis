#!/bin/bash
#SBATCH -J bcftools_view_MAF_pred
#SBATCH --time=0-10:00
#SBATCH --chdir=/work/hartmaca/bcftools_view_MAF_pred
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20

#output files
#SBATCH -o /work/%u/job_logs/%x-%j.out
#SBATCH -e /work/%u/job_logs/%x-%j.err

#loading modules
module load Miniforge3/24.3.0-0
source activate /data/laub-rna/hiwi/cajus/conda/bcftools

#path
VCF=/data/laub-rna/hiwi/cajus/thesis/data

#show version
bcftools --version

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
bcftools view -i 'MAF >= 0.01' $VCF/miRNA_positions_pred.vcf -O v -o miRNA_positions_pred_filtered_MAF0.01.vcf