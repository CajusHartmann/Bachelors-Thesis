#!/bin/bash
#SBATCH -J bcftools_filter_pred
#SBATCH --time=0-10:00
#SBATCH --chdir=/work/hartmaca/bcftools_pred
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
#command for filtering
bcftools view $VCF/1001genomes_snp-short-indel_only_ACGTN.vcf.gz -O v -o miRNA_positions_pred.vcf -R regions_pred.tsv