#!/bin/bash
#SBATCH -J VIVA_predicted_repeated
#SBATCH --time=0-05:00
#SBATCH --chdir=/work/hartmaca/VIVA_predicted_repeated
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G

#output files
#SBATCH -o /work/%u/job_logs/%x-%A-%a.out
#SBATCH -e /work/%u/job_logs/%x-%A-%a.err

#loading modules
module load Julia/1.4.1-linux-x86_64

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PATH=/data/laub-rna/hiwi/cajus/thesis/bin:$PATH
export JULIA_DEPOT_PATH=/data/laub-rna/hiwi/cajus/thesis/julia_package

#show versions
echo "Versions:"
julia --version
julia -e 'using Pkg; Pkg.status("VariantVisualization")'

FILE="predictions_names.csv"

NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE" | tr -d '\r')
echo "This job belongs to $NAME"

#command
VCF=/work/hartmaca/bcftools_pred
CSV=/data/laub-rna/hiwi/cajus/thesis/data/csv_all_predicted

viva -f "$VCF/miRNA_positions_pred.vcf" -l "$CSV/$NAME.csv" \
     -s pdf -o "/work/hartmaca/VIVA_predicted_repeated/$NAME" --heatmap genotype,read_depth -n \
     --avg_dp sample,variant