#!/bin/bash
#SBATCH -J bcftools_consensus
#SBATCH --time=0-02:00
#SBATCH --chdir=/work/hartmaca/bcftools_consensus
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20

#output files
#SBATCH -o /work/%u/job_logs/%x-%A-%a.out
#SBATCH -e /work/%u/job_logs/%x-%A-%a.err

#loading modules
module load Miniforge3/24.3.0-0
source activate /data/laub-rna/hiwi/cajus/conda/bcftools
module load GCC/12.2.0
module load SAMtools/1.17

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#show versions
bcftools --version
samtools --version

VCF=/data/laub-rna/hiwi/cajus/thesis/data
GENOME=/data/laub-rna/hiwi/cajus/arabid_genome
#read current line from bed file
read chr start end name strand <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" coordinates.bed | tr -d '\r')
echo "Processing $name ($chr:$start-$end:$strand)"
#modify start because bed files are 0 based but samtools wants 1 based
new_start=$((start+1))
outfile="${name}.fa"
#iterate over all accessions
for sample in $(cat accessions.txt); do
    #fasta header
    echo ">$sample" >> "$outfile"

    samtools faidx $GENOME/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa "$chr:$new_start-$end" | bcftools consensus -s "$sample" $VCF/miRNA_positions_filtered_MAF0.01.vcf.gz | grep -v "^>" >> "$outfile"
done