#!/bin/bash
#SBATCH -J RNAfold_predicted_repeated
#SBATCH --time=0-01:00
#SBATCH --chdir=/work/hartmaca/RNAfold_predicted_repeated
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G

#output files
#SBATCH -o /work/%u/job_logs/%x-%A-%a.out
#SBATCH -e /work/%u/job_logs/%x-%A-%a.err

#loading modules
module load Miniforge3/24.3.0-0
module load GCCcore/15.3.0
module load ImageMagick/7.1.2-28
source activate /data/laub-rna/hiwi/cajus/conda/viennarna

DATA=/data/laub-rna/hiwi/cajus/thesis/Haplotpye_Sequences/predicted_repeated
REL=/data/laub-rna/hiwi/cajus/thesis/Haplotpye_Sequences/relplot

#read current line from bed file
read chr start end name strand <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" $DATA/bcftools_consensus_pred_repeated/coordinates_pred.bed | tr -d '\r')
echo "Processing $name ($chr:$start-$end:$strand)"

#show version
RNAfold --version
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir $name
cd $name

#make headers unique (necessary because two haplotypes that have the frequency one have the same header (>1x))
awk 'NR==FNR {sub(/\r$/,""); if (/^>/) count[$0]++; next} {sub(/\r$/,""); if (/^>/) {if (count[$0] > 1) {seen[$0]++; print $0 "_" seen[$0]} else {print $0}; next} {print}}' \
$DATA/Results/$name/frequencies.fa $DATA/Results/$name/frequencies.fa | \
RNAfold -p --noLP > out.txt

for ss in *_ss.ps; do
    freq=${ss%_ss.ps} #extract frequency for naming
    $REL/relplot_nolegend.pl ${freq}_ss.ps ${freq}_dp.ps > ${freq}_rss.ps #add positional entropy to secondary structure
    magick -density 300 ${freq}_rss.ps ${freq}_rss.pdf #create pdf of ss
done