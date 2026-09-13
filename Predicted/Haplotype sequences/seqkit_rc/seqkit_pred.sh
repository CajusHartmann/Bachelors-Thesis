#!/bin/bash

#SBATCH -J seqkit_reverse_complement_pred_repeated
#SBATCH --time=0-05:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=5
#SBATCH --chdir=/work/hartmaca/seqkit_rc_pred_repeated

#output files
#SBATCH -o /work/%u/job_logs/%x-%j.out
#SBATCH -e /work/%u/job_logs/%x-%j.err

#loading modules
module load SeqKit/2.2.0

#extract the miRNAs on the - strand and create reverse complement
tr -d '\r' < /work/hartmaca/bcftools_consensus_pred_repeated/coordinates_pred.bed | awk '$5=="-"{print $4}' | 
while read name; do
    echo "$name"
    seqkit seq -r -p "/work/hartmaca/bcftools_consensus_pred_repeated/${name}.fa" -o "${name}_rc.fa"
done