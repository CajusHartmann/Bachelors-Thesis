#!/bin/bash

#SBATCH -J seqkit_reverse_complement
#SBATCH --time=0-05:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=5
#SBATCH --chdir=/work/hartmaca/seqkit_rc

#output files
#SBATCH -o /work/%u/job_logs/%x-%j.out
#SBATCH -e /work/%u/job_logs/%x-%j.err

#loading modules
module load SeqKit/2.2.0

#extract the miRNAs on the - strand and create reverse complement
tr -d '\r' < /work/hartmaca/bcftools_consensus/coordinates.bed | awk '$5=="-"{print $4}' | 
while read name; do
    echo "$name"
    seqkit seq -r -p "/work/hartmaca/bcftools_consensus/${name}.fa" -o "${name}_rc.fa"
done