#!/bin/bash
#SBATCH -J starseeker
#SBATCH --time=0-12:00
#SBATCH --chdir=/work/hartmaca/starseeker
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=cajus.hartmann@gmx.net
#SBATCH --mail-type=BEGIN,END,FAIL

#output files
#SBATCH -o /work/%u/job_logs/%x-%j.out
#SBATCH -e /work/%u/job_logs/%x-%j.err

#loading modules
module load Miniforge3/24.3.0-0
source activate /data/laub-rna/hiwi/cajus/conda/starseeker

python StarSeeker/Starseeker.py -m mature_filtered.fa -p hairpin_filtered.fa -l