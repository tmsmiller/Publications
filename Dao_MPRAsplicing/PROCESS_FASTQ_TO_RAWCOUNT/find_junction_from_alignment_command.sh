#!/usr/bin/bash
#SBATCH -n 2
#SBATCH --time=4-00:00:00
#SBATCH --mem=4G

report=$1
sample=$2
python3 find_junction_from_alignment.py $report $sample
