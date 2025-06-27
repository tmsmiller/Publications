#!/usr/bin/bash
#SBATCH -n 4
#SBATCH --mem=8GB
#SBATCH -t 1-00:00:00

# command to run bbmap for alignment
~/softwares/bbmap/bbmap.sh ref=$1 in=$2 out=$3 nodisk

# covert sam file to bam file
samtools view -S -b $3 > "$3".bam
# delete sam file to save space
#rm $3
