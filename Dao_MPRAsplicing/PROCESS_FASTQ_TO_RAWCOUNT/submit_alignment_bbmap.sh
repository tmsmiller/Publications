#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=4-00:00:00
#SBATCH --mem=4G

# script to use bbmap to align fastq files to reference sequence

for f in FASTQ_BY_reporters/*; do	
	reporter=$(echo ${f%.*} | cut -d '_' -f 4-5)
	ref_file=PTREseq_ref_seq/"$reporter".fa
	basename="$(basename $f .fastq)"
	out_sam=ALIGNMENT/"$basename".sam

	if [ -e $out_sam ]; then
		echo "$out_sam is done"
		continue
	else
		# limit submission to 50 jobs
		while [ $(squeue -u mkdao | wc -l) -gt 50 ]; do
			sleep 10s
		done
		# submit slurm jobs 
		sbatch alignment_command.sh $ref_file $f $out_sam
	fi
done
