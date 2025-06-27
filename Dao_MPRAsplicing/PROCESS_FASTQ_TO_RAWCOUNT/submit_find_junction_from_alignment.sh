#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=7-00:00:00
#SBATCH --mem=2G

	
for fn in REPORTS/*txt; do
	
	while [ $(squeue -u mkdao | wc -l) -gt 61 ]; do
		sleep 5s
	done

	f="${fn##*/}"
	sample=$(echo $f | cut -d '_' -f 1 | cut -d '-' -f 1-2)
	JOB_ID=$(sbatch --parsable find_junction_from_alignment_command.sh "$f" "$sample")
	echo "Submitting job $JOB_ID for $fn"
done
