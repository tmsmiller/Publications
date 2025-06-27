1. bbmerge.sh in1=HELA-1_S23_L001_R2_001.fastq.gz in2=HELA-1_S23_L001_R1_001.fastq.gz out=HELA-1_merged.fastq

2. python3 demultiplex_barcode.py HELA-1_merged.fastq HELA-1_demux.txt

3. python3 split_fastq_by_reporter.py HELA-1_demux.txt HELA-1_merged.fastq

4. python3 split_report_by_reporters.py HELA-1_demux.txt

5. sbatch submit_alignment_bbmap.sh
   rm slurm*

6. bash submit_find_junction_from_alignment.sh
   rm slurm*

7. python3 check_artifacts.py HELA-1

8. cat REPORTS/HELA-1* > HELA-1_report.txt

9. python3 combine_report_and_check_artifacts.py HELA-1

10. python3 combine_report_and_check_artifacts.py HELA-1
