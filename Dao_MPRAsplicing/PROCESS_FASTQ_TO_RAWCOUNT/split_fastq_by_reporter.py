# scripts to split report text file into separate files for each barcode
# e.g. for HELA-1 sample:
#	input = HELA-1_merged.fastq and HELA-1_report.txt 
#	output = HELA-1_ABAB_1.fastq, HELA-1_ABAB_2.fastq, HELA-1_ABAB_3.fastq... (6500 files corresponding to 6500 barcodes)

import os, sys
import pandas as pd
from Bio import SeqIO

def open_and_write(fn, mode, line):
	with open(fn, mode) as out:
		out.write(line)

reporter_ls = []

demux = sys.argv[1]
fastq = sys.argv[2]

sample = demux.replace('_demux.txt', '')

# read report as a dataframe
df = pd.read_csv(demux, sep = '\t', names = ['read', 'length', 'type', 'reporter']) 
# read corresponding fastq file
reader = SeqIO.parse(fastq, 'fastq') 

for idx, record in enumerate(reader):
	read_report = df.iloc[idx]
	read_type = read_report['type']
	if read_type != 'cannot_find_barcode' and read_type != 'RT_error' and read_type != 'missing_forward_primer':
		# write a separate fastq file for the 3 types of reads above
		reporter_ls = [i for i in read_report['reporter'].split(',')]
		for reporter in reporter_ls:
			with open(f'FASTQ_BY_reporters/{sample}_{reporter}_demux.fastq', 'a') as handle:
				SeqIO.write(record, handle, 'fastq')
	else:
		# write a separate fastq file for each barcode
		with open(f'FASTQ_BY_reporters/{sample}_error_demux.fastq', 'a') as handle:
			SeqIO.write(record, handle, 'fastq')


