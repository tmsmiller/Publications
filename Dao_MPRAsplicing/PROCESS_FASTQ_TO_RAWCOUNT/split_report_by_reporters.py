# scripts to split report text file into separate files for each barcode
# e.g. for HELA-1 sample:
#	input =  HELA-1_report.txt 
#	output = HELA-1_ABAB_1_report.txt, HELA-1_ABAB_2_report.txt, HELA-1_ABAB_3_report.txt... (6500 files corresponding to 6500 barcodes)

import os, sys
import pandas as pd
from Bio import SeqIO

def open_and_write(fn, mode, line):
	# function to open and write lines
	with open(fn, mode) as out:
		out.write(line)

reporter_ls = []

demux = sys.argv[1]
sample = demux.replace('_demux.txt', '')

with open(demux, 'r') as fn:
	for line in fn:
		spl = line.rstrip().split('\t')
		read_type = spl[2]
		reporter = spl[3]
		if read_type in ['cannot_find_barcode', 'RT_error', 'missing_forward_primer']:
			# write a separate report file for these 3 categories
			if 'error' not in reporter_ls:
				open_and_write(f'REPORTS/{sample}_error_report.txt', 'w', line)
				reporter_ls.append('error')
			else:
				open_and_write(f'REPORTS/{sample}_error_report.txt', 'a', line)
		else:
			# write a separate report file for each barcode
			reporters = reporter.split(',')
			for reporter in reporters:
				if reporter not in reporter_ls:	
					open_and_write(f'REPORTS/{sample}_{reporter}_report.txt', 'w', line)
					reporter_ls.append(reporter)
				else:
					open_and_write(f'REPORTS/{sample}_{reporter}_report.txt', 'a', line)
