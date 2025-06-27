import argparse
from Bio import SeqIO
from os.path import basename, splitext
from Bio.Seq import reverse_complement
import edlib

#######################################

def match(fastq_file, barcode_handle, file_format='fastq', forward = 'GCTGCCCGACAACCAC'):
	
	# read fastq file
	reader = SeqIO.index(fastq_file, file_format)

	read = []	
	read_len = []
	read_type = []
	reporters = []

	for record in reader: # iterate through all reads
		seq = str(reader[record].seq)  
		read.append(reader[record].id)
		read_len.append(len(seq))	
				
		forward_alignment = edlib.align(forward, seq[:50], mode = 'HW') # check if forward primer is present
		if forward_alignment['editDistance'] <= 2:  
			if 'CGGTCGGATGCATCCG' in seq: # check for palindromic inserttion
				read_type.append('RT_error') # I previously called palindromic insertation artifact as RT_error
				reporters.append('n/a')
			else:			
				dist_0 = []
				dist_1 = []

				for reporter, _barcode in barcodes.items():
					barcode = _barcode + 'GGTAC' # align barcode plus 4 downstream flanking nucleotide
					alignment = edlib.align(barcode, seq, mode = 'HW')

					# allow reads to be assigned to multiple barcodes
					if alignment['editDistance'] == 0: 
						dist_0.append(reporter)
					elif alignment['editDistance'] == 1: 
						dist_1.append(reporter)				
		
				if len(dist_0) == 0 and len(dist_1) == 0:
					best_read_type = 'cannot_find_barcode'
					best_reporter = 'n/a'
				elif len(dist_0) > 0: 
					# if found barcodes with perfect match, then only assign reads to those barcodes
					best_read_type = 'perfect_match'
					best_reporter = ','.join(dist_0)
				elif len(dist_1) > 0:
					# else if found barcodes with 1 mismatch, then assign reads to those barcodes
					best_read_type = 'one_mismatch'
					best_reporter = ','.join(dist_1)
				else:
					# deal with edge cases
					best_read_type = 'error'
					best_report = 'error'

				reporters.append(best_reporter)
				read_type.append(best_read_type)
		else:
			# for reads without the forward primer present (some samples were pooled with other experiments)
			read_type.append('missing_forward_primer')
			reporters.append('n/a')

	return read, read_type, read_len, reporters

################################################################################################

def parseArgs():
	parser = argparse.ArgumentParser(description='Demultiplex barcodes in fastq files. Input = merged fastq file and. Output = report text file.')
	parser.add_argument('merged', type=str, help='Merged fastq file')
	parser.add_argument('outfile', type=str, help='Output file')
	args = parser.parse_args()

	return args

#################################################################################################

if __name__== '__main__':	
	args = parseArgs()
	# read barcode seq
	barcodes = {}
	with open('data/barcode_indexed.txt', 'r') as f:
		for line in f:
			spl = line.rstrip().split('\t')
			barcodes[spl[0]] = spl[1]
	
	# assign reads to reporters by finding matched barcodes (allow 1 mismatch)
	read, read_type, read_len, reporters = match(args.merged, barcodes)

	# write a report file for all reads in the fastq file
	with open(args.outfile, 'w') as out:
		for name, l, t, re in zip(read, read_len, read_type, reporters):
			out.write(f'{name}\t{l}\t{t}\t{re}\n')


	



























