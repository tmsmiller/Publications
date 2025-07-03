# Script that use pre-trained spliceAI model to predict splice site probabilities
# Input: fasta file containing sequences of interest
# Output: 
#	text file, each row is a sequence in fasta file,  with 5 columns (sequence, highest donor position, highest donor probability, highest acceptor position, highest acceptor probability) 
#	(optional) spliceai prediction for each individual sequence with 4 columns (position, nucleotide, donor probability, acceptor probability)

from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
import pandas as pd
import os, sys
from Bio import SeqIO
import argparse

# Function to use spliceAI for splice prediction
# Input = sequence
# Output = 2 vectors, donor and acceptor probability for each position
def spliceai(input_sequence, models, context = 10000):
	x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
	y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

	acceptor_prob = y[0, :, 1]
	donor_prob = y[0, :, 2]
	
	if len(acceptor_prob) == 0 or len(donor_prob) == 0:
		print('Empty list')

	return donor_prob, acceptor_prob

def assign_site(df):
	# function to mark GT as donor site and AG as acceptor site
	for i in range(1, len(df) - 1):
		if df['seq'][i] == 'G' and df['seq'][i + 1] == 'T':
			df.loc[i - 1, 'label'] = 'donor'
		elif df['seq'][i] == 'A' and df['seq'][i + 1] == 'G':
			df.loc[i + 2, 'label'] = 'acceptor'
	return df

def get_best_sites(df):
	# retrieve the highest probability splice donor and splice acceptor
	best_donor = df.loc[df['donor']==df['donor'].max(), 'index'].values[0]
	donor_prob = df['donor'].max()
	
	df = df[df['index']>=best_donor]
	best_acceptor = df.loc[df['acceptor']==df['acceptor'].max(), 'index'].values[0]
	acceptor_prob = df['acceptor'].max()
	
	while best_donor > best_acceptor:
		if donor_prob > acceptor_prob:
			df = df[df['index']!=best_acceptor]
			best_acceptor = df.loc[df['acceptor']==df['acceptor'].max(), 'index'].values[0]
			acceptor_prob = df['acceptor'].max()
		elif donor_prob < acceptor_prob:
			df = df[df['index']!=best_donor]
			best_acceptor = df.loc[df['donor']==df['donor'].max(), 'index'].values[0]
			acceptor_prob = df['donor'].max()	

	return best_donor, donor_prob, best_acceptor, acceptor_prob

def parseArgs():
	parser = argparse.ArgumentParser(description='Predict donor and acceptor probabilities with spliceAI')
	parser.add_argument('fasta', type=str, help='Fasta file with sequences of interest')
	parser.add_argument('out', type=str, help='Output text files with best donor and acceptor for each sequence in input fasta')	
	parser.add_argument('--write_indiv', action='store_true', help='Write spliceAI output for individual sequences in fasta')
	args = parser.parse_args()

	return args

#-------------------------------------------------------------------------------------------------------#

if __name__== '__main__':

	# arguments	
	args = parseArgs()

	# retrieve pre-trained spliceAI model
	paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
	models = [load_model(resource_filename('spliceai', x)) for x in paths]

	# open output file
	output = open(args.out, 'w')
	output.write('seq\tbest_donor\tdonor_prob\tbest_acceptor\tacceptor_prob\n')

	# read fasta file
	reader = SeqIO.parse(args.fasta, 'fasta')

	for record in reader:
		plasmid = record.id
		input_sequence = str(record.seq)

		# predict donor and acceptor prob for each nucleotide in sequence
		donor_prob, acceptor_prob = spliceai(input_sequence, models)
		df = pd.DataFrame(data = {'seq': [i for i in input_sequence], 'donor': donor_prob, 'acceptor': acceptor_prob})
		
		# mark location GT and AG 
		df = assign_site(df.reset_index())

		# retrieve the highest probability donor and acceptor sites
		d, dp, a, ap = get_best_sites(df)
		output.write(f'{plasmid}\t{d}\t{dp}\t{a}\t{ap}\n')

		# write individual files if usage --write_indiv flag
		if args.write_indiv:
			df.to_csv(f'{plasmid}.spliceAI.txt', sep = '\t', index = False)
	
	output.close()
