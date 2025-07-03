# Script that use pre-trained Pangolin model to predict splice site probabilities
# Input: fasta file containing sequences of interest
# Output: 
#	text file, each row is a sequence in fasta file,  with 5 columns (sequence, highest donor position, highest donor probability, highest acceptor position, highest acceptor probability) 
#	(optional) spliceai prediction for each individual sequence with 4 columns (position, nucleotide, donor probability, acceptor probability)

from pkg_resources import resource_filename
from pangolin.model import *
import sys
import pandas as pd
from Bio import SeqIO
import argparse


# Change this to the desired models. The model that each number corresponds to is listed below.
model_nums = [2, 3]
# 0 = Heart, P(splice)
# 1 = Heart, usage
# 2 = Liver, P(splice)
# 3 = Liver, usage
# 4 = Brain, P(splice)
# 5 = Brain, usage
# 6 = Testis, P(splice)
# 7 = Testis, usage

# Change this to the desired sequences and strand for each sequence. If the sequence is N bases long, Pangolin will
# return scores for the middle N-10000 bases (so if you are interested in the score for a single site, the input should
# be: 5000 bases before the site, base at the site, 5000 bases after the site). Sequences < 10001 bases can be padded with 'N'.
def add_context(seq, context = 10001):
    added_seq = 'N'*(context//2) + seq + 'N'*(context//2)
    return added_seq

def one_hot_encode(seq, strand):
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')
    if strand == '+':
        seq = np.asarray(list(map(int, list(seq))))
    elif strand == '-':
        seq = np.asarray(list(map(int, list(seq[::-1]))))
        seq = (5 - seq) % 5  # Reverse complement
    return IN_MAP[seq.astype('int8')]

def assign_site(df):
	for i in range(1, len(df) - 1):
		if df['Seq'][i] == 'G' and df['Seq'][i + 1] == 'T':
			df.loc[i - 1, 'label'] = 'donor'
		elif df['Seq'][i] == 'A' and df['Seq'][i + 1] == 'G':
			df.loc[i + 2, 'label'] = 'acceptor'
	return df

def get_best_sites(df):
	donor_candidates = df[df['label']=='donor']
	acceptor_candidates = df[df['label']=='acceptor']

	best_donor = None
	best_acceptor = None
	max_donor_prob = 0
	max_acceptor_prob = 0

	for i, donor_row in donor_candidates.iterrows():
		downstream_acceptors = acceptor_candidates[acceptor_candidates['Position'] > donor_row['Position']]
		if downstream_acceptors.empty:
			continue
		else:
			best_downstream_acceptor = downstream_acceptors.loc[downstream_acceptors['liver_prob'].idxmax()]	
		
		if donor_row['liver_prob'] > max_donor_prob or best_downstream_acceptor['liver_prob'] > max_acceptor_prob:
			best_donor = donor_row['Position']
			best_acceptor = best_downstream_acceptor['Position']
			max_donor_prob = donor_row['liver_prob']
			max_acceptor_prob = best_downstream_acceptor['liver_prob']
		
	return best_donor, max_donor_prob, best_acceptor, max_acceptor_prob

def parseArgs():
	parser = argparse.ArgumentParser(description='Predict donor and acceptor probabilities with spliceAI')
	parser.add_argument('fasta', type=str, help='Fasta file with sequences of interest')
	parser.add_argument('out', type=str, help='Output text files with best donor and acceptor for each sequence in input fasta')	
	parser.add_argument('--write_indiv', action='store_true', help='Write Pangolin output for individual sequences in fasta')
	args = parser.parse_args()

	return args

#-------------------------------------------------------------------------------------------------------#

if __name__== '__main__':

	# arguments	
	args = parseArgs()
	
	# open output file
	output = open(args.out, 'w')
	output.write('seq\tbest_donor\tdonor_prob\tbest_acceptor\tacceptor_prob\n')

	# load models
	models = []
	for i in model_nums:
		for j in range(1, 6):
			model = Pangolin(L, W, AR)
			if torch.cuda.is_available():
			    model.cuda()
			    weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)))
			else:
			    weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)),
						 map_location=torch.device('cpu'))
			model.load_state_dict(weights)
			model.eval()
			models.append(model)

	# Get scores
	IN_MAP = np.asarray([[0, 0, 0, 0],
			     [1, 0, 0, 0],
			     [0, 1, 0, 0],
			     [0, 0, 1, 0],
			     [0, 0, 0, 1]])
	INDEX_MAP = {0:1, 1:2, 2:4, 3:5, 4:7, 5:8, 6:10, 7:11}
	
	# read fasta file
	reader = SeqIO.parse(args.fasta, 'fasta')

	for record in reader:

		plasmid = record.id
		seq = str(record.seq)
		df = pd.DataFrame({'Position': [i for i, n in enumerate(seq)], 'Seq': [i for i in seq]})

		seq_w_context = add_context(seq)
		seq_one_hot_encode = one_hot_encode(seq_w_context, '+').T
		seq_one_hot_encode = torch.from_numpy(np.expand_dims(seq_one_hot_encode, axis=0)).float()

		if torch.cuda.is_available():
			seq_one_hot_encode = seq_one_hot_encode.to(torch.device("cuda"))

		for j, (model_name, model_num) in enumerate(zip(['liver_prob', 'liver_usage'], model_nums)):
			score = []
			# Average across 5 models
			for model in models[5*j:5*j+5]:
				with torch.no_grad():
					score.append(model(seq_one_hot_encode)[0][INDEX_MAP[model_num],:].cpu().numpy())
		
			mean = np.mean(score, axis=0)
			df[model_name] = mean

		# mark location GT and AG 
		df = assign_site(df.reset_index(drop = True))

		# retrieve the highest probability donor and acceptor sites
		d, dp, a, ap = get_best_sites(df)
		output.write(f'{plasmid}\t{d}\t{dp}\t{a}\t{ap}\n')

		# write individual files if usage --write_indiv flag
		if args.write_indiv:
			df.to_csv(f'{plasmid}.Pangolin.txt', sep = '\t', index = False)
	
	output.close()
