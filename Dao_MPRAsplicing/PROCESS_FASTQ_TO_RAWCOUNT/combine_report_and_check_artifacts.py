# script to combine reports with check artifact output
# input:
#	name of sample (e.g. HELA-1)
#	must have the report file (e.g. HELA-1_report.txt) and artifact file (e.g. HELA-1_artifacts.txt) in the same directory
# output: raw count file (e.g. HELA-1_raw_count.txt)

import pandas as pd
import numpy as np
import sys


def get_correct_category(category, dist_RE_not_cloned, dist_from_cannot_find_barcode, nospacer_dist, dist_pum_chimera):
    if 'missing_forward_primer' in category or 'RT_error' in category:
        return category
    elif dist_RE_not_cloned == 0:
        return 'RE_not_cloned'
    elif 'cannot_find_barcode' in category:
        # align the sequence immediately upstream of the blank
        # if cannot find this sequence (dist >= 3), then blank mispriming
        if int(dist_from_cannot_find_barcode) >= 3 and int(dist_from_cannot_find_barcode) != -1:
            return 'blank_mispriming'
        else:
            return 'unidentifiable_barcode'
    elif int(nospacer_dist) <= 3 and int(nospacer_dist) != -1:
        return 'spacer_not_cloned'
    elif 'full_length' in category:
        return 'full_length'
    elif 'AG_junction_low_mismatch' in category:
        junction = category.split(',')[0].split('-')[1]
        return f'spliced-{junction}'
    elif 'AG_junction_high_mismatch' in category: # should the high mismatched stuffs be treated separate?
        #junction = category.split('-')[1]
        return category
    elif 'ambiguous' in category:
        # align the sequence immediately upstream of the spacer that chimerize with PUM
        # if cannot find this sequence (dist >= 9), then PUM-spacer chimera
        if int(dist_pum_chimera) > 4 and int(dist_pum_chimera) != -1 and int(dist_pum_chimera) != 9:
        #if dist_pum_chimera < 9 and dist_pum_chimera != -1:
            return 'PRE-SPACER_mispriming'
        else:
            return 'ambiguous_deletion'
    else:
        return 'PCR_chimera'
    
    
def get_reporter(reporters, distances, best_category):
    if best_category == 'full_length' or 'spliced' in best_category:
        reporters = reporters.split(',')
        distances = distances.split(',')
        if len(reporters) == 1:
            return reporters[0]
        else:
            lowest_dist = distances.index(min(distances))
            return reporters[lowest_dist]
    else:
        return 'n/a'

def process_report(sample):

	# read report
	df = pd.read_csv(
	    f'{sample}_report.txt', 
	    sep = '\t', 
	    names = ['read_name', 'length', 'read_type', 'reporters', 'cigar', 'distance', 'category', 'reporter']
	)
	# convert edit distance column to numeric data type and fill missing cells with -1
	df['distance'] = pd.to_numeric(df['distance'], errors = 'coerce').fillna(-1).astype(float)
	# in cases where sequencing reads are assigned to multiple reporters (see demultiplex_barcode.py script), 
	# only choose the reporter with the lowest edit distance
	df = df.loc[df.groupby('read_name')['distance'].idxmin()]

	# read artifact file
	artifacts = pd.read_csv(f'{sample}_artifacts.txt', sep = '\t')
	artifacts = artifacts.drop_duplicates()
	
	# merge report and artifact file
	df = df.merge(artifacts, on = 'read_name', how = 'outer')
	df = df.fillna(-1)
	
	# get best category
	df['best_category'] = df.apply(lambda g: get_correct_category(g['category'], g['dist_no_RE'], g['dist_barcode'], g['dist_no_spacer'], g['dist_pum_chim']), axis = 1)

	#df.to_csv(f'{sample}_with_best_category.txt', sep = '\t', index = False)
	df = df[df['best_category']!='missing_forward_primer']

	# group by barcode and count number of reads for each category
	count = df.groupby('reporter')['best_category'].value_counts().to_frame()
	count = count.rename(columns={'best_category': 'count'}).reset_index()
	count.to_csv(f'{sample}_raw_count.txt', sep = '\t', index = False)
	

sample = sys.argv[1]
process_report(sample)
