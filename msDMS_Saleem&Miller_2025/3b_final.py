import random
import numpy as np
from ReactivityProfile import ReactivityProfile as rprofile
import RNAStructureObjects as RNAtools

import matplotlib.pyplot as plot
import scipy.stats

import sys, os, math
import argparse
plot.rcParams['pdf.fonttype'] = 42
plot.rcParams['font.sans-serif'] = 'Arial'

parser = argparse.ArgumentParser()
parser.add_argument('--real', required=True, help='path to profile')
parser.add_argument('--random', required=True)
parser.add_argument('--output', required=True, help='output prefix')
parser.add_argument('--random2', required=True, help='2nd random dist')
parser.add_argument('--real2', required=True, help='2nd real dist')
args = parser.parse_args()


def contact_filter(protections,CT1, filter_distance=0):
    
    num_within_range = {}
    for i in protections:
        count = 0
        for j in protections:
            if i == j:
                continue
            if CT1.contactDistance(i, j) > filter_distance:
                
                continue
            else:
                count += 1
        num_within_range[i] = count
    return num_within_range

def parse_infile(input_f):
    nts = {}
    with open(input_f) as inp:
        for i,line in enumerate(inp):
            if i == 0:
                continue
            else:
                spl = line.split()
                if spl[2] == '0.0':
                    continue
                else:
                    nts[int(spl[0])] = int(spl[1])
    return nts

def violin_plot(name, labels, data, p):
    quants = []
    for i in data:
        quants.append([.25,.75])

    fig, ax = plot.subplots(nrows = 1, ncols = 1)
    ax.violinplot(data, showmedians=True)
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    #ax.set_ylim(0, 40)
    ax.set_ylabel('number of protections within 20A')
    ax.set_title('P = {}'.format(p))
    fig.tight_layout(pad = 2)
    fig.savefig(name + '.pdf')


random_counts = parse_infile(args.random)
counts = parse_infile(args.real)
random_counts2 = parse_infile(args.random2)
counts2 = parse_infile(args.real2)

x_labels = ['16s Observed', '23s real', '16s Random Distribution', '23s random']
data_toplot = [list(counts.values()), list(counts2.values()), list(random_counts.values()), list(random_counts2.values())]
u, p = scipy.stats.mannwhitneyu(list(counts.values()), list(random_counts.values()), alternative='two-sided')
u2, p2 = scipy.stats.mannwhitneyu(list(counts2.values()), list(random_counts2.values()), alternative='two-sided')
print(p)
print(p2)

violin_plot(args.output, x_labels, data_toplot, str(p) + '\t' + str(p2))