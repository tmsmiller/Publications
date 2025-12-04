
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
parser.add_argument('--profile', required=True, help='path to profile')
parser.add_argument('--ctfile', required=True)
parser.add_argument('--profile2', required=True, help='path to profile')
parser.add_argument('--ctfile2', required=True)
parser.add_argument('--output', required=True, help='output prefix')
parser.add_argument('--cutoff', type=int, default=3)
parser.add_argument('--filterdistance', type=int, default=10)
args = parser.parse_args()


def contact_filter(protections,CT1, filter_distance=0):
    
    num_within_range = {}
    for i in protections:
        count = 0
        for j in protections:
            if i == j:
                continue
            if CT1.contactDistance(i, j) > filter_distance:
                #print('filtered')
                continue
            else:
                count += 1
        num_within_range[i] = count
    return num_within_range

def violin_plot(name, labels, data, p):
    quants = []
    for i in data:
        quants.append([.25,.75])

    fig, ax = plot.subplots(nrows = 1, ncols = 1)
    ax.violinplot(data, showmedians=True, quantiles=quants)
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    #ax.set_ylim(0, 40)
    ax.set_ylabel('number of protections within CF 5')
    ax.set_title('P = {}'.format(p))
    fig.tight_layout(pad = 2)
    fig.savefig(name + '.pdf')


prof = rprofile(args.profile)
CT = RNAtools.CT(args.ctfile)
gs = []
prot = []
random_prot = []
for i,v in enumerate(prof.normprofile):
    if np.isfinite(v):
        gs.append(i+1)
    if v >= args.cutoff:
        prot.append(i +1)

random_prot = random.sample(gs, k=len(prot))
print(prot)
print(random_prot)
random_counts = contact_filter(random_prot, CT, args.filterdistance)
counts = contact_filter(prot, CT, args.filterdistance)

prof2 = rprofile(args.profile2)
CT2 = RNAtools.CT(args.ctfile2)
gs2 = []
prot2 = []
random_prot2 = []
for i,v in enumerate(prof2.normprofile):
    if np.isfinite(v):
        gs2.append(i+1)
    if v >= args.cutoff:
        prot2.append(i +1)

random_prot2 = random.sample(gs2, k=len(prot2))
print(prot2)
print(random_prot2)
random_counts2 = contact_filter(random_prot2, CT2, args.filterdistance)
counts2 = contact_filter(prot2, CT2, args.filterdistance)

x_labels = ['16s Observed', '23s real', '16s Random Distribution', '23s random']
data_toplot = [list(counts.values()), list(counts2.values()), list(random_counts.values()), list(random_counts2.values())]
u, p = scipy.stats.mannwhitneyu(list(counts.values()), list(random_counts.values()), alternative='two-sided')
u2, p2 = scipy.stats.mannwhitneyu(list(counts2.values()), list(random_counts2.values()), alternative='two-sided')
print(p)
print(p2)

violin_plot(args.output, x_labels, data_toplot, str(p) + '\t' + str(p2))

