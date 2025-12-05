import sys
from ReactivityProfile import ReactivityProfile
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import argparse
mpl.use('Agg')
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'

argparser = argparse.ArgumentParser()
argparser.add_argument('--prof1')
argparser.add_argument('--prof2')
argparser.add_argument('--out')
args = argparser.parse_args()

prof1 = ReactivityProfile(args.prof1)
prof2 = ReactivityProfile(args.prof2)

fig, ax = plt.subplots()

ax.step(prof1.nts, (prof1.rawprofile * 100), label=args.prof1, where='mid', linewidth=2, linestyle='-')
ax.step(prof1.nts, (prof2.rawprofile * 100), label=args.prof2, where='mid', linewidth=2)

#52 for ssu 72 for lsu
seq =  [str(i) for i in prof1.sequence]
#seq[72] = 'm7G'
seq[52] = 'm7G'

ax.set_xlabel('Nucleotide', fontsize=20)
ax.set_xticks(prof1.nts)
ax.set_xticklabels(seq)

#ax.set_xlim([65,82])
ax.set_xlim([44,64])


ax.set_ylabel('Mutation rate(%)', fontsize=20)
ax.set_ylim([0, 25])

ax.legend()
fig.tight_layout()

fig.savefig(args.out)

