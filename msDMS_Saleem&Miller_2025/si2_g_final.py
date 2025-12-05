import re
import sys
import argparse
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as pyplot
import numpy as np

from ReactivityProfile import ReactivityProfile as rprofile
import RNAStructureObjects as RNAtools

argparser = argparse.ArgumentParser()
argparser.add_argument('--n7profile')
argparser.add_argument('--ctfile')
argparser.add_argument('--out')
args = argparser.parse_args()

n7 = rprofile(args.n7profile)

bgprof = n7.backprofile
rawprof = n7.rawprofile
bgsub = n7.subprofile
normprof = n7.normprofile

react, seq = [],[]
x=0


while x < (len(n7.sequence)):
    seq.append(n7.sequence[x])
    react.append(n7.subprofile[x])
    x += 1


ct = RNAtools.CT(args.ctfile)
pairs = ct.pairedResidueList()

s,r = [],[]

# go through the RNA, only looking at Gs with defined reactivites (not NaN), check to see if the G measured is in the crystal structure as well
for i,v in enumerate(seq):
    if v == 'G' and np.isfinite(react[i]):
        r.append(react[i])
        if i in pairs:
            s.append(True)
        else:
            s.append(False)


# compute ROC curve and plot, along with AUC
tmp = roc_curve(s, r, drop_intermediate=False)
pyplot.plot(tmp[0], tmp[1], label='Basepairing: {:.2f}'.format(roc_auc_score(s,r)))
    

pyplot.plot([0, 1], [0, 1], linestyle="--")
pyplot.xlim([0.0, 1.0])
pyplot.ylim([0.0, 1.05])
pyplot.title('')
pyplot.legend()

pyplot.savefig(args.out)