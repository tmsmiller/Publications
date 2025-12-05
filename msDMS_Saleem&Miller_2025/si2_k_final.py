
import re
import sys
import argparse
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as pyplot
import numpy as np

from ReactivityProfile import ReactivityProfile as rprofile

argparser = argparse.ArgumentParser()
argparser.add_argument('--prof')
argparser.add_argument('--sasa')
argparser.add_argument('--out')
args = argparser.parse_args()
# read in the N7 reactivity file
n7 = rprofile(args.prof)

bgprof = n7.backprofile
rawprof = n7.rawprofile
bgsub = n7.subprofile
normprof = n7.normprofile

react, seq = [],[]
x=0


while x < (len(n7.sequence)):
    seq.append(n7.sequence[x])
    react.append(n7.normprofile[x])
    x += 1


# read in the sasa file output from pymol script
sasa = {}
with open(args.sasa) as inp:
    for line in inp:
        spl = line.split()
        sasa[int(spl[0])-1] = float(spl[1])
        
print(sasa)
# plot the data at several different cutoffs
for cut in (1,2):
    s,r = [],[]

    # go through the RNA, only looking at Gs with defined reactivites (not NaN), check to see if the G measured is in the crystal structure as well
    for i,v in enumerate(seq):
        if v == 'G' and np.isfinite(react[i]) and i in sasa:
            r.append(react[i])
            s.append(int(sasa[i] < cut))
    


    # compute ROC curve and plot, along with AUC
    tmp = roc_curve(s, r, drop_intermediate=False)
    pyplot.plot(tmp[0], tmp[1], label='{}: {:.2f}'.format(cut, roc_auc_score(s,r)))
    

pyplot.plot([0, 1], [0, 1], linestyle="--")
pyplot.xlim([0.0, 1.0])
pyplot.ylim([0.0, 1.05])
pyplot.title('')
pyplot.legend()

pyplot.savefig(args.out)


