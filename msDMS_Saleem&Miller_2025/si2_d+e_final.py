from json.tool import main
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys
import RNAtools2 as RNAtools
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
from ReactivityProfile import ReactivityProfile
plot.rcParams['pdf.fonttype'] = 42
plot.rcParams['font.sans-serif'] = 'Arial'
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('--prof1')
argparser.add_argument('--prof2')
argparser.add_argument('--ctfile')
argparser.add_argument('--out')
args = argparser.parse_args()
profobj = ReactivityProfile(args.prof1)
prof2obj = ReactivityProfile(args.prof2)
ctobj = RNAtools.CT(args.ctfile)
pairedlist = ctobj.pairedResidueList(False)

def getNtLists(profile, pairedlist, ntype):
    
        react = []
        ispaired = []
        for i,v in enumerate(profile.subprofile):               
            if v > -10 and profile.sequence[i] == ntype:
                react.append(v)
                ispaired.append( int((i+1) in pairedlist) )
            
        print(react)
        return react, ispaired


labels = ('A', 'C', 'G', 'U')
scores = []
scores2 = []

fig,ax = plot.subplots(1, 4, figsize=(8,2))

for i,nt in enumerate(labels):
    r, p = getNtLists(profobj, pairedlist, nt)
    scores.append(roc_auc_score(p,r))

for i,nt in enumerate(labels):
    r, p = getNtLists(prof2obj, pairedlist, nt)
    scores2.append(roc_auc_score(p,r))


fig,ax = plot.subplots()

barcount = 4
ind = np.arange(barcount)
width = 0.3

fig, ax = plot.subplots()

ax.bar(ind,scores, width, label='new')
ax.bar(ind+width , scores2, width, label='old')
ax.set_xticks(ind + width/2, labels)


fig.tight_layout()
plot.savefig(args.out)