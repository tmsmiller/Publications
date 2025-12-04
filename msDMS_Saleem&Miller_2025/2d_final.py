from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plot
import numpy as np
import argparse

from ReactivityProfile import ReactivityProfile as rprofile

argparser = argparse.ArgumentParser()
argparser.add_argument('--n7prof')
argparser.add_argumnet('--sasa')
args = argparser.parse_args()



n7 = rprofile(args.n7prof)

infile = args.sasa


# read in the sasa file output from pymol script
sasa = {}
with open(infile) as inp:
    for line in inp:
        spl = line.split()
        sasa[int(spl[0])] = float(spl[1])


prot_cut = 1.6
sasa_cut = 1

x = 0
protected = {}
cfprotected = {}
while x < len(n7.normprofile):
    ##loop through the normprofile and take any normalized values above what we consider "protected"
    if np.isfinite(n7.normprofile[x]) and n7.normprofile[x] >= prot_cut and n7.sequence[x] == 'G' and n7.nts[x] in sasa:
        protected[n7.nts[x]] = n7.normprofile[x]
    x += 1


positives = {}
negatives = {}
for i in sasa:
        ##if measured area is greater than our cutoff, it is unprotected and therefore a "negative"
        if sasa[i] >= sasa_cut:
            negatives[i] = sasa[i]
        ##if measured area is lower than our cutoff it is protected and therefore "positive"
        if sasa[i] < sasa_cut:
            positives[i] =  sasa[i]

TP = 0
FP = 0
TN = 0
FN = 0
for i in sasa:
    ##compare our dict of protected values against the positive and negatives at every position, tick up appropriate bin
    if i in protected and i in positives:
        TP += 1
    if i in protected and i in negatives:
        FP += 1
    if i in positives and i not in protected:
        FN += 1
    if i in negatives and i not in protected:
        TN += 1

print(' TP: {}\n FP: {}\n TN: {}\n FN: {}\n'.format(TP, FP, TN, FN))
print('Sensitivity: {}\n'.format(TP/(TP+FN)))
print('Specificity: {}\n'.format(TN/(TN+FP)))
print('PPV: {}\n'.format(TP/(TP+FP)))
print('NPV: {}\n'.format(TN/(TN+FN)))
print('Diagnostic Accuracy: {}'.format((TP + TN)/(TP + FP + FN + TN))) 
print('Protection cutoff: {}\nSASA cutoff:{}\n'.format(prot_cut,sasa_cut))
print(len(sasa))
print(TP+FP+TN+FN)

