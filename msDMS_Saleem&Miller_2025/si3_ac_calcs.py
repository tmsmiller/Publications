import argparse
import matplotlib.pyplot as plot
import numpy as np
import glob
from ReactivityProfile import ReactivityProfile as rprofile

argparser = argparse.ArgumentParser()
argparser.add_argument('--sasa')
argparser.add_argument('--subunit', choices=['16', '23'])
argparser.add_argument('--out')
args = argparser.parse_args()

prot_cut = 1.6

sasa_cut = 1

infile = args.sasa


out = open(args.out,'w')

adjustframe = 0
def adjust_frame(profile, difference):
    profile.nts = [x + difference for x in profile.nts]

#match end half of profile names
false_p = []
false_out = open('{}_false_negatives.txt'.format(sasa.split('.')[0]), 'w')
for file in glob.glob("*ec{}S_profile.txtga".format(args.subunit)):
    print(file)
    n7 = rprofile(file)
    adjust_frame(n7,adjustframe)

    # read in the sasa file output from pymol script
    sasa = {}
    with open(infile) as inp:
        for line in inp:
            spl = line.split()
            sasa[int(spl[0])] = float(spl[1])

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
            false_p.append(str(i))
        if i in negatives and i not in protected:
            TN += 1

    out.write(' TP: {}\n FP: {}\n TN: {}\n FN: {}\n'.format(TP, FP, TN, FN))
    out.write('Sensitivity: {}\n'.format(TP/(TP+FN)))
    out.write('Specificity: {}\n'.format(TN/(TN+FP)))
    out.write('PPV: {}\n'.format(TP/(TP+FP)))
    out.write('NPV: {}\n'.format(TN/(TN+FN)))
    out.write('Diagnostic Accuracy: {}\n'.format((TP + TN)/(TP + FP + FN + TN))) 
    out.write('Protection cutoff: {}\nSASA cutoff:{}\n'.format(prot_cut,sasa_cut))
    out.write(str(len(sasa)))
    out.write('\n')
    out.write(str(TP+FP+TN+FN))
    out.write('\n')

print(len(set(false_p)))
for nt in set(false_p):
    false_out.write(nt+' ')

out.close()

