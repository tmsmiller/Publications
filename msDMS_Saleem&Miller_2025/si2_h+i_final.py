import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plot
import numpy as np
import scipy.stats
import math
import itertools
from ReactivityProfile import ReactivityProfile as rprofile
import argparse

plot.rcParams['pdf.fonttype'] = 42
plot.rcParams['font.sans-serif'] = 'Arial'

argparser = argparse.ArgumentParser()
argparser.add_argument('--prof1', required=True)
argparser.add_argument('--prof2')
argparser.add_argument('--out', required=True)
args = argparser.parse_args()


# read in the N7 reactivity file
n7 = rprofile(args.prof1)
n72 = rprofile(args.prof2)



def reactivityByneighbors(self, add=None, GA=False, Triplet=False, pos=1, byNT=False, name=''):
    '''
    Takes a reactivity profile as an argument.
    Plots reactivity rates based on neighbors. 
    Pass a negative number to look upstream, positive for downstream. 
    Will group purines and pyrimidines by default.
    '''

    def adjacent_values(vals, q1, q3):
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value


    def set_axis_style(ax, labels):
        ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
        ax.set_xlim(0.25, len(labels) + 0.75)
        ax.set_xlabel('Sample name')

    if add is not None:
        self.subprofile = np.concatenate([self.subprofile, add.subprofile])
        self.sequence = np.concatenate([self.sequence, add.sequence])

    def stars(p):
    ##get stars for annotating plots using a p value
        if p < 0.000001:
            return '******'
        elif p < 0.00001:
            return '*****'
        elif p < 0.0001:
            return "****"
        elif (p < 0.001):
            return "***"
        elif (p < 0.01):
            return "**"
        elif (p < 0.05):
            return "*"
        else:
            return "-"

    if Triplet:
            nts = ['A','G','U','C']
            _triplets = {}
            for i in nts:
                for j in nts:
                    for k in nts:
                        _triplets[i + j + k] = []
                
            for n in nts:
                for i, v in enumerate(self.subprofile):
                    if np.isfinite(v) and self.sequence[i] == n:
                        pre = self.sequence[i - 1]
                        nxt = self.sequence[i + 1]
                        _triplets[pre + n + nxt].append(v)


    elif not GA:
            '''
            Data is ordered (A) purines, pyrimidines then (G) purines, pyrimidines, etc.
            
            '''


            nts = ['A','G','U','C']
            data = []
            labels = []
            for n in nts:
                A = []
                G = []
                T = []
                C = []
                purines = []
                pyrimidine = []
                for i,v in enumerate(self.subprofile):
                    if np.isfinite(v) and self.sequence[i] == n:
                        if byNT:
                            if self.sequence[i + pos] == 'A':
                                A.append(v)
                            elif self.sequence[i + pos] == 'G':
                                G.append(v)
                            elif self.sequence[i + pos] == 'U':
                                T.append(v)
                            else:
                                C.append(v)
                        else:
                            if self.sequence[i + pos] == 'A' or self.sequence[i + pos] == 'G':
                                purines.append(v)
                            else:
                                pyrimidine.append(v)
                if byNT:
                    data.append(A.copy())
                    data.append(G.copy())
                    data.append(T.copy())
                    data.append(C.copy())
                    labels.append(str(n) + ' A')
                    labels.append(str(n) + ' G')
                    labels.append(str(n) + ' T')
                    labels.append(str(n) + ' C')
                else:
                    data.append(purines.copy())
                    data.append(pyrimidine.copy())
                    labels.append(str(n) + ' Purines')
                    labels.append(str(n) + ' Pyrimadines')

    else:
        if byNT:
            A = []
            G = []
            T = []
            C = []

            for i, v in enumerate(self.subprofile):
                if np.isfinite(v):
                    if self.sequence[i + pos] == 'A':
                        A.append(v)
                    elif self.sequence[i + pos] == 'G':
                        G.append(v)
                    elif self.sequence[i + pos] == 'U':
                        T.append(v)
                    else:
                        C.append(v)
            labels = ['A','G','U','C']

            data = [A,G,T,C]


        else:
                purines = []
                pyrimidine = []
                print('GA')
                for i, v in enumerate(self.normprofile):
                    if np.isfinite(v):
                        if self.sequence[i + pos] == 'A' or self.sequence[i + pos] == 'G':
                            purines.append(v)
                        else:
                            pyrimidine.append(v)

                labels = ['purines', 'pyrimidines']
                data = [purines, pyrimidine]
    if not Triplet:
        if GA:
            fig, ax = plot.subplots(nrows = 2, ncols = 1, figsize=(6,8))
        elif not GA and not byNT:
            fig, ax = plot.subplots(nrows = 2, ncols = 1, figsize=(8,16))
        else:
            fig, ax = plot.subplots(nrows = 2, ncols = 1, figsize=(30,30))
        quants = []
        for i in data:
            quants.append([.25,.75])


        ax[0].hist(data[0], histtype = 'step', label='Purines')
        ax[0].hist(data[1], histtype = 'step', label='Pyrimidines')
        ax[0].legend()

        for i in data:
            quantiles = np.quantile(i, np.array([ 0.75]))
            print(quantiles)

        tabledata = [['Comparison', 'P value']]
        for i,j in itertools.combinations(data, 2):

                u_stat, p = scipy.stats.mannwhitneyu(i, j, alternative='two-sided')

                p_value = p
                print(p_value)
                print('{}   {}'.format(data.index(i), data.index(j)))

                #get max and min for putting the annotation in the right place
                y_max = np.max(np.concatenate((i, j)))
                y_min = np.min(np.concatenate((i, j)))

                tabledata.append([labels[data.index(i)] + ' vs. ' + labels[data.index(j)], str(p_value)])

            
        ax[1].table(cellText=tabledata, loc='center' )

        fig.tight_layout(pad = 2)
        fig.savefig(name + '.pdf')
    else:
        out = open('{}.txt'.format(name), 'w')
        out.write('triplet\tmean_rep1\tmedian_rep1\ttriplet_count_rep1\n')
        quants = []
        labels = []
        data = []
        for i in _triplets:
            print(i)
            try:
                quantiles = np.quantile(_triplets[i], np.array([0.5]))
                mean = np.mean(_triplets[i])
                out.write('{}\t{:.6}\t{:.6}\t{}\n'.format(i, mean, quantiles[0], len(_triplets[i])))
                labels.append(i)
                data.append(_triplets[i])
                quants.append([.25,.75])
            except:
                continue

        fig, ax = plot.subplots(nrows = 1, ncols = 1, figsize=(30,30))

        ax.violinplot(data, showmedians=True, quantiles=quants)

        ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
        ax.set_ylim(0,0.25)
        ax.set_ylabel('BG sub reactivity')
        fig.tight_layout(pad = 2)
        fig.savefig(name + '.pdf')


reactivityByneighbors(n7, add=n72, GA=True, byNT=False, Triplet=False, pos=1, name=args.out)
