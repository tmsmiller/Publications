import argparse
import matplotlib.pyplot as plot
import numpy as np
import random
import copy

plot.rcParams['pdf.fonttype'] = 42
plot.rcParams['font.sans-serif'] = 'Arial'

def read_prot_nts(fname):
    sizes = []
    with open(fname) as inp:
        for line in inp:
            sizes.append(int(line.strip()))
    return sizes

def read_primary_nts(fname):
    p_locus = []
    with open(fname) as inp:
        for line in inp:
            spl = line.split('-')
            for i in range(int(spl[0]), int(spl[1]) + 1):
                p_locus.append(i)
    return p_locus

def sample_fold(size_l,loci_nts, rna):
    rna_copy = copy.deepcopy(rna)
    length = len(rna_copy)
    percent_cover = []
    encapsulated = 0
    selected = {}
    for i in size_l:
        nts_in = []
        selected[i] = set()
        while len(selected[i]) == 0:
            prospective = set()
            start = random.choice(rna_copy)
            end = start + i
            if end > length:
                continue
            for j in range(start, end + 1):
                prospective.add(j)


            if prospective.issubset(set(rna_copy)):
                selected[i] = prospective
        for j in selected[i]:
            rna_copy.remove(j)
            if j in loci_nts:
                nts_in.append(j)
        cover = len(nts_in) / len(selected[i])
        percent_cover.append(cover)
        if cover == 1.0:
            encapsulated += 1

    return encapsulated

def histogram(en, label = '', title='',name='sfd_hist.pdf'):
    fig, ax = plot.subplots()
    ax.hist(en, histtype = 'bar', bins=[0,1,2,3,4,5,6,7,8,9,10],label=label, align='left')
    ax.axvline(x=np.percentile(en, [50]), label=label + " median " + str(np.percentile(en, [50])))
    ax.axvline(x=np.mean(en), label=label + " mean " + str(np.mean(en)))
    ax.axvline(x=6, label='Observed')

    
    ax.set_ylabel('Count')
    ax.set_xlabel('# of SFD encapsulated')
    ax.set_xticks(range(10))
    
    ax.set_title(title)
    ax.legend(loc='upper right')
    fig.tight_layout()
    plot.savefig(name)

if __name__ == '__main__':
    #1542 16s
    #2904 23s
    prot_nts = read_prot_nts('23s_selffoldingnts.txt')
    prim_nts = read_primary_nts('23s_primarynts.txt')
    #prot_nts = read_prot_nts('16s_selffoldingnts.txt')
    #prim_nts = read_primary_nts('16s_primarynts.txt')
    ssu = []
    lsu = []
    
    #for i in range(1, 1542 + 1):
    #    ssu.append(i)
    
    for i in range(1, 2904 + 1):
        ssu.append(i)

    #lsu_en = []
    ssu_en = []

    while len(ssu_en) < 10000:
        ssu_en.append(sample_fold(prot_nts, prim_nts, ssu))

    hit = []
    for i in ssu_en:
        if i == 6:
            hit.append(i)
    print(len(hit)/ len(ssu_en))
    histogram(ssu_en, title='23s', name='23s_10k_sfd_hist.pdf')
    
