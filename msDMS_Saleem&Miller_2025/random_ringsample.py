import argparse
import matplotlib.pyplot as plot
import numpy as np
import random
import copy

plot.rcParams['pdf.fonttype'] = 42
plot.rcParams['font.sans-serif'] = 'Arial'

def read_ring(input_rings):
    length = 0
    n1rings = []
    n7rings = []
    with open(input_rings) as ringin:
        
        for i,line in enumerate(ringin):
            if i > 1:
                spl = line.split()
                n1 = int(spl[0])
                n2 = int(spl[1])
                a = float(spl[9])
                #print(a)
                if a < 2:
                    continue
                else:
                    if n1 <= length and n2 <= length: #filter n1 rings
                        n1rings.append((n1,n2))
                    else:
                        n7rings.append((n1,n2))
            else:
                if i == 0:
                    length = int(line.split()[0])

    #print(n1rings)
    return length,n1rings, n7rings

def read_primary_nts(fname):
    p_locus = []
    with open(fname) as inp:
        for line in inp:
            spl = line.split('-')
            for i in range(int(spl[0]), int(spl[1]) + 1):
                p_locus.append(i)
    print(p_locus)
    return p_locus

def histogram(en, label = '', title='',name='sfd_hist.pdf'):
    fig, ax = plot.subplots()
    ax.hist(en, histtype = 'bar', bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],label=label, align='left')
    ax.axvline(x=np.percentile(en, [50]), label=label + " median " + str(np.percentile(en, [50])))
    ax.axvline(x=np.mean(en), label=label + " mean " + str(np.mean(en)))
    ax.axvline(x=.84, label='Observed')

    
    ax.set_ylabel('Count')
    ax.set_xlabel('% rings between tertiary sites')

    
    ax.set_title(title)
    ax.legend(loc='upper right')
    fig.tight_layout()
    plot.savefig(name)

def sample_fold(rings,loci_nts, rna):
    encapsulated = 0
    for i in rings:
        rna_copy = copy.deepcopy(rna)
        nt1 = random.choice(rna_copy)
        rna_copy.remove(nt1)
        if nt1 - 1 in rna_copy:
            rna_copy.remove(nt1 - 1)
        if nt1 + 1 in rna_copy:
            rna_copy.remove(nt1 + 1)
        nt2 = random.choice(rna_copy)
        #print(nt1, nt2)
        if nt1 in loci_nts and nt2 in loci_nts:
            encapsulated += 1
    return encapsulated / len(rings)

if __name__ == '__main__':
    # tpp_tertnts.txt p546_tertnts.txt
    loci = read_primary_nts('p546_tertnts.txt')

    # p546_r2r3_intersect_shared_corrbuffer.txt satTPP_intersect_shared_corrbuffer.txt
    rna_len,n1,n7 = read_ring('p546_r2r3_intersect_shared_corrbuffer.txt')

    rna_nts = []
    for i in range(1, rna_len + 1):
        rna_nts.append(i)

    samples = []
    while len(samples) < 10000:
        samples.append(sample_fold(n1, loci, rna_nts))
    print(samples)

    
    hit = []
    for i in samples:
        if i >= 0.83: # 
            hit.append(i)
    print(len(hit)/ len(samples))
    histogram(samples, title='P546', name='P546_10k_ringsample_hist.pdf')
   
