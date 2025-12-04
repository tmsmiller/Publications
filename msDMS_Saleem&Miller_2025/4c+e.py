from sklearn.metrics import roc_curve, roc_auc_score

import numpy as np

from ReactivityProfile import ReactivityProfile as rprofile
import RNAStructureObjects as RNAtools
import matplotlib as mpl
import matplotlib.pyplot as plot
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] ='Arial'

import sys, os, math
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('--pdbfile', required=True, help='path to crystal structure')
parser.add_argument('--chain', type=str, required=True, help='RNA chain id in crystal structure')
parser.add_argument('--inputRings', required=True, help='path to ringfile')
parser.add_argument('--adjustFrame', type=int, default=0, help='adjust frame of profile to match crystal structure if necessary')
parser.add_argument('--profile', required=True, help='path to profile')
parser.add_argument('--ringType', choices=('N1', 'N1N7', 'N7'), required=True, help='type of rings to plot')
parser.add_argument('--crystalStart', type=int, required=True, help='first position of the RNA chain')
parser.add_argument('--crystalEnd', type=int, required=True, help='last position of the RNA chain')
parser.add_argument('--outputName', required=True, help='Output prefix for pymol session')
parser.add_argument('--contactFilter', type=int, default=0, help='contact filter distance')
parser.add_argument('--concatenated', action='store_true', help='plot all types of rings')
parser.add_argument('--filterNegative', action='store_true', help='filter negative rings')
parser.add_argument('--ctFile', required=True, help='path to ct file')
parser.add_argument('--APCcutoff', type=int, default=2, help='filter rings with a score lower than or equal to this number')
parser.add_argument('--approx', action='store_true', help='draw ring to closest nt if nt is missing from crystal structure')
parser.add_argument('--renderall', action='store_true', help='render all the rings if a total concat file is passed')
parser.add_argument('--distanceFilter', type=int, default=99999, help='check average contact distance of things within filter(angstroms)')
args = parser.parse_args()



'''
tool for parsing ring data and drawing it onto a crystal structure

takes a profile, ringfile
takes a chain of interest for a crystal structure
takes first and last residue numbers of that chain to get things in frame
takes a frameshift integer
'''

cmd.load(args.pdbfile)

profile = rprofile(args.profile)

rings = []
residues = []

def adjust_frame(profile, ring_list, difference):
    profile.nts = [x + difference for x in profile.nts]
    for i in ring_list:
        
        i[0] += difference
        i[1] += difference


def parse_ringmap(cutoff, name='', filterneg=False):
    rings = []
    with open(name) as inp:
        for p, line in enumerate(inp):
            if p == 0 or p == 1:
                continue
            else:
                split = line.split()
                #split[2] for sig split[9] for alpha
                if float(split[9]) <= cutoff:
                    continue
                    
                if filterneg:
                    if split[3] == '-1':
                        
                        continue
                rings.append([int(split[0]), int(split[1]), float(split[9])])
    return rings

def contact_filter(inrings, filter_distance=0, ctfile=None):
    
    rings = []
    CT1 = None
    if ctfile is not None:
        CT1 = RNAtools.CT(ctfile)
    
    for i in inrings:
            #print(CT1.contactDistance(i[0], i[1]))
            #print(str(i[0]) + ' ' + str(i[1]))
            if CT1.contactDistance(i[0], i[1]) <= filter_distance:
                #print('filtered')
                continue
            else:
                
                i.append(CT1.contactDistance(i[0],i[1]))
                print("contact filter add {}".format(i))
                rings.append(i)

    return rings

def split_concat(concat_list, profile):
    last =  profile.nts[-1]
    n1_rings = []
    n1n7_rings = []
    n7_rings = []
    for i in concat_list:
        if i[0] > last and i[1] > last:
            i[0] -= last
            i[1] -= last
            n7_rings.append(i)
        elif i[0] < last and i[1] < last:
            n1_rings.append(i)
        else:
            n1n7_rings.append(i)
            if i[0] > last:
                i[0] -= last
            else:
                i[1] -= last
    print(len(n1_rings), len(n1n7_rings), len(n7_rings))
    return(n1_rings, n1n7_rings, n7_rings)

def closest(lst, K):
     
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]


def draw(ring_list, chain, crystalstart, crystalend, approx, ringtype, render=False, distancefilter=0):
    distances = []
    scores = []

    
    #get all the residues in a chain. some crystal structures are missing residues in the middle
    
    dchain = 'chain {}'.format(chain)
    cmd.iterate(dchain,'residues.append(int(resi))')
    fresidues = [*set(residues)]
    #print('residues {}'.format(fresidues))
    
    
    if ringtype == 'N1':
        colors = [(44,123,182), (44,123,182), (171,217,233), (255,255,255), (253,174,97), (215,25,28), (215,25,28)]
    else:
        colors = [(91, 218, 190), (91, 218, 190), (0, 95, 154), (255, 255, 255), (218, 91, 172), (83, 48, 95), (83, 48, 95)]

    #change to (20,100) for sig
    cutoffs = [-1e10, -5, -2, 0, 2, 5, 1e10]
    #alpha = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
    gradiant = [False, True, False, False, True, False]

    for v,i in enumerate(ring_list):
        linestart = i[0]
        startstr = 'start_{}'.format(v+1)
        endstr = 'end_{}'.format(v+1)
        lineend = i[1]

        
        if linestart not in fresidues:
            if approx:
                linestart = closest(fresidues, linestart)
            else:
                continue
        if lineend not in fresidues:
            if approx:
                lineend = closest(fresidues, lineend)
            else:
                continue

        if linestart  <= crystalstart or lineend <= crystalstart:
            
            continue
        elif linestart >= crystalend or lineend >= crystalend:
            
            continue

        else:
            distance = cmd.get_distance(atom1='chain {} and resi {} and name P'.format(chain, linestart), atom2='chain {} and resi {} and name P'.format(chain, lineend))
            
            if distance < distancefilter:
                try:
                    print(i)
                    scores.append(i[3] / distance)
                except:
                    print('FAIL at ring {}'.format(i))
            

            distances.append(distance)
            
            for ind in range(len(cutoffs) - 1):
                if cutoffs[ind] < i[2] <= cutoffs[ind + 1]:
                    if gradiant[ind]:
                        calccolor = colorGrad(i[2], colors[ind], colors[ind + 1], cutoffs[ind], cutoffs[ind + 1])
                    else:
                        calccolor = colors[ind]

            cmd.set_color("color {}".format(v+1),calccolor)

            cmd.select(startstr, 'chain {} and resi {} and name P'.format(chain, linestart))
            start = cmd.get_coords('chain {} and resi {} and name P'.format(chain, linestart), 1)
            
            cmd.pseudoatom('points', name='start_{}'.format(v+1), pos=tuple(start[0]), color="color {}".format(v+1))
        
            cmd.deselect()
        
            cmd.select(endstr, 'chain {} and resi {} and name P'.format(chain, lineend))
            end = cmd.get_coords('end_{}'.format(v+1), 1)
            cmd.pseudoatom('points', name='end_{}'.format(v+1), pos=tuple(end[0]), color="color {}".format(v+1))

            cmd.deselect()
            if render:
                print('rendering rings')
                cmd.bond('points and name start_{}'.format(v+1), 'points and name end_{}'.format(v+1))

    print("scores {}".format(scores))
    try:
        print("average score {}".format(sum(scores)/len(scores)))
    except:
        print('no rings in category')
    return distances


def filter_unprot_pos(ring_list, gaprofile):
    '''
    should be run after adjust_frame if applicable
    '''

    normvalues = {}
    for i, v in enumerate(gaprofile.nts):
        if np.isfinite(gaprofile.normprofile[i]):
            normvalues[v] =  gaprofile.normprofile[i]

    for i in ring_list:
        if i[0] in normvalues and normvalues[i[0]] <= 2:
            ring_list.remove(i)
        if i[1] in normvalues and normvalues[i[1]] <= 2:
            ring_list.remove(i)

    return ring_list


def filter_sig(ring_list, cutoff):
    for i in ring_list:
        if float(i[3]) < cutoff:
            ring_list.remove(i)
    return ring_list

def filter_unconcat(concat_list, unconcat_list):
    print('concat length {}'.format(len(concat_list)))
    num_matches = 0
    
    non_matches = []
    original = concat_list.copy()
    for i in unconcat_list:
        if i in original:
            num_matches += 1
            concat_list.remove(i)
            
        elif i not in original:
            non_matches.append(i)
    

    for i in non_matches:
        i[0] += 1
    for i in non_matches:
        if i in concat_list:
            num_matches += 1
    print(concat_list, '\n')
    print(non_matches)
    print('{} matches'.format(num_matches))

    print(len(unconcat_list))
    print("{} filtered length".format(len(concat_list)))
    return concat_list



def dist_histogram(dist_list, bins=10, title='',name='dist_hist.pdf'):
    fig, ax = plot.subplots()
    label = ''
    colors = ['blue', 'orange', 'green']
    for i, v in enumerate(dist_list):
        if len(v) <= 0:
            continue
        if i == 0:
            label = "N7"
        elif i == 1:
            label = 'N1'
        elif i == 2:
            label = "N7"
        
        ax.hist(v, histtype = 'step', bins = bins, label=label)
        ax.axvline(x=np.percentile(v, [50]), color=colors[i],label=label + " median " + str(np.percentile(v, [50])))
        ax.axvline(x=np.mean(v), color=colors[i],label=label + " mean " + str(np.mean(v)))

    
    ax.set_ylabel('Count')
    ax.set_xlabel('P-P Distance A')
    ax.set_xlim(0, 70)
    
    #ax.set_ylim(20,150)
    ax.set_title(title)
    ax.legend(loc='upper right')
    fig.tight_layout()
    plot.savefig(name)

def colorGrad(value, colorMin, colorMax, minValue, maxValue, log=False):
        """ return an interpolated rgb color value """
        
        if log:
            value = 10**-value
            minValue = 10**-minValue
            maxValue = 10**-maxValue

        if value > maxValue:
            return colorMax
        elif value < minValue:
            return colorMin
        else:
            v = value - min(maxValue, minValue)
            v /= abs(maxValue-minValue)
        
            col = []
            for i in range(3):
                col.append( v*(colorMax[i] - colorMin[i]) + colorMin[i] )
            
            return col
    

rings = parse_ringmap(args.APCcutoff, args.inputRings, args.filterNegative)
N1_rings, N1N7_rings, N7_rings = split_concat(rings, profile)

if args.concatenated:
    args.ringType= 'N7'
    N7_rings.extend(N1N7_rings)
    #N1_rings.extend(N7_rings)
    pass

rings_toplot = []
if args.ringType == 'N7':
    rings_toplot = N7_rings.copy()
elif args.ringType == 'N1N7':
    rings_toplot = N1N7_rings.copy()
else:
    rings_toplot = N1_rings.copy()

print('precf {}'.format(rings_toplot))
if args.contactFilter > 0:
    rings_toplot = contact_filter(rings_toplot, filter_distance=args.contactFilter, ctfile=args.ctFile)
    N1_rings =  contact_filter(N1_rings, filter_distance=args.contactFilter, ctfile=args.ctFile)
    N1N7_rings = contact_filter(N1N7_rings, filter_distance=args.contactFilter, ctfile=args.ctFile)
    N7_rings =  contact_filter(N7_rings, filter_distance=args.contactFilter, ctfile=args.ctFile)




adjust_frame(profile, rings_toplot, args.adjustFrame)
dist_data = []
dist_data1 = draw(rings_toplot, args.chain, args.crystalStart, args.crystalEnd, args.approx, args.ringType, render=True)
dist_data2 = draw(N1_rings, args.chain, args.crystalStart, args.crystalEnd, args.approx, args.ringType, render=False)
dist_data.append(dist_data1)
dist_data.append(dist_data2)

dist_histogram(dist_data, bins=np.linspace(0,50,10), title='{} {}cf filterneg={} rings'.format(args.inputRings, args.contactFilter, str(args.filterNegative)), name='{} {}cf filterneg={} rings.pdf'.format(args.inputRings, args.contactFilter, str(args.filterNegative)))
plot.close('all')
cmd.save(args.outputName + '.pse')