import numpy as np
from ReactivityProfile import ReactivityProfile as rprofile
import argparse
from pymol import cmd

def get_nts_profile(profile):
    nt_list = []
    prof = rprofile(profile)
    for i,v in enumerate(prof.normprofile):
        if np.isfinite(v) and v >= 1.6:
            nt_list.append(i + 1)
    return nt_list

def adjust_frame(number, listnts):
    adjustednts = []
    for i in listnts:
        adjustednts.append(i+number)
    return adjustednts



def get_nts(input_file):
    nts = []
    with open(input_file, 'r') as inp:
        for i,line in enumerate(inp):
            if i == 0 or i == 1:
                continue
            spl= line.split()
            nts.append(int(spl[0]))
    return nts

def get_nt_list(inputfile):
    with open(inputfile, 'r') as inp:
        for line in inp:
            nts = line.split()
    return nts


def color_nts(chain, nt_list):
    cmd.set_color('protected', (83, 48, 95))
    for G in nt_list:
        try:
            cmd.select('protected_nts', 'chain {} and resi {}'.format(chain, G), merge=1)
            cmd.color('protected', 'chain {} and resi {}'.format(chain, G))
        except:
            print('Nt {} not in structure'.format(G))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdbfile', required=True, help='path to crystal structure')
    parser.add_argument('--chain', type=str, required=True, help='RNA chain id in crystal structure')
    parser.add_argument('--input', required=True, help='path to infile')
    parser.add_argument('--outputName', required=True)
    parser.add_argument('--adjustframe', type=int)
    parser.add_argument('--profile', action="store_true")
    args = parser.parse_args()
    cmd.load(args.pdbfile)
    #to_color = get_nts(args.input)
    if not args.profile:
        to_color = get_nt_list(args.input)
    else:
        to_color = get_nts_profile(args.input)
    if args.adjustframe is not None:
        to_color = adjust_frame(args.adjustframe,to_color)

    color_nts(args.chain, to_color)

    cmd.save(args.outputName + '.pse')