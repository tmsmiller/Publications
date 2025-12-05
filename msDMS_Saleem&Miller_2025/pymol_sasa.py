
#__main__.pymol_argv = ['pymol','-qc']

import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--pdbfile')
parser.add_argument('--chain')
parser.add_argument('--output')

args = parser.parse_args()

radius = 4

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 3)
cmd.set('solvent_radius', radius)

cmd.load(args.pdbfile)  



rnumbers = []
cmd.iterate('chain {} and resn G and name N7'.format(args.chain), 'rnumbers.append(resi)')

out = open(args.output,'w')
for i in rnumbers:
    
    if cmd.count_atoms('chain {} and resi {} and name N7'.format(args.chain, i)) == 1:
        
        v = cmd.get_area('chain {} and resi {} and name N7'.format(args.chain, i))
        out.write('{} {}\n'.format(i,v))

out.close()


