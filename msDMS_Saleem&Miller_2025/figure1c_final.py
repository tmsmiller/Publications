from ReactivityProfile import ReactivityProfile as rprofile
import argparse
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'
from matplotlib import pyplot as plot

argparser = argparse.ArgumentParser()
argparser.add_argument('--sat')
argparser.add_argument('--unsat')
argparser.add_argument('--out')
args = argparser.parse_args()

satTPP = rprofile(args.sat)
unsatTPP = rprofile(args.unsat)

fig, ax = plot.subplots()


#manual bg sub
unsat_man_bg = np.nan_to_num(unsatTPP.rawprofile[18:120]) - np.nan_to_num(unsatTPP.backprofile[18:120])
sat_man_bg = np.nan_to_num(satTPP.rawprofile[18:120]) - np.nan_to_num(satTPP.backprofile[18:120])


ax.step(satTPP.nts[18:120], unsat_man_bg, label='Unsaturated', where='mid')
ax.step(satTPP.nts[18:120], sat_man_bg, label='Saturated', where='mid')



ax.set_title('N7 reactivities')
ax.legend(loc='upper right')
fig.tight_layout()
plot.savefig(args.out)