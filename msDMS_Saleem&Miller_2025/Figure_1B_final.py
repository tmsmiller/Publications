import argparse
from ReactivityProfile import ReactivityProfile as rprofile
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'
from matplotlib import pyplot as plot

argparser = argparse.ArgumentParser()
argparser.add_argument('--ssu_old')
argparser.add_argument('--ssu_new')
argparser.add_argument('--lsu_old')
argparser.add_argument('--lsu_new')
argparse.add_argument('--out')
args = argparser.parse_args()
small_old = rprofile(args.ssu_old)
large_old = rprofile(args.lsu_old)
small_new = rprofile(args.ssu_new)
large_new = rprofile(args.lsu_new)


old = [small_old.rawprofile[52], large_old.rawprofile[72]]
new = [small_new.rawprofile[52],large_new.rawprofile[72]]


xaxis = ['16s', '23s']
barcount = 2
ind = np.arange(barcount)
width = 0.3

fig, ax = plot.subplots()

ax.bar(ind,new, width, label='new')
ax.bar(ind+width , old, width, label='old')
ax.set_xticks(ind + width/2, ('16s', '23s'))


ax.set_title('m7g reactivities')
ax.legend(loc='upper right')
fig.tight_layout()
plot.savefig(args.out)