import pandas as pd
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'
from matplotlib import pyplot as plt

#df = pd.read_csv('2mmph9_Modified_16sm7_mutation_counts.txt',sep='\t',header=(0))
df = pd.read_csv('ph9_Modified_23sm7_mutation_counts.txt',sep='\t',header=(0))

#23s
m7g = df.iloc[72]
#16s
#m7g = df.iloc[52]


labels = 'GA', 'Other'
data = [m7g['GA'], m7g['GT'] + m7g['GC'] + m7g['complex_insertion'] + m7g['complex_deletion'] + m7g['G_multinuc_mismatch']]

fig, ax = plt.subplots()
ax.pie(data, labels=labels, autopct='%1.1f%%')

print(m7g['GA'])

fig.tight_layout()
plt.savefig('SIFigure_1C_23sfinalRT.pdf')