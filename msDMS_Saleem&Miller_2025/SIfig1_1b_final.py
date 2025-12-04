from ReactivityProfile import ReactivityProfile as rp
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'
from matplotlib import pyplot as plot
from matplotlib.ticker import AutoMinorLocator

def get_data(listofprofs, subunit):
    m7gs_mn = []
    bgs_mn = []
    if subunit == 'large':
        x = 72
    else:
        print('smallsub')
        x = 52
    for i in listofprofs:
        print(i.rawprofile)
        m7gs_mn.append(i.rawprofile[x])
        allreactivities = 0
        for j in i.rawprofile:
            if np.isfinite(j):
                allreactivities += j
        allreactivities -= i.rawprofile[x]
        bg = allreactivities / (len(i.rawprofile) - 1)
        bgs_mn.append(bg)
    return m7gs_mn, bgs_mn


Original_16s = rp('ph83_16s_16sm7_profile.txt')
half_mn_16s = rp('05mM_16sm7_profile.txt')
two_mn_16s =rp('2mM_16sm7_profile.txt')
two5_mn_16s = rp('25ph83vs2ph83_16sm7_profile.txt')
three_mn_16s = rp('3mM_16sm7_profile.txt')


ph75_16s = rp('ph75_16sm7_profile.txt')
ph8_16s = rp('ph8_16sm7_profile.txt')
ph85_16s = rp('ph85_16sm7_profile.txt')
finalcondition_16s = rp('2mmph9_16sm7_profile.txt')

Original_23s = rp('ph83_23sm7_profile.txt')
half_mn_23s = rp('05mM_23sm7_profile.txt')
two_mn_23s = rp('2mM_23sm7_profile.txt')
two5_mn_23s = rp('25mM_23sm7_profile.txt')
three_mn_23s = rp('3mM_23sm7_profile.txt')

ph75_23s = rp('ph75_23sm7_profile.txt')
ph8_23s = rp('ph8_23sm7_profile.txt')
ph85_23s = rp('ph85_16sm7_profile.txt')
finalcondition_23s = rp('2mmph9_23sm7_profile.txt')

ssu_manganese_m7gs, ssu_manganese_bgs = get_data([half_mn_16s, Original_16s, two_mn_16s, two5_mn_16s, three_mn_16s], '')
ssu_ph_m7gs, ssu_ph_bgs = get_data([ph75_16s, ph8_16s, Original_16s, ph85_16s], '')
ssu_final_m7g, ssu_final_bg = get_data([finalcondition_16s], '')

lsu_manganese_m7gs, lsu_manganese_bgs = get_data([half_mn_23s, Original_23s, two_mn_23s, two5_mn_23s, three_mn_23s], 'large')
lsu_ph_m7gs, lsu_ph_bgs = get_data([ph75_23s, ph8_23s, Original_23s, ph85_23s], 'large')
lsu_final_m7g, lsu_final_bg = get_data([finalcondition_23s], 'large')


manganese_m7gs = [np.mean(i) for i in zip(ssu_manganese_m7gs,lsu_manganese_m7gs)]
ph_m7gs = [np.mean(i) for i in zip(ssu_ph_m7gs,lsu_ph_m7gs)]
final_m7g = [np.mean(i) for i in zip(ssu_final_m7g,lsu_final_m7g)]


manganese_bgs = [np.mean(i) for i in zip(ssu_manganese_bgs,lsu_manganese_bgs)]
ph_bgs = [np.mean(i) for i in zip(ssu_ph_bgs,lsu_ph_bgs)]
final_bg = [np.mean(i) for i in zip(ssu_final_m7g,lsu_final_bg)]


#add in empty data so final condition can be included in plot
ph_m7gs.append(0.0)
ph_bgs.append(0.0)

labels_1 = ['0.5', '1', '2', '2.5', '3']
labels_2 = ['7.5', '8.0', '8.3', '8.5', '9']
labels_3 = ['2mM Mn pH 8.3', '1mM Mn pH 9', '2mM Mn pH 9']
fig,ax = plot.subplots(1, 4, figsize=(10,2))

#ax[2] = ax[0].twinx()
#share yaxis
ax[0].set_title('Manganese Concentration')
ax[0].scatter(labels_1, manganese_m7gs)
ax[0].plot(labels_1, manganese_m7gs)
#add bg data
ax[0].scatter(labels_1, manganese_bgs)
ax[0].scatter(labels_1, [0.0,0.0,final_m7g[0],0.0,0.0])
ax[0].plot(labels_1, manganese_bgs)
ax[0].set_yscale('log', base=10)


ax[1].set_title('Acidity')
ax[1].scatter(labels_2, ph_m7gs)
ax[1].plot(labels_2, ph_m7gs)
#add bg data
ax[1].scatter(labels_2, ph_bgs)
ax[1].plot(labels_2, ph_bgs)
ax[1].scatter(labels_2, [0.0, 0.0 ,0.0 ,0.0, final_m7g[0]])
ax[1].set_yscale('log', base=10)
#ax[3] = ax[1].twinx()




fig.tight_layout()
plot.savefig('SIFigure_1A_averages.pdf')