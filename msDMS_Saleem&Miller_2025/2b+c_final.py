from ReactivityProfile import ReactivityProfile as rprof
import matplotlib.pyplot as plot
import numpy as np
import argparse
plot.rcParams['pdf.fonttype'] = 42
plot.rcParams['font.sans-serif'] = 'Arial'


def reactivityByNt(self, resnums = None, nts=None, name = None):
        """ return a list of reactivities for a given set of nts (array), or nt type"""
        
        pro = self.profile(name)
        
        mask = np.isfinite(pro)
        #with np.errstate(invalid='ignore'):
        #    mask = mask & (pro > -0.3) & (pro < 4.0)

        try:
            ntit = iter(nts)
            ntmask = (self.sequence == next(ntit))
            for n in ntit:
                ntmask = ntmask | (self.sequence == n) 
            
            mask = mask & ntmask
        except TypeError:
            pass

        try:
            resnums = set(resnums)
            mask = mask & np.array([i in resnums for i in self.nts])
        except TypeError:
            pass

        return pro[mask]


def mutHistogram(self, name = None,name_2 = None, name_3=None, nts = None, resnums = None,
                     bins=100,axes = None, pool=False, writename='mutHist.pdf', **kwargs):
        

        
        write = False
        if axes is None:
            fig, axes = plot.subplots()
            write = True

        if pool:
            rxn = ((self.reactivityByNt(name='norm', nts='G', resnums = resnums)))
            rxn_2 = ((name.reactivityByNt(name='norm', nts='G', resnums = resnums)))
            rxn_3 = ((name_2.reactivityByNt(name='norm', nts='G', resnums = resnums)))
            rxn_4 = ((name_3.reactivityByNt(name='norm', nts='G', resnums = resnums)))
            
            bins=np.histogram(np.hstack((rxn,rxn_2,rxn_3)), bins=(np.arange(0.0, 10)))[1]
            
            axes.hist(rxn, histtype = 'step', bins = bins, label='InCell', **kwargs)
            
            axes.hist(rxn_2, histtype= 'step', bins=bins, label='Cell Free', **kwargs)
            
            axes.hist(rxn_3, histtype= 'step', bins=bins, label='Urea', **kwargs)
            axes.hist(rxn_4, histtype= 'step', bins=bins, label='Cell Free -Mg', **kwargs)
        
        else:
        
            rxn = ((self.reactivityByNt(name='norm', nts='G', resnums = resnums)))        
            axes.hist(rxn, histtype = 'step', bins = bins, label='In Cell', **kwargs)

            if name is not None:

                rxn_2 = ((name.reactivityByNt(name='norm', nts='G', resnums = resnums)))

                axes.hist(rxn_2, histtype= 'step', bins=bins, label='Cell Free', **kwargs)
        
            if name_2 is not None:

                rxn_3 = ((name_2.reactivityByNt(name='norm', nts='G', resnums = resnums)))

                axes.hist(rxn_3, histtype= 'step', bins=bins, label='Urea', **kwargs)
            if name_3 is not None:

                rxn_4 = ((name_3.reactivityByNt(name='norm', nts='G', resnums = resnums)))

                axes.hist(rxn_4, histtype= 'step', bins=bins, label='Cell Free -Mg', **kwargs)


        

        axes.legend(loc='upper left')
        axes.set_ylabel('Count')
        axes.set_xlabel('Norm N7 reactivity')
        axes.set_xlim(-3.0, 3.0)
        axes.set_ylim(0,350)
        axes.yaxis.set_ticks(np.arange(0, 375, 25))
        axes.set_title('23s')

        if write:
            plot.savefig(writename)




def winsorize(prof):
    for i, v in enumerate(prof):
        if prof[i] > 3:
            prof[i] = 3
        if prof[i] < -3:
            prof[i] = -3

    return prof


if __name__ == '__main__':
    ntorder = ('A','C','G','U')

    masknt = ('A', 'C', 'U')
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--incell')
    argparser.add_argument('--cf')
    argparser.add_argument('--cfnomg')
    argparser.add_argument('--urea')
    args = argparser.parse_args()
    #read in reactivity profiles
    incell = rprof(args.incell)
    cellfree = rprof(args.cf)
    cf_nomg = rprof(args.cfnomg)
    urea = rprof(args.urea)

    incell.normprofile = winsorize(incell.normprofile)
    cellfree.normprofile = winsorize(cellfree.normprofile)
    urea.normprofile = winsorize(urea.normprofile)
    cf_nomg.normprofile= winsorize(cf_nomg.normprofile)

    #pass in incell first, then denatured
    mutHistogram(incell,cellfree, urea,cf_nomg, bins=np.linspace(-3.0, 3.0, num=25),writename ='Fig2B_23s_InCellvsCFvsUrea_norm_top_08052024.pdf')