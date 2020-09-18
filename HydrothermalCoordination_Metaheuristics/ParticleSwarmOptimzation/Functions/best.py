# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 18:27:39 2020

@author: Ana Silva
"""
import numpy as np

def _swarm_best(mpop_init, mpop_ant, nper, pop, cpinit, cvinit, ccinit, cpant, cvant, ccant):
    #pop_init VS pop_ant (Local Best/bi)
    
    mpopbest = np.zeros([nper+1, pop])
    mcmgbest = np.zeros([1, pop])
    mphidrbest = np.zeros([nper, pop])
    mptermbest = np.zeros([nper, pop])
    mvturbbest = np.zeros([nper, pop])
    mvdispbest = np.zeros([nper, pop])
    mvsobrbest = np.zeros([nper, pop])
    for j in range(pop):
        if(mpop_init[nper, j] < mpop_ant[nper, j]):
            mpopbest[nper, j] = mpop_init[nper, j]
            mcmgbest[0, j] = ccinit.mcmg[0, j]
            for i in range(nper):
                mpopbest[i, j] = mpop_init[i, j]
                
                mphidrbest[i, j] = cpinit.phidr[i, j]
                mptermbest[i, j] = cpinit.pterm[i, j]
                mvturbbest[i, j] = cvinit.vturb[i, j]
                mvdispbest[i, j] = cvinit.vdisp[i, j]
                mvsobrbest[i, j] = cvinit.vsobr[i, j]
        else:
            mpopbest[nper, j] = mpop_ant[nper, j]
            mcmgbest[0, j] = ccant.mcmg[0, j]
            for i in range(nper):
                mpopbest[i, j] = mpop_ant[i, j]
                mphidrbest[i, j] = cpant.phidr[i, j]
                mptermbest[i, j] = cpant.pterm[i, j]
                mvturbbest[i, j] = cvant.vturb[i, j]
                mvdispbest[i, j] = cvant.vdisp[i, j]
                mvsobrbest[i, j] = cvant.vsobr[i, j]
    #Best of best(Global Best/bg)
    auxcmin = np.power(10, 9)
   
    for j in range(pop):
        if(mpopbest[nper, j] <  auxcmin):
            auxcmin = mpopbest[nper, j]
            nbg = j
            ncmgbg = mcmgbest[0, j]
            ncostg = mpopbest[nper, j]
    
    class _BestLocal:
        def __init__(self, popbest, cmgbest, phidrbest, ptermbest, vturbbest, vdispbest, vsobrbest):
            self.popbest = popbest
            self.cmgbest = cmgbest
            self.phidrbest = phidrbest
            self.ptermbest = ptermbest
            self.vturbbest = vturbbest
            self.vdispbest = vdispbest
            self.vsobrbest = vsobrbest
    
    class _BestGlobal:
        def __init__(self, nbglobal, ncmgbglobal, ncostbglobal):
            self.nbglobal = nbglobal
            self.ncmgbglobal = ncmgbglobal
            self.ncostbglobal = ncostbglobal
    c_bi =_BestLocal(mpopbest, mcmgbest, mphidrbest, mptermbest, mvturbbest, mvdispbest, mvsobrbest)
    c_bg = _BestGlobal(nbg, ncmgbg, ncostg)
    
    return c_bi, c_bg
    
