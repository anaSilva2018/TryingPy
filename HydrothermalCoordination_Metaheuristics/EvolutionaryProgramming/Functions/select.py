# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np

def _cSelect(pop, nper, mmutate, ccost, cvol, cpot):
    
    mcost = ccost.mcost
    mcaux = np.zeros([1, pop])
    mvaux = np.zeros([nper, pop])
    mphaux = np.zeros([nper, pop])
    mptaux = np.zeros([nper, pop])
    mpsaux = np.zeros([nper, pop])
    mvsaux = np.zeros([nper, pop])
    mvtaux = np.zeros([nper, pop])
    mvdispaux = np.zeros([nper, pop])
    mvhaux = np.zeros([nper, pop])
    mvdaux = np.zeros([nper, pop])
    mcmgaux = np.zeros([1, pop])
    
    mcsel = np.zeros([1, pop])
    mvsel = np.zeros([nper, pop])
    mphsel = np.zeros([nper, pop])
    mptsel = np.zeros([nper, pop])
    mpssel = np.zeros([nper, pop])
    mvssel = np.zeros([nper, pop])
    mvtsel = np.zeros([nper, pop])
    mvdispsel = np.zeros([nper, pop])
    mvhsel = np.zeros([nper, pop])
    mvdsel = np.zeros([nper, pop])
    mcmgsel = np.zeros([1, pop])
    
    for k in range(pop):
        auxcmin = np.power(10, 9)
        for j in range(2*pop):
            if(mcost[nper, j] < auxcmin):
                auxcmin = mcost[nper, j]
                locmin = j
        best = locmin 
        for i in range(nper):
            mvaux[i, k] = mmutate[i, best]
            mphaux[i, k] = cpot.phidr[i, best]
            mptaux[i, k] = cpot.pterm[i, best]
            mpsaux[i, k] = cpot.psobr[i, best]
            mvsaux[i, k] = cvol.vsobr[i, best]
            mvtaux[i, k] = cvol.vturb[i, best]
            mvdispaux[i, k] = cvol.vdisp[i, best]
            mvhaux[i, k] = cvol.vherd[i, best]
            mvdaux[i, k] = cvol.vdesc[i, best]
            mcmgaux[0, k] = ccost.mcmg[0, best]
        mcaux[0, k] = mcost[nper, best]
        mcost[nper, best] = np.power(10, 9)
    for j in range(pop):
        for i in range(nper):
            mvsel[i, j] = mvaux[i, j]
            mphsel[i, j] = mphaux[i, j]
            mptsel[i, j] = mptaux[i, j]
            mpssel[i, j] = mpsaux[i, j]
            mvssel[i, j] = mvsaux[i, j]
            mvtsel[i, j] = mvtaux[i, j]
            mvdispsel[i, j] = mvdispaux[i, j]
            mvdsel[i, j] = mvdaux[i, j]
            mvhsel[i, j] = mvhaux[i, j]
            mcmgsel[0, j] = mcmgaux[0, j]
        mcsel[0, j] = mcaux[0, j]
    
    class _PopSelect:
        def __init__(self, mcaux, mvaux, mcsel, mvsel, mcmg):
            self.mcaux = mcaux
            self.mvaux = mvaux
            self.mcsel = mcsel
            self.mvsel = mvsel
            self.mcmg = mcmg
    class _cpotbest:
        def __init__(self, phidr, pterm, psobr):
            self.phidr = phidr
            self.pterm = pterm
            self.psobr = psobr
    class _cvolbest:
        def __init__(self, vsobr, vturb, vdisp, vherd, vdesc):
            self.vsobr = vsobr
            self.vturb = vturb
            self.vdisp = vdisp
            self.vherd = vherd
            self.vdesc = vdesc
    cpotb = _cpotbest(mphsel, mptsel, mpssel)
    cvolb = _cvolbest(mvssel, mvtsel, mvdispsel, mvhsel, mvdsel)
    pselec = _PopSelect(mcaux, mvaux, mcsel, mvsel, mcmgsel)
    return pselec, cpotb, cvolb
