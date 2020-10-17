# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np

def _cselect(cpar, nper, mpopd, cvol, cpop, ccost, cwicmd):
    mnpop = mpopd
    pop = cpar.pop
    mcaux = np.zeros([1, pop])
    mvaux = np.zeros([nper+2, pop])
    mphaux = np.zeros([nper, pop])
    mptaux = np.zeros([nper, pop])
    mpsaux = np.zeros([nper, pop])
    mvsaux = np.zeros([nper, pop])
    mvtaux = np.zeros([nper, pop])
    mvdispaux = np.zeros([nper, pop])
    mvhaux = np.zeros([nper, pop])
    mvdaux = np.zeros([nper, pop])
    mcmgaux = np.zeros([1, pop])
    mwiaux = np.zeros([nper+1, pop])
    mwcaux = np.zeros([nper+1, pop])
    mwmaux = np.zeros([nper+1, pop])

    mpopsel = np.zeros([nper+2, pop])
    mphsel = np.zeros([nper, pop])
    mptsel = np.zeros([nper, pop])
    mpssel = np.zeros([nper, pop])
    mvssel = np.zeros([nper, pop])
    mvtsel = np.zeros([nper, pop])
    mvdispsel = np.zeros([nper, pop])
    mvhsel = np.zeros([nper, pop])
    mvdsel = np.zeros([nper, pop])
    mcmgsel = np.zeros([1, pop])
    mwisel = np.zeros([nper+1, pop])
    mwcsel = np.zeros([nper+1, pop])
    mwmsel = np.zeros([nper+1, pop])

    for k in range(pop):
        auxcmin = np.power(10, 9)
        for j in range(2*pop):
            if(mnpop[nper+1, j] < auxcmin):
                auxcmin = mnpop[nper+1, j]
                locmin = j
        best = locmin
        for i in range(nper):
            mvaux[i, k] = mnpop[i, best]
            mphaux[i, k] = cpop.phidr[i, best]
            mptaux[i, k] = cpop.pterm[i, best]
            mpsaux[i, k] = cpop.psobr[i, best]
            mvsaux[i, k] = cvol.vsobr[i, best]
            mvtaux[i, k] = cvol.vturb[i, best]
            mvdispaux[i, k] = cvol.vdisp[i, best]
            mvhaux[i, k] = cvol.vherd[i, best]
            mvdaux[i, k] = cvol.vdesc[i, best]
            mcmgaux[0, k] = ccost.mcmg[0, best]
            mwiaux[i, k] = cwicmd.mwi[i, best]
            mwcaux[i, k] = cwicmd.mwc[i, best]
            mwmaux[i, k] = cwicmd.mwm[i, best]
        mvaux[nper, k] = mnpop[nper, best]
        mwiaux[nper, k] = cwicmd.mwi[nper, best]
        mwcaux[nper, k] = cwicmd.mwc[nper, best]
        mwmaux[nper, k] = cwicmd.mwm[nper, best]
        mcaux[0, k] = mnpop[nper+1, best]
        mnpop[nper+1, best] = np.power(10, 9)
   
    for j in range(pop):
        for i in range(nper):
            mpopsel[i, j] = mvaux[i, j]
            mphsel[i, j] = mphaux[i, j]
            mptsel[i, j] = mptaux[i, j]
            mpssel[i, j] = mpsaux[i, j]
            mvssel[i, j] = mvsaux[i, j]
            mvtsel[i, j] = mvtaux[i, j]
            mvdispsel[i, j] = mvdispaux[i, j]
            mvdsel[i, j] = mvdaux[i, j]
            mvhsel[i, j] = mvhaux[i, j]
            mcmgsel[0, j] = mcmgaux[0, j]
            mwisel[i, j] = mwiaux[i, j]
            mwcsel[i, j] = mwcaux[i, j]
            mwmsel[i, j] = mwmaux[i, j]
        mpopsel[nper, j] = mvaux[nper, j]
        mpopsel[nper+1, j] = mcaux[0, j]
        mwisel[nper, j] = mwiaux[nper, j]
        mwcsel[nper, j] = mwcaux[nper, j]
        mwmsel[nper, j] = mwmaux[nper, j]
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
    class _cWinit:
        def __init__(self, mwi, mwc, mwm):
            self.mwi = mwi
            self.mwc = mwc
            self.mwm = mwm
    cpotb = _cpotbest(mphsel, mptsel, mpssel)
    cvolb = _cvolbest(mvssel, mvtsel, mvdispsel, mvhsel, mvdsel)
    cwicmb = _cWinit(mwisel, mwcsel, mwmsel)
    return mpopsel, mcmgsel, cvolb, cpotb, cwicmb
