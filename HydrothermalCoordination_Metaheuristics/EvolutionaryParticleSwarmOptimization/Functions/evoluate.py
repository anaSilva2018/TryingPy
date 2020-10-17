# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np

def _pop_evaluate(pop, nper, chydth, nvolinic, mafl, mmutate, mload, valdef):
    khidr = (chydth.vturbmax*(10/3))/3600
    vherd = np.zeros([nper, pop])
    vdisp = np.zeros([nper, pop])
    vturb = np.zeros([nper, pop])
    vsobr = np.zeros([nper, pop])
    vdesc = np.zeros([nper, pop])
    psobr = np.zeros([nper, pop])
    phidr = np.zeros([nper, pop])
    pterm = np.zeros([nper, pop])
    mcmg = np.zeros([1, pop])
    mpen = np.zeros([nper, pop])
    mval = np.zeros([nper, pop])
    mcost = np.zeros([nper+1, pop])
    mpop = np.zeros([nper+2, pop])
    for i in range(nper):
        for j in range(pop):
            mpop[i, j] = mmutate[i, j]
            mpop[nper, j] = mmutate[nper, j]
            if(i == 0):
                vherd[i, j] = nvolinic
            else:
                vherd[i, j] = vsobr[i-1, j]
            vdisp[i, j] = vherd[i, j]+mafl[i, 0]
            if(vdisp[i, j] > chydth.vdispmax):
                vdisp[i, j] = chydth.vdispmax
                vdesc[i, j] = vdisp[i, j]-chydth.vdispmax
            vturb[i, j] = mmutate[i, j]*vdisp[i, j]
            vsobr[i, j] = vdisp[i, j]-vturb[i, j]
            psobr[i, j] = (10/3)*(vsobr[i, j]/3600)
            phidr[i, j] = (10/3)*(vturb[i, j]/3600)
            pterm[i, j] = mload[i, 0]- phidr[i, j]
    for i in range(nper):
        for j in range(pop):
            if(phidr[i, j] > khidr):
                mpen[i, j] = (phidr[i, j]-khidr)+mpen[i, j]
            if(phidr[i, j] < 0):
                mpen[i, j] = (phidr[i, j])+mpen[i, j]
            if(pterm[i, j] > chydth.ptermmax):
                mpen[i, j] = (pterm[i, j]-chydth.ptermmax)+mpen[i, j]
            if(pterm[i, j] < 0):
                mpen[i, j] = pterm[i, j]+mpen[i, j]
            if(psobr[nper-1, j] > 0 and np.all(pterm[:, j] <= chydth.ptermmax) and np.all(phidr[:, j] <= khidr)):
                mcmg[0, j] = valdef
                mval[nper-1, j] = psobr[nper-1, j]
            mcost[i, j] = 2000+100*pterm[i, j] + 1.5*np.power(pterm[i, j], 2) + chydth.penaliz*np.power(mpen[i, j], 2) - mcmg[0, j]*mval[i, j]
            mcost[nper, j] = mcost[i, j]+mcost[nper, j]
            mpop[nper+1, j] = mcost[nper, j]
    class _PopVol:
        def __init__(self, vherd, vdisp, vturb, vsobr, vdesc):
            self.vherd = vherd
            self.vdisp = vdisp
            self.vturb = vturb
            self.vsobr = vsobr
            self.vdesc = vdesc
    class _PopPot:
        def __init__(self, psobr, phidr, pterm):
            self.psobr = psobr
            self.phidr = phidr
            self.pterm = pterm
    class _Ccost:
        def __init__(self, mcmg, mpen, mval, mcost):
            self.mcmg = mcmg
            self.mpen = mpen
            self.mval = mval
            self.mcost = mcost
    p_vol = _PopVol(vherd, vdisp, vturb, vsobr, vdesc)
    p_pot = _PopPot(psobr, phidr, pterm)
    p_cost = _Ccost(mcmg, mpen, mval, mcost)
    return p_vol, p_pot, p_cost, mpop
