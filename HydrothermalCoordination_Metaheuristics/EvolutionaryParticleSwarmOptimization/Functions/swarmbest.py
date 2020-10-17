# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np

def _swarm_best(cpar, nper, mpopsel, mcmgsel, cpotb, cvolb):
    plocalbest = np.zeros([nper+2, cpar.pop])
    pglbest = np.zeros([nper+2, 1])
    cmgbest = 0
    phbest = np.zeros([nper, 1])
    ptbest = np.zeros([nper, 1])
    psbest = np.zeros([nper, 1])
    vsbest = np.zeros([nper, 1])
    vtbest = np.zeros([nper, 1])
    vdispbest = np.zeros([nper, 1])
    vhbest = np.zeros([nper, 1])
    vdbest = np.zeros([nper, 1])
    for j in range(cpar.pop):
        for i in range(nper+2):
            plocalbest[i, j] = mpopsel[i, j]
            pglbest[i, 0] = mpopsel[i, 0]
        for i in range(nper):
            cmgbest = mcmgsel[0, 0]
            phbest[i, 0] = cpotb.phidr[i, 0]
            ptbest[i, 0] = cpotb.pterm[i, 0]
            psbest[i, 0] = cpotb.psobr[i, 0]
            vsbest[i, 0] = cvolb.vsobr[i, 0]
            vtbest[i, 0] = cvolb.vturb[i, 0]
            vdispbest[i, 0] = cvolb.vdisp[i, 0]
            vhbest[i, 0] = cvolb.vherd[i, 0]
            vdbest[i, 0] = cvolb.vdesc[i, 0]
    class _cpopbest:
        def __init__(self, mbg, mbi, mcmg, ncost):
            self.mbg = mbg
            self.mbi = mbi
            self.mcmg = mcmg
            self.ncost = ncost
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
    cpbest = _cpopbest(pglbest, plocalbest, cmgbest, pglbest[nper+1, 0])
    cptbest = _cpotbest(phbest, ptbest, psbest)
    cvbest = _cvolbest(vsbest, vtbest, vdispbest, vhbest, vdbest)
    return cpbest, cptbest, cvbest
