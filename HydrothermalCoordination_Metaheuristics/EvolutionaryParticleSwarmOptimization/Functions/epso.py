# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 00:18:09 2020

@author: Sofia
"""
import numpy as np
from Functions.selection import _cselect

from Functions.duplicate import _popdupl
from Functions.mutate import _popmut
from Functions.evoluate import _pop_evaluate

from Functions.swarmbest import _swarm_best


def _dupmut(cwicm, mpopi, cpar, nper):
    mwid = _popdupl(cpar, cwicm.mwi, nper+1, nper+1)
    mwcd = _popdupl(cpar, cwicm.mwc, nper+1, nper+1)
    mwmd = _popdupl(cpar, cwicm.mwm, nper+1, nper+1)
    mwim = _popmut(cpar, nper, mwid, nper+1, nper+1)
    mwcm = _popmut(cpar, nper, mwcd, nper+1, nper+1)
    mwmm = _popmut(cpar, nper, mwmd, nper+1, nper+1)
    popid = _popdupl(cpar, mpopi, nper+2, nper+1)
    popim = _popmut(cpar, nper, popid, nper+2, nper+1)

    class _cWinit:
        def __init__(self, mwi, mwc, mwm):
            self.mwi = mwi
            self.mwc = mwc
            self.mwm = mwm
    cwicf = _cWinit(mwim, mwcm, mwmm)
    popf = popim
    return cwicf, popf

def _swarminit(cpar, nper, chydth, volinic, mdata_afl, mdata_load, valdef, cswarm):
    mpop = np.zeros([nper+2, cpar.pop])
    mpopant = np.zeros([nper+2, cpar.pop])
    mpopnew = np.zeros([nper+2, cpar.pop])
    mwi = np.zeros([nper+1, cpar.pop])
    mwc = np.zeros([nper+1, cpar.pop])
    mwm = np.zeros([nper+1, cpar.pop])
    mwin = np.zeros([nper+1, cpar.pop])
    mwcn = np.zeros([nper+1, cpar.pop])
    mwmn = np.zeros([nper+1, cpar.pop])
    
    class _cWinit:
        def __init__(self, mwi, mwc, mwm):
            self.mwi = mwi
            self.mwc = mwc
            self.mwm = mwm
    for i in range(nper):
        for j in range(cpar.pop):
            mpop[i, j] = np.random.uniform(0, 1)
            mpop[nper, j] = cpar.sigma
            mwi[i, j] = cswarm._wi
            mwc[i, j] = cswarm._wc
            mwm[i, j] = cswarm._wm
    
    cwicm = _cWinit(mwi, mwc, mwm)
    [cwicmd, mpopd] = _dupmut(cwicm, mpop, cpar, nper)
    [pvol, ppot, pcost, popsel] = _pop_evaluate(2*cpar.pop, nper, chydth, volinic, mdata_afl, mpopd, mdata_load, valdef)  
  
    #popant and popnew
    for i in range(nper+2):
        for j in range(cpar.pop, 2*cpar.pop, 1):
            mpopnew[i, j-cpar.pop] = popsel[i, j-cpar.pop]
            mpopant[i, j-cpar.pop] = popsel[i, j]
    #popbest, bg e bi
    [mselp, mselcmg, cselv, cselpot, cselcw] = _cselect(cpar, nper, popsel, pvol, ppot, pcost, cwicmd)
    [popb, potb, volb] = _swarm_best(cpar, nper, mselp, mselcmg, cselpot, cselv)
    #mwi, mwm, mwc
    for i in range(nper+1):
        for j in range(cpar.pop):
            mwin[i, j] = cwicmd.mwi[i, j]
            mwcn[i, j] = cwicmd.mwc[i, j]
            mwmn[i, j] = cwicmd.mwm[i, j]
    ncwicm = _cWinit(mwin, mwcn, mwmn)
    return ncwicm, popb, potb, volb, mpopant, mpopnew

def _updatebest(cpopb, cpotb, cvolb, mpopsel, cmgs, cpots, cvols, nper, cpar):
    phidrb = cpotb.phidr
    ptermb = cpotb.pterm
    psobrbt = cpotb.psobr
    mb_bi = cpopb.mbi
    mb_bg = cpopb.mbg
    mcmgb = cpopb.mcmg
    ncostb = cpopb.ncost
    bvturb = cvolb.vturb
    bvsobr = cvolb.vsobr
    bvherd = cvolb.vherd
    bvdesc = cvolb.vdesc
    bvdisp = cvolb.vdisp
    
    for j in range(cpar.pop):
        if(mpopsel[nper+1, j] < mb_bg[nper+1, 0]):
            mb_bg[nper+1, 0] = mpopsel[nper+1, j]
            mb_bg[nper, 0] = mpopsel[nper, j]
            ncostb = mpopsel[nper+1, j]
            mcmgb = cmgs[0, j]
            for i in range(nper):
                mb_bg[i, 0] = mpopsel[i, j]
                phidrb[i, 0] = cpots.phidr[i, j]
                ptermb[i, 0] = cpots.pterm[i, j]
                psobrbt[i, 0] = cpots.psobr[i, j]
                bvturb[i, 0] = cvols.vturb[i, j]
                bvsobr[i, 0] = cvols.vsobr[i, j]
                bvherd[i, 0] = cvols.vdesc[i, j]
                bvdesc[i, 0] = cvols.vdesc[i, j]
                bvdisp[i, 0] = cvols.vdisp[i, j]
        if(mpopsel[nper+1, j] < mb_bi[nper+1, j]):
            for i in range(nper+2):
                mb_bi[i, j] = mpopsel[i, j]
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
    cpbestn = _cpopbest(mb_bg, mb_bi, mcmgb, ncostb)
    cptbestn = _cpotbest(phidrb, ptermb, psobrbt)
    cvbestn = _cvolbest(bvsobr, bvturb, bvdisp, bvherd, bvdesc)
    return cpbestn, cptbestn, cvbestn

def _move_updatebest(cpar, nper, mpopn, mpopant, cwicm, cpopi, cpoti, cvoli, valdef, mload, mafl, nvolinic, chydth):
    mbi = cpopi.mbi
    mbg = cpopi.mbg
    ger_pbest = np.zeros([nper, cpar.germax])
    ger_phidr = np.zeros([nper, cpar.germax])
    ger_pterm = np.zeros([nper, cpar.germax])
    ger_psobr = np.zeros([nper, cpar.germax])
    ger_vturb = np.zeros([nper, cpar.germax])
    ger_vsobr = np.zeros([nper, cpar.germax])
    ger_vdesc = np.zeros([nper, cpar.germax])
    ger_vherd = np.zeros([nper, cpar.germax])    
    ger_vdisp = np.zeros([nper, cpar.germax])
    ger_ncmg = np.zeros([1, cpar.germax])
    ger_ncost = np.zeros([1, cpar.germax])
    for ger in range(cpar.germax):
        [mcwif, mPopul] = _dupmut(cwicm, mpopn, cpar, nper)
        auxnewpop = np.zeros([nper+2, 2*cpar.pop])
        #MOVE
        for i in range(nper):
            for j in range(cpar.pop):
                for k in range(2*cpar.pop):
                    auxnewpop[i, k] = mcwif.mwi[i, k]*(mPopul[i, k] - mpopant[i, j])
                    auxnewpop[i, k] = auxnewpop[i, k] + mcwif.mwm[i, k]*(mbi[i, j] - mPopul[i, k])
                    auxnewpop[i, k] = auxnewpop[i, k] + mcwif.mwc[i, k]*(mbg[i, 0] - mPopul[i, k]) 
                    auxnewpop[i, k] = auxnewpop[i, k] + mPopul[i, k]
                    if(auxnewpop[i, k] > 1):
                        auxnewpop[i, k] = 1
                    elif(auxnewpop[i, k] < 0):
                        auxnewpop[i, k] = 0
                    auxnewpop[nper, k] = mPopul[nper, k]
                    auxnewpop[nper+1, k] = mPopul[nper+1, k]
        #Update popant
        mpopant = mpopn
        #Update mPopul            
        mPopul = auxnewpop
        [vnew, pnew, cnew, mPopul] = _pop_evaluate(2*cpar.pop, nper, chydth, nvolinic, mafl, mPopul, mload, valdef)
        #SELECT and update mpopn and cwim
        [mpopn, cmgsel, volsel, potsel, cwicn] = _cselect(cpar, nper, mPopul, vnew, pnew, cnew, mcwif)
        #Update Best
        [cpopi, cpoti, cvoli] = _updatebest(cpopi, cpoti, cvoli, mpopn, cmgsel, potsel, volsel, nper, cpar)
        mbi = cpopi.mbi
        mbg = cpopi.mbg
        ger_ncmg[0, ger] = cpopi.mcmg
        ger_ncost[0, ger] = cpopi.ncost
        for i in range(nper):
            ger_pbest[i, ger] = cpopi.mbg[i, 0]
            ger_psobr[i, ger] = cpoti.psobr[i, 0]
            ger_phidr[i, ger] = cpoti.phidr[i, 0]
            ger_pterm[i, ger] = cpoti.pterm[i, 0]
            ger_vturb[i, ger] = cvoli.vturb[i, 0]
            ger_vsobr[i, ger] = cvoli.vsobr[i, 0]
            ger_vdesc[i, ger] = cvoli.vdesc[i, 0]
            ger_vherd[i, ger] = cvoli.vherd[i, 0]    
            ger_vdisp[i, ger] = cvoli.vdisp[i, 0]
    class _cgerpbest:
        def __init__(self, mpbest, mcmg, mcost):
            self.mpbest = mpbest
            self.mcmg = mcmg
            self.mcost = mcost
    class _cgerpot:
        def __init__(self, phidr, pterm, psobr):
            self.phidr = phidr
            self.pterm = pterm
            self.psobr = psobr
    class _cgervol:
        def __init__(self, vturb, vsobr, vdesc, vherd, vdisp):
            self.vturb = vturb
            self.vsobr = vsobr
            self.vdesc = vdesc
            self.vherd = vherd
            self.vdisp = vdisp     
    cgerpop = _cgerpbest(ger_pbest, ger_ncmg, ger_ncost)    
    cgerpot = _cgerpot(ger_phidr, ger_pterm, ger_psobr)
    cgervol = _cgervol(ger_vturb, ger_vsobr, ger_vdesc, ger_vherd, ger_vdisp)
    return cgerpop, cgerpot, cgervol
    