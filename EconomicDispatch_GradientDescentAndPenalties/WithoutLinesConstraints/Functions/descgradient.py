# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 21:25:04 2020

@author: Ana Silva
"""
import numpy as np

def _calctcost(data_gen, pgval, cinit, cpar, mpen):
    nger = cinit.ngen
    mems = 0
    mcprod = 0
    mcost = 0
    for i in range(nger):
        mems = pgval[i, 0]*data_gen[i, 5]+mems
        mcprod = data_gen[i, 2]+pgval[i, 0]*data_gen[i, 3]+mcprod
        mcprod = np.power(pgval[i, 0], 2)*data_gen[i, 4]+mcprod
    mcost = cpar.nlamb*mcprod+(1-cpar.nlamb)*mems+np.sum(mpen)
    return mcost, mems, mcprod
def _calcgrad(cpar, cgen, cinitl):
    auxtload = cinitl.newload
    mpgv = cinitl.mpgval
    mconst = cinitl.mpgconst
    nger = cinitl.ngen
    mpen = np.zeros([2*nger + 1, 1])
    mgrad = np.zeros([nger, 1])
    k = 1
    cmin = np.power(10,9)
    while(auxtload != cpar.pload):
        #print(f"************k = {k}***************")
        #gradient and penalties values
        for i in range(nger):
            #Prodution Cost
            mgrad[i, 0] = cpar.nlamb*(cgen[i, 3] + 2*cgen[i, 4]*mpgv[i, 0])
            mgrad[i, 0] = (1-cpar.nlamb)*cgen[i, 5] + mgrad[i, 0]
            if(mpgv[i, 0] < mconst[i, 0]):
                #Pmin penalty
                mpen[i, 0] = k*np.power((mpgv[i, 0] - mconst[i, 0]), 2)
                mgrad[i, 0] = mgrad[i, 0] - 2*k*(mconst[i, 0] - mpgv[i, 0])
            if(mpgv[i, 0] > mconst[i+nger, 0]):
                #Pmax penalty
                mpen[i+nger, 0] = k*np.power((mpgv[i, 0] - mconst[i+nger, 0]), 2)
                mgrad[i, 0] = mgrad[i, 0] + 2*k*(mpgv[i, 0] - mconst[i+nger, 0])
            if(auxtload != mconst[2*nger, 0]):
                #Pload penalty
                mpen[2*nger, 0] = k*np.power((auxtload - mconst[2* nger, 0]), 2)
                mgrad[i, 0] = mgrad[i, 0] + 2*k*(auxtload-mconst[2*nger, 0])
        totalpen = np.sum(mpen)
        #maximum gradient
        posmax = np.argmax(mgrad, axis=0)
        #totalcost
        [auxcost, auxems, auxcprod] = _calctcost(cgen, mpgv, cinitl, cpar, mpen)
        #update step; update pg and load values
        #print(f"\t\t Generator {posmax+1} / Load : {auxtload} MW")
        auxstep = cpar.learnrate
        auxnext = mpgv[posmax, 0] - auxstep*mgrad[posmax, 0]
        auxtload = 0
        for i in range(nger):
            if(i != posmax):
                auxtload = auxtload + mpgv[i, 0]
        while(auxnext + auxtload < cpar.pload or auxnext < mconst[posmax, 0] or auxnext > mconst[posmax+nger, 0]):
            auxstep = auxstep*0.5
            auxnext = mpgv[posmax, 0] - auxstep*mgrad[posmax, 0]
        auxtload = auxtload + auxnext
        #print(f"Newvalue->{auxnext} TotalLoad->{auxtload} Step->{auxstep}\n")
        mpgv[posmax, 0] = auxnext
        if(auxcost < cmin):
            cmin = auxcost
            cbest = cmin
            loadbest = auxtload
            mvalb = mpgv
            penbest = totalpen
            mpenbest = mpen
            mgradbest = mgrad
            posbest = posmax
            cprod = auxcprod
            cems = auxems
            kbest = k   
        k = k+1
    class _Rcost:
        def __init__(self, ntotalcost, ntotalpen, mpen, cemis, cprod):
            self.ntotalcost = ntotalcost
            self.ntotalpen = ntotalpen
            self.mpen = mpen
            self.cemis = cemis
            self.cprod = cprod
    class _Rval:
        def __init__(self, mpval, nload, mgrad, nmax, nk):
            self.mpval = mpval
            self.nload = nload
            self.mgrad = mgrad
            self.nmax = nmax
            self.nk = nk
    crcost = _Rcost(cbest, penbest, mpenbest, cems, cprod)
    crval = _Rval(mvalb, loadbest, mgradbest, posbest, kbest)
    print(f"\nTotalcost: {cbest:.2f} € -> Generation Cost {cprod:.2f} € / Emissions {cems:.2f}")
    return crcost, crval
