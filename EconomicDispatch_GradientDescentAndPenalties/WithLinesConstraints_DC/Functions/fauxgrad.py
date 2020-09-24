# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 23:52:05 2020

@author: Ana Silva
"""
import numpy as np
#Constraints
def _mconst(data_gen, cpar, data_lines, clines, cbus):
    nger = cbus.tngen
    nlines = clines.nlines
    mconst = np.zeros([2*nger+nlines+1, 1])
    for i in range(nger):
        mconst[i, 0] = data_gen[i, 6]
        mconst[i+nger, 0] = data_gen[i, 7]
        mconst[2*nger, 0] = cpar.pload
    for i in range(nlines):
        mconst[2*nger+1+i, 0] = data_lines[i, 4]*0.85
    return mconst
#Pgvalues
def _mpgval(data_gen, cbus):
    nger = cbus.tngen
    pgval = np.zeros([nger, 1])
    for i in range(nger):
        pgval[i, 0] = data_gen[i, 7]
    return pgval
#Total_Load
def _vtload(pgval, cbus):
    nger = cbus.tngen
    tload = 0
    for i in range(nger):
        tload = tload + pgval[i, 0]
    return tload
#Pinj
def _mpinj(pgval, data_bus, cbus, data_gen):
    nP = cbus.nP
    mbusP = cbus.mbusP
    nger = cbus.tngen
    pinj = np.zeros([nP, 1])
    for i in range(nP):
        pos = int(mbusP[0, i])-1
        if(data_bus[pos, 2] == 999):
            pinj[i, 0] = -data_bus[pos, 3]
        else:
            pinj[i, 0] = data_bus[pos, 2] - data_bus[pos, 3]
        for j in range(nger):
            if(pos+1 == data_gen[j, 1]):
                pinj[i, 0] = pinj[i, 0]+pgval[j, 0]
    return pinj
#Pij
def _mpij(pinj, clines, cbus):
    msens = clines.msens1
    nlines = clines.nlines
    nP = cbus.nP
    mpij = np.zeros([nlines, 1])
    for i in range(nlines):
        for j in range(nP):
            mpij[i, 0] = mpij[i, 0]+msens[i, j]*pinj[j, 0]
    return mpij
#Gradient_Matrix
def _mgrad(pgval, tload, mpij, mconst, cbus, clines, data_gen, cpar, k):
    nger = cbus.tngen
    nlines = clines.nlines
    ngp = cbus.ngsens
    gP = cbus.mgsens
    msens = clines.msens1
    mgrad = np.zeros([nger, 1])
    for i in range(nger):
        mgrad[i, 0] = cpar.nlamb*(data_gen[i, 3]+2*data_gen[i, 4]* pgval[i, 0])
        mgrad[i, 0] = mgrad[i, 0]+(1-cpar.nlamb)*data_gen[i, 5]
        #Pgmin
        if(pgval[i, 0] < mconst[i, 0]):
            mgrad[i, 0] = -2*k*(mconst[i, 0]-pgval[i, 0])+mgrad[i, 0]
        #Pgmax
        if(pgval[i, 0] > mconst[i+nger, 0]):
            mgrad[i, 0] = 2*k*(pgval[i, 0]-mconst[i+nger, 0])+mgrad[i, 0]
        #Pload
        if(tload != mconst[2*nger, 0]):
            mgrad[i, 0] = 2*k*(tload-mconst[2*nger, 0])+mgrad[i, 0]
    for i in range(ngp):
        gerj = int(gP[0, i])-1
        for j in range(nlines):
            #Pij
            if(abs(np.absolute(mpij[j, 0])) - mconst[2*nger+1+j, 0] > 0):
                mgrad[gerj, 0] = 2*k*(abs(np.absolute(mpij[j, 0]))-mconst[2*nger+1+j, 0])*msens[j, gerj]+mgrad[gerj, 0]
    return mgrad
#Penalties
def _mpen(pgval, tload, mpij, mconst, cbus, clines, k):
    nger = cbus.tngen
    nlines = clines.nlines
    ngp = cbus.ngsens
    mpen = np.zeros([nlines+2*nger+1, 1])
    for i in range(nger):
        #Pgmin
        if(pgval[i, 0] < mconst[i, 0]):
            mpen[i, 0] = np.power(k*(mconst[i, 0]-pgval[i, 0]), 2)
        #Pgmax
        if(pgval[i, 0] > mconst[i+nger, 0]):
            mpen[i+nger, 0] = np.power(k*(pgval[i, 0]-mconst[i+nger, 0]), 2)
        #Pload
        if(tload != mconst[2*nger, 0]):
            mpen[2*nger, 0] = np.power(k*(tload-mconst[2*nger, 0]), 2)
    for i in range(ngp):
        for j in range(nlines):
            #Pij
            if(abs(np.absolute(mpij[j, 0])) - mconst[2*nger+1+j, 0] > 0):
                mpen[2*nger+1+j, 0] = np.power(k*(abs(np.absolute(mpij[j, 0])) - mconst[2*nger+1+j, 0]), 2)
    return mpen
#Cost
def _calctcost(data_gen, pgval, cbus, cpar, mpen):
    nger = cbus.tngen
    mems = 0
    mcprod = 0
    mcost = 0
    for i in range(nger):
        mems = pgval[i, 0]*data_gen[i, 5]+mems
        mcprod = data_gen[i, 2]+pgval[i, 0]*data_gen[i, 3]+mcprod
        mcprod = np.power(pgval[i, 0], 2)*data_gen[i, 4]+mcprod
    mcost = cpar.nlamb*mcprod+(1-cpar.nlamb)*mems+np.sum(mpen)
    return mcost, mems, mcprod
#GradientDescent_DC
def _fgradpenext(cpar, cbus, clines, data_gen, data_bus, data_lines, auxstop, b1):
    cmin = np.power(10, 9)
    nger = cbus.tngen
    totalpen = 0
    auxnext = 0
    auxload = 0
    penbest = np.power(10, 9)
    k = 1
    mconst = _mconst(data_gen, cpar, data_lines, clines, cbus)
    mpval = _mpgval(data_gen, cbus)
    auxload = _vtload(mpval, cbus)
    mpinj = _mpinj(mpval, data_bus, cbus, data_gen)
    mpij = _mpij(mpinj, clines, cbus)
    
    while(abs(auxload-cpar.pload) > auxstop):
        #print(f"************k = {k}***************")
        #initial step
        auxstep = cpar.learnrate
        #maximum gradient
        mgrad = _mgrad(mpval, auxload, mpij, mconst, cbus, clines, data_gen, cpar, k)
        posmax = int(np.argmax(mgrad, axis=0))
        #total penalty
        mpen = _mpen(mpval, auxload, mpij, mconst, cbus, clines, k)
        totalpen = np.sum(mpen)
        #totalcost
        [auxcost, auxems, auxcprod] = _calctcost(data_gen, mpval, cbus, cpar, mpen)
        #print(f"\t\t Generator {posmax+1} / Load : {auxload} MW")
        #update pgval and load
        auxnext = mpval[posmax, 0] - auxstep*mgrad[posmax, 0]
        auxload = 0
        for i in range(nger):
            if(i != posmax):
                auxload = auxload + mpval[i, 0]
        while(auxnext + auxload < cpar.pload or auxnext < mconst[posmax, 0] or auxnext > mconst[posmax+nger, 0]):
            auxstep = auxstep*b1
            auxnext = mpval[posmax, 0] - auxstep*mgrad[posmax, 0]
        auxload = auxload + auxnext
       # print(f"Newvalue->{auxnext} TotalLoad->{auxload} Step->{auxstep}\n")
        mpval[posmax, 0] = auxnext
        #update Pinj and Pij
        mpinj = _mpinj(mpval, data_bus, cbus, data_gen)
        mpij = _mpij(mpinj, clines, cbus)
        #update best
        if(totalpen < penbest):
            cmin = auxcost
            cbest = cmin
            loadbest = auxload
            mvalb = mpval
            mpbest = mpinj
            mpijbest = mpij
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
        def __init__(self, mpinj, mpij, mpval, nload, mgrad, nmax, nk):
            self.mpinj = mpinj
            self.mpij = mpij
            self.mpval = mpval
            self.nload = nload
            self.mgrad = mgrad
            self.nmax = nmax
            self.nk = nk
    crcost = _Rcost(cbest, penbest, mpenbest, cems, cprod)
    crval = _Rval(mpbest, mpijbest, mvalb, loadbest, mgradbest, posbest, kbest)
    print(f"\nTotalcost: {cbest:.2f} € -> Generation Cost {cprod:.2f} € / Emissions {cems:.2f}")
    return crcost, crval
    
