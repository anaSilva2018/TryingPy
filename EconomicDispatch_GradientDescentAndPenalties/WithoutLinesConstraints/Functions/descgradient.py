# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 21:25:04 2020

@author: Ana Silva
"""
import numpy as np
def _calcgrad(cpar, cgen, cinitl):
    auxtload = cinitl.newload
    mpgv = cinitl.mpgval
    mconst = cinitl.mpgconst
    nger = cinitl.ngen
    mpen = np.zeros([2*nger + 1, 1])
    mgrad = np.zeros([nger, 1])
    k = 1
    while(auxtload != cpar.pload):
        print(f"************k = {k}***************")
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
        #maximum gradient
        posmax = np.argmax(mgrad, axis=0)
        #update step; update pg and load values
        print(f"\t\t Generator {posmax+1} / Load : {auxtload} MW")
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
        
        print(f"Newvalue->{auxnext} TotalLoad->{auxtload} Step->{auxstep}\n")
        mpgv[posmax, 0] = auxnext
        k = k+1   
    
    totalpen = np.sum(mpen)
    class _Cgrad:
        def __init__(self, mpgen, mgradient, mpenalis, totpen, totload):
            self.mpgen = mpgen
            self.mgradient = mgradient
            self.mpenalis = mpenalis
            self.totpen = totpen
            self.totload = totload
    cauxg = _Cgrad(mpgv, mgrad, mpen, totalpen, auxtload)
    return cauxg
    
