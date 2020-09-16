# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:55:58 2020

@author: Ana Silva
"""
import numpy as np

def _cSelect(pop, nper, mmutate, mcost):
    mcaux = np.zeros([1, pop])
    mvaux = np.zeros([nper, pop])
    mcsel = np.zeros([1, 2*pop])
    mvsel = np.zeros([nper, 2*pop])
    newpop = np.zeros([nper, pop])
   
    for _l in range(2*pop):
        for _m in range(nper):
            mvsel[_m, _l] = mmutate[_m, _l]
            mcsel[0, _l] = mcost[nper, _l]
    for k in range(pop):
        auxcmin = np.power(10, 9)
        for j in range(2*pop):
            if(mcsel[0, j] < auxcmin):
                auxcmin = mcsel[0, j]
                locmin = j
        auxbest = locmin 
        for i in range(nper):
            mvaux[i, k] = mvsel[i, auxbest]
        mcaux[0, k] = mcsel[0, auxbest]
        mcsel[0, auxbest] = np.power(10, 9)
    for j in range(pop):
        for i in range(nper):
            mvsel[i, j] = mvaux[i, j]
            newpop[i, j] = mvsel[i, j]
        mcsel[0, j] = mcaux[0, j]
    for i in range(2*pop):
        mcost[nper, i] = mcsel[0, i]
    class _PopSelect:
        def __init__(self, mcaux, mvaux, mcostm, mvsel, newpop, best):
            self.mcaux = mcaux
            self.mvaux = mvaux
            self.mcostm = mcostm
            self.mvsel = mvsel
            self.newpop = newpop
            self.best = best
    pselec = _PopSelect(mcaux, mvaux, mcost, mvsel, newpop, auxbest)
    return pselec
    
