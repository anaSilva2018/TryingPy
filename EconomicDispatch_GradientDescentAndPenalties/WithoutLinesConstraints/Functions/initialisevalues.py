# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 22:55:31 2020

@author: Ana Silva
"""
import numpy as np
def _iniatilise(cgen, cpar):
    auxtload = 0
    ngen = np.size(cgen, axis=0)
    pgconst = np.zeros([2*ngen + 1, 1])
    pgval = np.zeros([ngen, 1])
    text_const = ""
    text_Pg = ""
    for i in range(ngen):
        text_const = text_const + "_Pg" + str(int(i+1))+"min"
        text_Pg = text_Pg + "_Pg" + str(int(i+1))
    for i in range(ngen):
        pgconst[i, 0] = cgen[i, 6]
        text_const = text_const + "_Pg" + str(int(i+1))+"max"
        pgconst[i+ngen, 0] = cgen[i, 7]
        pgval[i, 0] = cgen[i, 7]
        auxtload = auxtload + pgval[i, 0]
        pgconst[2*ngen, 0] = cpar.pload
    text_const = text_const + "_Ptotaload"
    class _Cinit:
        def __init__(self, mpgval, mpgconst, newload, ngen):
            self.mpgval = mpgval
            self.mpgconst = mpgconst
            self.newload = newload
            self.ngen = ngen
    class _TextInit:
        def __init__(self, tcst, tpg):
            self.tcst = tcst
            self.tpg = tpg
    cauxinit = _Cinit(pgval, pgconst, auxtload, ngen)
    cauxt = _TextInit(text_const, text_Pg)
    return cauxinit, cauxt
    
