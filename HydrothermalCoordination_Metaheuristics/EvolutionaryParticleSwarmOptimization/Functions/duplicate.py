# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 23:34:32 2020

@author: Sofia
"""
import numpy as np

def _popdupl(cpar, mpop, indxa, indxb):
    mdupl = np.zeros([indxa, 2*cpar.pop])
    for i in range(indxb):
        for j in range(cpar.pop):
            mdupl[i, j] = mpop[i, j]
            mdupl[i, j+cpar.pop] = mpop[i, j]
    return mdupl
