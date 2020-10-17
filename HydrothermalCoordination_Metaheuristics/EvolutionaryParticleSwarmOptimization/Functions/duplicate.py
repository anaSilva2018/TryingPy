# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np

def _popdupl(cpar, mpop, indxa, indxb):
    mdupl = np.zeros([indxa, 2*cpar.pop])
    for i in range(indxb):
        for j in range(cpar.pop):
            mdupl[i, j] = mpop[i, j]
            mdupl[i, j+cpar.pop] = mpop[i, j]
    return mdupl
