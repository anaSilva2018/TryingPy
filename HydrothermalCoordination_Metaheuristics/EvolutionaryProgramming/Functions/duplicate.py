# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:22:22 2020

@author: Ana Silva
"""
import numpy as np

def _pop_duplicate(mpop, pop, nper):
    mduplicate = np.zeros([nper, 2*pop])
    for i in range(nper):
        for j in range(pop):
            mduplicate[i, j] = mpop[i, j]
            mduplicate[i, j+pop] = mduplicate[i, j]
    return mduplicate
    
