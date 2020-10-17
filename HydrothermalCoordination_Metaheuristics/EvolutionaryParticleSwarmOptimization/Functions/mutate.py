# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 23:45:55 2020

@author: Sofia
"""
import numpy as np

def _popmut(cpar, nper, mdupl, indxa, indxb):
    mmut = np.zeros([indxa, 2*cpar.pop])
    munif = np.zeros([indxb, 2*cpar.pop])
    mgauss = np.zeros([indxb, 2*cpar.pop])
    for i in range(indxb):
        for j in range(cpar.pop):
            mmut[i, j] = mdupl[i, j]
    for j in range(cpar.pop, 2*cpar.pop, 1):
        munif[nper, j] = np.random.uniform(0, 1, 1)
        mgauss[nper, j] = (np.power((-2*np.log(munif[nper, j])), 0.5))*np.cos(2*np.pi* munif[nper, j])
        mmut[nper, j] = mdupl[nper, j] + cpar.tau*mgauss[nper, j]
        for i in range(nper):
            munif[i, j] = np.random.uniform(0, 1, 1)
            mgauss[i, j] = (np.power((-2*np.log(munif[i, j])), 0.5))*np.cos(2*np.pi* munif[i, j])
            mmut[i, j] = mdupl[i, j] + mmut[nper, j]*mgauss[i, j]
            if mmut[i, j] > 1:
                mmut[i, j] = 1
            elif mmut[i, j] < 0:
                mmut[i, j] = 0
    return mmut
