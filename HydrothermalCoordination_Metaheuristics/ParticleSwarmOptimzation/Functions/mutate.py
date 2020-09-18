# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:29:00 2020

@author: Ana Silva
"""
import numpy as np

def _pop_mutate(pop, nper, mpop, sigma):
    mmutate = np.zeros([nper, pop])
    for i in range(nper):
        for k in range(pop):
            mmutate[i, k] = mpop[i, k]
        for j in range(pop):
            auxunif = np.random.uniform(0, 1)
            auxgauss = (np.power((-2*np.log(auxunif)), 0.5))*np.cos(2*np.pi*auxunif)
            mmutate[i, j] = mmutate[i, j] + sigma*auxgauss
            if(mmutate[i, j] > 1):
                mmutate[i, j] = 1
            elif(mmutate[i, j] < 0):
                mmutate[i, j] = 0
    return mmutate
    
