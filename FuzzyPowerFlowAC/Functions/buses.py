# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 12:45:51 2020

@author: Ana Silva
"""
import numpy as np

def iniciateBuses(data_bus,bus):
    nQ = 0
    nP = 0
    busP = np.zeros([1,1])
    busQ = np.zeros([1,1])
    for i in range(bus):
        if(data_bus[i,1] != 1):
        #busPQ
            if(data_bus[i,1] == 2): 
                nP = nP+1
                busP.resize(1,nP)
                busP[0,nP-1] = data_bus[i,0]
                
                nQ = nQ+1
                busQ.resize(1,nQ)
                busQ[0,nQ-1] = data_bus[i,0]
        #busPV
            elif(data_bus[i,1] == 3):
                nP = nP+1
                busP.resize(1,nP)
                busP[0,nP-1] = data_bus[i,0]
        #busRef
        else:
            busref = int(data_bus[i,0])
    nzcent = np.size(busP,1)+np.size(busQ,1)
    
    return busP, busQ, busref, nP, nQ, nzcent
    
