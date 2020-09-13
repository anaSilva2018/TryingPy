# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 12:45:51 2020

@author: Ana Silva
"""
import numpy as np

def iniciateBuses(data_bus,bus):
    nP = 0
    busP = np.zeros([1,1])
    text_theta = ""
    text_Pinj = ""
    for i in range(bus):
        if(data_bus[i,1] != 1):
        #busPQ
            nP = nP+1
            busP.resize(1,nP)
            busP[0,nP-1] = data_bus[i,0]
            
            text_theta = text_theta + "_theta" + str(int(data_bus[i,0]))
            text_Pinj = text_Pinj + "_P" + str(int(data_bus[i,0]))

        else:
            busref = int(data_bus[i,0])
    nzcent = np.size(busP,1)
    
    text_theta = "Rows:" + text_theta 
    text_Pinj = "Rows:" + text_Pinj 

    return busP, busref, nP, nzcent, text_theta, text_Pinj
