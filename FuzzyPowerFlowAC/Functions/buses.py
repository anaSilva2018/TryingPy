# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np
def _iniciateBuses(data_bus):
    bus = np.size(data_bus, axis=0)
    nQ = 0
    nP = 0
    busP = np.zeros([1, 1])
    busQ = np.zeros([1, 1])
    for i in range(bus):
        if(data_bus[i, 1] != 1):
        #busPQ
            if(data_bus[i, 1] == 2):
                nP = nP+1
                busP.resize(1, nP)
                busP[0, nP-1] = data_bus[i, 0]
                nQ = nQ+1
                busQ.resize(1, nQ)
                busQ[0, nQ-1] = data_bus[i, 0]
        #busPV
            elif(data_bus[i, 1] == 3):
                nP = nP+1
                busP.resize(1, nP)
                busP[0, nP-1] = data_bus[i, 0]
        #busRef
        else:
            busref = int(data_bus[i, 0])
    nzcent = np.size(busP, 1)+np.size(busQ, 1)
    class _Cbus:
        def __init__(self, mbusP, mbusQ, nbusref, nP, nQ, nz, tbus):
            self.mbusP = mbusP
            self.mbusQ = mbusQ
            self.nbusref = nbusref
            self.nP = nP
            self.nQ = nQ
            self.nz = nz
            self.tbus = tbus
    cauxbus = _Cbus(busP, busQ, busref, nP, nQ, nzcent, bus)
    return cauxbus
    
