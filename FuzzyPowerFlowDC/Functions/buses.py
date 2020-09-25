# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np
def _iniciateBuses(data_bus):
    nP = 0
    busP = np.zeros([1, 1])
    text_theta = ""
    text_Pinj = ""
    bus = np.size(data_bus, axis=0)
    for i in range(bus):
        if(data_bus[i, 1] != 1):
        #busP
            nP = nP+1
            busP.resize(1, nP)
            busP[0, nP-1] = data_bus[i, 0]
            text_theta = text_theta + "_theta" + str(int(data_bus[i, 0]))
            text_Pinj = text_Pinj + "_P" + str(int(data_bus[i, 0]))
        else:
            busref = int(data_bus[i, 0])
    nzcent = np.size(busP, 1)
    text_theta = "Rows:" + text_theta
    text_Pinj = "Rows:" + text_Pinj
    class _Cbus:
        def __init__(self, mbusp, nbusref, nP, nz, ttheta, tpinj, tbus):
            self.mbusp = mbusp
            self.nbusref = nbusref
            self.nP = nP
            self.nz = nz
            self.ttheta = ttheta
            self.tpinj = tpinj
            self.tbus = tbus
    cauxb = _Cbus(busP, busref, nP, nzcent, text_theta, text_Pinj, bus)
    return cauxb
    
