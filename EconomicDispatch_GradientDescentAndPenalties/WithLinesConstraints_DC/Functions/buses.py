# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np

def _iniciateBuses(data_bus, data_gen, data_lines):
    bus = np.size(data_bus, axis=0)
    ngen = np.size(data_gen, axis=0)
    nlines = np.size(data_lines, axis=0)
    genP = np.zeros([1, 1])
    mger = 0
    nP = 0
    busP = np.zeros([1, 1])
    text_Pg = ""
    text_Pmin = ""
    text_Pmax = ""
    text_Pij = ""
    text_const = ""
    #Pinj
    for i in range(bus):
        if(data_bus[i, 1] != 1):
            nP = nP+1
            busP.resize(1, nP)
            busP[0, nP-1] = data_bus[i, 0]
        else:
            busref = int(data_bus[i, 0])
    for i in range(ngen):
        text_Pmin = text_Pmin + "_Pgmin" + str(int(i+1))
        text_Pmax = text_Pmax + "_Pgmax" + str(int(i+1))
        text_Pg = text_Pg + "_Pg" + str(int(i+1))    
    #Gen_Msens
    for i in range(nlines):
        posx = int(data_lines[i, 0])
        posy = int(data_lines[i, 1])
        text_Pij = text_Pij + "_P" + str(int(posx))+str(int(posy))
        for j in range(ngen):
            if(posx == int(data_gen[j, 1]) or posy == int(data_gen[j, 1])):
                if(int(data_gen[j, 1]) != busref):
                    if(np.all(genP != data_gen[j, 0])):
                        mger = mger+1
                        genP.resize(1, mger)
                        genP[0, mger-1] = data_gen[j, 0]
    text_const = text_const + text_Pmin
    text_const = text_const + text_Pmax
    text_const = text_const + "_Ptotaload"
    text_const = text_const + text_Pij
    class _Cbuses:
        def __init__(self, mbusP, nbusref, nP, nbus, ngsens, mgsens, tngen):
            self.mbusP = mbusP
            self.nbusref = nbusref
            self.nP = nP
            self.nbus = nbus
            self.ngsens = ngsens
            self.mgsens = mgsens
            self.tngen = tngen
    class _Ctext:
        def __init__(self, tpg, tconst, tpij):
            self.tpg = tpg
            self.tconst = tconst
            self.tpij = tpij
    auxcbus = _Cbuses(busP, busref, nP, bus, mger, genP, ngen)
    ctxt = _Ctext(text_Pg, text_const, text_Pij)
    return auxcbus, ctxt
