# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np
def _mZdif(data_difuse, cbus, Sb):
    bus = cbus.tbus
    nP = cbus.nP
    busP = cbus.mbusp
    zdif = np.zeros([nP, 3])
    zcent = np.zeros([nP, 1])
    auxDeltaZ = np.zeros([nP, 3])
    for i in range(bus):
        for j in range(nP):
            if(data_difuse[i, 0] == busP[0, j]):
                aux1 = np.matrix([data_difuse[i, 1], data_difuse[i, 2], data_difuse[i, 3]])
                aux2 = np.matrix([data_difuse[i, 9], data_difuse[i, 8], data_difuse[i, 7]])
                paux = np.subtract(aux1, aux2)/Sb
                zdif[j, :] = paux
                zcent[j, 0] = zdif[j, 1]
                auxDeltaZ[j, :] = zdif[j, :] - zcent[j, 0]
    class _CZdif:
        def __init__(self, mzdif, mzcent, mdeltaz):
            self.mzdif = mzdif
            self.mzcent = mzcent
            self.mdeltaz = mdeltaz
    cdif = _CZdif(zdif, zcent, auxDeltaZ)
    return  cdif
def _matrixA(clines, cbus, data_lines):
    Zaux = clines.mzauxlin
    X_lines = clines.mxlin
    lines = clines.nlines
    busref = cbus.nbusref
    nP = cbus.nP
    mSens = np.zeros([lines, nP])
    mposSens = np.zeros([lines, nP])
    mnegSens = np.zeros([lines, nP])
    posx = 0
    posy = 0
    textauxP = ""
    for i in range(lines):
        posx = int(data_lines[i, 0])
        posy = int(data_lines[i, 1])
        textauxP = textauxP + "_P" + str(int(posx)) +  str(int(posy))
        for j in range(nP):
            if(posx != busref and posy != busref):
                mSens[i, j] = (Zaux[posx-2, j] - Zaux[posy-2, j])/X_lines[posx-1, posy-1]
            elif(posx == busref):
                mSens[i, j] = (0-Zaux[posy-2, j])/X_lines[posx-1, posy-1]
            elif(posy == busref):
                mSens[i, j] = (Zaux[posx-2, j]-0)/X_lines[posx-1, posy-1]
        for j in range(nP):
            if(mSens[i, j] < 0):
                mnegSens[i, j] = mSens[i, j]
            else:
                mposSens[i, j] = mSens[i, j]
    textauxP = "Rows:" + textauxP
    class _CSens:
        def __init__(self, msens, msensn, msensp, txtpij):
            self.msens = msens
            self.msensn = msensn
            self.msensp = msensp
            self.txtpij = txtpij
    csens = _CSens(mSens, mnegSens, mposSens, textauxP)
    return csens
def _mXdif(zaux, cdif, nP):
    auxDeltaZ = cdif.mdeltaz
    zcent = cdif.mzcent
    auxDeltaX = np.zeros([nP, 3])
    #Radians
    auxDeltaX = (np.matmul(zaux, auxDeltaZ)*180)/np.pi
    auxXcent = (np.matmul(zaux, zcent)*180)/np.pi
    auxXdif = auxDeltaX + auxXcent
    class _CXdif:
        def __init__(self, mdeltax, mxcent, mxdif):
            self.mdeltax = mdeltax
            self.mxcent = mxcent
            self.mxdif = mxdif
    cxdif = _CXdif(auxDeltaX, auxXcent, auxXdif)
    return cxdif
def _Pdifuse(zdif, csens, nP, nlines):
    mnegSens = csens.msensn
    mposSens = csens.msensp
    zmin = np.zeros([nP, 1])
    zcent = np.zeros([nP, 1])
    zmax = np.zeros([nP, 1])
    pmin = np.zeros([nlines, 1])
    pcent = np.zeros([nlines, 1])
    pmax = np.zeros([nlines, 1])
    zmin[:, 0] = zdif[:, 0]
    zcent[:, 0] = zdif[:, 1]
    zmax[:, 0] = zdif[:, 2]
    
    pmin = np.matmul(mposSens, zmin) + np.matmul(mnegSens, zmax)
    pcent = np.matmul(mposSens, zcent) + np.matmul(mnegSens, zcent)
    pmax = np.matmul(mposSens, zmax) + np.matmul(mnegSens, zmin)
    pdif = np.concatenate((pmin, pcent, pmax), axis=1)
    return pdif
    
