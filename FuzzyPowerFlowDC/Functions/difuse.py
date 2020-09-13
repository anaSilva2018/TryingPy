# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:57:03 2020

@author: Ana Silva
"""
import numpy as np


def mZdif(data_bus, data_difuse, bus, nP, busP, Sb):
    zdif = np.zeros([nP,3])
    zcent = np.zeros([nP,1])
    auxDeltaZ = np.zeros([nP,3])
    for i in range (bus):
        for j in range (nP):
            if(data_difuse[i,0] == busP[0,j]):
                aux1 = np.matrix([data_difuse[i,1], data_difuse[i,2], data_difuse[i,3]])
                aux2 = np.matrix([data_difuse[i,9], data_difuse[i,8], data_difuse[i,7]])
                paux = np.subtract(aux1,aux2)/Sb
                zdif[j,:] = paux
                zcent[j,0] = zdif[j,1]
                auxDeltaZ[j,:] = zdif[j,:] - zcent[j,0]
    
    return  zdif, zcent, auxDeltaZ


def matrixA(X_lines, Zaux, busref, lines, nP, data_lines):
    mSens = np.zeros([lines,nP])
    mposSens = np.zeros([lines,nP])
    mnegSens = np.zeros([lines,nP])
    posx = 0
    posy = 0
    textauxP = ""
    
    for i in range(lines):
        posx = int(data_lines[i,0])
        posy = int(data_lines[i,1])
        
        textauxP = textauxP + "_P" + str(int(posx)) +  str(int(posy))
        for j in range(nP):
            if(posx != busref and posy != busref):
               mSens[i,j] = (Zaux[posx-2,j] - Zaux[posy-2,j])/X_lines[posx-1,posy-1]
            elif(posx == busref):
               mSens[i,j] = (0-Zaux[posy-2,j])/X_lines[posx-1,posy-1]
            elif(posy == busref):
               mSens[i,j] = (Zaux[posx-2,j]-0)/X_lines[posx-1,posy-1]
        
        for j in range(nP):
            if(mSens[i,j] <0 ):
                mnegSens[i,j] = mSens[i,j]
            else: 
                mposSens[i,j] = mSens[i,j]
    textauxP = "Rows:" + textauxP     
    return mSens, mnegSens, mposSens, textauxP

def mXdif(zaux, auxDeltaZ, zcent, nP):
    auxDeltaX = np.zeros([nP, 3])
    #Radians
    auxDeltaX = (np.matmul(zaux, auxDeltaZ)*180)/np.pi
    auxXcent = (np.matmul(zaux,zcent)*180)/np.pi
    auxXdif = auxDeltaX + auxXcent
    
    return auxDeltaX, auxXcent, auxXdif


def Pdifuse(zdif, mnegSens, mposSens, nP, nlines):
    zmin = np.zeros([nP,1])
    zcent = np.zeros([nP,1])
    zmax = np.zeros([nP,1])
    
    pmin = np.zeros([nlines,1])
    pcent = np.zeros([nlines,1])
    pmax = np.zeros([nlines,1])


    zmin[:,0] = zdif[:,0]
    zcent[:,0] = zdif[:,1]
    zmax[:,0] = zdif[:,2]
    
    pmin = np.matmul(mposSens, zmin) + np.matmul(mnegSens, zmax)
    pcent = np.matmul(mposSens, zcent) + np.matmul(mnegSens, zcent)
    pmax = np.matmul(mposSens, zmax) + np.matmul(mnegSens, zmin)

    pdif = np.concatenate((pmin, pcent, pmax), axis = 1)
    
    return pdif
    
    
