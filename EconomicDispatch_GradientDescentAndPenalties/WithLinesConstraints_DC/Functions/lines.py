# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 12:48:15 2020

@author: Ana Silva
"""
import numpy as np

def _initlines(bus, data_lines, busref, nP):
    lines = np.size(data_lines, axis=0)
    def _R(bus, data_lines, lines):
        r1_lines = np.zeros([bus, bus])
        for i in range(lines):
            x = int(data_lines[i, 0])-1
            y = int(data_lines[i, 1])-1
            r1_lines[x, y] = data_lines[i, 2]
            r1_lines[y, x] = r1_lines[x, y]
        return r1_lines
    def _X(bus, data_lines, lines):
        x1_lines = np.zeros([bus, bus])
        for i in range(lines):
            x = int(data_lines[i, 0])-1
            y = int(data_lines[i, 1])-1
            x1_lines[x, y] = data_lines[i, 3]
            x1_lines[y, x] = x1_lines[x, y]
        return x1_lines
    def _Z(bus, data_lines, lines):
        z1_lines = np.zeros([bus, bus], dtype=complex)
        for i in range(lines):
            x = int(data_lines[i, 0])-1
            y = int(data_lines[i, 1])-1
            a = data_lines[i, 2]
            b = data_lines[i, 3]
            z1_lines[x, y] = complex(a, b)
            z1_lines[y, x] = z1_lines[x, y]
        return z1_lines
    def _mZaux(bus, X_lines, busref):
        auxB = np.zeros([bus, bus])
        b1_lines = np.zeros([bus-1, bus-1])
        for i in range(bus):
            for j in range(bus):
                if(X_lines[i, j] != 0):
                    auxB[i, j] = -1/(X_lines[i, j])
                    auxB[i, i] = auxB[i, i]+1/(X_lines[i, j])
        b1_lines = np.delete(auxB, busref-1, 0)
        b1_lines = np.delete(b1_lines, busref-1, 1)
        zaux = np.linalg.inv(b1_lines)
        return zaux
    def _matrixA(x_lines, zaux, busref, nP, data_lines, lines):
        msens = np.zeros([lines, nP])
        posx = 0
        posy = 0
        for i in range(lines):  
            posx = int(data_lines[i, 0])
            posy = int(data_lines[i, 1])
            for j in range(nP):
                if(posx != busref and posy != busref):
                    msens[i, j] = (zaux[posx-2, j] - zaux[posy-2, j])/x_lines[posx-1, posy-1]
                elif(posx == busref):
                    msens[i, j] = (0-zaux[posy-2, j])/x_lines[posx-1, posy-1]
                elif(posy == busref):
                    msens[i, j] = (zaux[posx-2, j]-0)/x_lines[posx-1, posy-1]
        return msens
    mr1 = _R(bus, data_lines, lines)
    mx1 = _X(bus, data_lines, lines)
    mz1 = _Z(bus, data_lines, lines)
    m1zaux = _mZaux(bus, mx1, busref)
    m1sens = _matrixA(mx1, m1zaux, busref, nP, data_lines, lines)
    class _Clines:
        def __init__(self, mr, mx, mz, mzaux, msens1, nlines):
            self.mr = mr
            self.mx = mx
            self.mz = mz
            self.mzaux = mzaux
            self.msens1 = msens1
            self.nlines = nlines
    cauxl = _Clines(mr1, mx1, mz1, m1zaux, m1sens, lines)
    return cauxl
    
