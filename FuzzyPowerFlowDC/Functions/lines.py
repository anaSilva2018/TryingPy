# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np

def _R(cbus, data_lines):
    lines = np.size(data_lines, axis=0)
    bus = cbus.tbus
    r1_lines = np.zeros([bus, bus])
    for i in range(lines):
        x = int(data_lines[i, 0])-1
        y = int(data_lines[i, 1])-1
        r1_lines[x, y] = data_lines[i, 2]
        r1_lines[y, x] = r1_lines[x, y]
    return r1_lines
def _X(cbus, data_lines):
    lines = np.size(data_lines, axis=0)
    bus = cbus.tbus
    x1_lines = np.zeros([bus, bus])
    for i in range(lines):
        x = int(data_lines[i, 0])-1
        y = int(data_lines[i, 1])-1
        x1_lines[x, y] = data_lines[i, 3]
        x1_lines[y, x] = x1_lines[x, y]
    return x1_lines
def _Z(cbus, data_lines):
    lines = np.size(data_lines, axis=0)
    bus = cbus.tbus
    z1_lines = np.zeros([bus, bus], dtype=complex)
    for i in range(lines):
        x = int(data_lines[i, 0])-1
        y = int(data_lines[i, 1])-1
        a = data_lines[i, 2]
        b = data_lines[i, 3]
        z1_lines[x, y] = complex(a, b)
        z1_lines[y, x] = z1_lines[x, y]
    return z1_lines
def _mZaux(cbus, X_lines):
    bus = cbus.tbus
    busref = cbus.nbusref
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
def _flines(cbus, data_lines):
    mr = _R(cbus, data_lines)
    mx = _X(cbus, data_lines)
    mz = _Z(cbus, data_lines)
    mzaux = _mZaux(cbus, mx)
    lines = np.size(data_lines, axis = 0)
    class _Clines:
        def __init__(self, mrlin, mxlin, mzlin, mzauxlin, nlines):
            self.mrlin = mrlin
            self.mxlin = mxlin
            self.mzlin = mzlin
            self.mzauxlin = mzauxlin
            self.nlines = nlines
    clines = _Clines(mr, mx, mz, mzaux, lines)
    return clines
    
