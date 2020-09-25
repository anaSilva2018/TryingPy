# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np
def _R(cbus, data_lines):
    bus = cbus.tbus
    lines = np.size(data_lines, axis=0)
    r1_lines = np.zeros([bus, bus])
    for i in range(lines):
        x = int(data_lines[i, 0])-1
        y = int(data_lines[i, 1])-1
        r1_lines[x, y] = data_lines[i, 2]
        r1_lines[y, x] = r1_lines[x, y]
    return r1_lines

def _X(cbus, data_lines):
    bus = cbus.tbus
    lines = np.size(data_lines, axis=0)
    x1_lines = np.zeros([bus, bus])
    for i in range(lines):
        x = int(data_lines[i, 0])-1
        y = int(data_lines[i, 1])-1
        x1_lines[x, y] = data_lines[i, 3]
        x1_lines[y, x] = x1_lines[x, y]
    return x1_lines

def _Z(cbus, data_lines):
    bus = cbus.tbus
    lines = np.size(data_lines, axis=0)
    z1_lines = np.zeros([bus, bus], dtype=complex)
    for i in range(lines):
        x = int(data_lines[i, 0])-1
        y = int(data_lines[i, 1])-1
        a = data_lines[i, 2]
        b = data_lines[i, 3]
        z1_lines[x, y] = complex(a, b)
        z1_lines[y, x] = z1_lines[x, y]
    return z1_lines

def _Y(cbus, z_lines, data_lines):
    bus = cbus.tbus
    lines = np.size(data_lines, axis=0)
    y1_lines = np.zeros([bus, bus], dtype=complex)
    for i in range(lines):
        x = int(data_lines[i, 0])-1
        y = int(data_lines[i, 1])-1
        y1_lines[x, y] = 1/(z_lines[x, y])
        y1_lines[y, x] = 1/(z_lines[y, x])
    return y1_lines

def _flines(cbus, data_lines):
    mr = _R(cbus, data_lines)
    mx = _X(cbus, data_lines)
    mz = _Z(cbus, data_lines)
    my = _Y(cbus, mz, data_lines)
    lines = np.size(data_lines, axis=0)
    class _Clines:
        def __init__(self, mrlin, mxlin, mzlin, mylin, nlines):
            self.mrlin = mrlin
            self.mxlin = mxlin
            self.mzlin = mzlin
            self.mylin = mylin
            self.nlines = nlines
    clines = _Clines(mr, mx, mz, my, lines)
    return clines
    
