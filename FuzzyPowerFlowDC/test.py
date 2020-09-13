# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 13:05:10 2020

@author: Ana Silva
"""

import numpy as np

#Sb=100MW
Sb = 100

#                     Bus Type Pg(MW) Pc(MW) Qg(Mvar) Qc(Mvar)      Vn(kV) 
data_bus = np.matrix([[1, 1, 999,    20,    999,     5,              150    ],
                      [2, 2, 0,      60,    0,       25,             150    ],
                      [3, 2, 0,      40,    0,       15,             150    ],
                      [4, 3, 999,    30,    999,     15,             150    ],
                      [5, 2, 0,      100,   0,       40,             150    ]])

#NOTA: 1->REF; 2->PQ;3->PV

#                      Busi Busj R(p.u.) X(p.u.) Smax(MVA)            
data_lines = np.matrix([[1, 2, 0.03,    0.08,    62   ],
                        [1, 3, 0.05,    0.16,    70   ],
                        [1, 5, 0.09,    0.32,    45   ], 
                        [2, 3, 0.05,    0.08,   100   ],
                        [3, 4, 0.04,    0.10,    77   ],
                        [4, 5, 0.04,    0.10,    97   ]])

#                        Bus Pga(MW) Pgb(MW) Pgc(MW) Qga(Mvar) Qgb(Mvar) Qgc(Mvar) Pca(MW) Pcb(MW) Pcc(MW) Qca(Mvar) Qcb(Mvar) Qcc(Mvar)             
data_difuse = np.matrix([[1, 999,   999,    999,    999,     999,       999,       16,    20,     23,        3,      5,       7],
                         [2, 0,     0,      0,      0,       0,         0,         55,    60,     66,        20,     25,      30],
                         [3, 0,     0,      0,      0,       0,         0,         37,    40,     45,        12,     15,      17], 
                         [4, 130,   135,    140,    999,     999,       999,       23,    30,     32,        12,     15,      18],
                         [5, 0,     0,      0,      0,       0,         0,         95,    100,    106,       35,     40,      46]])


#teta e P->PQ e PV
#v e Q->PQ

bus = np.size(data_bus, axis = 0)
nlines = np.size(data_lines, axis = 0)

print(f" *****Fuzzy Power Flow DC****** \n Study Case: {bus} buses, {nlines} lines and Sb = {Sb} MW .\n ")
from Functions.buses import iniciateBuses as fbuses

[mbusP, busref, nP, nZ, textx, textz] = fbuses(data_bus, bus)


from Functions.lines import R as fR, X as fX, Z as fZ, mZaux as fZaux
r_lines = fR(nlines, bus, data_lines)
x_lines = fX(nlines, bus, data_lines)
z_lines = fZ(nlines, bus, data_lines)
zaux_lines = fZaux(nlines, bus, x_lines, busref)

from Functions.difuse import mZdif, matrixA, mXdif, Pdifuse

[zdifu, zcent, deltaz] = mZdif(data_bus, data_difuse, bus, nP, mbusP, Sb)

[sens, negSens, posSens,textp] = matrixA(x_lines, zaux_lines, busref, nlines, nP, data_lines)

[deltax, xcent, xdifu] = mXdif(zaux_lines, deltaz, zcent, nP)

pdifu = Pdifuse(zdifu, negSens, posSens, nP, nlines)

class CLines:
    def __init__(self, r, x, z, zaux, tbus, tlines):
        self.r = r
        self.x = x
        self.z = z
        self.zaux = zaux
        self.tbus = tbus
        self.tlines = tlines

Mlines = CLines(r = r_lines, x = x_lines, z = z_lines, zaux = zaux_lines, tbus = bus, tlines = nlines)

class Difuse:
    def __init__(self, zdif, xdifu, pdif):
        
        self.zdif = zdif
        self.xdifu = xdifu
        self.pdif = pdif
      
MZfinal = Difuse(zdif = zdifu, xdifu = xdifu, pdif = pdifu)        

np.savetxt('Outputs/zdif.csv', MZfinal.zdif, delimiter = ",", header = '[Zmin, Zcent, Zmax]', footer = textz)
np.savetxt('Outputs/xdif.csv', MZfinal.xdifu, delimiter = ",", header = '[Xmin, Xcent, Xmax]', footer = textx)
np.savetxt('Outputs/pijdif.csv', MZfinal.pdif, delimiter = ",", header = '[Plmin, Plcent, Plmax]', footer =textp)


print("Open folder -> Outputs <-")

