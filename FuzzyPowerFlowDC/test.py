# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""

import numpy as np
from Functions.buses import _iniciateBuses
from Functions.lines import _flines
from Functions.difuse import _mZdif, _matrixA, _mXdif, _Pdifuse
#Sb=100MW
Sb = 100
#                     Bus Type Pg(MW) Pc(MW) Qg(Mvar) Qc(Mvar)      Vn(kV)
data_bus = np.matrix([[1, 1, 999, 20, 999, 5, 150],
                      [2, 2, 0, 60, 0, 25, 150],
                      [3, 2, 0, 40, 0, 15, 150],
                      [4, 3, 999, 30, 999, 15, 150],
                      [5, 2, 0, 100, 0, 40, 150]])
#NOTA: 1->REF; 2->PQ;3->PV
#                      Busi Busj R(p.u.) X(p.u.) Smax(MVA)
data_lines = np.matrix([[1, 2, 0.03, 0.08, 62],
                        [1, 3, 0.05, 0.16, 70],
                        [1, 5, 0.09, 0.32, 45],
                        [2, 3, 0.05, 0.08, 100],
                        [3, 4, 0.04, 0.10, 77],
                        [4, 5, 0.04, 0.10, 97]])
#                        Bus Pga(MW) Pgb(MW) Pgc(MW) Qga(Mvar) Qgb(Mvar) Qgc(Mvar) Pca(MW) Pcb(MW) Pcc(MW) Qca(Mvar) Qcb(Mvar) Qcc(Mvar)      
data_difuse = np.matrix([[1, 999, 999, 999, 999, 999, 999, 16, 20, 23, 3, 5, 7],
                         [2, 0, 0, 0, 0, 0, 0, 55, 60, 66, 20, 25, 30],
                         [3, 0, 0, 0, 0, 0, 0, 37, 40, 45, 12, 15, 17],
                         [4, 130, 135, 140, 999, 999, 999, 23, 30, 32, 12, 15, 18],
                         [5, 0, 0, 0, 0, 0, 0, 95, 100, 106, 35, 40, 46]])
#teta e P->PQ e PV
#v e Q->PQ
Mbus = _iniciateBuses(data_bus)
Mlines = _flines(Mbus, data_lines)
print(f" *****Fuzzy Power Flow DC****** \n Study Case: {Mbus.tbus} buses, {Mlines.nlines} lines and Sb = {Sb} MW .\n ")
Mdif = _mZdif(data_difuse, Mbus, Sb)
Msens = _matrixA(Mlines, Mbus, data_lines)
Mxdif = _mXdif(Mlines.mzauxlin, Mdif, Mbus.nP)
pdifu = _Pdifuse(Mdif.mzdif, Msens, Mbus.nP, Mlines.nlines)
np.savetxt('Outputs/zdif.csv', Mdif.mzdif, delimiter=",", header='[Zmin, Zcent, Zmax]', footer=Mbus.tpinj)
np.savetxt('Outputs/xdif.csv', Mxdif.mxdif, delimiter=",", header='[Xmin, Xcent, Xmax]', footer=Mbus.ttheta)
np.savetxt('Outputs/pijdif.csv', pdifu, delimiter=",", header='[Plmin, Plcent, Plmax]', footer=Msens.txtpij)
print("Open folder -> Outputs <-")
del(Mbus, Mlines, Mdif, Msens, Mxdif, pdifu)
del(data_bus, data_difuse, data_lines, Sb)
