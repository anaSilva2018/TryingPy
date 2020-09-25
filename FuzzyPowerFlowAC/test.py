# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np
from Functions.buses import _iniciateBuses
from Functions.lines import _flines
from Functions.difuse import _auxZiniciar, _auxDeltaX, _Xdifuse, _auxdeltaPQ, _PQdifuse, _Lossdifuse
#Sb=100MW
Sb = 100

#                    Bus Type Pg(MW) Pc(MW) Qg(Mvar) Qc(Mvar) V(p.u)     Vn(kV) teta(º)    Pinj(p.u.)   Qinj(p.u.)
data_bus = np.matrix([[1, 1, 999,    20,    999,     5,        1.01,      150,     0.0000,    1.0023,   0.3652],
                      [2, 2, 0,      60,    0,       25,       0.9797,    150,    -2.0600,   -0.6000,  -0.2500],
                      [3, 2, 0,      40,    0,       15,       0.9906,    150,    -1.8700,   -0.4000,  -0.1500],
                      [4, 3, 999,    30,    999,     15,       1.02,      150,    -0.7300,    1.0500,   0.5757],
                      [5, 2, 0,      100,   0,       40,       0.9539,    150,    -4.3900,   -1.0000,  -0.4000]])

#NOTA: 1->REF; 2->PQ;3->PV

#                      Busi Busj R(p.u.) X(p.u.) Smax(MVA)  Pij(p.u.)  Pji(p.u.) Qij(p.u.)  Qji(p.u.)           
data_lines = np.matrix([[1, 2, 0.03,    0.08,    62,       0.5190,   -0.5099,  0.1960,    -0.1723],
                        [1, 3, 0.05,    0.16,    70,       0.2210,   -0.2188,  0.0570,    -0.0484],
                        [1, 5, 0.09,    0.32,    45,       0.2619,   -0.2550,  0.1122,    -0.0870], 
                        [2, 3, 0.05,    0.08,   100,      -0.0900,    0.0908, -0.0780,     0.0788],
                        [3, 4, 0.04,    0.10,    77,      -0.2720,    0.2760, -0.1805,     0.1910],
                        [4, 5, 0.04,    0.10,    97,       0.7740,   -0.7450,  0.3844,    -0.3130]])

#                        Bus Pga(MW) Pgb(MW) Pgc(MW) Qga(Mvar) Qgb(Mvar) Qgc(Mvar) Pca(MW) Pcb(MW) Pcc(MW) Qca(Mvar) Qcb(Mvar) Qcc(Mvar)             
data_difuse = np.matrix([[1, 999,   999,    999,    999,     999,       999,       16,    20,     23,        3,      5,       7],
                         [2, 0,     0,      0,      0,       0,         0,         55,    60,     66,        20,     25,      30],
                         [3, 0,     0,      0,      0,       0,         0,         37,    40,     45,        12,     15,      17], 
                         [4, 130,   135,    140,    999,     999,       999,       23,    30,     32,        12,     15,      18],
                         [5, 0,     0,      0,      0,       0,         0,         95,    100,    106,       35,     40,      46]])


                     # teta2(º) teta3(º)   teta4(º)   teta5(º)   V2(p.u)     V3(p.u.)     V5(p.u.)
data_Jacob = np.matrix([[19.39,  -8.70,     0.00,      0.00,    8.92,      -5.53,     0.00],  # P2
                        [-8.74,  23.02,    -8.64,      0.00,   -5.53,      10.34,     0.00],  # P3
                        [ 0.00,  -8.78,    17.36,     -8.58,    0.00,      -3.34,    -2.95],  # P4
                        [ 0.00,   0.00,    -8.16,     10.88,    0.00,       0.00,     3.02],  # P5
                        [-9.94,   5.48,     0.00,      0.00,   19.29,      -8.79,     0.00],  # Q2
                        [ 5.42, -11.04,     3.66,      0.00,   -8.92,      22.93,     0.00],  # Q3
                        [ 0.00,   0.00,     3.88,     -4.88,    0.00,       0.00,    10.57]]) # Q5

#teta e P->PQ e PV
#v e Q->PQ
Mbus = _iniciateBuses(data_bus)
Mlines = _flines(Mbus, data_lines)
print(f" *****Fuzzy Power Flow AC****** \n Study Case: {Mbus.tbus} buses, {Mlines.nlines} lines and Sb = {Sb} MW .\n ")
Mzdif = _auxZiniciar(Mbus, data_difuse, data_bus, Sb)
Mdeltax = _auxDeltaX(data_Jacob, Mzdif.mdeltaz, Mbus.nz)
Mxdif =  _Xdifuse(data_bus, Mzdif.mdeltaz, Mbus)
Mdeltapq = _auxdeltaPQ(data_lines, data_bus, Mlines, Mbus, Mdeltax.mjacobinv)
Mpqdif = _PQdifuse(Mdeltapq, Mbus.nz, Mlines.nlines, Mdeltax)
Mlosses = _Lossdifuse(Mdeltax, Mbus.nz, Mlines.nlines, Mdeltapq, data_lines)

np.savetxt('Outputs/Zdif.csv', Mzdif.mzdif, delimiter=",", header='[Zmin, Zcent, Zmax]', footer=Mzdif.txtz)
np.savetxt('Outputs/Xdif.csv', Mxdif.mxdif, delimiter=",", header='[Xmin, Xcent, Xmax]', footer=Mxdif.txtx)
np.savetxt('Outputs/Pijdif.csv', Mpqdif.mpdif, delimiter=",", header='[Plmin, Plcent, Plmax]', footer=Mdeltapq.txtP)
np.savetxt('Outputs/Qijdif.csv', Mpqdif.mqdif, delimiter=",", header='[Qlmin, Qlcent, Qlmax]', footer=Mdeltapq.txtQ)
np.savetxt('Outputs/Plossesdif.csv', Mlosses.mpdif, delimiter=",", header='[plmin, plcent, plmax]', footer=Mdeltapq.txtp)
np.savetxt('Outputs/Qlossesdif.csv', Mlosses.mqdif, delimiter=",", header='[qlmin, qlcent, qlmax]', footer=Mdeltapq.txtq)
print("Open folder -> Outputs <-\n")
