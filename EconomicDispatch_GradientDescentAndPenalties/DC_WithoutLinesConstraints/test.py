# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 14:55:37 2020

@author: Ana Silva
"""

import numpy as np
from Functions.descgradient import _CalcGrad
from Functions.initialisevalues import _Iniatilise
#                    Gen_i  Bus_i  a(€/h)       b(€/MWh)       c(€/(MW)^2.h)   Emissions(kgCO2/MWh)     Pgmin(MW)   Pgmax(MW)       
data_gen = np.matrix([[1, 1, 800, 20, 0.06, 80, 30, 90],
                      [2, 1, 800, 16, 0.1, 200, 20, 60],
                      [3, 4, 800, 20, 0.08, 175, 40, 80], 
                      [4, 4, 800, 22, 0.11, 60, 30, 70]])
#                      Bus_i Bus_j  Smax(MVA)   
data_lines = np.matrix([[1, 2, 62],
                        [1, 3, 70],
                        [1, 5, 45], 
                        [2, 3, 100],
                        [3, 4, 77],
                        [4, 5, 97]])
class _DataPar:
    def __init__(self, pload, nlamb, learnrate):
        self.pload = pload
        self.nlamb = nlamb
        self.learnrate = learnrate

MPar = _DataPar(250, 0.95, 0.005)
print("******Economical Dispatch*******")
print(f"Case Study-> Load Power: {MPar.pload} MW Lambda: {MPar.nlamb} LearningRate: {MPar.learnrate}\n")
[auxcinit, ctext] = _Iniatilise(data_gen, MPar)
cGrad = _CalcGrad(MPar, data_gen, auxcinit)
print(f"PgValues(MW)\n {cGrad.mpgen}\n")
print(f"GradientValues\n {cGrad.mgradient}\n")
print(f"TOTAL PENALTY: {cGrad.totpen} €\n")    
np.savetxt('Outputs/Pgvalues.csv', cGrad.mpgen, delimiter=",", header="Load_Dispatch_(MW)", footer="Rows:"+ctext.tpg)
np.savetxt('Outputs/MGradient.csv', cGrad.mgradient, delimiter=",", header="Gradient_Matrix", footer="Rows:"+ctext.tpg)
np.savetxt('Outputs/Totalpenalty.csv', cGrad.mpenalis, delimiter=",", header="Penalties(€)", footer="Rows:"+ctext.tcst)
del(ctext.tpg, ctext.tcst, ctext)
del(cGrad.mpgen, cGrad.mgradient, cGrad.mpenalis, cGrad.totpen, cGrad)
del(auxcinit.mpgval, auxcinit.mpgconst, auxcinit.newload, auxcinit.ngen, auxcinit)
del(MPar.nlamb, MPar.pload, MPar.learnrate)
del(MPar, data_gen, data_lines)
