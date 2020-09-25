# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""

import numpy as np
from Functions.descgradient import _calcgrad
from Functions.initialisevalues import _iniatilise
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
[auxcinit, ctext] = _iniatilise(data_gen, MPar)
[ccost, cval] = _calcgrad(MPar, data_gen, auxcinit)
print(f"PgValues(MW)\n {cval.mpval}\n")
print(f"GradientValues\n {cval.mgrad}\n")
print(f"TOTAL PENALTY: {ccost.ntotalpen} €\n")

text_file = open("Outputs/totalcost.txt", "w")
text_file.write(f"Totalcost: {ccost.ntotalcost:.2f} €/ GenerationCost {ccost.cprod:.2f} € /Emissions {ccost.cemis:.2f} €/ Penalties {ccost.ntotalpen:.2f} €")
text_file.close()
np.savetxt('Outputs/Pgvalues.csv', cval.mpval, delimiter=",", header="Load_Dispatch_(MW)", footer="Rows:"+ctext.tpg)
np.savetxt('Outputs/MGradient.csv', cval.mgrad, delimiter=",", header="Gradient_Matrix", footer="Rows:"+ctext.tpg)
np.savetxt('Outputs/Totalpenalty.csv', ccost.mpen, delimiter=",", header="Penalties(€)", footer="Rows:"+ctext.tcst)

del(ctext.tpg, ctext.tcst, ctext)
del(cval.mpval, cval.mgrad, ccost.mpen, ccost.ntotalpen, ccost, cval)
del(auxcinit.mpgval, auxcinit.mpgconst, auxcinit.newload, auxcinit.ngen, auxcinit)
del(MPar.nlamb, MPar.pload, MPar.learnrate)
del(MPar, text_file, data_gen, data_lines)
