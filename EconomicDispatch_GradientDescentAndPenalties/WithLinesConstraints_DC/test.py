# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 23:58:40 2020

@author: Ana Silva
"""
import numpy as np
from Functions.buses import _iniciateBuses
from Functions.lines import _initlines
from Functions.fauxgrad import  _fgradpenext

#                    Gen_i  Bus_i  a(€/h)       b(€/MWh)       c(€/(MW)^2.h)   Emissions(kgCO2/MWh)     Pgmin(MW)   Pgmax(MW)       
data_gen = np.matrix([[1, 1, 800, 20, 0.06, 80, 30, 90],
                      [2, 1, 800, 16, 0.1, 200, 20, 60],
                      [3, 4, 800, 20, 0.08, 175, 40, 80], 
                      [4, 4, 800, 22, 0.11, 60, 30, 70]])
#                      Busi Busj R(p.u.) X(p.u.) Smax(MVA)        
data_lines = np.matrix([[1, 2, 0.03, 0.08, 62],
                        [1, 3, 0.05, 0.16, 70],
                        [1, 5, 0.09, 0.32, 45],
                        [2, 3, 0.05, 0.08, 100],
                        [3, 4, 0.04, 0.10, 77],
                        [4, 5, 0.04, 0.10, 97]])
#                     Bus Type Pg(MW) Pc(MW) Qg(Mvar) Qc(Mvar)      Vn(kV) 
data_bus = np.matrix([[1, 1, 999, 20, 999, 5, 150],
                      [2, 2, 0, 60, 0, 25, 150],
                      [3, 2, 0, 40, 0, 15, 150],
                      [4, 3, 999, 30, 999, 15, 150],
                      [5, 2, 0, 100, 0, 40, 150]])
class _DataPar:
    def __init__(self, pload, nlamb, learnrate, nsb):
        self.pload = pload
        self.nlamb = nlamb
        self.learnrate = learnrate
        self.nsb = nsb
MPar = _DataPar(250, 0.95, 0.75, 100)
print("******Economical Dispatch*******")
print(f"Case Study-> Load Power: {MPar.pload} MW Lambda: {MPar.nlamb} LearningRate: {MPar.learnrate}\n")
[Mbus, Mtxt] = _iniciateBuses(data_bus, data_gen, data_lines)
Mlines = _initlines(Mbus.nbus, data_lines, Mbus.nbusref, Mbus.nP)
[Ccost, Cval] = _fgradpenext(MPar, Mbus, Mlines, data_gen, data_bus, data_lines, 0.0001, 0.000002)

print(f"PgValues(MW)\n {Cval.mpval}\n")
print(f"PijValues(MW)\n {Cval.mpij}\n")
print(f"GradientValues\n {Cval.mgrad}\n")
print(f"TOTAL PENALTY: {Ccost.ntotalpen} €\n")  
text_file = open("Outputs/totalcost.txt", "w")
text_file.write(f"Totalcost: {Ccost.ntotalcost:.2f} €/ GenerationCost {Ccost.cprod:.2f} € /Emissions {Ccost.cemis:.2f} €/ Penalties {Ccost.ntotalpen:.2f} €")
text_file.close()
np.savetxt('Outputs/Pgvalues.csv', Cval.mpval, delimiter=",", header="Load_Dispatch_(MW)", footer="Rows:"+Mtxt.tpg)
np.savetxt('Outputs/Pijvalues.csv', Cval.mpij, delimiter=",", header="Power_FLow_(MW)", footer="Rows:"+Mtxt.tpij)
np.savetxt('Outputs/MGradient.csv', Cval.mgrad, delimiter=",", header="Gradient_Matrix", footer="Rows:"+Mtxt.tpg)
np.savetxt('Outputs/Totalpenalty.csv', Ccost.mpen, delimiter=",", header="Penalties(€)", footer="Rows:"+Mtxt.tconst)
del(Cval, Ccost)
del(Mbus, Mlines, MPar.pload, MPar.nlamb, MPar.nsb, MPar.learnrate)
del(MPar, Mtxt, text_file, data_bus, data_lines, data_gen)
