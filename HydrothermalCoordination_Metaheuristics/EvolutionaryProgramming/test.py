# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:08:44 2020

@author: Ana Silva
"""
import os
import numpy as np

from Functions.duplicate import _pop_duplicate
from Functions.mutate import _pop_mutate
from Functions.evaluate import _pop_evaluate
from Functions.select import _cSelect


class _DataPop:
    def __init__(self, sigma, pop, germax, perload, pervolin, perafl):
        self.sigma = sigma
        self.pop = pop
        self.germax = germax
        self.perload = perload
        self.pervolin = pervolin
        self.perafl = perafl

class _DataHydTh:
    def __init__(self, vdispmax, vturbmax, vinicial, penaliz, ptermmax):
        self.vdispmax = vdispmax
        self.vturbmax = vturbmax
        self.vinicial = vinicial
        self.penaliz = penaliz
        self.ptermmax = ptermmax

MPar = _DataPop(0.005, 20, 1000, 1, 1, 1)
MHydTh = _DataHydTh(150000, 80000, 40000, 1000, 80)
mdata_afl = np.matrix([[35000], [60000], [20000], [10000], [50000], [10000]])*MPar.perafl
mdata_load = np.matrix([[70], [80], [130], [50], [70], [110]])*MPar.perload
volinic = MPar.pervolin* MHydTh.vinicial
nPeriod = np.size(mdata_afl, 0)
pop_inicial = np.random.uniform(0, 1, (nPeriod, MPar.pop))
auxnewpop = np.zeros([nPeriod, MPar.pop])

print("-------->> Hydro-Thermal Coordination<<----- ")
print(f"Sigma= {MPar.sigma}, Population= {MPar.pop}, Generations= {MPar.germax}, load= {MPar.perload}p.u., Start Volume= {MPar.pervolin}p.u., Affluence= {MPar.perafl}p.u.")

ger_best = np.zeros([1, MPar.germax])
ger_bestcmg = np.zeros([1, MPar.germax])
ger_bestcost = np.zeros([1, MPar.germax])
ger_Phidr = np.zeros([nPeriod, MPar.germax])
ger_Pterm = np.zeros([nPeriod, MPar.germax])
ger_vturb = np.zeros([nPeriod, MPar.germax])
ger_vsobr = np.zeros([nPeriod, MPar.germax])

for ncmg in range(230, 410, 10):
    print("\n Marginal cost "+str(ncmg)+"€/MW")  
    for ger in range (MPar.germax):
        if(ger == 0):
            mPopula = pop_inicial
        else:
            mPopula = auxnewpop
            auxdup = _pop_duplicate(mPopula, MPar.pop, nPeriod)
            auxmut = _pop_mutate(MPar.pop, nPeriod, auxdup, MPar.sigma)
            [auxvol, auxpot, auxcost] = _pop_evaluate(nPeriod, MPar, MHydTh, volinic, mdata_afl, auxmut, mdata_load, ncmg)
            auxselect = _cSelect(MPar.pop, nPeriod, auxmut, auxcost.mcost)
            auxbest = auxselect.best
            auxnewpop = auxselect.newpop
        
            ger_best[0, ger] = auxbest
            ger_bestcost[0, ger] = auxselect.mcostm[nPeriod, auxbest]
            ger_bestcmg[0, ger] = auxcost.mcmg[0, auxbest]
            ger_Phidr[:, ger] = auxpot.phidr[:, auxbest]
            ger_Pterm[:, ger] = auxpot.pterm[:, auxbest]
            ger_vturb[:, ger] = auxvol.vturb[:, auxbest]
            ger_vsobr[:, ger] = auxvol.vsobr[:, auxbest]
    
    print("Exporting data...") 
    if(not os.path.exists("Outputs/"+str(ncmg))):
        os.mkdir("Outputs/" +str(ncmg))     
    np.savetxt('Outputs/'+str(ncmg)+'/nbest.csv', ger_best, delimiter=",", header="Position of the best individual", footer="NºGenerations="+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/bestcost.csv', ger_bestcost, delimiter=",", header="Best cost found (€)", footer="NºGenerations="+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/bestcmg.csv', ger_bestcmg, delimiter=",", header="Marginal cost of the best individual (€/MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/Phidr.csv', ger_Phidr, delimiter=",", header="Hydro Power of the best individual (MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/Pterm.csv', ger_Pterm, delimiter=",", header="Thermal Power of the best individual (MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/vturb.csv', ger_vturb, delimiter=",", header="Pumped Volume of the best individual (m3)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/vsobr.csv', ger_vsobr, delimiter=",", header="Remaining Volume of the best individual (m3)", footer="NºGenerations= "+str(MPar.germax))
        


del(ger_best, ger_bestcmg, ger_bestcost, ger_Phidr, ger_Pterm, ger_vturb, ger_vsobr)
del(mPopula, auxbest, auxcost, auxdup, auxmut, auxpot, auxselect, auxvol)
del(mdata_afl, mdata_load, pop_inicial, auxnewpop, volinic, ncmg, nPeriod, ger)
del(MPar.germax, MPar.perafl, MPar.perload, MPar.pervolin, MPar.pop, MPar.sigma)
del(MHydTh.penaliz, MHydTh.ptermmax, MHydTh.vdispmax, MHydTh.vinicial, MHydTh.vturbmax)
del(MPar, MHydTh)
print("END")
