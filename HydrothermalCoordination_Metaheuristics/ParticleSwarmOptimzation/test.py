# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:08:44 2020

@author: Ana Silva
"""
import os
import numpy as np

from Functions.moveandupdate import  _move_updatebest
from Functions.evaluate import _pop_evaluate
from Functions.mutate import _pop_mutate
from Functions.best import _swarm_best

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
        
class _DataSwarm:
    def __init__(self, _wi, _wm, _wc):
        self._wi = _wi
        self._wm = _wm
        self._wc = _wc
     
MPar = _DataPop(0.005, 20, 1000, 1, 1, 1)
MHydTh = _DataHydTh(150000, 80000, 40000, 1000, 80)
MSwarm = _DataSwarm(0.5, 1, 1) 
mdata_afl = np.matrix([[35000], [60000], [20000], [10000], [50000], [10000]])*MPar.perafl
mdata_load = np.matrix([[70], [80], [130], [50], [70], [110]])*MPar.perload
volinic = MPar.pervolin* MHydTh.vinicial
nPeriod = np.size(mdata_afl, 0)
print("-------->> Hydro-Thermal Coordination<<----- ")
print(f"Sigma= {MPar.sigma}, Population= {MPar.pop}, Generations= {MPar.germax}, load= {MPar.perload}p.u., Start Volume= {MPar.pervolin}p.u., Affluence= {MPar.perafl}p.u.")
print(f"wi = {MSwarm._wi}, wm = {MSwarm._wm}, wc = {MSwarm._wc}")

for ncmg in range(230, 410, 10):
    print(f"\nMarginal Cost {ncmg} €/MW")
    pop_initial = np.zeros([nPeriod+1, MPar.pop])
    pop_ant = np.zeros([nPeriod+1, MPar.pop])
    for i in range(nPeriod):
        for j in range(MPar.pop):
            pop_initial[i, j] = np.random.uniform(0, 1, 1)
    [mvaux_init, mpaux_init, mcaux_init, mPopula]=_pop_evaluate(nPeriod, MPar, MHydTh, volinic, mdata_afl, pop_initial, mdata_load, ncmg)
    auxpant = _pop_mutate(MPar.pop, nPeriod, mPopula, MPar.sigma)
    for i in range(nPeriod):
        for j in range(MPar.pop):
            pop_ant[i, j] = auxpant[i, j]
    [mvaux_ant, mpaux_ant, mcaux_ant, fpop_ant]=_pop_evaluate(nPeriod, MPar, MHydTh, volinic, mdata_afl, pop_ant, mdata_load, ncmg)
    [caux_bi, caux_bg]= _swarm_best(mPopula, fpop_ant, nPeriod, MPar.pop, mpaux_init, mvaux_init, mcaux_init, mpaux_ant, mvaux_ant, mcaux_ant)
    #move and update
    cger =_move_updatebest(caux_bi, caux_bg, mPopula, MPar.pop, nPeriod, MPar.germax, MSwarm, fpop_ant, MPar, MHydTh, volinic, mdata_afl, mdata_load, ncmg)
    print("Exporting data...") 
    if(not os.path.exists("Outputs/"+str(ncmg))):
        os.mkdir("Outputs/" +str(ncmg))     
    np.savetxt('Outputs/'+str(ncmg)+'/nbest.csv', cger.gbest, delimiter=",", header="Position of the best individual", footer="NºGenerations="+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/bestcost.csv', cger.gcost, delimiter=",", header="Best cost found (€)", footer="NºGenerations="+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/bestcmg.csv', cger.gcmg, delimiter=",", header="Marginal cost of the best individual (€/MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/Phidr.csv', cger.gphidr, delimiter=",", header="Hydro Power of the best individual (MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/Pterm.csv', cger.gpterm, delimiter=",", header="Thermal Power of the best individual (MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/vturb.csv', cger.gvturb, delimiter=",", header="Pumped Volume of the best individual (m3)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/vsobr.csv', cger.gvsobr, delimiter=",", header="Remaining Volume of the best individual (m3)", footer="NºGenerations= "+str(MPar.germax))
    

del(cger.gbest, cger.gcost, cger.gcmg, cger.gphidr, cger.gpterm, cger.gvturb, cger.gvsobr, cger)
del(caux_bi, caux_bg, fpop_ant, mpaux_ant, mvaux_ant, mcaux_ant)
del(auxpant, mPopula, mpaux_init, mvaux_init, mcaux_init)
del(pop_initial, pop_ant, ncmg, mdata_afl, mdata_load, volinic, nPeriod)
del(MPar.germax, MPar.perafl, MPar.perload, MPar.pervolin, MPar.pop, MPar.sigma)
del(MHydTh.penaliz, MHydTh.ptermmax, MHydTh.vdispmax, MHydTh.vinicial, MHydTh.vturbmax)
del(MSwarm._wc, MSwarm._wm, MSwarm._wi)
del(MPar, MHydTh, MSwarm)
print("END")
