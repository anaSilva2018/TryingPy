# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import os
import numpy as np

from Functions.epso import  _swarminit, _move_updatebest
class _DataPop:
    def __init__(self, sigma, pop, germax, perload, pervolin, perafl, tau):
        self.sigma = sigma
        self.pop = pop
        self.germax = germax
        self.perload = perload
        self.pervolin = pervolin
        self.perafl = perafl
        self.tau = tau
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
        
MPar = _DataPop(0.01, 20, 1000, 1, 1, 1, 0.001)
MHydTh = _DataHydTh(150000, 80000, 40000, 1000, 80)
MSwarm = _DataSwarm(0.8, 0.4, 0.4) 
mdata_afl = np.matrix([[35000], [60000], [20000], [10000], [50000], [10000]])*MPar.perafl
mdata_load = np.matrix([[70], [80], [130], [50], [70], [110]])*MPar.perload
volinic = MPar.pervolin* MHydTh.vinicial
nPeriod = np.size(mdata_afl, 0)
print("-------->> Hydro-Thermal Coordination<<----- ")
print(f"Sigma= {MPar.sigma}, Tau = {MPar.tau}, Population= {MPar.pop}, Generations= {MPar.germax}, load= {MPar.perload}p.u., Start Volume= {MPar.pervolin}p.u., Affluence= {MPar.perafl}p.u.")
print(f"wi = {MSwarm._wi}, wm = {MSwarm._wm}, wc = {MSwarm._wc}")

for ncmg in range(230, 410, 10):
    print(f"\nMarginal Cost {ncmg} €/MW")
    [mcwic, bpop, bpot, bvol, mpop_ant, mpopi] = _swarminit(MPar, nPeriod, MHydTh, volinic, mdata_afl, mdata_load, ncmg, MSwarm)
    [cgerp, cptger, cvger] = _move_updatebest(MPar, nPeriod, mpopi, mpop_ant, mcwic, bpop, bpot, bvol, ncmg, mdata_load, mdata_afl, volinic, MHydTh)
   
    print("Exporting data...") 
    if(not os.path.exists("Outputs/"+str(ncmg))):
        os.mkdir("Outputs/" +str(ncmg))     
    np.savetxt('Outputs/'+str(ncmg)+'/bestcost.csv', cgerp.mcost, delimiter=",", header="Best cost found (€)", footer="NºGenerations="+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/bestcmg.csv', cgerp.mcmg, delimiter=",", header="Marginal cost of the best individual (€/MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/Phidr.csv', cptger.phidr, delimiter=",", header="Hydro Power of the best individual (MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/Pterm.csv', cptger.pterm, delimiter=",", header="Thermal Power of the best individual (MW)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/vturb.csv', cvger.vturb, delimiter=",", header="Pumped Volume of the best individual (m3)", footer="NºGenerations= "+str(MPar.germax))
    np.savetxt('Outputs/'+str(ncmg)+'/vsobr.csv', cvger.vsobr, delimiter=",", header="Remaining Volume of the best individual (m3)", footer="NºGenerations= "+str(MPar.germax))

del(cgerp, cptger, cvger, ncmg)
del(mcwic, bpop, bpot, bvol, mpop_ant, mpopi)
del(nPeriod, volinic, mdata_afl, mdata_load)
del(MSwarm._wi, MSwarm._wc, MSwarm._wm, MSwarm)
del(MHydTh.penaliz, MHydTh.ptermmax, MHydTh.vdispmax, MHydTh.vinicial, MHydTh.vturbmax, MHydTh)
del(MPar.germax, MPar.perafl, MPar.perload, MPar.pervolin, MPar.pop, MPar.sigma, MPar.tau, MPar)

print("END")
