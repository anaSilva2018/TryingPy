# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:08:44 2020

@author: Ana Silva
"""
#import os
import numpy as np
from Functions.pso import  _PSO

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

_PSO(MPar, MHydTh, MSwarm, mdata_afl, mdata_load, volinic, nPeriod)

del(mdata_afl, mdata_load, volinic, nPeriod)
del(MPar.germax, MPar.perafl, MPar.perload, MPar.pervolin, MPar.pop, MPar.sigma)
del(MHydTh.penaliz, MHydTh.ptermmax, MHydTh.vdispmax, MHydTh.vinicial, MHydTh.vturbmax)
del(MSwarm._wc, MSwarm._wm, MSwarm._wi)
del(MPar, MHydTh, MSwarm)
print("END")
