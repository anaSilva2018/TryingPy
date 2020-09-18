# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 23:51:11 2020

@author: Ana Silva
"""
import numpy as np
from Functions.evaluate import _pop_evaluate
def _move_updatebest(cbi, cbg, mpop, pop, nper, ngermax, cSwarm, mpop_ant, cpar, chydth, nvolinic, mafl, mload,valdef):
   
    cbi_pop = cbi.popbest
    cbi_cmg = cbi.cmgbest
    cbi_phidr = cbi.phidrbest
    cbi_pterm = cbi.ptermbest
    cbi_vdisp = cbi.vdispbest
    cbi_vturb = cbi.vturbbest
    cbi_vsobr = cbi.vsobrbest
    cbg_nbg = cbg.nbglobal
    cbg_bgcmg = cbg.ncmgbglobal
    cbf_bgcost = cbg.ncostbglobal
    ger_best = np.zeros([1, cpar.germax])
    ger_bestcmg = np.zeros([1, cpar.germax])
    ger_bestcost = np.zeros([1, cpar.germax])
    ger_Phidr = np.zeros([nper, cpar.germax])
    ger_Pterm = np.zeros([nper, cpar.germax])
    ger_vturb = np.zeros([nper, cpar.germax])
    ger_vsobr = np.zeros([nper, cpar.germax])
    
    ger_best[0, 0] = cbg_nbg
    ger_bestcost[0, 0] = cbf_bgcost
    ger_bestcmg[0, 0] = cbg_bgcmg
    ger_Phidr[:, 0] = cbi_phidr[:, cbg_nbg]
    ger_Pterm[:, 0] = cbi_pterm[:, cbg_nbg]
    ger_vturb[:, 0] = cbi_vturb[:, cbg_nbg]
    ger_vsobr[:, 0] = cbi_vsobr[:, cbg_nbg]
    
    mPopul = mpop
    
    for ger in range(0,ngermax-1, 1):
        auxnewpop = np.zeros([nper+1, pop])
        for j in range(pop):
            for i in range(nper):
                nDEC = (1/(ger+2))
                auxnewpop[i, j] = nDEC*cSwarm._wi*(mPopul[i, j]-mpop_ant[i, j])
                auxnewpop[i, j] = auxnewpop[i, j]+ np.random.uniform(0, 1, 1)*cSwarm._wm*(cbi_pop[i, j]-mPopul[i, j])
                auxnewpop[i, j] = auxnewpop[i, j]+ np.random.uniform(0, 1, 1)*cSwarm._wc*(cbi_pop[i, cbg_nbg]-mPopul[i, j]) 
                if(auxnewpop[i, j] > 1):
                    auxnewpop[i, j] = 1
                elif(auxnewpop[i, j] < 0):
                    auxnewpop[i, j] = 0
        mPopul  = auxnewpop
        [mvaux_new, mpaux_new, mcaux_new,mPopul]=_pop_evaluate(nper, cpar, chydth, nvolinic, mafl, mPopul, mload, valdef)
        for j in range(pop):
            if(mPopul[nper, j] < cbi_pop[nper, j]):
                cbi_pop[nper, j] = mPopul[nper, j]
                cbi_cmg[0, j] =  mcaux_new.mcmg[0, j]
                for i in range(nper):
                    cbi_pop[i, j] = mPopul[i, j]
                    cbi_phidr[i, j] = mpaux_new.phidr[i, j]
                    cbi_pterm[i, j] = mpaux_new.pterm[i, j]
                    cbi_vturb[i, j] = mvaux_new.vturb[i, j]
                    cbi_vdisp[i, j] = mvaux_new.vdisp[i, j]
                    cbi_vsobr[i, j] = mvaux_new.vsobr[i, j]        
        for j in range(pop):
            if(mPopul[nper, j] < cbf_bgcost):
                cbg_nbg  = j
                cbg_bgcmg = cbi_cmg[0, j]
                cbf_bgcost = cbi_pop[nper, j]
        ger_best[0, ger+1] = cbg_nbg
        ger_bestcost[0, ger+1] = cbf_bgcost
        ger_bestcmg[0, ger+1] = cbg_bgcmg
        ger_Phidr[:, ger+1] = cbi_phidr[:, cbg_nbg]
        ger_Pterm[:, ger+1] = cbi_pterm[:, cbg_nbg]
        ger_vturb[:, ger+1] = cbi_vturb[:, cbg_nbg]
        ger_vsobr[:, ger+1] = cbi_vsobr[:, cbg_nbg]
              
    class _cGerb:
        def __init__(self, gbest, gcost, gcmg, gphidr, gpterm, gvturb, gvsobr):
            self.gbest =  gbest
            self.gcost =  gcost
            self.gcmg =  gcmg
            self.gphidr = gphidr
            self.gpterm = gpterm
            self.gvturb = gvturb
            self.gvsobr = gvsobr
   
    cnew_ger = _cGerb(ger_best, ger_bestcost, ger_bestcmg, ger_Phidr, ger_Pterm, ger_vturb, ger_vsobr)
    
    return  cnew_ger
    
