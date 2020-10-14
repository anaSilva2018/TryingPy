 -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np
def _pop_evaluate(nper, cinic, chydth, nvolinic, mafl, mmutate, mload, valdef):
    pop = cinic.pop
    khidr = (chydth.vturbmax*(10/3))/3600
    vherd = np.zeros([nper, pop])
    vdisp = np.zeros([nper, pop])
    vturb = np.zeros([nper, pop])
    vsobr = np.zeros([nper, pop])
    vdesc = np.zeros([nper, pop])
    psobr = np.zeros([nper, pop])
    phidr = np.zeros([nper, pop])
    pterm = np.zeros([nper, pop])
    mcmg = np.zeros([1, pop])
    mpen = np.zeros([nper, pop])
    mval = np.zeros([nper, pop])
    mcost = np.zeros([nper+1, pop])
    mpop = np.zeros([nper+1, pop])
    for i in range(nper):
        for j in range(pop):
            mpop[i, j] = mmutate[i, j]
            if(i == 0):
                vherd[i, j] = nvolinic
            else:
                vherd[i, j] = vsobr[i-1, j]
            vdisp[i, j] = vherd[i, j]+mafl[i, 0]
            if(vdisp[i, j] > chydth.vdispmax):
                vdisp[i, j] = chydth.vdispmax
                vdesc[i, j] = vdisp[i, j]-chydth.vdispmax
            vturb[i, j] = mmutate[i, j]*vdisp[i, j]
            vsobr[i, j] = vdisp[i, j]-vturb[i, j]
            psobr[i, j] = (10/3)*(vsobr[i, j]/3600)
            phidr[i, j] = (10/3)*(vturb[i, j]/3600)
            pterm[i, j] = mload[i, 0]- phidr[i, j]
    for i in range(nper):
        for j in range(pop):
            if(phidr[i, j] > khidr):
                mpen[i, j] = (phidr[i, j]-khidr)+mpen[i, j]
            if(phidr[i, j] < 0):
                mpen[i, j] = (phidr[i, j])+mpen[i, j]
            if(pterm[i, j] > chydth.ptermmax):
                mpen[i, j] = (pterm[i, j]-chydth.ptermmax)+mpen[i, j]
            if(pterm[i, j] < 0):
                mpen[i, j] = pterm[i, j]+mpen[i, j]
            if(psobr[nper-1, j] > 0 and np.all(pterm[:, j] <= chydth.ptermmax) and np.all(phidr[:, j] <= khidr)):
                mcmg[0, j] = valdef
                mval[nper-1, j] = psobr[nper-1, j]
            mcost[i, j] = 2000+100*pterm[i, j]+1.5*np.power(pterm[i, j], 2)+chydth.penaliz*np.power(mpen[i, j], 2)-mcmg[0, j]*mval[i, j]
            mcost[nper, j] = mcost[i, j]+mcost[nper, j]
            mpop[nper, j] = mcost[nper, j]
            
    class _PopVol:
        def __init__(self, vherd, vdisp, vturb, vsobr, vdesc):
            self.vherd = vherd
            self.vdisp = vdisp
            self.vturb = vturb
            self.vsobr = vsobr
            self.vdesc = vdesc
            
    class _PopPot:
        def __init__(self, psobr, phidr, pterm):
            self.psobr = psobr 
            self.phidr = phidr
            self.pterm = pterm
            
    class _Ccost:
        def __init__(self, mcmg, mpen, mval, mcost):
            self.mcmg = mcmg
            self.mpen = mpen
            self.mval = mval
            self.mcost = mcost
            
    p_vol = _PopVol(vherd, vdisp, vturb, vsobr, vdesc)
    p_pot = _PopPot(psobr, phidr, pterm)
    p_cost = _Ccost(mcmg, mpen, mval, mcost)  
    return p_vol, p_pot, p_cost, mpop        

def _move_updatebest(cbi, cbg, mpop, pop, nper, ngermax, cSwarm, mpop_ant, cpar, chydth, nvolinic, mafl, mload, valdef):
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
    
    w = cSwarm._wm+cSwarm._wc
    k = 2/np.abs(2-w-np.sqrt(np.power(w, 2)-4*w))
    
    mPopul = mpop
    for ger in range(0, ngermax-1, 1):
        auxnewpop = np.zeros([nper+1, pop])
        for j in range(pop):
            for i in range(nper):
                nDEC = (1/(ger+2))
                auxnewpop[i, j] = nDEC*cSwarm._wi*(mPopul[i, j]-mpop_ant[i, j])
                auxnewpop[i, j] = auxnewpop[i, j]+ np.random.uniform(0, 1, 1)*cSwarm._wm*(cbi_pop[i, j]-mPopul[i, j])
                auxnewpop[i, j] = auxnewpop[i, j]+ np.random.uniform(0, 1, 1)*cSwarm._wc*(cbi_pop[i, cbg_nbg]-mPopul[i, j]) 
                auxnewpop[i, j] = auxnewpop[i, j]*k
                auxnewpop[i, j] = auxnewpop[i, j]+mPopul[i, j]
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
    
