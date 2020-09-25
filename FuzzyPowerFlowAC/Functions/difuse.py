# -*- coding: utf-8 -*-
"""
@author: Ana Silva
"""
import numpy as np
def _auxZiniciar(Mbus, data_difuse, data_bus, Sb):
    bus = Mbus.tbus
    nP = Mbus.nP
    nQ = Mbus.nQ
    busP = Mbus.mbusP
    busQ = Mbus.mbusQ
    nzcent = Mbus.nz
    zdif = np.zeros([nzcent, 3])
    zcent = np.zeros([nzcent, 1])
    auxdelta = np.zeros([nzcent, 3])
    textP = ""
    textQ = ""
    for i in range(bus):
        for j in range(nP):
            if(busP[0, j] == data_difuse[i, 0]):
                aux1 = np.matrix([data_difuse[i, 1], data_difuse[i, 2], data_difuse[i, 3]])
                aux2 = np.matrix([data_difuse[i, 9], data_difuse[i, 8], data_difuse[i, 7]])
                paux = np.subtract(aux1, aux2)/Sb
                zdif[j, :] = paux
                textP = textP+"_P"+str(int(busP[0, j]))
            if(busP[0, j] == data_bus[i, 0]):
                zcent[j, 0] = data_bus[i, 9]
        for j in range(nQ):
            if(busQ[0, j] == data_difuse[i, 0]):
                aux1 = np.matrix([data_difuse[i, 4], data_difuse[i, 5], data_difuse[i, 6]])
                aux2 = np.matrix([data_difuse[i, 12], data_difuse[i, 11], data_difuse[i, 10]])
                qaux = np.subtract(aux1, aux2)/Sb
                zdif[j+nP, :] = qaux
                textQ = textQ+"_Q" + str(int(busQ[0, j]))
            if(busQ[0, j] == data_bus[i, 0]):
                zcent[j+nP, 0] = data_bus[i, 10]
    for i in range(nzcent):
        auxdelta[i, :] = zdif[i, :] - zcent[i, 0]
    textZ = "Rows" + textP + textQ
    class _CZdif:
        def __init__(self, mzcent, mzdif, mdeltaz, txtz):
            self.mzcent = mzcent
            self.mzdif = mzdif
            self.mdeltaz = mdeltaz
            self.txtz = txtz
    cauxzdif = _CZdif(zcent, zdif, auxdelta, textZ)
    return cauxzdif

def _auxDeltaX(data_Jacob, deltaZ, nzcent):
    jacob_inv = np.linalg.inv((data_Jacob))
    jacob_inv_pos = np.zeros([nzcent, nzcent])
    jacob_inv_neg = np.zeros([nzcent, nzcent])
    for i in range(nzcent):
        for j in range(nzcent):
            if (jacob_inv[i, j] < 0):
                jacob_inv_neg[i, j] = jacob_inv[i, j]
            else:
                jacob_inv_pos[i, j] = jacob_inv[i, j]
    deltazmin = np.zeros([nzcent, 1])
    deltazmax = np.zeros([nzcent, 1])
    deltazcent = np.zeros([nzcent, 1])
    for k in range(nzcent):
        deltazmin[k, 0] = deltaZ[k, 0]
        deltazcent[k, 0] = deltaZ[k, 1]
        deltazmax[k, 0] = deltaZ[k, 2]
    deltaxmin = np.matmul(jacob_inv_neg, deltazmax) + np.matmul(jacob_inv_pos, deltazmin)
    deltaxcent = np.matmul(jacob_inv_neg, deltazcent) + np.matmul(jacob_inv_pos, deltazcent)
    deltaxmax = np.matmul(jacob_inv_neg, deltazmin) + np.matmul(jacob_inv_pos, deltazmax)
    deltax = np.concatenate((deltaxmin, deltaxcent, deltaxmax), axis=1)
    class _Deltax:
        def __init__(self, mdeltax, mjacobinv, mdzmin, mdzcent, mdzmax):
            self.mdeltax = mdeltax
            self.mjacobinv = mjacobinv
            self.mdzmin = mdzmin
            self.mdzcent = mdzcent
            self.mdzmax = mdzmax
    cdeltax = _Deltax(deltax, jacob_inv, deltazmin, deltazcent, deltazmax)
    return cdeltax

def _Xdifuse(data_bus, deltaX, cbus):
    nzcent = cbus.nz
    busP = cbus.mbusP
    busQ = cbus.mbusQ
    nP = cbus.nP
    nQ = cbus.nQ
    bus = cbus.tbus
    auxXctr = np.zeros([nzcent, 1])
    auxXdif = np.zeros([nzcent, 3])
    auxt = 0
    text_theta = ""
    text_vi = ""
    for i in range(nP):
        for j in range(bus):
            if(data_bus[j, 0] == int(busP[0, i])):
                auxXctr[auxt, 0] = data_bus[j, 8]
                auxt = auxt + 1
                text_theta = text_theta + "_theta" + str(int(busP[0, i]))
    for i in range(nQ):
        for j in range(bus):
            if(data_bus[j, 0] == int(busQ[0, i])):
                auxXctr[auxt, 0] = data_bus[j, 6]
                auxt = auxt + 1
                text_vi = text_vi + "_V" + str(int(busQ[0, i]))
    auxXdif = deltaX + auxXctr
    text_X = "Rows:" + text_theta + text_vi
    class _Cxdif:
        def __init__(self, mxctr, mxdif, txtx):
            self.mxctr = mxctr
            self.mxdif = mxdif
            self.txtx = txtx
    cxdif = _Cxdif(auxXctr, auxXdif, text_X)
    return cxdif

def _auxdeltaPQ(data_lines, data_bus, clines, cbus, inv_Jacob):
    lines = clines.nlines
    y_lines = clines.mylin
    nzcent = cbus.nz
    busP = cbus.mbusP
    busQ = cbus.mbusQ
    nP = cbus.nP
    nQ = cbus.nQ
    bus = cbus.tbus
    posx = 0
    posy = 0
    vx = 0
    vy = 0
    tetax = 0
    tetay = 0
    matrixG = np.zeros([bus, bus])
    matrixB = np.zeros([bus, bus])
    matrixderP = np.zeros([2*lines, nzcent])
    matrixderQ = np.zeros([2*lines, nzcent])
    matrixPlinescent = np.zeros([2*lines, 1])
    matrixQlinescent = np.zeros([2*lines, 1])
    posf = int(0)
    tetaik = 0
    textauxP = ""
    textauxQ = ""
    textauxp = ""
    textauxq = ""
    for i in range(lines):
        posx = int(data_lines[i, 0]) - 1
        posy = int(data_lines[i, 1]) - 1
        matrixG[posx, posy] = np.real(y_lines[posx, posy])
        matrixG[posy, posx] = np.real(y_lines[posy, posx])
        matrixB[posx, posy] = np.imag(y_lines[posx, posy])
        matrixB[posy, posx] = np.imag(y_lines[posy, posx])
        matrixPlinescent[i, 0] = data_lines[i, 5]
        matrixPlinescent[i+lines, 0] = data_lines[i, 6]
        matrixQlinescent[i, 0] = data_lines[i, 7]
        matrixQlinescent[i+lines, 0] = data_lines[i, 8]
        textauxP = textauxP + "_P" + str(int(posx+1)) +  str(int(posy+1))
        textauxQ = textauxQ + "_Q" + str(int(posx+1)) +  str(int(posy+1))
        textauxp = textauxp + "_P" + str(int(posx+1)) +  str(int(posy+1))
        textauxq = textauxq + "_Q" + str(int(posx+1)) +  str(int(posy+1))
        for j in range(bus):
            if(data_bus[j, 0] == posx+1):
                vx = data_bus[j, 6]
                tetax = (data_bus[j, 8]*np.pi)/180
            elif(data_bus[j, 0] == posy+1):
                vy = data_bus[j, 6]
                tetay = (data_bus[j, 8]*np.pi)/180
        tetaik = tetax - tetay
        tetaki = tetay - tetax
        for k in range(nP):
            if(int(busP[0, k]) == posx+1):
                posf = k
                #dPik/Dtetai->tetaik
                matrixderP[i, posf] = vx*vy*(-matrixG[posx, posy]*np.sin(tetaik) + matrixB[posx, posy]*np.cos(tetaik))
                #dPki/Dtetak->tetaki
                matrixderP[i+lines, posf] = -vx*vy*(-matrixG[posy, posx]*np.sin(tetaki) + matrixB[posy, posx]*np.cos(tetaki))
                #dQik/Dtetai->tetaik
                matrixderQ[i, posf] = vx*vy*(matrixG[posx, posy]*np.cos(tetaik) + matrixB[posx, posy]*np.sin(tetaik))
                #dQki/Dtetak->tetaki
                matrixderQ[i+lines, posf] = -vx*vy*(matrixG[posy, posx]*np.cos(tetaki) + matrixB[posy, posx]*np.sin(tetaki))
            elif(int(busP[0, k]) == posy+1):
                posf = k
                 #dPik/Dtetak->tetaik
                matrixderP[i, posf] = -vx*vy*(-matrixG[posx, posy]*np.sin(tetaik) + matrixB[posx, posy]*np.cos(tetaik))
                 #dPki/Dtetai->tetaki
                matrixderP[i+lines, posf] = vx*vy*(-matrixG[posy, posx]*np.sin(tetaki) + matrixB[posy, posx]*np.cos(tetaki))
                #dQik/Dtetak->tetaik
                matrixderQ[i, posf] = -vx*vy*(matrixG[posx, posy]*np.cos(tetaik) + matrixB[posx, posy]*np.sin(tetaik))
                #dQki/Dtetai->tetaki
                matrixderQ[i+lines, posf] = vx*vy*(matrixG[posy, posx]*np.cos(tetaki) + matrixB[posy, posx]*np.sin(tetaki))
        for k in range(nQ):
            if(int(busQ[0, k]) == posx+1):
                posf = k+nP
                #dPik/Dvi->tetaik
                matrixderP[i, posf] = -2*matrixG[posx, posy]*vx + vy*(matrixG[posx, posy]*np.cos(tetaik) + matrixB[posx, posy]*np.sin(tetaik))
                #dPki/Dvk->tetaki
                matrixderP[i+lines, posf] = vy*(matrixG[posy, posx]*np.cos(tetaki) + matrixB[posy, posx]*np.sin(tetaki))
                #dQik/Dvi->tetaik
                matrixderQ[i, posf] = 2*matrixB[posx, posy]*vx + vy*(matrixG[posx, posy]*np.sin(tetaik) - matrixB[posx, posy]*np.cos(tetaik))
                #dQki/Dvk->tetaki
                matrixderQ[i+lines, posf] = vy*(matrixG[posy, posx]*np.sin(tetaki) - matrixB[posy, posx]*np.cos(tetaki))
            elif(int(busQ[0, k]) == posy+1):
                posf = k+nP
                #dPik/Dvk->tetaik
                matrixderP[i, posf] = vx*(matrixG[posx, posy]*np.cos(tetaik) + matrixB[posx, posy]*np.sin(tetaik))
                #dPki/Dvi->tetaki
                matrixderP[i+lines, posf] = -2*matrixG[posy, posx]*vy + vx*(matrixG[posy, posx]*np.cos(tetaki) + matrixB[posy, posx]*np.sin(tetaki))
                #dQik/Dvk->tetaik
                matrixderQ[i, posf] = vx*(matrixG[posx, posy]*np.sin(tetaik) - matrixB[posx, posy]*np.cos(tetaik))
                #dQki/Dvi->tetaki
                matrixderQ[i+lines, posf] = 2*matrixB[posy, posx]*vy + vx*(matrixG[posy, posx]*np.sin(tetaki) - matrixB[posy, posx]*np.cos(tetaki))
        matrixSensP = np.matmul(matrixderP, inv_Jacob)
        matrixSensQ = np.matmul(matrixderQ, inv_Jacob)
    textauxP = "Rows:_Pik e_Pki_|Pik->" + textauxP
    textauxQ = "Rows:_Qik e_Qki_|Qik->" + textauxQ
    textauxp = "Rows:_Pik->" + textauxp
    textauxq = "Rows:_Qik->" + textauxq
    class _CDeltaPQ:
        def __init__(self, mPlinecent, msensp, mQlinecent, msensq, txtP, txtQ, txtp, txtq):
            self.mPlinecent = mPlinecent
            self.msensp = msensp
            self.mQlinecent = mQlinecent
            self.msensq = msensq
            self.txtP = txtP
            self.txtQ = txtQ
            self.txtp = txtp
            self.txtq = txtq
    cdeltapq = _CDeltaPQ(matrixPlinescent, matrixSensP, matrixQlinescent, matrixSensQ, textauxP, textauxQ, textauxp, textauxq)
    return cdeltapq

def _PQdifuse(cpq, nzcent, lines, cdx):
    pcent = cpq.mPlinecent
    sensP = cpq.msensp
    qcent = cpq.mQlinecent
    sensQ = cpq.msensq
    deltazmin = cdx.mdzmin
    deltazcent = cdx.mdzcent
    deltazmax = cdx.mdzmax
    sensPpos = np.zeros([2*lines, nzcent])
    sensPneg = np.zeros([2*lines, nzcent])
    sensQpos = np.zeros([2*lines, nzcent])
    sensQneg = np.zeros([2*lines, nzcent])
    for i in range(2*lines):
        for j in range(nzcent):
            if(sensP[i, j] < 0):
                sensPneg[i, j] = sensP[i, j]
            else:
                sensPpos[i, j] = sensP[i, j]
            if(sensQ[i,j] < 0):
                sensQneg[i, j] = sensQ[i, j]
            else:
                sensQpos[i, j] = sensQ[i, j]           
    deltaPmin = np.matmul(sensPpos, deltazmin) + np.matmul(sensPneg, deltazmax)
    deltaPcent = np.matmul(sensPpos, deltazcent) + np.matmul(sensPneg, deltazcent)
    deltaPmax = np.matmul(sensPpos, deltazmax) + np.matmul(sensPneg, deltazmin)
    deltaQmin = np.matmul(sensQpos, deltazmin)+ np.matmul(sensQneg, deltazmax)
    deltaQcent = np.matmul(sensQpos, deltazcent)+ np.matmul(sensQneg, deltazcent)
    deltaQmax = np.matmul(sensQpos, deltazmax)+ np.matmul(sensQneg, deltazmin)
    deltaP = np.concatenate((deltaPmin, deltaPcent, deltaPmax), axis=1)
    deltaQ = np.concatenate((deltaQmin, deltaQcent, deltaQmax), axis=1)
    pdif = deltaP + pcent
    qdif = deltaQ + qcent
    class _Cpqdif:
        def __init__(self, mdP, mpdif, mdQ, mqdif):
            self.mdP = mdP 
            self.mpdif = mpdif
            self.mdQ = mdQ
            self.mqdif = mqdif
    cpqdif = _Cpqdif(deltaP, pdif, deltaQ, qdif)
    return cpqdif
            
def _Lossdifuse(cdx, nzcent, lines, cpq, data_lines):
    deltazmin = cdx.mdzmin
    deltazcent = cdx.mdzcent
    deltazmax = cdx.mdzmax
    sensP = cpq.msensp
    sensQ = cpq.msensq
    mSenspij = np.zeros([lines, nzcent])
    mSensqij = np.zeros([lines, nzcent])
    mposSenspij = np.zeros([lines, nzcent])
    mnegSenspij = np.zeros([lines, nzcent])
    mposSensqij = np.zeros([lines, nzcent])
    mnegSensqij = np.zeros([lines, nzcent])
    mpcent = np.zeros([lines, 1])
    mqcent = np.zeros([lines, 1])
    for i in range(lines):
        for j in range(nzcent):
            mSenspij[i, j] = sensP[i, j] + sensP[i+lines, j]
            mSensqij[i, j] = sensQ[i, j] + sensQ[i+lines, j]
        for j in range(nzcent):
            if(mSenspij[i, j] < 0):
                mnegSenspij[i, j] = mSenspij[i, j]
            else:
                mposSenspij[i, j] = mSenspij[i, j]
            if(mSensqij[i, j] < 0):
                mnegSensqij[i, j] = mSensqij[i, j]
            else:
                mposSensqij[i, j] = mSensqij[i, j]
        mpcent[i, 0]= data_lines[i, 5]+ data_lines[i, 6]
        mqcent[i, 0]= data_lines[i, 7]+ data_lines[i, 8]
        
    deltapmin = np.matmul(mposSenspij, deltazmin)+np.matmul(mnegSenspij, deltazmax)
    deltapcent = np.matmul(mposSenspij, deltazcent)+np.matmul(mnegSenspij, deltazcent)
    deltapmax = np.matmul(mposSenspij, deltazmax)+np.matmul(mnegSenspij, deltazmin)
    deltaqmin = np.matmul(mposSensqij, deltazmin) + np.matmul(mnegSensqij, deltazmax)
    deltaqcent = np.matmul(mposSensqij, deltazcent) + np.matmul(mnegSensqij, deltazcent)
    deltaqmax = np.matmul(mposSensqij, deltazmax) + np.matmul(mnegSensqij, deltazmin)
    deltap = np.concatenate((deltapmin, deltapcent, deltapmax), axis=1)
    deltaq = np.concatenate((deltaqmin, deltaqcent, deltaqmax), axis=1)
    matrixp = deltap + mpcent
    matrixq = deltaq + mqcent
    class _Closses:
        def __init__(self, mdp, mdq, mpdif, mqdif):
            self.mdp = mdp
            self.mdq = mdq
            self.mpdif = mpdif
            self.mqdif = mqdif
    closses = _Closses(deltap, deltaq, matrixp, matrixq)
    return closses
    
