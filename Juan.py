import random
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import complex as cp

N = 100
Centroides = [1+1j,1-1j,-1+1j,-1-1j]
epsilon = 5.e-1  
nb_centres = 4
SNR = 10
        
def Estime_h(snr,N):
    m1 = complex(1,1)
    m2 = complex(-1,1)
    m3 = complex(-1,-1)
    m4 = complex(1,-1)
    mk = [m1,m2,m3,m4]
    R = random.random()
    theta = np.pi*random.random()/2.
    H = R*np.exp(complex(0,theta))
    mi = np.array([complex(0,0) for i in range(N)])
    bruit = np.array([complex(0,0) for i in range(N)])
    for i in range(N): #bruits gaussiens
        mi[i] = mk[random.randint(0,3)]
        bruit[i]= random.gauss(0,1) + complex(0,1)*random.gauss(0,1)
    #ajustement du snr
    Ps = np.absolute(H)**2*(mi.dot(mi)/N)
    Pb = bruit.dot(bruit)/N
    alpha = np.sqrt(Ps/Pb/snr)
    Data = H*mi+alpha*bruit    

    def associe_centroide(Point,centroides):
        n0 = len(centroides)
        indice_centroide_PP = 0
        dist = np.linalg.norm(Point-centroides[0])**2
        for i in range(n0):
            if dist > np.linalg.norm(Point-centroides[i])**2:
                dist = np.linalg.norm(Point-centroides[i])**2
                indice_centroide_PP = i
        return indice_centroide_PP
            
    def rempli_T(data,centroides):
        t = np.zeros((N,nb_centres))
        for i in range(N):
            for j in range(nb_centres):
                if j == associe_centroide(data[i],centroides):
                    t[i][j] = 1
                    break
        return t
    
    def L(data,h):
        L = cp(0,0)
        centr = h*np.array(Centroides)
        TT = rempli_T(A,centr)    
        for k in range(N):
            for j in range(nb_centres):
                L += TT[k][j]*np.linalg.norm(data[k] - centr[j])**2
        return 0.5*L
        
    def jac_L(data,r,phi):
        jac = np.zeros(2)
        centr = r*np.exp(1j*phi)*np.array(Centroides)
        
        T = rempli_T(data,centr)
        for k in range(N):
            for j in range(nb_centres):
                if (T[k][j] ==1):   
                    X = data[k]
                    jac[0] += 2*r + np.real(X)*(np.sin(phi)*np.imag(Centroides[j]) - np.cos(phi)*np.real(Centroides[j])) -np.imag(X)*(np.cos(phi)*np.imag(Centroides[j])+np.sin(phi)*np.real(Centroides[j]))
                    
                    jac[1] += r*(np.sin(phi)*np.real(Centroides[j])*np.real(X) + np.cos(phi)*np.imag(Centroides[j])*np.real(X) + np.sin(phi)*np.imag(Centroides[j])*np.imag(X) - np.cos(phi)*np.real(Centroides[j])*np.imag(X))
                    break
        return 0.5*jac 
    
    def descente_grad(H,data,h_0,err,pas,maxiter=200):
        step = pas
        iterations = 0
    
        r_old = np.linalg.norm(h_0)
        phi_old = np.arctan(np.imag(h_0)/np.real(h_0))
        if np.real(h_0) < 0:
            phi_old += np.pi
    
        while (np.linalg.norm(jac_L(data,r_old,phi_old)) > err and iterations < maxiter):
            jac = jac_L(data,r_old,phi_old)
            r_old = r_old - step*jac[0]
            phi_old = phi_old - step*jac[1]
            #print(np.linalg.norm(jac_L(data,r_old,phi_old)))
            iterations += 1
            
        h = r_old*np.exp(1j*phi_old)
        err = np.linalg.norm(h-H)
        
        hh = h
        for i in range(1,4):
            if err > np.linalg.norm((1j)**i*hh - H):
                h = (1j)**i*hh #erreur à pi/2 pres
                err = np.linalg.norm(h - H)
                
        return [h,iterations] # h_new correspond à 1/H

    h_0 = random.random()*np.exp(1j*2*np.pi*random.random())
    h = descente_grad(H,Data,h_0,epsilon,0.01)[0]
    return 100*np.linalg.norm(h-H)/np.linalg.norm(H)