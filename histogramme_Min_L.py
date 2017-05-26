import random
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import complex as cp

N = 100
sigma = .3
nb_centres = 4
epsilon = 1.e-1
Centroides = [1+1j,1-1j,-1+1j,-1-1j]  
hist = []        

for i in range(600):
    print(i*1./6,'%')
    A = (1+1j)*np.zeros(N)
    B = (1+1j)*np.zeros(N)
    
    # generation des données EMETTEUR
    for k in range(N):
        x = np.sign(random.gauss(0,1))
        y = np.sign(random.gauss(0,1))
        A[k] = cp(x,y)
        B[k] = cp(x,y)
    
    # genération du H: on suppose le module de H inférieur à 1
    
    #pour sigma = 0.3 99% des tirages sont inférieurs à 0.9
    H = random.random()*np.exp(1j*2*np.pi*random.random())
    norme = np.linalg.norm(H)
    
    #generation des données: à peu pres réparties sur un carré avec bruit
    #la transformation par H
    A = H*A
    B = H*B
    #le rajout du bruit gaussien
    for k in range(N):
        A[k] += cp(random.gauss(0,sigma),random.gauss(0,sigma))
    
    #je rajoute au signal reçu du bruit gaussien un grand nombre de fois pour eliminer le bruit -> il faut d'abord avoir une estimation de la variance du bruit b qui s est ajouté 
    #alpha = 10
    #sigma_tilde = alpha*sigma #disons que l'on se trompe d'un facteur alpha
    
    #Preconditionnement
    # for i in range(100):
    #     for k in range(N):
    #         A[k] += cp(random.gauss(0,sigma),random.gauss(0,sigma)) + B[k]
    # for k in range(N):
    #     A[k] = A[k]/100
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
    h = descente_grad(H,A,h_0,epsilon,0.01)[0]
    hist += [100*np.linalg.norm(h-H)/np.linalg.norm(H)]

plt.subplot(2,1,1)
eff, val, patches = plt.hist(histogramme2, range = (0, 10), bins = 30, edgecolor = 'black', normed=False)
mu = np.mean(histogramme2)
sigma = np.std(histogramme2)
max = np.max(eff)
plt.plot([mu,mu],[0,max],color = 'r', ls=':', linewidth = 3, label = '$ \sigma = 0.01$')
plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, max), xytext=(mu+0.2, max-2), color="red")
plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, max*3/4), xytext=(mu+sigma+0.2, max*3/4), color="green")
#plt.plot([mu-sigma,mu-sigma],[0,0.75*max], color = 'g', ls=':', linewidth = 3)
plt.plot([mu+sigma,mu+sigma],[0,0.75*max], color = 'g', ls=':', linewidth = 3)
plt.ylabel('Frequency of occurency')
plt.title('Minimization of cost function $\mathcal{L}(h)$ Algorithm')

plt.legend(loc = 1)

plt.subplot(2,1,2)
eff2, val2, patches2 = plt.hist(histogramme, range = (0, 10), bins = 30, edgecolor = 'black', normed=False)
mu2 = np.mean(histogramme)
sigma2 = np.std(histogramme)
max2 = np.max(eff2)
plt.plot([mu2,mu2],[0,max2],color = 'r', ls=':', linewidth = 3, label = '$ \sigma = 0.1$')
plt.annotate(r'$\mu_2$ = '+str(round(mu2, 2))+'%', xy=(mu2, max2), xytext=(mu2+0.2, max2-5), color="red")
plt.annotate(r'$\sigma_2$ = '+str(round(sigma, 2))+'%', xy=(mu2+sigma, max2*3/4), xytext=(mu2+sigma2+0.2, max2*3/4), color="green")
plt.plot([mu2-sigma2,mu2-sigma2],[0,0.75*max2], color = 'g', ls=':', linewidth = 3)
plt.plot([mu2+sigma2,mu2+sigma2],[0,0.75*max2], color = 'g', ls=':', linewidth = 3)
plt.legend(loc = 1)
plt.xlabel('Error on $H$ (%)')
plt.ylabel('Frequency of occurency')

plt.savefig('Min_L(h)_histogram_2.png')