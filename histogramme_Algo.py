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

#plt.scatter([np.real(A[k]) for k in range(N)],[np.imag(A[k]) for k in range(N)])        
##
hist = []
##
for i in range(200):
    print(i*0.5,'%')
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
    alpha = 10
    sigma_tilde = alpha*sigma #disons que l'on se trompe d'un facteur alpha
    
    #Preconditionnement
    for i in range(1000):
        for k in range(N):
            A[k] += cp(random.gauss(0,sigma_tilde),random.gauss(0,sigma_tilde)) + B[k]
    for k in range(N):
        A[k] = A[k]/1000
        
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
        
    def jac_L(data,r,phi,T):
        jac = np.zeros(2)
        
        #T = rempli_T(data,Centroides)
        for k in range(N):
            for j in range(nb_centres):
                if (T[k][j] ==1):
                    r_k = np.linalg.norm(data[k])
                    phi_k = np.arctan(np.imag(data[k])/np.real(data[k]))
                    if np.real(data[k]) < 0:
                        phi_k += np.pi
                        
                    jac[0] += r_k*(np.cos(phi+phi_k)*(r*r_k*np.cos(phi+phi_k) - np.real(Centroides[j])) + np.sin(phi+phi_k)*(r*r_k*np.sin(phi+phi_k) - np.imag(Centroides[j])))
                    jac[1] += r*r_k*(np.sin(phi+phi_k)*np.real(Centroides[j]) - np.cos(phi+phi_k)*np.imag(Centroides[j]))
                    break
        return jac
        
    def descente_grad(H,data,h_0,err,pas,maxiter = 200):
        step = pas
        T = rempli_T(data,Centroides)
        iterations = 0
    
        r_old = np.linalg.norm(h_0)
        phi_old = np.arctan(np.imag(h_0)/np.real(h_0))
        if np.real(h_0) < 0:
            phi_old += np.pi
        
        jjac = jac_L(data,r_old,phi_old,T)
        r_new = r_old - step*jjac[0]
        phi_new = phi_old - step*jjac[1]
        
        while (np.linalg.norm(jac_L(data,r_new,phi_new,T)) > err and iterations < maxiter):
            r_old = r_new
            phi_old = phi_new
            
            jac = jac_L(data,r_old,phi_old,T)
            r_new = r_old - step*jac[0]
            phi_new = phi_old - step*jac[1]
            
            iterations += 1
        if (iterations == maxiter):
            print('Maxiter reached')
        h_new = r_new*np.exp(1j*phi_new)
        h = 1/h_new
        err = np.linalg.norm(h-H)
        
        hh = h
        for i in range(1,4):
            if err > np.linalg.norm((1j)**i*hh - H):
                h = (1j)**i*hh #erreur à pi/2 pres
                err = np.linalg.norm(h - H)
                
        return [h,iterations] # h_new correspond à 1/H
    
    h = descente_grad(H,A,random.random()*np.exp(1j*2*np.pi*random.random()),epsilon,0.01)[0]
    hist += [100*np.linalg.norm(h-H)/np.linalg.norm(H)]

eff, val, patches = plt.hist(histogramme, range = (0, 10), bins = 30, edgecolor = 'black', normed=False)
mu = np.mean(histogramme)
sigma = np.std(histogramme)
max = np.max(eff)
plt.plot([mu,mu],[0,max],color = 'r', ls=':', linewidth = 3)
plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, max), xytext=(mu+0.2, max), color="red")
plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, max*3/4), xytext=(mu+sigma+0.2, max*3/4), color="green")
plt.plot([mu-sigma,mu-sigma],[0,0.75*max], color = 'g', ls=':', linewidth = 3)
plt.plot([mu+sigma,mu+sigma],[0,0.75*max], color = 'g', ls=':', linewidth = 3)
plt.xlabel('Error on $H$ (%)')
plt.ylabel('Frequency of occurency')
plt.title('Minimization of cost function $\mathcal{L}(h)$ Algorithm')
plt.savefig('Min_L(h)_histogram_2.png')