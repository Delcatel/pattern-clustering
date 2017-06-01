from scipy.cluster import *
import numpy as np
import matplotlib.pyplot as plt
import math
from random import random
from numpy import argmax, argmin, mean, std



def H(nbpoints,nbtests,algo,nom_algo):
    #on veut des graphes de H en fct du snr. algo doit renvoyer l'erreur sur H.
    L = [5,10,20,30,50]
    H = []
    for snr in L:
        erreur = 0
        for k in range(nbtests):
            if (k+1)%50==0:
                print("Test " + str(k+1) + " sur " + str(nbtests))
            erreur+=algo(nbpoints,snr)
        H.append(erreur/nbtests)
    plt.plot(L,H)
    #plt.show()
    plt.xlabel("SNR")
    plt.axis([0,50,0,4])
    plt.ylabel("Erreur relative sur H")
    #plt.title(nom_algo + "Algorithm")
    
        





def histogram_erreur_relative(nbTests, dispIntervals,algo_test,nbpoints,snr,nom_algo,absc):
    #l'algo_test doit renvoyer une erreur relative en %.. L'algo ne prend en parametre que le nb de points et le snr
    #nom_algo est un string
    #absc est de type (0, 3) et definit les abscicsses
    Data = []
    for i in range(1, nbTests+1):
        if(i%dispIntervals==0):
            print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
        Data.append(algo_test(nbpoints,snr))
    
    eff, val, patches = plt.hist(Data, range = absc, bins = 25, edgecolor = 'black', normed=True)
    #plt.xlabel("Erreur sur H (%)")
    #plt.ylabel("Densité d'occurrence sur "+str(nbTests)+" tests")
    #plt.title(nom_algo + " Algorithm")
    
    mu = mean(Data)
    sigma = std(Data)
    print("mu = "+str(mu))
    print("std = "+str(sigma))
    l = max(eff)
    plt.plot([mu, mu], [0, l], color="red", ls="--",label = 'SNR = ' + str(snr) + ', N = ' + str(nbpoints))
    plt.legend(loc=1,prop={'size':8})
    plt.plot([mu-sigma, mu-sigma], [0, l*3/4], color="green", ls=":")
    plt.plot([mu+sigma, mu+sigma], [0, l*3/4], color="green", ls=":")
    plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, l), xytext=(mu+0.2, 4*l/5), color="red")
    plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=((mu+sigma), l*3/4), xytext=(mu+4*sigma/5+0.2, l*3/5), color="green")
    
    #plt.savefig("histogramOba"+str(nbTests)+".png")
    #plt.show()
    
def histo_erreur_binaire(nbTests, dispIntervals,algo,transmission,nbpoints,snr,nom_algo,absc):
    #algo doit renvoyer H_estime et prendre H et data en parametres.. 
    #transmission correspond a la fct transmission dans transmission_bits
    #nom_algo est un string
    #absc est de type (0, 3) et definit les abscisses
    L = [0]*nbpoints
    for k in range(nbpoints):
        r = random()
        if r > 0.5:
            L[k]+=1
    Data = []
    for i in range(1, nbTests+1):
        if(i%dispIntervals==0):
            print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
        Data.append(transmission(L,snr,algo))
    
    eff, val, patches = plt.hist(Data, range = absc, bins = 25, edgecolor = 'black', normed=True)
    #plt.xlabel("Erreur sur le train binaire (%)")
    #plt.ylabel("Densité d'occurrence sur "+str(nbTests)+" tests")
    #plt.title(nom_algo + " Algorithm")
    
    mu = mean(Data)
    sigma = std(Data)
    print("mu = "+str(mu))
    print("std = "+str(sigma))
    l = max(eff)
    plt.plot([mu, mu], [0, l], color="red", ls="--",label = 'SNR = ' + str(snr) + ', N = ' + str(nbpoints))
    plt.legend(loc=1,prop={'size':8})
    plt.plot([mu-sigma, mu-sigma], [0, l*3/4], color="green", ls=":")
    plt.plot([mu+sigma, mu+sigma], [0, l*3/4], color="green", ls=":")
    plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, l), xytext=(mu+0.2, 4*l/5), color="red")
    plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, l*3/4), xytext=(mu+4*sigma/5+0.2, l*3/5), color="green")
    
    #plt.savefig("histogramOba"+str(nbTests)+".png")
    #plt.show()


#l'idee etant de reperer des "catastrophes" (aka foirage total) sur de nombreux tirages
# eps est le seuil d'ecart relatif en % pour decider d'appeler cela une catastrophe
# on renvoie le pourcentage de catastrophe, et la liste des ecarts associes
def tests_intensifs(eps,nbtests,nbpoints,snr,algo_test):
    cata = 0
    l = []
    for k in range(nbtests):
        ecart = algo_test(nbpoints,snr)
        if ecart > eps:
            cata +=1
            l.append(ecart)
    return (l,cata/nbtests*100)
    

