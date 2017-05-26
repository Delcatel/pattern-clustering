from scipy.cluster import *
import numpy as np
import matplotlib.pyplot as plt
import math
from random import random
from numpy import argmax, argmin, mean, std








def histogram_erreur_relative(nbTests, dispIntervals,algo_test,nbpoints,snr,nomalgo,absc):
    #l'algo_test doit renvoyer une erreur relative en %.. L'algo ne prend en parametre que le nb de points et le snr
    #nomalgo est un string
    #absc est de type (0, 3) et definit les abscicsses
    Data = []
    for i in range(1, nbTests+1):
        if(i%dispIntervals==0):
            print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
        Data.append(algo_test(nbpoints,snr))
    
    eff, val, patches = plt.hist(Data, range = absc, bins = 25, edgecolor = 'black', normed=True)
    plt.xlabel("Erreur sur H (en %)")
    plt.ylabel("Proportion d'occurrence sur "+str(nbTests)+" tests")
    plt.title(nomalgo + " Algorithm")
    
    mu = mean(Data)
    sigma = std(Data)
    print("mu = "+str(mu))
    print("std = "+str(sigma))
    l = max(eff)
    plt.plot([mu, mu], [0, l], color="red", ls="--")
    plt.plot([mu-sigma, mu-sigma], [0, l*3/4], color="green", ls=":")
    plt.plot([mu+sigma, mu+sigma], [0, l*3/4], color="green", ls=":")
    plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, l), xytext=(mu+0.2, l), color="red")
    plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, l*3/4), xytext=(mu+sigma+0.2, l*3/4), color="green")
    
    #plt.savefig("histogramOba"+str(nbTests)+".png")
    plt.show()


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
    

