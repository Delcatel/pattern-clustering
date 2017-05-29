from scipy.cluster import *
import numpy as np
import matplotlib.pyplot as plt
import math
from random import random
from GenereCarre import GenerationCarre2
from analyse import*

def viterbi(data,H):
    H_estime=0
    N= len(data)
    for k in range(N):
        H_estime+= data[k]**4
    H_estime= (-1/4/N)*H_estime
    H_estime = H_estime**(1/4)
    #on va lever le doute sur la rotation
    H_est = [H_estime,complex(0,1)*H_estime,-H_estime,-complex(0,1)*H_estime]
    ecart = np.absolute(H-H_estime)
    indice = 0
    for i in range(1,4):
        if np.absolute(H_est[i]-H) < ecart:
            ecart = np.absolute(H_est[i]-H)
            indice = i
    return H_est[indice]
    
    
def viterbi_test(nbpoints,snr):
    (R,theta,data) = GenerationCarre2(nbpoints,snr)
    H = R*np.exp(complex(0,theta))
    H_estime = viterbi(data,H)
    return np.absolute((H-H_estime)/H)*100


        
def main():
    N=500
    snr = 8
    (R,theta,data) = GenerationCarre2(N,snr)
    H = R*np.exp(complex(0,theta))
    plt.scatter(np.real(data),np.imag(data))
    plt.show()
    H_estime = viterbi(data,H)
    print(H)
    print(H_estime)



histogram_erreur_relative(1000,1,viterbi_test,100,10,"Viterbi",(0, 5))
#(l,taux_cata) = tests_intensifs(5,100,10000,5,viterbi_test)
#print(l)
#print(taux_cata)