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


        



"""plt.subplot(3,1,1)
histogram_erreur_relative(1000,50,viterbi_test,1000,5,"Viterbi",(0, 5))
plt.title("Viterbi algorithm")
plt.subplot(3,1,2)
histogram_erreur_relative(1000,50,viterbi_test,1000,10,"Viterbi",(0, 5))
plt.ylabel("Densité d'occurrence sur 1000 tests")
plt.subplot(3,1,3)
histogram_erreur_relative(1000,50,viterbi_test,1000,30,"Viterbi",(0, 5))
plt.xlabel("Erreur sur H (%)")

plt.show()"""