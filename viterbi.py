from scipy.cluster import *
import numpy as np
import matplotlib.pyplot as plt
import math
from random import random
from generator import *
from oba import *
from GenereCarre import GenerationCarre2

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
    return np.absolute((H-H_estime)/H)
    

def tests_intensifs(eps,nbtests,nbpoints,snr):
    erreurs = 0
    for k in range(nbtests):
        ecart = viterbi_test(nbpoints,snr)
        if ecart > eps:
            erreurs +=1
    return erreurs/nbtests*100
    
    
    
def histogramOba(nbTests, dispIntervals):
	Data = []
	for i in range(1, nbTests+1):
		if(i%dispIntervals==0):
			print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
		Data.append(obaTest(500, 0.5, 128, "auto", False))

	eff, val, patches = plt.hist(Data, range = (0, 5), bins = 25, edgecolor = 'black', normed=True)
	plt.xlabel("Erreur sur H (en %)")
	plt.ylabel("Proportion d'occurrence sur "+str(nbTests)+" tests")
	plt.title("Overlapping Barycentric Algorithm")

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

	plt.savefig("histogramOba"+str(nbTests)+".png")
	plt.show()
	
	
def histoviterbi(nbpoints,snr,nbTests,dispIntervals):
    Data = []
    for i in range(1, nbTests+1):
        if(i%dispIntervals==0):
            print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
        Data.append(viterbi_test(nbpoints,snr))
    eff, val, patches = plt.hist(Data, range = (0, 1), bins = 25, edgecolor = 'black', normed=True)
    plt.xlabel("Erreur sur H (en %)")
    plt.ylabel("Proportion d'occurrence sur "+str(nbTests)+" tests")
    plt.title("Viterbi Algorithm")
    
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


histoviterbi(1000,4,1000,1)