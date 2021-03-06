from GenereCarre import*
from random import *
import numpy as np
import matplotlib.pyplot as plt
import math
from viterbi import viterbi
from analyse import*



#la fct prend une liste de bits émise par l'émetteur (on suppose le nbre de bits pair)
#et renvoie l'erreur binaire sur la liste reçue en %
def transmission(L_bits,snr,algo):
    N = int(len(L_bits)/2)
    #on convertit en complexes
    m1 = complex(1,1)
    m2 = complex(-1,1)
    m3 = complex(-1,-1)
    m4 = complex(1,-1)
    R = random()
    theta = np.pi*random()/2.
    H = R*np.exp(complex(0,theta))
    data = np.array([complex(0,0) for i in range(N)])
    bruit = np.array([complex(0,0) for i in range(N)])
    for i in range(N): 
        if L_bits[2*i]==0:
            if L_bits[2*i+1] == 0:
                data[i] = m2
            else:
                data[i] = m1
        else:
            if L_bits[2*i+1] == 0:
                data[i] = m3
            else:
                data[i] = m4
        bruit[i]= np.random.normal(0,1) + complex(0,1)*np.random.normal(0,1)
    #ajustement du snr
    Ps = np.absolute(H)**2*(data.dot(np.conj(data))/N)
    Pb = bruit.dot(np.conj(bruit))/N
    alpha = np.sqrt(Ps/Pb/snr)
    data = H*data + alpha*bruit
    H_estime = algo(data)
    #print("erreur relative sur H en % = " + str(np.absolute((H-H_estime)/H)*100))
    data = data/H_estime
    L_bis = [0]*2*N
    # on reconvertit en binaire
    for i in range(N):
        
        if data[i].real < 0:
            if data[i].imag < 0:
                L_bis[2*i] = 1
                L_bis[2*i+1] = 0
            else:
                L_bis[2*i] = 0
                L_bis[2*i+1] = 0
            
        else:
            if data[i].imag < 0:
                L_bis[2*i] = 1
                L_bis[2*i+1] = 1
            else:
                L_bis[2*i] = 0
                L_bis[2*i+1] = 1
    N = len(L_bits)
    erreurs = 0
    for k in range(N):
        if L_bits[k]!=L_bis[k]:
            erreurs+=1
    return erreurs/N*100

def transmission4bits(L_bits,snr,algo):
    N = int(len(L_bits)//4)
    #on convertit en complexes
    s = [-3+3j, -1+3j, 1+3j, 3+3j, -3+1j, -1+1j, 1+1j, 3+1j, -3-1j, -1-1j, 1-1j, 3-1j, -3-3j, -1-3j, 1-3j, 3-3j]
    bits = ["0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"]
    R = random()
    theta = np.pi*random()/2.
    H = R*np.exp(complex(0,theta))
    data = np.array([complex(0,0) for i in range(N)])
    bruit = np.array([complex(0,0) for i in range(N)])
    for i in range(N):
        data[i] = s[bits.index(L_bits[str(4*i)]+L_bits[str(4*i)+1]+L_bits[str(4*i)+2]+L_bits[str(4*i+3)])]
        bruit[i]= np.random.normal(0,1) + complex(0,1)*np.random.normal(0,1)
    #ajustement du snr
    Ps = np.absolute(H)**2*(data.dot(np.conj(data))/N)
    Pb = bruit.dot(np.conj(bruit))/N
    alpha = np.sqrt(Ps/Pb/snr)
    data = H*data + alpha*bruit
    H_estime = algo(data[0:4*N])
    #print("erreur relative sur H en % = " + str(np.absolute((H-H_estime)/H)*100))
    data = data/H_estime
    L_bis = [0]*(4*N)
    # on reconvertit en binaire
    for i in range(N):
        L_bis[4*i:4*i+3] =bits[np.argmin(np.array([abs(s[k]-data[i]) for k in range(16)]))]

    N = len(L_bits)
    erreurs = 0
    for k in range(N):
        if L_bits[k]!=L_bis[k]:
            erreurs+=1
    return erreurs/N*100
    
    


    



"""plt.subplot(3,1,1)
histo_erreur_binaire(1000,50,viterbi,transmission,100,10,"Viterbi",(0, 5))
plt.title("Viterbi algorithm")
plt.axis([0,5,0,5])
plt.subplot(3,1,2)
histo_erreur_binaire(1000,50,viterbi,transmission,1000,10,"Viterbi",(0, 5))
plt.axis([0,5,0,5])
plt.ylabel("Densité d'occurrence sur 1000 tests")
plt.subplot(3,1,3)
histo_erreur_binaire(1000,50,viterbi,transmission,10000,10,"Viterbi",(0, 5))
plt.xlabel("Erreur sur le train binaire (%)")
plt.axis([0,5,0,5])

plt.show()"""


            