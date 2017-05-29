"""
Created on Thu May 18 12:26:26 2017

@author: dekeyser
"""

from random import *
import numpy as np
import matplotlib.pyplot as plt
import math

#R et theta correspondent au module et a l'argument du facteur de deformation de canal
#data est la donnee des centroides sous la forme d'une liste de complexes et C est la donnee des clusters

#Genere des points autour des coins d'un carre tire aleatoirement. eps est entre 0 et 1. 
# Plus il est petit, plus les points sont proches des coins du carre.
#R et theta sont renvoyes, on peut s'en servir pour evaluer l'efficacite de nos algorithmes
#NE PLUS UTILISER CETTE GENERATION
def GenerationCarre(N,eps):
	m1 = complex(1,1)
	m2 = complex(-1,1)
	m3 = complex(-1,-1)
	m4 = complex(1,-1)
	mk = [m1,m2,m3,m4]
	R = random()
	theta = np.pi*random()/2.
	H = R*np.exp(complex(0,theta))
	mi = np.array([complex(0,0) for i in range(N)])
	for i in range(N):
		gaussienne = np.random.normal(0,eps,1)[0]
		phi = 2*np.pi*random()
		mi[i] = H*(mk[randint(0,3)]+gaussienne*np.exp(complex(0,phi)))
#on a eparpille au hasard des points autour des sommets du carre theorique et maintenant on effectue la transformation du canal.

	return (R,theta,mi)

#on a un pilotage de generation en snr et non plus en eps
# UTILISER CELLE CI
def GenerationCarre2(N,snr):
	m1 = complex(1,1)
	m2 = complex(-1,1)
	m3 = complex(-1,-1)
	m4 = complex(1,-1)
	mk = [m1,m2,m3,m4]
	R = random()
	theta = np.pi*random()/2.
	H = R*np.exp(complex(0,theta))
	mi = np.array([complex(0,0) for i in range(N)])
	bruit = np.array([complex(0,0) for i in range(N)])
	for i in range(N): #bruits gaussiens
		mi[i] = mk[randint(0,3)]
		bruit[i]= np.random.normal(0,1) + complex(0,1)*np.random.normal(0,1)
	#ajustement du snr
	Ps = np.absolute(H)**2*(mi.dot(np.conj(mi))/N)
	Pb = bruit.dot(np.conj(bruit))/N
	alpha = np.sqrt(Ps/Pb/snr)
	return (R,theta,H*mi+alpha*bruit)
	
		
	
	#PAS CELLE LA
def Generation4bits(N,eps):
	entiers=[-3,-1,1,3]
	mk=[[complex(i,j) for i in entiers] for j in entiers]
	mk=mk[0]+mk[1]+mk[2]+mk[3]
	R=random()
	theta = np.pi*random()/2.
	H = R*np.exp(complex(0,theta))
	mi = np.array([complex(0,0) for i in range(N)])
	for i in range(N):
		gaussienne = np.random.normal(0,eps,1)[0]
		phi = 2*np.pi*random()
		mi[i] = H*(mk[randint(0,15)]+gaussienne*np.exp(complex(0,phi)))
	return (R,theta,mi)


#PAS CELLE LA	
def Affichage2bits():
	N=500
	eps=0.3
	(R,theta,data)=GenerationCarre(N,eps)
	print('R=',R,'theta=',theta)
	plt.scatter(np.real(data),np.imag(data))
	
#MAIS CELLE CI	
def Affichage2bits2(N,snr):
	(R,theta,data)= GenerationCarre2(N,snr)
	print('R=',R,'theta=',theta)
	plt.scatter(np.real(data),np.imag(data))
	plt.show()
	

	
def Affichage4bits():
	N=500
	eps=0.3
	(R,theta,data)=Generation4bits(N,eps)
	print('R=',R,'theta=',theta)
	plt.scatter(np.real(data),np.imag(data))


#Affichage2bits()
#Affichage2bits2(500,10)