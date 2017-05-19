# -*- coding: utf-8 -*-
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
def GenerationCarre(N,eps):
	m1=complex(1,1)
	m2=complex(-1,1)
	m3=complex(-1,-1)
	m4=complex(1,-1)
	mk=[m1,m2,m3,m4]
	R=random()
	theta=np.pi*random()/2.
	mi=np.array([complex(0,0) for i in range(N)])
	for i in range(N):
		k=randint(0,3)
		e=np.random.normal(0,eps,1)[0]
		phi=2*np.pi*random()
		mi[i]=e*np.exp(complex(0,phi))+mk[k]
#on a eparpille au hasard des points autour des sommets du carre theorique et maintenant on effectue la transformation du canal.
	mi=R*np.exp(complex(0,theta))*mi
	return (R,theta,mi)
	
def main():
	N=500
	eps=0.3
	(R,theta,data)=GenerationCarre(N,eps)
	print('R=',R,'theta=',theta)
	plt.scatter(np.real(data),np.imag(data))
	plt.show()
	
# main()