# -*- coding: utf-8 -*-
"""
Created on Wed May 17 23:33:51 2017

"""
from random import *
import numpy as np
import matplotlib.pyplot as plt
import math
from cmath import phase
from GenereCarre import *

#L'idee : a chaque iteration de K-means on retient les barycentres les plus eloignes, ce qui definit
#une diagonale du carre, les sommets duquel definissent les nouveaux centroides
#Implementation sur un carre
###############################################RECOPIAGE#######################################################################
###############################################################################################################################

def associe_centroide(Point,centroides):
    n0 = len(centroides)
    indice_centroide_PP = 0
    dist = (Point[0]-centroides[0][0])**2 + (Point[1]-centroides[0][1])**2
    for i in range(n0):
        if dist > (Point[0]-centroides[i][0])**2 + (Point[1]-centroides[i][1])**2:
            dist = (Point[0]-centroides[i][0])**2 + (Point[1]-centroides[i][1])**2
            indice_centroide_PP = i
    return indice_centroide_PP
    
def clusterization(data,centroides):
    n = len(data)
    classes = [[]]*len(centroides)
    for i in range(n):
        Point = data[i]
        k = associe_centroide(Point,centroides)
    
        #for j in range(len(centroides)):
         #   s=0
          #  for l in range(len(classes[j])):
           #     if classes[j][l][0] == Point[0] and classes[j][l][1] == Point[1]:
            #        s=s+1
            #if s>0:
             #   break
            #else:
                
        classes[k] = classes[k]+[Point]   
    return classes
    
def barycentre_cluster(cluster):
    l = len(cluster)
    barycentre = np.zeros(2)
    for k in range(l):
        barycentre[0] = barycentre[0] + cluster[k][0]
        barycentre[1] = barycentre[1] + cluster[k][1]
    return (1./l)*barycentre
    


def ecart_centroides(centroides_new,centroides_old):
    ecart = 0
    for k in range(len(centroides_new)):
        ecart += (centroides_new[k][0]-centroides_old[k][0])**2 + (centroides_new[k][1]-centroides_old[k][1])**2
    return ecart

#un k-means modifie pour initialiser aleatoirement
def K_means(data,epsilon):
    n0=4
    # choix des centroides initiaux aleatoire
    init=[randint(0,len(data)-1) for i in range(4)]
    centroides_old = np.zeros((n0,2))
    centroides = np.zeros((n0,2))
    for i in range(n0):
        centroides[i] = data[init [i]] # on prend des points aleatoires des donnees
    clusters = [[[0,0]]]*n0
    barycentres = [[0,0]]*n0              
    compteur=0
    while (compteur<100 and ecart_centroides(centroides,centroides_old) >= epsilon):   
        compteur+=1
        for i in range(n0):
            if len(clusters[i])==0:
                barycentres[i]=np.copy(centroides[i])
            else:
                barycentres[i] = barycentre_cluster(clusters[i])
        centroides_old = np.copy(centroides)
        centroides = np.copy(barycentres)
        clusters = clusterization(data,centroides)
    #print(compteur)
    return [clusters,centroides]

#####################################MES FONCTIONS############################################################################
##############################################################################################################################
			
def K_means_diag_carre(data,epsilon):
    init=[randint(0,len(data)-1) for i in range(4)]
    n0=4
    centroides_old = np.zeros((n0,2))
    centroides = np.zeros((n0,2))
    for i in range(n0):
        centroides[i] = data[init [i]]
    barycentres = [[0,0]]*n0              
    clusters = clusterization(data,centroides)
    compteur=0
    while (compteur<100 and ecart_centroides(centroides,centroides_old) > epsilon):
        compteur+=1
        for i in range(n0):
            if len(clusters[i])==0:
                barycentres[i]=np.copy(centroides[i])
            else:
                barycentres[i] = barycentre_cluster(clusters[i])
        centroides_old = np.copy(centroides)
        (max,imax,jmax)=((barycentres[0][0]-barycentres[1][0])**2+(barycentres[0][1]-barycentres[1][1])**2,0,1)
        for i in range(n0):
            for j in range(n0):
                a=(barycentres[i][0]-barycentres[j][0])**2+(barycentres[i][1]-barycentres[j][1])**2
                if (a>max):
                    max=a
                    imax=i
                    jmax=j
        centroides[0]=np.copy(barycentres[imax])
        centroides[2]=np.copy(barycentres[jmax])
        a=complex(centroides[0][0],centroides[0][1])
        b=complex(centroides[2][0],centroides[2][1])
        omega=(a+b)/2.
        c=omega+complex(0,1)*(b-omega)
        d=omega+complex(0,1)*(a-omega)
        centroides[1][0]=d.real
        centroides[1][1]=d.imag
        centroides[3][0]=c.real
        centroides[3][1]=c.imag
        clusters = clusterization(data,centroides)
    #print(compteur)
    return [clusters,centroides]

#renvoie la somme des distances de chaque point a  son centroide
def ErreurK_means(N,eps):
	(R,theta,mi)=GenerationCarre(N,eps)
	data=[[x.real,x.imag] for x in mi]
	[clusters,centroides] = K_means(data,0.00001)
	s=0
	for k in range(len(clusters)):
		for x in clusters[k]:
			s+= (x[0]-centroides[k][0])**2 + (x[1]-centroides[k][1])**2
	return s
	
def ErreurK_means_diag_carre(N,eps):
	(R,theta,mi)=GenerationCarre(N,eps)
	data=[[x.real,x.imag] for x in mi]
	[clusters,centroides] = K_means_diag_carre(data,0.00001)
	s=0
	for k in range(len(clusters)):
		for x in clusters[k]:
			s+= (x[0]-centroides[k][0])**2 + (x[1]-centroides[k][1])**2
	return s

def ComparaisonK_meansK_means_diag_carre(N,eps):
	n=100
	s1=0
	s2=0
	for i in range(n):
		s1+=ErreurK_means(N,eps)
		s2+=ErreurK_means_diag_carre(N,eps)
	return('erreur k_means',s1*(1./n),'erreur k_means_diag_carre',s2*(1./n))

#renvoie R et theta compris entre 0 et pi/2
#data est sous forme de liste des coordonnees	
def Resolution(data):
	[clusters,centroides] = K_means(data,0.00001)
	(max,imax,jmax)=((centroides[0][0]-centroides[1][0])**2+(centroides[0][1]-centroides[1][1])**2,0,1)
	for i in range(4):
		for j in range(4):
			a=(centroides[i][0]-centroides[j][0])**2+(centroides[i][1]-centroides[j][1])**2
			if (a>max):
				max=a
				imax=i
				jmax=j
	print(max)
	R=np.sqrt(max/8.)
	print(R)
	p=0
	while (p==imax or p==jmax):
		p=p+1
	a=complex(centroides[imax][0], centroides[imax][1])
	b=complex(centroides[p][0], centroides[p][1])
	theta = phase(a-b)
	theta=theta%(np.pi/2.)
	return(R,theta)
###################################################################################################################################
###################################################################################################################################	

#Pour visualiser le resultat de k_means ou de k_means_diag_carre en fonction de la ligne decommentee
def test1():
	N=500
	eps=0.3
	(R,theta,mi)=GenerationCarre(N,eps)
	data=[[x.real,x.imag] for x in mi]
	[clusters,centroides] = K_means(data,0.00001)
	#[clusters,centroides] = K_means_diag_carre(data,0.0001)
	colors = ['b', 'c', 'y', 'm', 'r']
	for k in range(4):
		X=[x[0] for x in clusters[k]]
		Y=[x[1] for x in clusters[k]]
		plt.scatter(X,Y,color=colors[k])
		plt.scatter(centroides[k][0],centroides[k][1],color='k')
	plt.show()
	
#Permet de tester k_means ou k_means_diag_carre plusieurs fois sur les memes points et d'observer la dependance a l'initialisation,
#qui est aleatoire pour les 2
def test2(data):
	#[clusters,centroides] = K_means(data,0.00001)
	[clusters,centroides] = K_means_diag_carre(data,0.0001)
	colors = ['b', 'c', 'y', 'm', 'r']
	for k in range(4):
		X=[x[0] for x in clusters[k]]
		Y=[x[1] for x in clusters[k]]
		plt.scatter(X,Y,color=colors[k])
		plt.scatter(centroides[k][0],centroides[k][1],color='k')
	plt.show()
	
#compare le H approche au vrai H
def test3():
	N=500
	eps=0.3
	(R,theta,mi)=GenerationCarre(N,eps)
	data=[[x.real,x.imag] for x in mi]
	(Rapp,thetaapp)=Resolution(data)
	print(Rapp,thetaapp)
	H=R*np.exp(np.complex(0,theta))
	Happ=Rapp*np.exp(np.complex(0,thetaapp))
	plt.scatter(H.real,H.imag, color='r')
	plt.scatter(Happ.real,Happ.imag)
	plt.xlim(-1.,1.)
	plt.ylim(-1.,1.)
	plt.show()

N=100
eps=0.3
print(ComparaisonK_meansK_means_diag_carre(N,eps))
#(R,theta,mi)=GenerationCarre(N,eps)
#data=[[x.real,x.imag] for x in mi]
#test2(data)
#test1()
#test3()

