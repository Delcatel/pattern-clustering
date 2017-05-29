# -*- coding: utf-8 -*-
"""
Created on Wed May 17 23:33:51 2017

"""
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

#un k-means avec initialisation aleatoire
def K_means(nbclusters,data,epsilon):
    # choix des centroides initiaux aleatoire
    init=[randint(0,len(data)-1) for i in range(nbclusters)]
    centroides_old = np.zeros((nbclusters,2))
    centroides = np.zeros((nbclusters,2))
    for i in range(nbclusters):
        centroides[i] = np.copy(data[init[i]]) # on prend des points aleatoires des donnees
    clusters = clusterization(data,centroides)
    barycentres = [[0,0]]*nbclusters       
    compteur=0
    
    while (compteur<100 and ecart_centroides(centroides,centroides_old) >= epsilon):   
        compteur+=1
        for i in range(nbclusters):
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
			
def K_means_diag_carre_2bits(data,epsilon):
    init=[randint(0,len(data)-1) for i in range(4)]
    n0=4
    centroides_old = np.zeros((n0,2))
    centroides = np.zeros((n0,2))
    for i in range(n0):
        centroides[i] = data[init[i]]
    barycentres = [[0,0]]*n0              
    clusters = clusterization(data,centroides)
    compteur=0
    while (compteur<100 and ecart_centroides(centroides,centroides_old) > epsilon):
        compteur+=1
        for i in range(n0):
            if len(clusters[i])==0:
                barycentres[i]=np.copy(centroides[i])
            else:
                barycentres[i] = np.copy(barycentre_cluster(clusters[i]))
        centroides_old = np.copy(centroides)
        (dmax,imax,jmax)=((barycentres[0][0]-barycentres[1][0])**2+(barycentres[0][1]-barycentres[1][1])**2,0,1)
        for i in range(n0):
            for j in range(n0):
                a=(barycentres[i][0]-barycentres[j][0])**2+(barycentres[i][1]-barycentres[j][1])**2
                if (a>dmax):
                    dmax=a
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
				
#initialisation aleatoire
"""def K_means_diag_carre_4bits(data,epsilon):
    n0=16
    init=[randint(0,len(data)-1) for i in range(n0)]
    
    centroides_old = np.zeros((n0,2))
    centroides = np.zeros((n0,2))
    for i in range(n0):
        centroides[i] = data[init [i]]
    barycentres = [[0,0]]*n0              
    clusters = clusterization(data,centroides)
    compteur=0
    while (compteur<10000 and ecart_centroides(centroides,centroides_old) > epsilon):
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
        b=complex(centroides[1][0],centroides[1][1])
        omega=(a+b)/2.
        c=omega+complex(0,1)*(b-omega)
        d=omega+complex(0,1)*(a-omega)
        #on redefinit tous les centroides a partir de a, b, c, d
        p=0 #le premier centroide a definir est le numero 0
        pointcourant=complex(0,0)
        t=0.25
        for c1 in [a,b,c,d]:
            for c2 in [a,b,c,d]:
                pointcourant=t*c1+(1-t)*c2
                centroides[p][0]=pointcourant.real
                centroides[p][1]=pointcourant.imag
                p+=1
        clusters = clusterization(data,centroides)
    #print(compteur)
    return [clusters,centroides]"""
				
#L'initialisation n'est pas aleatoire: on prend les deux points de data les plus eloignes entre eux et on les choisit en premier
def K_means_diag_carre_4bits(data,epsilon):
    n0=16
    init=[randint(0,len(data)-1) for k in range (n0)]
    (dmax,imax,jmax)=((data[0][0]-data[1][0])**2+(data[0][1]-data[1][1])**2,0,1)
    for i in range(len(data)):
        for j in range(len(data)):
            a=(data[i][0]-data[j][0])**2+(data[i][1]-data[j][1])**2
            if (a>dmax):
                dmax=a
                imax=i
                jmax=j
    init[0]=imax
    init[1]=jmax
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
                barycentres[i] = np.copy(barycentre_cluster(clusters[i]))
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
        centroides[1]=np.copy(barycentres[jmax])
        a=complex(centroides[0][0],centroides[0][1])
        b=complex(centroides[1][0],centroides[1][1])
        omega=(a+b)/2.
        c=omega+complex(0,1)*(b-omega)
        d=omega+complex(0,1)*(a-omega)
        #on redefinit tous les centroides a partir de a, b, c, d
        p=0 #le premier centroide a definir est le numero 0
        pointcourant=complex(0,0)
        t=1./3.
        for c1 in [a,b,c,d]:
            for c2 in [a,b,c,d]:
                pointcourant=t*c1+(1-t)*c2
                centroides[p][0]=pointcourant.real
                centroides[p][1]=pointcourant.imag
                p+=1
        clusters = clusterization(data,centroides)
    #print(compteur)
    return [clusters,centroides]

#renvoie R et theta compris entre 0 et pi/2
#data est sous forme de liste des coordonnees	
def Resolution_K_means_diag_carre_2bits(data):
	n0=4
	[clusters,centroides] = K_means_diag_carre_2bits(data,0.00001)
	(dmax,imax,jmax)=((centroides[0][0]-centroides[1][0])**2+(centroides[0][1]-centroides[1][1])**2,0,1)
	for i in range(n0):
		for j in range(n0):
			a=(centroides[i][0]-centroides[j][0])**2+(centroides[i][1]-centroides[j][1])**2
			if (a>dmax):
				dmax=a
				imax=i
				jmax=j
	R=np.sqrt(dmax/8.)
	
	a=complex(centroides[imax][0], centroides[imax][1])
	b=complex(centroides[jmax][0], centroides[jmax][1])
	theta = phase(a-b)
	return(R,theta-0.25*np.pi)
	
def Resolution_K_means_diag_carre_4bits(data):
	n0=16
	[clusters,centroides] = K_means_diag_carre_4bits(data,0.00001)
	(dmax,imax,jmax)=((centroides[0][0]-centroides[1][0])**2+(centroides[0][1]-centroides[1][1])**2,0,1)
	for i in range(n0):
		for j in range(n0):
			a=(centroides[i][0]-centroides[j][0])**2+(centroides[i][1]-centroides[j][1])**2
			if (a>dmax):
				dmax=a
				imax=i
				jmax=j
	R=np.sqrt(dmax/72.)
	
	a=complex(centroides[imax][0], centroides[imax][1])
	b=complex(centroides[jmax][0], centroides[jmax][1])
	theta = phase(a-b)
	return(R,theta-0.25*np.pi)
	
#Ne fonctionne qu'avec n0=4 ou 16
def Resolution_K_means(n0,data):
	[clusters,centroides] = K_means(n0,data,0.00001)
	(dmax,imax,jmax)=((centroides[0][0]-centroides[1][0])**2+(centroides[0][1]-centroides[1][1])**2,0,1)
	for i in range(n0):
		for j in range(n0):
			a=(centroides[i][0]-centroides[j][0])**2+(centroides[i][1]-centroides[j][1])**2
			if (a>dmax):
				dmax=a
				imax=i
				jmax=j
	if n0==4:
		R=np.sqrt(dmax/8.)
	if n0==16:
		R=np.sqrt(dmax/72.)

	a=complex(centroides[imax][0], centroides[imax][1])
	b=complex(centroides[jmax][0], centroides[jmax][1])
	theta = phase(a-b)
	return(R,theta-0.25*np.pi)
###################################################################################################################################
###################################################################################################################################	

#Pour visualiser le resultat de k_means ou de k_means_diag_carre_2bits en fonction de la ligne decommentee
def test1_2bits(SNR):
	N=500
	(R,theta,mi)=GenerationCarre2(N,SNR)
	data=[[x.real,x.imag] for x in mi]
	#[clusters,centroides] = K_means_diag_carre_2bits(data,0.00001)
	[clusters,centroides] = K_means(4,data,0.00001)
	colors = ['b', 'c', 'y', 'm', 'r']
	for k in range(4):
		X=[x[0] for x in clusters[k]]
		Y=[x[1] for x in clusters[k]]
		plt.scatter(X,Y,color=colors[k])
		plt.scatter(centroides[k][0],centroides[k][1],color='k')
	plt.show()

def test1_4bits(SNR):
	N=500
	(R,theta,mi)=Generation4bits2(N,SNR)
	data=[[x.real,x.imag] for x in mi]
	[clusters,centroides] = K_means_diag_carre_4bits(data,0.00001)
	#[clusters,centroides] = K_means(16,data,0.00001)
	colors = [(random(),random(),random()) for k in range(16)]  # en RGB avec 3 nb entre 0 et 1
	for k in range(16):
		X=[x[0] for x in clusters[k]]
		Y=[x[1] for x in clusters[k]]
		plt.scatter(X,Y,color=colors[k])
		plt.scatter(centroides[k][0],centroides[k][1],color='k')
	plt.show()
	
#Permet de tester k_means ou k_means_diag_carre plusieurs fois sur les memes points et d'observer la dependance a l'initialisation,
#qui est aleatoire pour les 2
def test2(data):
	[clusters,centroides] = K_means(4,data,0.00001)
	#[clusters,centroides] = K_means_diag_carre(data,0.0001)
	colors = ['b', 'c', 'y', 'm', 'r']
	for k in range(4):
		X=[x[0] for x in clusters[k]]
		Y=[x[1] for x in clusters[k]]
		plt.scatter(X,Y,color=colors[k])
		plt.scatter(centroides[k][0],centroides[k][1],color='k')
	plt.show()
	
def norme(z):
	return z.real**2+z.imag**2
#compare le H approche au vrai H avec l'erreur relative en %
#n0=4 ou 16
def erreur_K_means_2bits(N,SNR):
	n0=4
	(R,theta,mi)=GenerationCarre2(N,SNR)
	data=[[x.real,x.imag] for x in mi]
	(Rapp,thetaapp)=Resolution_K_means(n0,data)
	H=R*np.exp(np.complex(0,theta))
	Happ=Rapp*np.exp(np.complex(0,thetaapp))
	return (100*min([norme(Happ*np.exp(np.complex(0,k*np.pi*0.5))-H)/norme(H) for k in range(4)]))
	
def erreur_K_means_4bits(N,SNR):
	n0=16
	(R,theta,mi)=GenerationCarre2(N,SNR)
	data=[[x.real,x.imag] for x in mi]
	(Rapp,thetaapp)=Resolution_K_means(n0,data)
	H=R*np.exp(np.complex(0,theta))
	Happ=Rapp*np.exp(np.complex(0,thetaapp))
	return (100*min([norme(Happ*np.exp(np.complex(0,k*np.pi*0.5))-H)/norme(H) for k in range(4)]))
	
def erreur_K_means_diag_carre_2bits(N,SNR):
	(R,theta,mi)=GenerationCarre2(N,SNR)
	data=[[x.real,x.imag] for x in mi]
	(Rapp,thetaapp)=Resolution_K_means_diag_carre_2bits(data)
	H=R*np.exp(np.complex(0,theta))
	Happ=Rapp*np.exp(np.complex(0,thetaapp))
	return (100*min([norme(Happ*np.exp(np.complex(0,k*np.pi*0.5))-H)/norme(H) for k in range(4)]))
	
def erreur_K_means_diag_carre_4bits(N,SNR):
	(R,theta,mi)=Generation4bits2(N,SNR)
	data=[[x.real,x.imag] for x in mi]
	(Rapp,thetaapp)=Resolution_K_means_diag_carre_4bits(data)
	H=R*np.exp(np.complex(0,theta))
	Happ=Rapp*np.exp(np.complex(0,thetaapp))
	return (100*min([norme(Happ*np.exp(np.complex(0,k*np.pi*0.5))-H)/norme(H) for k in range(4)]))

def erreur_moy_et_nb_erreurs(fonction_erreur,N,SNR):
	Nb_echantillons=500
	SNR=10
	s=0
	a=0
	nbechec=0
	for i in range(Nb_echantillons):
		a=fonction_erreur(N,SNR)
		if (a>5):
			#print(a)
			nbechec+=1
		s+=a
	return([s*(1./N),nbechec])
	
#print(erreur_moy_et_nb_erreurs(erreur_K_means_2bits,500,10))
	
fonctions_erreurs_2bits=[erreur_K_means_2bits,erreur_K_means_diag_carre_2bits]
fonctions_erreurs_4bits=[erreur_K_means_4bits,erreur_K_means_diag_carre_4bits]
valeurs_SNR=[10,20,30,40,50,60]
valeurs_N=[50,100,500,1000]