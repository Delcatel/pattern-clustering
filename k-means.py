import random
import numpy as np
import matplotlib.pyplot as plt

##
N = 100
A = np.zeros((N,2))
eps = 1.e-2
nb_centres = 3

# generation des données: répartition uniforme sur [0,1]x[0,1]
for k in range(N):
    x = random.random()
    y = random.random()
    A[k] = [x,y]

##
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
    classes = [[[0,0]]]*len(centroides)
    for i in range(n):
        Point = data[i]
        k = associe_centroide(Point,centroides)
        #print(classes[k])
        classes[k]+= [Point]
        #print(classes[k])
    return classes
    
def barycentre_cluster(cluster):
    l = len(cluster)
    barycentre = (1./l)*np.array([sum(cluster[0:l-1][0]),sum(cluster[0:l-1][1])])
    return barycentre
    
def ecart_centroides(centroides_new,centroides_old):
    ecart = 0
    for k in range(len(centroides_new)):
        ecart += (centroides_new[k][0]-centroides_old[k][0])**2 + (centroides_new[k][1]-centroides_old[k][1])**2
    return ecart

def K_means(data,n0,epsilon):
    
    # n0 nombre de centroïdes 
    # data l ensemble des points
    # n nombre de points
    # epsilon marge à partir de laquelle on considère le système fixe
    
    # choix des centroides initiaux
    centroides_old = np.zeros((n0,2))
    centroides = np.zeros((n0,2))
    for i in range(n0):
        centroides[i] = data[i] # on prend les 4 premiers points des données
    clusters = [[[0,0]]]*n0
    barycentres = [[0,0]]*n0
    # tant que les centroides sont trop mobiles, on refait ce qui suit 
    while ( ecart_centroides(centroides,centroides_old) >= epsilon):   
        print(ecart_centroides(centroides,centroides_old) >= epsilon)
    
        # clusterization
        clusters = clusterization(data,centroides)
        print(len(clusters[0]))
        # calcul des barycentres pour chaque cluster
        for i in range(n0):
            barycentres[i] = barycentre_cluster(clusters[i])
        centroides_old = centroides
        centroides = barycentres
    return clusters

K0 = K_means(A,nb_centres,eps)[0]
l0 = len(K0)

X0 = np.zeros(l0)
Y0 = np.zeros(l0)
for i in range(l0):
    X0[i],Y0[i] = K0[i][0],K0[i][1]
    
K1 = K_means(A,nb_centres,eps)[1]
l1 = len(K1)

X1 = np.zeros(l1)
Y1 = np.zeros(l1)
for i in range(l1):
    X1[i],Y1[i] = K1[i][0],K1[i][1]

K2 = K_means(A,nb_centres,eps)[2]
l2 = len(K2)

X2 = np.zeros(l2)
Y2 = np.zeros(l2)
for i in range(l2):
    X2[i],Y2[i] = K2[i][0],K2[i][1]
  
print(K0)
print(K1)    
#plt.scatter(X0,Y0, color = 'b')
#plt.scatter(X1,Y1, color = 'g')
#plt.scatter(X2,Y2, color = 'r')
plt.show()