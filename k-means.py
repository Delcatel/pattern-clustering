import random
import numpy as np
import matplotlib.pyplot as plt

##
N = 100
A = np.zeros((N,2))
eps = 1.e-2
nb_centres = 4

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

def K_means(data,n0,epsilon):
    
    # n0 nombre de centroïdes 
    # data l ensemble des points
    # n nombre de points
    # epsilon marge à partir de laquelle on considère le système fixe
    
    # choix des centroides initiaux
    centroides_old = np.zeros((n0,2))
    centroides = np.zeros((n0,2))
    for i in range(n0):
        centroides[i] = data[i] # on prend les n0 premiers points des données
    clusters = [[[0,0]]]*n0
    barycentres = [[0,0]]*n0              
    # tant que les centroides sont trop mobiles, on refait ce qui suit 
    #clusters = clusterization(data,centroides)
    while ( ecart_centroides(centroides,centroides_old) >= epsilon):   
    #print(ecart_centroides(centroides,centroides_old))

        # clusterization
        clusters = clusterization(data,centroides)
        #print(len(clusters[0]))
        # calcul des barycentres pour chaque cluster
        for i in range(n0):
            barycentres[i] = np.copy(barycentre_cluster(clusters[i]))
        centroides_old = np.copy(centroides)
        centroides = np.copy(barycentres)
    return [clusters,centroides]

K0 = K_means(A,nb_centres,eps)[0][0]
l0 = len(K0)

X0 = np.zeros(l0)
Y0 = np.zeros(l0)
for i in range(l0):
    X0[i],Y0[i] = K0[i][0],K0[i][1]

    
K1 = K_means(A,nb_centres,eps)[0][1]
l1 = len(K1)

X1 = np.zeros(l1)
Y1 = np.zeros(l1)
for i in range(l1):
    X1[i],Y1[i] = K1[i][0],K1[i][1]

K2 = K_means(A,nb_centres,eps)[0][2]
l2 = len(K2)

X2 = np.zeros(l2)
Y2 = np.zeros(l2)
for i in range(l2):
    X2[i],Y2[i] = K2[i][0],K2[i][1]
    
K3 = K_means(A,nb_centres,eps)[0][3]
l3 = len(K3)

X3 = np.zeros(l3)
Y3 = np.zeros(l3)
for i in range(l3):
    X3[i],Y3[i] = K3[i][0],K3[i][1]    
  
centroides = K_means(A,nb_centres,eps)[1]
#print(K0)
#print(K1)    
plt.scatter([centroides[0][0],centroides[1][0],centroides[2][0],centroides[3][0]],[centroides[0][1],centroides[1][1],centroides[2][1],centroides[3][1]],marker =  '+')
plt.scatter(X0,Y0, color = 'b')
plt.scatter(X1,Y1, color = 'g')
plt.scatter(X2,Y2, color = 'r')
plt.scatter(X3,Y3, color = 'purple')
plt.show()
