import sys
sys.setrecursionlimit(10000)
import pylab as plt
from scipy.optimize import fmin
import math
import random
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


##
N = 2000
A = np.zeros((N,2))
B = np.zeros((N,2))

sigmax = 0.03
sigmay = 0.03

nb_centres = 16
eps = 1.e-3
  
#def creeligne(P0,P1):
#    t = np.linspace(0,1,num=NL)
#    x=P0[0]+t*(P1[0]-P0[0])
#    y=P0[1]+t*(P1[1]-P0[1])
#    ligne=np.array([x,y])
#    #plt.scatter(ligne[0,:],ligne[1,:])
#    return ligne
def dist(P0,P1):
    return np.sqrt((P1[0]-P0[0])**2 +(P1[1]-P0[1])**2)

def transform(X,theta,tx=0,ty=0): #rotation de theta et translation de [tx,ty]
    a=math.cos(theta)
    b=math.sin(theta)
    R=np.mat([[a, -b], [b,a]])
    Y=R*X
    Y[0,:] = Y[0,:]+tx
    Y[1,:] = Y[1,:]+ty
    return Y

def error_fun(p,X,Y): #la fonction qu'on minimise
    return np.linalg.norm(transform(X,p[0],p[1],p[2])-Y)

def match(X,Y,p0): #optimisation
    p_opt = fmin(error_fun, p0, args=(X,Y),retall=True)
    return p_opt[0],p_opt[1]


def plot_transform(A,p): #affichage après application de transform
    Afound = transform(A,p[0],p[1],p[2])
    plt.scatter(Afound[0,:],Afound[1,:])

def rotate(origin, point, angle):

    ox, oy = origin[0],origin[1]
    px, py = point[0],point[1]

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return [qx,qy]

# generation des données: répartition uniforme sur [0,1]x[0,1]
for k in range(N):
    x = 2*random.random()-1
    y = 2*random.random()-1
    A[k] = [x,y]

#generation des données: à peu pres réparties sur un carré avec bruit
indice = 0
nb = int(np.sqrt(nb_centres))
for k in range(nb):
    for j in range(nb):
        for i in range(int(N/nb_centres)):
            indice = i+int(N/nb_centres)*(k*nb+j)
            x = random.gauss((float(k)+0.5)/nb,sigmax)
            y = random.gauss((float(j)+0.5)/nb,sigmay)
            #B[indice]=[x,y]
            B[indice]= rotate([0,0],[x,y],math.radians(20))
           

plt.scatter([B[i][0] for i in range(N)],[B[i][1] for i in range(N)],color = 'r')
plt.axis("equal")
plt.show()

kmeans = KMeans(nb_centres).fit(B)

EPS = 1./nb_centres
plt.scatter([kmeans.cluster_centers_[k][0] for k in range(nb_centres)],[kmeans.cluster_centers_[k][1] for k in range(nb_centres)],color = 'r')
plt.axis("equal")
plt.show()
#plt.scatter([rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[0] for k in range(nb_centres)],[rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[1] for k in range(nb_centres)],color = 'r')


    
xr = [kmeans.cluster_centers_[k][0] for k in range(nb_centres)]
yr = [kmeans.cluster_centers_[k][1] for k in range(nb_centres)]
Ar = np.array([xr,yr]) #Ar est constituée des barycentres

# ds : distance minimale entre 2 barycentres
ds = dist(kmeans.cluster_centers_[0],kmeans.cluster_centers_[1])
for i in range(nb_centres):
    for j in range(i):
        if ds > dist(kmeans.cluster_centers_[i],kmeans.cluster_centers_[j]):
            ds = dist(kmeans.cluster_centers_[i],kmeans.cluster_centers_[j])

# grille ' a la main ' (16 points) centrée en [0,0]
if nb_centres == 16:
    x=[-3*ds/2,-3*ds/2,-3*ds/2,-3*ds/2,-ds/2,-ds/2,-ds/2,-ds/2,ds/2,ds/2,ds/2,ds/2,3*ds/2,3*ds/2,3*ds/2,3*ds/2]
    y=[-3*ds/2,-ds/2,ds/2,3*ds/2,-3*ds/2,-ds/2,ds/2,3*ds/2,-3*ds/2,-ds/2,ds/2,3*ds/2,-3*ds/2,-ds/2,ds/2,3*ds/2]      
# carré centré en [0,0]
if nb_centres == 4:
    x = [-ds/2,-ds/2,ds/2,ds/2]
    y = [-ds/2,ds/2,-ds/2,ds/2]
A = np.array([x,y]) #A est la grille qu'on veut superposer à Ar
p0 = np.random.rand(3) # parametres initiaux
p_opt, allsol = match(A,Ar,p0) # optimisation
#plt.scatter(Ar[0,:],Ar[1,:],linewidth=2)
#plt.scatter(A[0,:],A[1,:],linewidth=2)

#affichage
for k,i in enumerate(range(0,len(allsol),int(len(allsol)/15)+1)):
    plt.subplot(4,4,k+1)
    plt.scatter(Ar[0,:],Ar[1,:])
    plot_transform(A,allsol[i])
    plt.axis('equal')
    #plt.show()
#plt.scatter([P0[0],P1[0],P2[0],P3[0]],[P0[1],P1[1],P2[1],P3[1]])
plt.axis("equal")