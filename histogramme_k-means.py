import sys
sys.setrecursionlimit(10000)
import pylab as plt
from scipy.optimize import fmin
from math import *
import random
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


##
N = 1000
A = np.zeros((N,2))
B = np.zeros((N,2))

sigmax = 0.03
sigmay = 0.03

nb_centres = 4
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

def transform(X,theta,tx=0,ty=0,s=1): #rotation de theta et translation de [tx,ty]
    a=math.cos(theta)
    b=math.sin(theta)
    R=np.mat([[a*s, -b*s], [b*s,a*s]])
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
def GenerationCarre2(N,snr):
	m1 = complex(1,1)
	m2 = complex(-1,1)
	m3 = complex(-1,-1)
	m4 = complex(1,-1)
	mk = [m1,m2,m3,m4]
	R = random.random()
	theta = np.pi*random.random()/2.
	H = R*np.exp(complex(0,theta))
	mi = np.array([complex(0,0) for i in range(N)])
	bruit = np.array([complex(0,0) for i in range(N)])
	for i in range(N): #bruits gaussiens
		mi[i] = mk[random.randint(0,3)]
		bruit[i]= np.random.normal(0,1) + complex(0,1)*np.random.normal(0,1)
	#ajustement du snr
	Ps = np.absolute(H)**2*(mi.dot(mi)/N)
	Pb = bruit.dot(bruit)/N
	alpha = np.sqrt(Ps/Pb/snr)
	return (R,theta,H*mi+alpha*bruit)
hist = []  # SNR = 20
hist2=[] # SNR = 30
hist3 = [] #SNR = 50
for n in range(200):   
    print(n*1./6,'%')
    (R,theta,data)= GenerationCarre2(N,10)
    #print('R=',R,'theta=',theta)
    #plt.scatter(np.real(data),np.imag(data))
    #plt.axis("equal")
    #plt.show()
    for i in range(N):
        B[i][0] = np.real(data)[i]
    for i in range(N):
        B[i][1] = np.imag(data)[i]    
        
    #indice = 0
    #nb = int(np.sqrt(nb_centres))
    #for k in range(nb):
    #    for j in range(nb):
    #        for i in range(int(N/nb_centres)):
    #            indice = i+int(N/nb_centres)*(k*nb+j)
    #            x = random.gauss((float(k)+0.5)/nb,sigmax)
    #            y = random.gauss((float(j)+0.5)/nb,sigmay)
    #           #B[indice]=[x,y]
    #            B[indice]= rotate([0,0],[x,y],math.radians(20))
    
    
    
    #plt.scatter([B[i][0] for i in range(N)],[B[i][1] for i in range(N)],color = 'r')
    #plt.axis("equal")
    #plt.show()
    
    kmeans = KMeans(nb_centres).fit(B)
    
    EPS = 1./nb_centres
    #plt.scatter([kmeans.cluster_centers_[k][0] for k in range(nb_centres)],[kmeans.cluster_centers_[k][1] for k in range(nb_centres)],color = 'r')
    #plt.axis("equal")
    #plt.show()
    #plt.scatter([rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[0] for k in range(nb_centres)],[rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[1] for k in range(nb_centres)],color = 'r')
    
    
        
    xr = [kmeans.cluster_centers_[k][0] for k in range(nb_centres)]
    yr = [kmeans.cluster_centers_[k][1] for k in range(nb_centres)]
    Ar = np.array([xr,yr]) #Ar est constituée des barycentres
    
    # ds : distance minimale entre 2 barycentres
    ds = dist(kmeans.cluster_centers_[0],kmeans.cluster_centers_[1])
    #print('ds=',ds)
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
    #print('ds = ',ds)
    #affichage
    #for k,i in enumerate(range(0,len(allsol),int(len(allsol)/15)+1)):
    #    plt.subplot(4,4,k+1)
    #    plt.scatter(Ar[0,:],Ar[1,:])
    #    plot_transform(A,allsol[i])
    #    plt.axis('equal')
        #plt.show()
    #plt.scatter([P0[0],P1[0],P2[0],P3[0]],[P0[1],P1[1],P2[1],P3[1]])
    RelEquiv = [theta,theta + pi/2, theta +pi, theta +3*pi/2, theta + 2*pi, theta - pi/2, theta -pi, theta -3*pi/2, theta - 2*pi]
    angle = p_opt[0]
    mini = theta
    for i in range(9):
        if abs(RelEquiv[i]-p_opt[0]) < abs(mini-p_opt[0]):
            mini= RelEquiv[i]
           
    #mini = theta a pi/2 pres        
    #angle =  p_opt[0] + mini - theta        
    #plt.axis("equal")
    #print('angle = ',angle,' theta = ', theta)
    #print(ds/2*np.exp(complex(0,angle)))
    #print(R*np.exp(complex(0,mini)))
    erreur = 100*abs(ds/2*np.exp(complex(0,angle))-R*np.exp(complex(0,mini)))/abs(R*np.exp(complex(0,mini)))
    #print('erreur= ', erreur)
    if(erreur<30):
        hist +=[erreur]
for n in range(200):   
    print(n*1./6,'%')
    (R,theta,data)= GenerationCarre2(N,30)
    #print('R=',R,'theta=',theta)
    #plt.scatter(np.real(data),np.imag(data))
    #plt.axis("equal")
    #plt.show()
    for i in range(N):
        B[i][0] = np.real(data)[i]
    for i in range(N):
        B[i][1] = np.imag(data)[i]    
        
    #indice = 0
    #nb = int(np.sqrt(nb_centres))
    #for k in range(nb):
    #    for j in range(nb):
    #        for i in range(int(N/nb_centres)):
    #            indice = i+int(N/nb_centres)*(k*nb+j)
    #            x = random.gauss((float(k)+0.5)/nb,sigmax)
    #            y = random.gauss((float(j)+0.5)/nb,sigmay)
    #           #B[indice]=[x,y]
    #            B[indice]= rotate([0,0],[x,y],math.radians(20))
    
    
    
    #plt.scatter([B[i][0] for i in range(N)],[B[i][1] for i in range(N)],color = 'r')
    #plt.axis("equal")
    #plt.show()
    
    kmeans = KMeans(nb_centres).fit(B)
    
    EPS = 1./nb_centres
    #plt.scatter([kmeans.cluster_centers_[k][0] for k in range(nb_centres)],[kmeans.cluster_centers_[k][1] for k in range(nb_centres)],color = 'r')
    #plt.axis("equal")
    #plt.show()
    #plt.scatter([rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[0] for k in range(nb_centres)],[rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[1] for k in range(nb_centres)],color = 'r')
    
    
        
    xr = [kmeans.cluster_centers_[k][0] for k in range(nb_centres)]
    yr = [kmeans.cluster_centers_[k][1] for k in range(nb_centres)]
    Ar = np.array([xr,yr]) #Ar est constituée des barycentres
    
    # ds : distance minimale entre 2 barycentres
    ds = dist(kmeans.cluster_centers_[0],kmeans.cluster_centers_[1])
    #print('ds=',ds)
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
    #print('ds = ',ds)
    #affichage
    #for k,i in enumerate(range(0,len(allsol),int(len(allsol)/15)+1)):
    #    plt.subplot(4,4,k+1)
    #    plt.scatter(Ar[0,:],Ar[1,:])
    #    plot_transform(A,allsol[i])
    #    plt.axis('equal')
        #plt.show()
    #plt.scatter([P0[0],P1[0],P2[0],P3[0]],[P0[1],P1[1],P2[1],P3[1]])
    RelEquiv = [theta,theta + pi/2, theta +pi, theta +3*pi/2, theta + 2*pi, theta - pi/2, theta -pi, theta -3*pi/2, theta - 2*pi]
    angle = p_opt[0]
    mini = theta
    for i in range(9):
        if abs(RelEquiv[i]-p_opt[0]) < abs(mini-p_opt[0]):
            mini= RelEquiv[i]
           
    #mini = theta a pi/2 pres        
    #angle =  p_opt[0] + mini - theta        
    #plt.axis("equal")
    #print('angle = ',angle,' theta = ', theta)
    #print(ds/2*np.exp(complex(0,angle)))
    #print(R*np.exp(complex(0,mini)))
    erreur = 100*abs(ds/2*np.exp(complex(0,angle))-R*np.exp(complex(0,mini)))/abs(R*np.exp(complex(0,mini)))
    #print('erreur= ', erreur)
    if(erreur<30):
        hist2 +=[erreur]
for n in range(200):   
    print(n*1./6,'%')
    (R,theta,data)= GenerationCarre2(N,30)
    #print('R=',R,'theta=',theta)
    #plt.scatter(np.real(data),np.imag(data))
    #plt.axis("equal")
    #plt.show()
    for i in range(N):
        B[i][0] = np.real(data)[i]
    for i in range(N):
        B[i][1] = np.imag(data)[i]    
        
    #indice = 0
    #nb = int(np.sqrt(nb_centres))
    #for k in range(nb):
    #    for j in range(nb):
    #        for i in range(int(N/nb_centres)):
    #            indice = i+int(N/nb_centres)*(k*nb+j)
    #            x = random.gauss((float(k)+0.5)/nb,sigmax)
    #            y = random.gauss((float(j)+0.5)/nb,sigmay)
    #           #B[indice]=[x,y]
    #            B[indice]= rotate([0,0],[x,y],math.radians(20))
    
    
    
    #plt.scatter([B[i][0] for i in range(N)],[B[i][1] for i in range(N)],color = 'r')
    #plt.axis("equal")
    #plt.show()
    
    kmeans = KMeans(nb_centres).fit(B)
    
    EPS = 1./nb_centres
    #plt.scatter([kmeans.cluster_centers_[k][0] for k in range(nb_centres)],[kmeans.cluster_centers_[k][1] for k in range(nb_centres)],color = 'r')
    #plt.axis("equal")
    #plt.show()
    #plt.scatter([rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[0] for k in range(nb_centres)],[rotate([0,0],kmeans.cluster_centers_[k],math.radians(20))[1] for k in range(nb_centres)],color = 'r')
    
    
        
    xr = [kmeans.cluster_centers_[k][0] for k in range(nb_centres)]
    yr = [kmeans.cluster_centers_[k][1] for k in range(nb_centres)]
    Ar = np.array([xr,yr]) #Ar est constituée des barycentres
    
    # ds : distance minimale entre 2 barycentres
    ds = dist(kmeans.cluster_centers_[0],kmeans.cluster_centers_[1])
    #print('ds=',ds)
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
    #print('ds = ',ds)
    #affichage
    #for k,i in enumerate(range(0,len(allsol),int(len(allsol)/15)+1)):
    #    plt.subplot(4,4,k+1)
    #    plt.scatter(Ar[0,:],Ar[1,:])
    #    plot_transform(A,allsol[i])
    #    plt.axis('equal')
        #plt.show()
    #plt.scatter([P0[0],P1[0],P2[0],P3[0]],[P0[1],P1[1],P2[1],P3[1]])
    RelEquiv = [theta,theta + pi/2, theta +pi, theta +3*pi/2, theta + 2*pi, theta - pi/2, theta -pi, theta -3*pi/2, theta - 2*pi]
    angle = p_opt[0]
    mini = theta
    for i in range(9):
        if abs(RelEquiv[i]-p_opt[0]) < abs(mini-p_opt[0]):
            mini= RelEquiv[i]
           
    #mini = theta a pi/2 pres        
    #angle =  p_opt[0] + mini - theta        
    #plt.axis("equal")
    #print('angle = ',angle,' theta = ', theta)
    #print(ds/2*np.exp(complex(0,angle)))
    #print(R*np.exp(complex(0,mini)))
    erreur = 100*abs(ds/2*np.exp(complex(0,angle))-R*np.exp(complex(0,mini)))/abs(R*np.exp(complex(0,mini)))
    #print('erreur= ', erreur)
    if(erreur<30):
        hist3 +=[erreur]     
        
plt.figure()        
plt.subplot(3,1,1)
eff, val, patches = plt.hist(hist, range = (0, 10), bins = 30, edgecolor = 'black', normed=False)
mu = np.mean(hist)
sigma = np.std(hist)
max = np.max(eff)
plt.plot([mu,mu],[0,max],color = 'r', ls=':', linewidth = 3, label = '$ SNR = 10$')
plt.annotate(r'$\mu_1$ = '+str(round(mu, 2))+'%', xy=(mu, max), xytext=(mu+0.2, max-2), color="red")
plt.annotate(r'$\sigma_1$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, max*3/4), xytext=(mu+sigma+0.2, max*3/4), color="green")
#plt.plot([mu-sigma,mu-sigma],[0,0.75*max], color = 'g', ls=':', linewidth = 3)
plt.plot([mu+sigma,mu+sigma],[0,0.75*max], color = 'g', ls=':', linewidth = 3)
#plt.ylabel('Frequency of occurency')

plt.title('k-means algorithm')

plt.legend(loc = 1)

plt.subplot(3,1,2)
eff2, val2, patches2 = plt.hist(hist2, range = (0, 10), bins = 30, edgecolor = 'black', normed=False)
mu2 = np.mean(hist2)
sigma2 = np.std(hist2)
max2 = np.max(eff2)
plt.plot([mu2,mu2],[0,max2],color = 'r', ls=':', linewidth = 3, label = '$ SNR = 30$')
plt.annotate(r'$\mu_2$ = '+str(round(mu2, 2))+'%', xy=(mu2, max2), xytext=(mu2+0.2, max2-5), color="red")
plt.annotate(r'$\sigma_2$ = '+str(round(sigma2, 2))+'%', xy=(mu2+sigma2, max2*3/4), xytext=(mu2+sigma2+0.2, max2*3/4), color="green")
#plt.plot([mu2-sigma2,mu2-sigma2],[0,0.75*max2], color = 'g', ls=':', linewidth = 3)
plt.plot([mu2+sigma2,mu2+sigma2],[0,0.75*max2], color = 'g', ls=':', linewidth = 3)
plt.legend(loc = 1)
#plt.xlabel('Error on $H$ (%)')
plt.ylabel('Frequency of occurency')
plt.subplot(3,1,3)
eff3, val3, patches3 = plt.hist(hist3, range = (0, 10), bins = 30, edgecolor = 'black', normed=False)
mu3 = np.mean(hist3)
sigma3 = np.std(hist3)
max3 = np.max(eff3)
plt.plot([mu3,mu3],[0,max3],color = 'r', ls=':', linewidth = 3, label = '$ SNR = 50$')
plt.annotate(r'$\mu_3$ = '+str(round(mu3, 2))+'%', xy=(mu3, max3), xytext=(mu3+0.2, max3-5), color="red")
plt.annotate(r'$\sigma_3$ = '+str(round(sigma3, 2))+'%', xy=(mu2+sigma3, max3*3/4), xytext=(mu3+sigma3+0.2, max3*3/4), color="green")
#plt.plot([mu3-sigma3,mu3-sigma3],[0,0.75*max3], color = 'g', ls=':', linewidth = 3)
plt.plot([mu3+sigma3,mu3+sigma3],[0,0.75*max3], color = 'g', ls=':', linewidth = 3)
plt.legend(loc = 1)
plt.xlabel('Error on $H$ (%)')
#plt.ylabel('Frequency of occurency')
plt.savefig('kmeans_histogram_2.png')
#p[0] : theta
      