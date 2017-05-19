from cmath import *
from math import floor
from numpy import argmax, argmin, mean
import matplotlib.pyplot as plt


def graphics(set, color='black', linewidth=3, marker='o'):
	plt.scatter([z.real for z in set], [z.imag for z in set], color=color, linewidth=linewidth, marker=marker)

def rectToComplex(set):
	return [complex(a, b) for a, b in set]

def ComplexToRect(set):
	return [(z.real, z.imag) for z in set]


def oba(set, nb4Sect=4, surfSect=2):
	sectDeg = pi/(2*nb4Sect)
	nbPerSect = [0 for i in range (4*nb4Sect)]

	for z in set:
		nbPerSect[int(floor((pi-phase(z))/sectDeg))]+=1

	choice = argmax([sum([sum([nbPerSect[(i+nb4Sect*k+l)%(4*nb4Sect)] for l in range(surfSect)]) for k in range(4)]) for i in range(nb4Sect)])

	phi = pi-(choice+surfSect/2)*pi/2/nb4Sect

	refs = [exp(complex(0, phi+k*pi/2)) for k in range(4)]
	for z in refs:
		z0 = exp(complex(0, sectDeg*surfSect/2))
		A = abs(z)
		plt.plot([0, (z*z0/A).real], [0, (z*z0/A).imag], color="black", ls=":")
		plt.plot([0, (z/z0/A).real], [0, (z/z0/A).imag], color="black", ls=":")

	clusters = [[], [], [], []]
	for z in set:
		clusters[int(floor(2*(phase(z)-phi)/pi+0.5))%4].append(z)

	centroides = [0, 0, 0, 0]
	print("Centroides avant correction :")
	for c in range(4):
		centroides[c] = complex(mean([z.real for z in clusters[c]]), mean([z.imag for z in clusters[c]]))
		print('R=',abs(centroides[c]),'theta=',phase(centroides[c])) 
    
    
	return clusters, centroides
	

def correction(centroides):
	R = mean([abs(centroides[c]) for c in range(4)])
	phis = [phase(centroides[c]) for c in range(4)]
	m = argmin(phis)
	theta = mean([phase(centroides[(m+k)%4])-pi/2*k for k in range(4)])
	square = [R*exp(complex(0, theta+pi/2*k)) for k in range(4)]

	return square

def blackHole(set):
	pass