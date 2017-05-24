from cmath import *
from math import floor, sqrt
from numpy import argmax, argmin, mean, std
import matplotlib.pyplot as plt


def graphics(set, color='black', linewidth=3, marker='o'):
	plt.scatter([z.real for z in set], [z.imag for z in set], color=color, linewidth=linewidth, marker=marker)

def rectToComplex(set):
	return [complex(a, b) for a, b in set]

def ComplexToRect(set):
	return [(z.real, z.imag) for z in set]

	
def oba(set, nb4Sect=32, surfSect=1, threshold=0.6827):
	sectDeg = pi/(2*nb4Sect)
	value = 0
	surfSect-=1

	while value < len(set)*threshold:
		surfSect+=1
		nbPerSect = [0 for i in range (4*nb4Sect)]
		for z in set:
			nbPerSect[int(floor((pi-phase(z))/sectDeg))]+=1

		indice = argmax([sum([sum([nbPerSect[(i+nb4Sect*k+l)%(4*nb4Sect)] for l in range(surfSect)]) for k in range(4)]) for i in range(nb4Sect)])
		value = sum([sum([nbPerSect[(indice+nb4Sect*k+l)%(4*nb4Sect)] for l in range(surfSect)]) for k in range(4)])

	phi = pi-(indice+surfSect/2)*pi/2/nb4Sect
	deviationAngle = sectDeg*surfSect/2
	clusters = [[], [], [], []]
	for z in set:
		clusters[int(floor(2*(phase(z)-phi)/pi+0.5))%4].append(z)

	centroides = [0, 0, 0, 0]
	for c in range(4):
			centroides[c] = complex(mean([z.real for z in clusters[c]]), mean([z.imag for z in clusters[c]]))    
    
	return clusters, centroides, phi, deviationAngle

def oda (set, nb4Sect=64, surfSect=1, threshold=0.6827):
	""" Algo not based on barycentres but standard deviation """
	sectDeg = pi/(2*nb4Sect)
	value = 0
	surfSect-=1

	while value < len(set)*threshold:
		surfSect+=1
		nbPerSect = [0 for i in range (4*nb4Sect)]
		for z in set:
			nbPerSect[int(floor((pi-phase(z))/sectDeg))]+=1

		indice = argmax([sum([sum([nbPerSect[(i+nb4Sect*k+l)%(4*nb4Sect)] for l in range(surfSect)]) for k in range(4)]) for i in range(nb4Sect)])
		value = sum([sum([nbPerSect[(indice+nb4Sect*k+l)%(4*nb4Sect)] for l in range(surfSect)]) for k in range(4)])

	theta = pi-(indice+surfSect/2)*pi/2/nb4Sect
	R = mean([abs(z) for z in set])
	deviationAngle = sectDeg*surfSect/2

	square = [rect(R, theta+pi/2*k) for k in range(4)]
    
	return R/sqrt(2), (theta-pi/4)%(pi/2), square, deviationAngle

def correction(centroides):
	R = mean([abs(centroides[c]) for c in range(4)])
	phis = [phase(centroides[c]) for c in range(4)]
	m = argmin(phis)
	angles = [phase(centroides[(m+k)%4])-pi/2*k for k in range(4)]
	if angles[-1] < -2*pi:
		angles[-1]+=2*pi
	theta = mean(angles)
	square = [rect(R, theta+pi/2*k) for k in range(4)]

	return R/sqrt(2), (theta-pi/4)%(pi/2), square

def blackHole(set):
	pass