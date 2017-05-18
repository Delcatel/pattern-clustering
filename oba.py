from cmath import *
from math import floor
from numpy import argmax, mean
import matplotlib.pyplot as plt


def graphics(set, color):
	plt.scatter([z.real for z in set], [z.imag for z in set], color=color)

def rectToComplex(set):
	return [complex(a, b) for a, b in set]

def ComplexToRect(set):
	return [(z.real, z.imag) for z in set]


def oba(set, nb4Sect=4, surfSect=2, nbIter=3):
	sectDeg = pi/(2*nb4Sect)
	nbPerSect = [0 for i in range (4*nb4Sect)]

	for z in set:
		nbPerSect[int(floor((pi-phase(z))/sectDeg))]+=1

	choice = argmax([sum([sum([nbPerSect[(i+nb4Sect*k+l)%(4*nb4Sect)] for l in range(surfSect)]) for k in range(4)]) for i in range(nb4Sect)])

	phi = pi-(choice+surfSect/2)*pi/2/nb4Sect

	refs = [exp(complex(0, phi+k*pi/2)) for k in range(4)]
	graphics(refs, "black")

	clusters = [[], [], [], []]
	for z in set:
		clusters[int(floor(2*(phase(z)-phi)/pi+0.5))%4].append(z)

	centroides = [0, 0, 0, 0]
	for c in range(4):
		centroides[c] = complex(mean([z.real for z in clusters[c]]), mean([z.imag for z in clusters[c]]))
		print(centroides[c]) 
    
    
	# setRect = ComplexToRect(set)
	# centroidesRect = ComplexToRect(centroides)
	colors=["blue", "green", "purple", "orange"]

	for i in range(4):
		graphics(clusters[i], colors[i])
	graphics(centroides, "red")

	plt.show()

