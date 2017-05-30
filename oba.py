# -*- coding: utf-8 -*-

from cmath import *
from math import floor, sqrt
from numpy import argmax, argmin, mean, std, zeros
import matplotlib.pyplot as plt
from operator import itemgetter


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

def correction2bits(centroides):
	R = mean([abs(centroides[c]) for c in range(4)])
	phis = [phase(centroides[c]) for c in range(4)]
	m = argmin(phis)
	angles = [phase(centroides[(m+k)%4])-pi/2*k for k in range(4)]
	if angles[-1] < -2*pi:
		angles[-1]+=2*pi
	theta = mean(angles)
	square = [rect(R, theta+pi/2*k) for k in range(4)]

	return R/sqrt(2), (theta-pi/4)%(pi/2), square

def correction4bits(centroides):
	points = sorted([polar(c) for c in centroides], key=itemgetter(0))
	Rint = 
	Rext = mean([pt[0] for pt in points[]])
	phis = [phase(centroides[c]) for c in range(4)]
	m = argmin(phis)
	angles = [phase(centroides[(m+k)%4])-pi/2*k for k in range(4)]
	if angles[-1] < -2*pi:
		angles[-1]+=2*pi
	theta = mean(angles)
	square = [rect(R, theta+pi/2*k) for k in range(4)]

	return R/sqrt(2), (theta-pi/4)%(pi/2), square

class point():
	def __init__(self, z, num):
		self.z = z
		self.index = num
		self.weight = 1

def fusion(points, i, j):
	for index0, pt in enumerate(points):
		if pt.index == i:
			i0 = index0
		if pt.index == j:
			j0 = index0
			break

	points[i0].z = (points[i0].weight*points[i0].z+points[j0].weight*points[j0].z)/(points[i0].weight+points[j0].weight)
	points[i0].weight = points[i0].weight+points[j0].weight
	del(points[j0])

	return i0

def bha(set, nbClusters):
	""" Black Hole Algorithm """
	nbPoints = len(set)
	points = [point(z, num) for num, z in enumerate(set)]
	sortedPoints = sorted([(i, j, abs(points[i].z-points[j].z)) for i in range(nbPoints) for j in range(i+1, nbPoints)], key=itemgetter(2)) 
	# i < j est un invariant utilisé dans fusion !

	for k in range(nbPoints-nbClusters):
		i, j = sortedPoints[0][0:2]
		i0 = fusion(points, i, j)
		sortedPoints = [x for x in sortedPoints if x[0]!=i and x[1]!=i and x[0]!=j and x[1]!=j]
		for pt in [pt for pt in points if pt.index != i]:
			sortedPoints.append((min(i, pt.index), max(i, pt.index), abs(points[i0].z-pt.z)))
		sortedPoints = sorted(sortedPoints, key=itemgetter(2))
		# print([pt.index for pt in points])

	return [pt.z for pt in points]