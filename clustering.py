# -*- coding: utf-8 -*-

from scipy.cluster import *
from random import random
# https://docs.scipy.org/doc/scipy/reference/cluster.html
from generator import *
from oba import *
from GenereCarre import GenerationCarre

def square_kmeans(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform, k=4):
	return vq.kmeans(square_data(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform), k)

def obaTest(N, eps):
	(R, theta, set)=GenerationCarre(N,eps)
	print("Paramètres de génération :")
	print('R=',R,'theta=',theta)

	clusters, centroides = oba(set, 32, 8)
	colors=["blue", "green", "purple", "orange"]
	for i in range(4):
		graphics(clusters[i], colors[i])
	graphics(centroides, "red")

	square = correction(centroides)
	graphics(square, "yellow", 0.5)

	plt.show()

obaTest(500, 0.3)