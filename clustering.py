from scipy.cluster import *
# https://docs.scipy.org/doc/scipy/reference/cluster.html
from generator import *
from oba import *

def square_kmeans(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform, k=4):
	return vq.kmeans(square_data(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform), k)

def obaTest():
	set = []
	for a, b in [(-1, 1), (-1, -1), (1, 1), (1, -1)]:
		set.extend(square_data(a-0.5, a+0.5, b-0.5, b+0.5, 100))

	# oba(rectToComplex(set))
	oba(rectToComplex(set), 8, 3)

obaTest()