from scipy.cluster import *
# https://docs.scipy.org/doc/scipy/reference/cluster.html
from generator import *

def square_kmeans(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform, k=4):
	return vq.kmeans(square_data(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform), k)