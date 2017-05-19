from scipy.stats import *
from oba import rectToComplex
# https://docs.scipy.org/doc/scipy/reference/stats.html

def square_data(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform):
	return zip(loi.rvs(loc=x1, scale=x2-x1, size=num), loi.rvs(loc=y1, scale=y2-y1, size=num))

# import matplotlib.pyplot as plt
# l = square_data(0, 1, 0, 1, 10)
# 
# plt.scatter([a for a, b in l], [b for a, b in l])
# plt.show()

def genere_square_data(N):
	set = []
	for a, b in [(-1, 1), (-1, -1), (1, 1), (1, -1)]:
		set.extend(square_data(a-0.5, a+0.5, b-0.5, b+0.5, N))

	phaseAleatoire = 2*pi*random()
	return [exp(complex(0, phaseAleatoire))*z for z in rectToComplex(set)]