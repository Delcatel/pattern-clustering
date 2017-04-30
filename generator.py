from scipy.stats import *
# https://docs.scipy.org/doc/scipy/reference/stats.html

def square_data(x1=0, x2=1, y1=0, y2=1, num=100, loi=uniform):
	return zip(loi.rvs(loc=x1, scale=x2-x1, size=num), loi.rvs(loc=y1, scale=y2-y1, size=num))