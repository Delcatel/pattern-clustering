from matplotlib.pyplot import *
from numpy import array
from clustering import *

def muOda():	# N=500, 1000 tests, eps = 0.5 -> snr = 5
	thresholdOda = array([0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.68])
	muOda = array([3.43, 3.41, 3.44, 3.56, 3.74, 4.05, 4.49])
	sigmaOda = array([1.05, 1.02, 1.05, 1.09, 1.12, 1.27, 1.6])

	#for k in threshold:
	#	a, b = histogramOda(1000, 50, k)
	#	muOda.append(a)
	#	sigmaOda.append(b)

	plot(thresholdOda, muOda, color="red")
	plot(thresholdOda, muOda+sigmaOda, color="green")
	plot(thresholdOda, muOda-sigmaOda, color="green")
	xlim(0, 0.7)
	ylim(0, 7)
	ylabel("Erreur relative sur H (en %)")
	xlabel("seuil de clusterisation")
	text(0.3, 4, r'$\mu$', color="red")
	text(0.3, 5.2, r'$\mu+\sigma$', color ="green")
	text(0.3, 2.2, r'$\mu-\sigma$', color ="green")
	savefig("images/errorOfThresholdOda.png")
	show()

def muOba():	# N=500, 1000 tests, eps = 0.5 -> snr = 5
	thresholdOba = []#array([0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.68])
	muOba = []#array([1.49, 1.51, 1.5, 1.46, 1.42, 1.44, 1.42])
	sigmaOba = []#array([0.82, 0.83, 0.81, 0.78, 0.73, 0.76])

	plot(thresholdOba, muOba, color="red")
	plot(thresholdOba, muOba+sigmaOba, color="green")
	plot(thresholdOba, muOba-sigmaOba, color="green")
	ylabel("Erreur relative sur H")
	xlabel("seuil de clusterisation (en %)")
	xlim(0, 0.7)
	ylim(0, 7)
	savefig("images/errorOfThresholdOba.png")
	show()

muOda()
