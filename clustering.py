# -*- coding: utf-8 -*-

from scipy.cluster import *
from analyse import *
from random import random
# https://docs.scipy.org/doc/scipy/reference/cluster.html
from oba import *
from GenereCarre import GenerationCarre, GenerationCarre2, Generation4bits, Generation4bits2

def verbosity(verbose, error, err_threshold, R, theta, R2, theta2):	# O(1)
	if verbose==True or (verbose=="auto" and error >= 5):
		print("Paramètres de génération :")
		print('R = '+str(R)+'; theta = '+str(theta))
		print("Paramètres estimés à la réception :")
		print('R = '+str(R2)+'; theta = '+str(theta2))
		print('Erreur sur H : '+str(round(error, 1))+' %\n')

def obaTest(nbPoints, snr, precision=32, verbose="auto", display="auto", clusteringThreshold=0.1):	# O(nbPoints)
	(R, theta, set) = GenerationCarre2(nbPoints,snr)	# O(nbPoints)
	clusters, centroides, phi, deviationAngle = oba(set, precision, 1, clusteringThreshold)	# O(nbPoints + precision^2)
	R2, theta2, square = correction2bits(centroides)	# O(1)
	error = min([100*abs(rect(R2, theta2+pi/2*k)-rect(R, theta))/R for k in range(4)])

	verbosity(verbose, error, 5, R, theta, R2, theta2)	# O(1)
	
	if display==True or (display=="auto" and error >= 5):	# O(nbPoints)
		colors=["blue", "green", "purple", "orange"]
		for i in range(4):
			graphics(clusters[i], colors[i])
		graphics(centroides, "red")
		phaseRefs = [phi+k*pi/2 for k in range(4)]
		for p in phaseRefs:
			plt.plot([0, 1.5*R2*cos(p+deviationAngle)], [0, 1.5*R2*sin(p+deviationAngle)], color="black", ls=":")
			plt.plot([0, 1.5*R2*cos(p-deviationAngle)], [0, 1.5*R2*sin(p-deviationAngle)], color="black", ls=":")

		graphics(square, "yellow", 0.5)
		for i in range(4):
			plt.plot([square[i].real, square[(i+1)%4].real], [square[i].imag, square[(i+1)%4].imag], color="yellow", ls="--")

		plt.savefig("images/oba"+str(nbPoints)+"_"+str(snr)+".png")
		plt.show()

	return error

def histogramOba(nbTests, dispIntervals, clusteringThreshold=0.1):
	Data = []
	for i in range(1, nbTests+1):
		if(i%dispIntervals==0):
			print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
		Data.append(obaTest(500, 8, 128, "auto", False, clusteringThreshold))

	eff, val, patches = plt.hist(Data, range = (0, 5), bins = 25, edgecolor = 'black', normed=True)
	plt.xlabel("Erreur sur H (en %)")
	plt.ylabel("Proportion d'occurrence sur "+str(nbTests)+" tests")
	plt.title("Overlapping Barycentric Algorithm")

	mu = mean(Data)
	sigma = std(Data)
	print("mu = "+str(mu))
	print("std = "+str(sigma))
	l = max(eff)
	plt.plot([mu, mu], [0, l], color="red", ls="--")
	plt.plot([mu-sigma, mu-sigma], [0, l*3/4], color="green", ls=":")
	plt.plot([mu+sigma, mu+sigma], [0, l*3/4], color="green", ls=":")
	plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, l), xytext=(mu+0.2, l), color="red")
	plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, l*3/4), xytext=(mu+sigma+0.2, l*3/4), color="green")

	plt.savefig("images/histogramOba"+str(nbTests)+"_"+str(clusteringThreshold)+".png")
	# plt.show()
	return mu, sigma

def odaTest(nbPoints, snr, precision=64, verbose="auto", display="auto", clusteringThreshold=0.1):
	(R, theta, set) = GenerationCarre2(nbPoints, snr)
	R2, theta2, square, deviationAngle = oda(set, precision, 1, clusteringThreshold)
	error = min([100*abs(rect(R2, theta2+pi/2*k)-rect(R, theta))/R for k in range(4)])

	verbosity(verbose, error, 10, R, theta, R2, theta2)

	if display==True or (display=="auto" and error >= 10):
		graphics(set, "blue")
		graphics(square, "red")
		z0 = rect(1, deviationAngle)
		for z in square:
			plt.plot([0, 1.5*(z*z0).real], [0, 1.5*(z*z0).imag], color="black", ls=":")
			plt.plot([0, 1.5*(z/z0).real], [0, 1.5*(z/z0).imag], color="black", ls=":")
		for i in range(4):
			plt.plot([square[i].real, square[(i+1)%4].real], [square[i].imag, square[(i+1)%4].imag], color="red", ls="--")

		plt.savefig("images/oda"+str(nbPoints)+"_"+str(snr)+".png")
		plt.show()

	return error

def histogramOda(nbTests, dispIntervals, clusteringThreshold=0.1):
	Data = []
	for i in range(1, nbTests+1):
		if(i%dispIntervals==0):
			print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
		Data.append(odaTest(500, 0.5, 128, "auto", False, clusteringThreshold))

	eff, val, patches = plt.hist(Data, range = (0, 10), bins = 50, edgecolor = 'black', normed=True)
	plt.xlabel("Erreur sur H (en %)")
	plt.ylabel("Proportion d'occurrence sur "+str(nbTests)+" tests")
	plt.title("Overlapping Deviation Algorithm")

	mu = mean(Data)
	sigma = std(Data)
	print("mu = "+str(mu))
	print("std = "+str(sigma))
	l = max(eff)
	plt.plot([mu, mu], [0, l], color="red", ls="--")
	plt.plot([mu-sigma, mu-sigma], [0, l*3/4], color="green", ls=":")
	plt.plot([mu+sigma, mu+sigma], [0, l*3/4], color="green", ls=":")
	plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, l), xytext=(mu+0.2, l), color="red")
	plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, l*3/4), xytext=(mu+sigma+0.2, l*3/4), color="green")

	plt.savefig("images/histogramOda"+str(nbTests)+"_"+str(clusteringThreshold)+".png")
	# plt.show()

	return mu, sigma

def bhaTest2bits(nbPoints, snr, verbose="auto", display="auto", imageName="bha2bits"):	# O(nbPoints^3 ln(nbPoints))
	(R, theta, set) = GenerationCarre2(nbPoints,snr)	# O(nbPoints)
	graphics(set, "blue")

	centroides = bha(set, 4)	# O(nbPoints^3 ln(nbPoints))
	R2, theta2, square = correction2bits(centroides)
	error = min([100*abs(rect(R2, theta2+pi/2*k)-rect(R, theta))/R for k in range(4)])

	verbosity(verbose, error, 10, R, theta, R2, theta2)

	if display==True or (display=="auto" and error >= 10):
		graphics(set, "blue")
		graphics(centroides, "red")
		graphics(square, "yellow", 0.5)
		for i in range(4):
			plt.plot([square[i].real, square[(i+1)%4].real], [square[i].imag, square[(i+1)%4].imag], color="yellow", ls="--")

		plt.savefig("images/"+imageName+str(nbPoints)+"_"+str(snr)+".png")
		plt.show()

	return error

def bhaTest4bits(nbPoints, snr, verbose="auto", display="auto", imageName="bha4bits"):	# O(nbPoints^3 ln(nbPoints))
	(R, theta, set) = Generation4bits2(nbPoints,snr)	# O(nbPoints)

	centroides = bha(set, 16)	# O(nbPoints^3 ln(nbPoints))
	R2, theta2, square = correction4bits(centroides)	# O(1)
	error = min([100*abs(rect(R2, theta2+pi/2*k)-rect(R, theta))/R for k in range(4)])

	verbosity(verbose, error, 5, R, theta, R2, theta2)	# O(1)

	if display==True or (display=="auto" and error >= 5):	# O(nbPoints)
		graphics(set, "blue")
		graphics(centroides, "red")
		graphics(square, "yellow", 0.5)
		for i in range(4):
			plt.plot([square[4*i].real, square[4*i+3].real], [square[4*i].imag, square[4*i+3].imag], color="yellow", ls="--")
			plt.plot([square[i].real, square[i+12].real], [square[i].imag, square[i+12].imag], color="yellow", ls="--")
		plt.savefig("images/"+imageName+str(nbPoints)+"_"+str(snr)+".png")
		plt.show()

	return error

def histogramBha4bits(nbTests, dispIntervals, nbPoints):
	Data = []
	for i in range(1, nbTests+1):
		if(i%dispIntervals==0):
			print("Test numéro "+str(i)+" sur "+str(nbTests)+"\n")
		Data.append(bhaTest4bits(nbPoints, 0.3, "auto", False, "bha"+str(i)))

	eff, val, patches = plt.hist(Data, range = (0, 10), bins = 50, edgecolor = 'black', normed=True)
	plt.xlabel("Erreur sur H (en %)")
	plt.ylabel("Proportion d'occurrence sur "+str(nbTests)+" tests")
	plt.title("Overlapping Deviation Algorithm")

	mu = mean(Data)
	sigma = std(Data)
	print("mu = "+str(mu))
	print("std = "+str(sigma))
	l = max(eff)
	plt.plot([mu, mu], [0, l], color="red", ls="--")
	plt.plot([mu-sigma, mu-sigma], [0, l*3/4], color="green", ls=":")
	plt.plot([mu+sigma, mu+sigma], [0, l*3/4], color="green", ls=":")
	plt.annotate(r'$\mu$ = '+str(round(mu, 2))+'%', xy=(mu, l), xytext=(mu+0.2, l), color="red")
	plt.annotate(r'$\sigma$ = '+str(round(sigma, 2))+'%', xy=(mu+sigma, l*3/4), xytext=(mu+sigma+0.2, l*3/4), color="green")

	plt.savefig("images/histogramBha"+str(nbTests)+"_"+str(nbPoints)+".png")
	plt.show()

	return mu, sigma

def Test(nbTests, dispIntervals, nbPoints, algo, lendata, geometrie, nomHisto, nomAlgo):
	pass

# obaTest(500, 10, 128, True, True)
# odaTest(500, 10, 128, True, True)
# bhaTest2bits(500, 10, 128, True)
# bhaTest4bits(500, 20, 128, True)

def tripleHistoH(nbTests, nbPoints listSnr, nomAlgo, intervalHist=(0, 5)):
	n = len(listSnr)

	for k in listSnr
		plt.subplot(n,1,k+1)
		histogram_erreur_relative(1000,50,viterbi_test,1000,listSnr[k], nomAlgo,intervalAlgo)

	plt.title("Viterbi algorithm")
	plt.ylabel("Densité d'occurrence sur {} tests de {} points".format(nbTests, nbPoints))
	plt.xlabel("Erreur sur H (%)")

	savefig("images/tripleHisto_{}_{}_{}.png".format(nomAlgo, nbTests, nbPoints))
	plt.show()