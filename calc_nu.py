from __future__ import division
import numpy as np
import os, sys, re, urllib2
from beam import beam
from time import sleep,asctime
from loadtxt import loadtxt
from constants import pdgData
"""
This code uses the particle decays in a beamline, duplicates the particles by a certain number of times, force them to decay instanteously in G4Beamline, and filter out only the neutrinos from the beam loss files.
Modification date: 06/09/15;
				   06/10/15:  Bug fixes. G4BLjob content needed a bit modification.
"""

__all__ = ['calcFlux2Body']

def calcFlux2Body(beamMat,detSize):
	"""Given a pion/kaon beam matrix, calculate the resulting neutrino flux from a 2-body decay at the detector.
	"""

	# Find the decay pattern:
	parentPDGid = beamMat[0,7]
	parentDecay = pdgData[parentPDGid]['decay']
	nuFlux = []

	#
	count = beamMat.shape[0]
	beamP = np.sqrt(beamMat[:,3]**2+beamMat[:,4]**2+beamMat[:,5]**2)
	beamE = np.sqrt(beamP**2+pdgData[parentPDGid]['mass']**2)
	beamGamma = beamE/pdgData[parentPDGid]['mass']
	beamBeta = np.sqrt(1-1/beamGamma**2)
	beamDirection = np.zeros([count,3])
	for jj in range(3): beamDirection[:,jj] = beamMat[:,jj+3]/beamP
	detDirection = (np.array(detSize)/np.sqrt(detSize[0]**2+detSize[1]**2+detSize[2]**2)).reshape([3,1])
	# Angle between the parent particle direction and the det.
	theta = np.arccos(np.dot(beamDirection,detDirection).reshape([count,]))

	for jj in range(len(parentDecay)):
		if len(parentDecay[jj]) > 5: continue
		if 14 not in map(np.abs,parentDecay[jj][1:]) and 12 not in map(np.abs,parentDecay[jj][1:]): continue
		for kk in range(1,3):
			if np.abs(parentDecay[jj][kk]) != 14 and np.abs(parentDecay[jj][kk]) != 12:
				# Find the children pdgid, including the meson and neutrino.
				childPDGid = parentDecay[jj][kk]
				children = parentDecay[jj][1:-2]
				children.remove(childPDGid)
				nuPDGid = children[0]
		# Store the neutrinos in a N-by-3 matrix: 3=energy,flux,pdgid
		thisnuFlux = np.zeros([beamMat.shape[0],3])

		# Flux formula, can be found in Donega's thesis.
		nuE_rest = (pdgData[parentPDGid]['mass']**2-pdgData[childPDGid]['mass']**2)/(2*pdgData[parentPDGid]['mass'])
		nuE_lab = beamGamma*nuE_rest*(1+beamBeta*(beamBeta-np.cos(theta))/(beamBeta*np.cos(theta)-1))
		probAtDet = 1/(4*np.pi)*(detSize[0]*detSize[1]/detSize[2]**2)*(1-beamBeta**2)/(beamBeta*np.cos(theta)-1)**2

		# Need to take the decay branch ratio into account:
		probAtDet *= parentDecay[jj][-2]

		thisnuFlux[:,0] = nuE_lab
		thisnuFlux[:,1] = probAtDet
		thisnuFlux[:,2] = nuPDGid
		nuFlux.append(thisnuFlux)

	nuFlux = np.vstack(nuFlux)
	return nuFlux

def main():
	# Detector Size.
	# Assuming a 12-by-14 m far detector for DUNE..
	detSize = [12e3,14e3,1300e6]

	allNeutrinos = []
	for jj in sys.argv[1:]:
		os.system('echo "Now processing %s..." >> %s'%(jj,logFileName))
		allNeutrinos.append(calcFlux2Body(loadtxt_fast(jj,3)))
	allNeutrinos = np.vstack(allNeutrinos)
	np.savetxt('./neutrinos_%s'%("_".join(asctime().split())),allNeutrinos)

if __name__ == '__main__':
	main()