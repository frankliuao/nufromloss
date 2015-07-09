from __future__ import division
import numpy as np
import sys
from beam import beam
from loadtxt import loadtxt
from constants import pdgData
"""
This code is to obtain the particle decays in a beamline.
The argument(s) should be the beam loss file(s) from the G4Beamline simulations.
"""

__all__ = ['finddecayed','group_by_PDG']

# Find particles that decayed in the beamline:
def finddecayed(lossBeam):
	lossBeamCopy = np.copy(lossBeam)
	lossBeamCopySorted = lossBeamCopy[np.argsort(lossBeamCopy[:,8]),:]
	numberLoss = lossBeamCopy.shape[0]
	bool_keep = np.zeros(numberLoss)     # Boolean vector that marks which particle to keep
	iterator1 = iterator2 = 0
	while iterator2 <= numberLoss-2:
		iterator2 += 1
		# If there are no more than 2 particles for this eventID, skip;
		if lossBeamCopySorted[iterator1,8] != lossBeamCopySorted[iterator2,8]:
			if iterator2 - iterator1 >= 2:
				# iterator2 has covered all the particles with the same eventID, get them:
				thisID = lossBeamCopySorted[iterator1:iterator2,:]
				bool_keep[np.nonzero(np.in1d(thisID[:,9],thisID[:,10]))[0]+iterator1] = 1
			iterator1 = iterator2
	if iterator1 != numberLoss-1:   # This means at least the last two eventID are the same
		thisID = lossBeamCopySorted[iterator1:,:]
		bool_keep[np.nonzero(np.in1d(thisID[:,9],thisID[:,10]))[0]+iterator1] = 1

	# All the decayed particles in this lossBeam:
	allParents = lossBeamCopySorted[bool_keep.astype(bool)]

	return allParents

def group_by_PDG(beamMat):
	"""Returns a list, with each element being the beam matrix for each PDGid found in beamMat.

	len(return) = len(beamMat[:,7])
	"""
	beamMatSorted = beamMat[np.argsort(beamMat[:,7]),:]
	beamList = []
	lastIndex = 0
	for ii in range(1,beamMat.shape[0]):
		if beamMatSorted[ii,7] != beamMatSorted[ii-1,7]:
			thisPDGbeam = beamMatSorted[lastIndex:ii,:]
			thisPDGbeam = thisPDGbeam[np.argsort(thisPDGbeam[:,8]),:]
			beamList.append(thisPDGbeam)
	return beamList

def main():
	# Summarize all the loss files in G4BL's ASCII format
	particleLossFiles = []
	for ii in range(1,len(sys.argv)):
		particleLossFiles.append(loadtxt(sys.argv[ii],3))

	# Loop thru all the loss files.
	parentsInFiles = []
	for ii in range(len(sys.argv)):
		parentsInFiles.append(finddecayed(particleLossFiles[ii]))

	# stack the beam
	allDecays = np.row_stack(parentsInFiles)

	# Cut the beam by their loss location.
	# allDecays = allDecays[allDecays[:,2]>=58059,:]      # Straight section starts at 58059.

	# For each PDGid, save the beam
	allDecaysGrouped = group_by_PDG(allDecays)
	for ii in allDecaysGrouped:
		thisBeam = beam()
		thisBeam.loadBeam(ii)
		thisBeam.writeBeam('./allDecays_PDG%i.beam'%(ii[0,7]))

if __name__ == '__main__':
	main()