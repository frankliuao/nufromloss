"""
Calculate the neutrino flux at a far detector.


"""

from __future__ import division
import numpy as np
from loadtxt import loadtxt
from calc_nu import calcFlux2Body
import analyzeDecay
from time import sleep,asctime

__all__ = ['setcut','setdetector','getnuflux']

_detectorSize = np.zeros(3)
_cutFrom = -1

def setdetector(detSize):
    """Set the detector size

    detSize should be a list with three elements,
    or an numpy ndarray of the shape [3,]
    """
    try:
        detSizeCp = np.array(detSize)
        if detSizeCp.shape != (3,):
            print "detSize should have 3 elements."
            return
        else:
            global _detectorSize
            _detectorSize = detSizeCp

    except Exception as e:
        print e
        print "Can not convert detSize to an numpy ndarray"

def setcut(cutFromZ):
    """Set the longitudianl z for the cut

    beam before this z will be discarded.
    """
    global _cutFrom
    _cutFrom = cutFromZ

def getnuflux(lossbeam):
    """Calculate the flux at the FD, directly based on the loss beam in the beamline

    setdetector and setcut must have already been used before getnuflux.
    """
    if _cutFrom == -1 or len(np.nonzero(_detectorSize)[0]) == 0:
        print "setdetector and setcut first."
        return None
    parents = analyzeDecay.finddecayed(lossBeam=lossbeam)
    parents = parents[parents[:,2]>=_cutFrom,:]
    parentsGroup = analyzeDecay.group_by_PDG(parents)
    neutrinoGroup = {}
    for ii in parentsGroup:
        neutrinoGroup[int(ii[0,7])] = calcFlux2Body(ii,_detectorSize)
    return neutrinoGroup