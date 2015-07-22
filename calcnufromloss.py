#coding=ASCII
"""Calculate the neutrino flux at a far detector from a G4Beamline loss file.

Install the package and add the path to the PYTHONPATH environment variable.
In the Python code, do

import nufromloss

nufromloss.set_detector(det_size), where det_size is a list = [width,
height, distance] of the detector

nufromloss.set_startfrom(z), where z is a float = G4Beamline longitudinal z to
start calculating contributions to the nu flux. Particles before this z will
be discarded.

nufromloss.get_flux(lossfile)
Parameter:
    lossfile: str, the path to the G4BL loss file.
Return: A dictionary. Each key is the PDGid of the parent particle, and the
value is a N-by-3 numpy, columns = nu energy, flux, nu PDGid

"""

from __future__ import division
import numpy as np
from time import sleep,asctime

from src.loadtxt import loadtxt
from src.flux_kinematics import flux_2body

__all__ = ['set_cut', 'set_detector', "get_flux", "find_decay", "group_by_PDG"]

_det_size = np.zeros(3)
_cut_from_z = -1


def set_detector(det_size):
    """Set the detector size

    detSize should be a list with three elements,
    or an numpy ndarray of the shape [3,]
    """
    try:
        det_size_cp = np.array(det_size)
        if det_size_cp.shape != (3,):
            print "detSize should have 3 elements."
            return
        else:
            global _det_size
            _det_size = det_size_cp

    except Exception as e:
        print e
        print "Can not convert detSize to an numpy ndarray"


def set_cut(z):
    """Set the longitudianl z for the cut

    beam before this z will be discarded.
    """
    global _cut_from_z
    _cut_from_z = z


def get_flux(lossfile):
    """Calculate the flux at the FD, directly based on the loss beam in the beamline

    set_detector and set_cut must have already been used before getnuflux.
    """
    if _cutFrom == -1 or len(np.nonzero(_detectorSize)[0]) == 0:
        print "setdetector and setcut first."
        return None
    parents = find_decay(lossfile)
    parents = parents[parents[:,2]>=_cut_from_z, :]
    parents_list = group_by_PDG(parents)
    neutrino_group = {}
    for ii in parents_list:
        neutrino_group[int(ii[0, 7])] = flux_2body(ii, _det_size)
    return neutrino_group


def find_decay(lossfile):
    """Find the particles that decayed in a G4Beamline simulation.

    At least one daughter of the decay needs to be recorded in the loss file
    when doing the G4BL simulation.
    """
    lossbeam = loadtxt(lossfile)
    lossbeam_sorted = lossbeam[np.argsort(lossbeam[:, 8]), :]
    number_of_loss = lossbeam.shape[0]

    # Boolean vector that marks which particle to keep
    bool_keep = np.zeros(number_of_loss)

    # Iterating through the beam is much faster than filtering out each event.
    iterator1 = iterator2 = 0
    while iterator2 <= number_of_loss-2:
        iterator2 += 1

        # If there are no more than 2 particles for this eventID, skip;
        if lossbeam_sorted[iterator1, 8] != lossbeam_sorted[iterator2, 8]:
            if iterator2 - iterator1 >= 2:
                # iterator2 - iterator1 covers all the particles with the
                # same eventID, get them:
                this_event = lossbeam_sorted[iterator1:iterator2, :]
                bool_keep[np.nonzero(np.in1d(this_event[:, 9],
                                             this_event[:, 10]))[0] +
                          iterator1] = 1
            iterator1 = iterator2

    if iterator1 != number_of_loss - 1:
        # This means at least the last two eventID are the same
        this_event = lossbeam_sorted[iterator1:, :]
        bool_keep[np.nonzero(np.in1d(thisID[:, 9],
                                     thisID[:,10]))[0] +
                  iterator1] = 1

    # All the decayed particles in this lossbeam:
    all_parents = lossbeam_sorted[bool_keep.astype(bool), :]

    return allParents


def group_by_PDG(beam_array):
    """Group the beam by the particle PDGid.

    ------
    Parameter:
        beam_array: numpy ndarray, N-by-12, standard G4BL ASCII beam format;

    ------
    Return:
    beam_list, a list, [beam_array for PDGid 1, for PDGid 2, ...]

    """
    beam_array_sorted = beam_array[np.argsort(beam_array[:, 7]), :]
    beam_list = []
    iter_index = 0
    for ii in range(1, beam_array.shape[0]):
        if beam_array_sorted[ii, 7] != beam_array_sorted[ii-1, 7]:
            this_PDG = beam_array_sorted[iter_index:ii, :]
            this_PDG = this_PDG[np.argsort(this_PDG[:, 8]), :]
            beam_list.append(this_PDG)
    last_PDG = beam_array_sorted[iter_index:, :]
    last_PDG = last_PDG[np.argsort(last_PDG[:, 8]), :]
    beam_list.append(last_PDG)

    return beam_list