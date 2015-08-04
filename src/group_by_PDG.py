#coding=ASCII
"""Group the particles in the beam array by the PID into a list.

Given a G4Beamline ASCII beam array (numpy ndarray, N-by-12), group the
particles by their PID numbers.

------
Parameters:
    beam_array: numpy ndarray (N-by-12), a G4Beamline standard ASCII beam
    array.

------
Return:
    beam_list: a list. Each list member is an individual beam array (
    N-by-12), in which all particles have the same PID.

"""

from __future__ import division
import numpy as np

def group_by_PDG(beam_array):
    """Group the beam by the particle PDGid.

    ------
    Parameter:
        beam_array: numpy ndarray (N-by-12), standard G4BL ASCII beam format;

    ------
    Return:
    beam_list, a list, [beam array for PDGid 1, for PDGid 2, ...]

    """
    beam_array_sorted = beam_array[np.argsort(beam_array[:, 7]), :]
    beam_list = []
    iter_index = 0
    for ii in range(1, beam_array.shape[0]):
        if beam_array_sorted[ii, 7] != beam_array_sorted[ii-1, 7]:
            this_PDG = beam_array_sorted[iter_index:ii, :]
            this_PDG = this_PDG[np.argsort(this_PDG[:, 8]), :]
            beam_list.append(this_PDG)
            iter_index = ii
    last_PDG = beam_array_sorted[iter_index:, :]
    last_PDG = last_PDG[np.argsort(last_PDG[:, 8]), :]
    beam_list.append(last_PDG)

    return beam_list