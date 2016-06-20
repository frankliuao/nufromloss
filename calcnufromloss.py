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

from src.loadtxt import loadtxt
from src.flux_kinematics import flux_2body, flux_mu
from src.group_by_PDG import group_by_PDG

__all__ = ['set_cut', 'set_detector', "get_flux", "find_decay"]

_det_size = np.zeros(3)
_cut_from_z = None


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
    if _cut_from_z == None or len(np.nonzero(_det_size)[0]) == 0:
        print "setdetector and setcut first."
        return None
    parents_list = find_decay(lossfile)
    parents_list_cut = [ii_parent[ii_parent[:, 2]>=_cut_from_z, :]
                        for ii_parent in parents_list]
    parents_list = parents_list_cut
    neutrino_group = {}
    for ii_parent in parents_list:
        # for each type of parent particle, get its decay products
        if len(ii_parent) == 0: continue
        neutrino_group[int(ii_parent[0, 7])] = []
        flux_2body_tmp = flux_2body(ii_parent, _det_size)
        if flux_2body_tmp is not None:
            neutrino_group[int(ii_parent[0, 7])].append(flux_2body_tmp)
        flux_mu_tmp = flux_mu(ii_parent, _det_size)
        if flux_mu_tmp is not None:
            neutrino_group[int(ii_parent[0, 7])].append(flux_mu_tmp)
        if len(neutrino_group[int(ii_parent[0, 7])]) > 0:
            neutrino_group[int(ii_parent[0, 7])] = \
                np.row_stack(neutrino_group[int(ii_parent[0, 7])])
        else:
            neutrino_group[int(ii_parent[0, 7])] = np.array([])
    return neutrino_group


def find_decay(lossfile):
    """Find the particles that decayed in a G4Beamline simulation.

    At least one daughter of the decay needs to be kept in the loss file
    when doing the G4BL simulation.

    ------
    Parameters:
    lossfile, str, the path of the G4BL loss file

    ------
    Returns:
    all_parents_list, a list of numpy ndarrys (N-by-12), each list member is
    a kind of the decayed parent particles at the points of the decay,
    in G4Beamline ASCII beam format.

    """
    lossbeam = loadtxt(lossfile, 3)
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
        bool_keep[np.nonzero(np.in1d(this_event[:, 9],
                                     this_event[:,10]))[0] +
                  iterator1] = 1

    # All the decayed particles in this lossbeam:
    all_parents = lossbeam_sorted[bool_keep.astype(bool), :]
    all_parents_list = group_by_PDG(all_parents)

    return all_parents_list
