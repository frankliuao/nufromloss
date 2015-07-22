# coding=ASCII
"""Calculate neutrino flux from a G4Beamline beam array.

Modification date: 07/22/15
"""
from __future__ import division
import numpy as np
import os
import sys
from time import asctime

from constants import pdgData
from loadtxt import loadtxt

__all__ = ["flux_2body"]


def flux_2body(beam_array, det_size):
    """Given a pion/kaon beam array, calculate the resulting neutrino flux
    from a 2-body decay at the detector.
    """

    # Find the decay pattern:
    parentPDGid = beam_array[0,7]
    parent_decay_pattern = pdgData[parentPDGid]["decay"]
    tot_nu_flux = []

    #
    particle_quantity = beam_array.shape[0]
    beamP = np.sqrt(beam_array[:,3]**2 + beam_array[:,4]**2 +
                    beam_array[:,5]**2)
    beamE = np.sqrt(beamP**2 + pdgData[parentPDGid]["mass"]**2)
    beam_gamma = beamE / pdgData[parentPDGid]["mass"]
    beam_beta = np.sqrt(1 - 1/beam_gamma**2)
    beam_direction = np.zeros( [particle_quantity, 3] )
    for jj in range(3):
        beam_direction[:,jj] = beam_array[:, jj+3] / beamP
    '''det_direction is approximated to be [0, 0, 1]
    det_direction = (np.array(det_size) /
                     np.sqrt(det_size[0]**2 +
                             det_size[1]**2 +
                             det_size[2]**2)).reshape([3, 1])'''
    det_direction = np.array([0, 0, 1]).reshape([3, 1])
    # Angle between the parent particle direction and the det.
    theta = np.arccos(np.dot(beam_direction,
                              det_direction).reshape([particle_quantity,]))

    for jj in range(len(parent_decay_pattern)):
        if len(parent_decay_pattern[jj]) > 5:
            continue
        if 14 not in map(np.abs, parent_decay_pattern[jj][1:]) \
           and \
           12 not in map(np.abs, parent_decay_pattern[jj][1:]):
            continue
        for kk in range(1, 3):
            # Assumes (requires) at least one non-neutrino daughter of the
            # decay is kept in the loss file.
            if np.abs(parent_decay_pattern[jj][kk]) != 14 and \
               np.abs(parent_decay_pattern[jj][kk]) != 12:
                # Find the children pdgid, including the meson and neutrino.
                childPDGid = parent_decay_pattern[jj][kk]
                children = parent_decay_pattern[jj][1:-2]
                children.remove(childPDGid)
                nuPDGid = children[0]
        # Store the neutrinos in a N-by-3 matrix: columns = energy,flux,pdgid
        this_nu_flux = np.zeros([beam_array.shape[0],3])

        # Flux formula, can be found in Donega's thesis.
        nuE_rest = (pdgData[parentPDGid]["mass"]**2-
                    pdgData[childPDGid]["mass"]**2) / \
                   (2*pdgData[parentPDGid]["mass"])
        nuE_lab = beam_gamma * nuE_rest * \
                  (1+beam_beta*(beam_beta-np.cos(theta)) /
                   (beam_beta*np.cos(theta)-1))
        prob_hit_det = 1 / (4*np.pi) * \
                       (det_size[0]*det_size[1]/det_size[2]**2) *\
                       (1-beam_beta**2)/(beam_beta*np.cos(theta)-1)**2

        # Need to take the decay branch ratio into account:
        prob_hit_det *= parent_decay_pattern[jj][-2]

        # Also the particle weights:
        prob_hit_det *= beam_array[:, 11]

        this_nu_flux[:,0] = nuE_lab
        this_nu_flux[:,1] = prob_hit_det
        this_nu_flux[:,2] = nuPDGid
        tot_nu_flux.append(this_nu_flux)

    if len(tot_nu_flux) >= 1:
        tot_nu_flux = np.vstack(tot_nu_flux)

    return tot_nu_flux


def main(det_size):
    """Calculate the neutrino flux at a given detector from 2-body decays.

    ------
    Parameter:
        det_size: list, [width, height, distance] of the detector.

    ------
    Return:
        tot_nu_flux: N-by-3 numpy ndarray, columns = nu energy, flux, nu PDGid

    """

    all_nu_files = []
    for jj in sys.argv[1:]:
        os.system('echo "Now processing %s..."'%jj)
        all_nu_files.append(flux_2body(loadtxt(jj, 3), det_size))
    all_nu_files = np.vstack(all_nu_files)
    np.savetxt('./neutrinos_%s'%("_".join(asctime().split())), all_nu_files)