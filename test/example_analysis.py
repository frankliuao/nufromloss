#coding=ASCII
"""Test the nufromloss module using a loss file example.
"""
import numpy as np

import nufromloss

nufromloss.set_detector([14e3, 12e3, 1300e6])
nufromloss.set_cut(40e3)
nu_by_pid = nufromloss.get_flux('example_loss.beam')
for ii in nu_by_pid.keys():
    np.savetxt('./nu_from%i.txt'%ii, nu_by_pid[ii])
