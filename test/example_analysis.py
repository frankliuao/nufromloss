import numpy as np

from nufromloss import calcnufromloss
from nufromloss.loadtxt import loadtxt

example_lossbeam = loadtxt('./example_loss.beam',3)
calcnufromloss.setdetector([14e3,12e3,1300e6])
calcnufromloss.setcut(40e3)
nu_by_pid = calcnufromloss.getnuflux(example_lossbeam) 
for ii in nu_by_pid.keys():
    np.savetxt('./nu_from%i.txt'%ii,nu_by_pid[ii])
