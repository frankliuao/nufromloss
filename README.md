# nufromloss
A Python package that generates the neutrino flux at a given FAR detector from the G4Beamline tracking result.

Requires:
    1. The detector at which the neutrino flux is desired should be far enough from the particle decay points, that the size of the detector is much smaller compared to the distance between the detector and the decay points.
    2. In the G4Beamline tracking, one is allowed to select the particle kinds to be kept. One is also allowed to save the particle loss in an ASCII beam file. The loss file will record every particle loss incurred by particle hitting the apertures, or by decaying. In order to calculate the nu flux contribution from particle kind PIDa, at least ONE of the offspring particle kinds need to be kept in the loss file. e.g. If the neutrino flux from the pi+ decay to mu+ and nu-mu is needed, at least one of mu+ and nu-mu need to be present in the loss file. The loss file does not need to be processed further after G4BL writes it.
    3. A longitudinal coordinate is needed to filter any decays before that coordinate.

To use:
Add the folder to your PYTHONPATH.
In your Python code, import nufromloss.
Set the detector information by:
    nufromloss.setdetector(detector_size), where detector_size is a list, with the first and second elements being the width and height of the detector, and the thrid element is the longitudinal coordinate of the detector. The units are mm. e.g. nufromloss.setdetector([12e3,14e3,1300e6])
    nufromloss.setcut(cut_before_z), where cut_before_z is a float. All the particles that decayed before this point won't be used for neutrino production.
    nu_from_pid = nufromloss.getnuflux(lossbeam), where lossbeam is a N-by-12 G4Beamline loss beam array. You can obtain the lossbeam by using nufromloss.loadtxt('YOUR_LOSS_FILE',3). 
    nu_from_pid will be a dictionary. Each key of nu_from_pid is one PID of the particles that decayed in the lossfile. If this PID does not produce any neutrinos at the FD, the value of its key will be 0. Otherwise, the value will be a numpy ndarray with the first column being the neutrino energies, the second being the neutrino weights, and the third being the neutrino PIDs.

!! An example case is included in the ./example/ folder.
