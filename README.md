# nufromloss
A Python package that generates the neutrino flux at a given FAR detector from the G4Beamline tracking result.

Requires:
    1. The detector at which the neutrino flux is desired should be far enough from the particle decay points, that the size of the detector is much smaller compared to the distance between the detector and the decay points.
    2. In the G4Beamline tracking, one is allowed to select the particle kinds to be kept. One is also allowed to save the particle loss in an ASCII beam file. The loss file will record every particle loss incurred by particle hitting the apertures, or by decaying. In order to calculate the nu flux contribution from particle kind PIDa, at least ONE of the offspring particle kinds need to be kept in the loss file. e.g. If the neutrino flux from the pi+ decay to mu+ and nu-mu is needed, at least one of mu+ and nu-mu need to be present in the loss file. The loss file does not need to be processed further after G4BL writes it.
    3. 

To use:
Add the path to this folder to your PYTHONPATH.
In your Python code, do 

	import nufromloss

Set the detector information by:
	nufromloss.set_detector(detector_size)

where detector_size = list, [width, height, distance] of the detector. The units are mm. e.g. for DUNE's far detector, use nufromloss.setdetector([12e3,14e3,1300e6])

To select the decays beyond a certain longitudinal coordinate, do:
	nufromloss.set_cut(z), where z is a float. All the particles that decayed before this z will be discarded.

Get all the neutrino flux from a loss file by:
	nu_from_pid = nufromloss.get_flux(lossfile)
nu_from_pid is a dictionary. Each key (int) of nu_from_pid is PID of a kind of parent particle, while the value (numpy ndarray, columns = nu energy, flux, nu PDGid) gives the neutrino flux from all the decayed particles in the loss file.

######
!! A run example is included in the ./test/ folder.
In that folder, an example loss beam file is given. The Python script calls the nufromloss module to calculate all the neutrino flux genereated at a far detector.

The example uses a G4Beamline loss beam, which includes at least one of the daughters in any decay pattern. If instead of using a loss beam file, one wants to use a "parent beam" only, the following tricks can be done. 
Option 1:
    In G4Beamline, set the decay lifetime of the particle you are intestigating to a very small value, so they decay right away when they are put in a simulation. Keep at least one child particle. Then use the loss file with both parent + child(ren) as the argument of get_flux
Option 2:
    Trick the Python code. The nufromloss expects at least one daughter in the decay is kept for each particle decay, so if one can add the children mannually, e.g.
    lossile before adding children:
	    0 0 0 0 0 0 0 211 1 0 0 1
    lossfile after adding:
	    0 0 0 0 0 0 0 211 1 0 0 1
            0 0 0 0 0 0 0 -13 1 1 0 1
	    This is because the code looks for an event's daughters to confirm it really decayed. Therefore, if you provide a daughter manually (the x,y,z,px,py,pz don't need to be accurate, it's just showing the child exists), and set the event number the same as the parent, but change the parentID of the child to the parent's trackID, the parent will be treated as "decayed" 
######
