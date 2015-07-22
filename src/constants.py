import urllib2

def get_pdg_data():
	""" Get the pdg data from Ao Liu (frankliuao)'s website, get the pdgData dict.
	To get all the decay patterns of a pi+, use pdgData[211]["decay"], the returned data is a list, e.g.
	 [['munu', -13, 14, 0.999877, 0.999877], ['enu', -11, 12, 0.000123, 1]]
	"""
	pdgDataFile = urllib2.urlopen('http://frankliuao.com/downloads/Codes/pdg_py.dict')
	pdgDataFileCont = pdgDataFile.read()
	pdgDataFile.close()
	exec(pdgDataFileCont)
	return pdgData

pdgData = get_pdg_data()
