import numpy as np
import re
def loadtxt(filename, skiprows=0, delimiter=" ", comment="#"):
    """Works like the loadtxt function in numpy, but much faster.
    """
    def iter_func():
        with open(filename,'r') as inputFile:
            for ii in range(skiprows):
                next(inputFile)
            for line in inputFile:
                if line[0:len(comment)] == comment:
                    next(inputFile)
                line = line.rstrip()
                line = re.sub('%s+'%delimiter, delimiter, line)
                line = line.split(delimiter)
                for column in line:
                    yield float(column)
            loadtxt.numCol = len(line)
    data = np.fromiter(iter_func(), dtype=float)
    data = data.reshape(-1, loadtxt.numCol)
    return data