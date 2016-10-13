__author__ = 'gallegos'

from pyfits import getdata
import numpy as np
import os, subprocess
from math import *

#Open mask file
data_cube, header_data_cube = getdata(fits, 0, header=True)
data_cube = data_cube[0]
umasked = np.where(data_cube==0)

#Open sky file
sky = np.loadtxt('skyspec.txt').T
z = sky[0].astype('int')
wav = sky[1].astype('float')
flux = sky[2].astype('float')
sigma = np.std(flux)
