import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LinearSegmentedColormap
from math import *
from pyfits import getdata, PrimaryHDU
import h5py
import numpy as np
import os, sys
from math import sqrt
import astropy as ap
from astropy.io import fits
from pyfits import getdata, PrimaryHDU, getheader
# from pylab import *
import scipy.ndimage as ndimage


cmap='gray_r'
cubename = 'DATACUBE-HDFS-1.35-PROPVAR.fit'


hdu_list = fits.open(cubename)
head = getheader(cubename)
img = np.nanmean(hdu_list[0].data), 0)

cubename = 'DATACUBEFINALuser_20141021T055145_212058ad.fits'


plt.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, extent=[xmin, xmax, ymax, ymin])