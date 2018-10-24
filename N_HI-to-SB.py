#!/usr/bin/env python
__author__ = 'nefinia'
import numpy as np
from pyfits import getdata, PrimaryHDU
cubename = '../../EAGLE/Column_density_base_grid_512_npref_12_uvb_ss.ts0000.fits'
cube = getdata(cubename, 0)
cubelog = np.log10(cube)
low = cubelog<=17.5
mid = (cubelog>17.5) & (cubelog<=20.)
high = cubelog>20.
SB=2.4
cube[low] = SB*(cubelog[low]*.1145-1.7737)
cube[mid] = SB*(cubelog[mid]*.309-5.16)
cube[high] = SB
cube[cube<0] = 0
hdu = PrimaryHDU()
#hdu.data = cube
#hdu.writeto('../../EAGLE/SBcube.fits')
pix2arcsec = 0.1937778538085323
hdu.data = cube * pix2arcsec**2
hdu.writeto('../../EAGLE/Fluxcube.fits')