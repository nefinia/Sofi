#!/usr/bin/env python
import numpy as np
from pyfits import getdata
from math import pi
import os, subprocess
from tools_sofi import distance  # , hubble
import astropy.io.fits as fits
from astropy import wcs

def l2pix(l):
    """Convert wavelength to pixels"""
    l0 = 4750.
    pixsize = 1.25  # 1 pixel = 1.25 Angstrom
    return (l - l0) / pixsize + 1


def lya_ob(z):
    """Convert lya to observed wavelength"""
    lya_rest = 1215.67
    return lya_rest * (1 + z)

centerfix = True
theia = False
fitcat = 'HDFS'#'HDFS'

if theia:
    folder = '/scratch/gallegos/MUSE/'
else:
    folder = '../../'

cat = '%s/%s/cats/laes.fits'%(folder, fitcat)
data = getdata(cat, 1)
if fitcat=='HDFS':
    conf = data['confidence']
    good = np.where(conf >= 4)[0]
    sconf = data['Sofi_confid'][good]
    ids = data['id'][good]
    ra = data['ra'][good]
    dec = data['dec'][good]
    m300 = data['m300'][good]
    xs = data['X'][good]
    ys = data['Y'][good]
    redshift = data['z'][good]
    zs = l2pix(lya_ob(redshift))
    fs = data['LYALPHA_FLUX'][good]

if fitcat=='UDF':
    conf = data['confidence']
    good = np.where(conf >= 1)[0]
    sconf = data['Sofi_confid'][good]
    ids = data['ID'][good]
    ra = data['RA'][good]
    dec = data['DEC'][good]
    redshift = data['Z_MUSE'][good]
    wavelength = data['LYALPHA_LBDA_OBS'][good]
    fs = data['LYALPHA_FLUX'][good]

pairs = []

if fitcat=='HDFS':
    filename = '%s/%s/DATACUBE-HDFS-1.35-PROPVAR.csub.fits'%(folder, fitcat)
    data_cube, header_data_cube = getdata(filename, 0, header=True)

if fitcat=='UDF':
    filename="%s/%s/DATACUBEFINALuser_20141021T055145_212058ad.fits"%(folder, fitcat)
    hdulist = fits.open(filename)
    data_cube, header_data_cube = getdata(filename, 0, header=True)
    w = wcs.WCS(header_data_cube, hdulist)
    xs, ys, zs = w.all_world2pix(ra, dec, redshift, 1)
    zs = l2pix(wavelength)

if centerfix:
    xf = data['x_fixed'][good]
    yf = data['y_fixed'][good]
    zf = data['z_fixed'][good]
    fix = data['fix'][good]
    xs[fix == 1] = xf[fix == 1]
    ys[fix == 1] = yf[fix == 1]
    zs[fix == 1] = zf[fix == 1]

if 1:
    for id1, x1, y1, z1, ra1, dec1, red1, sc1, f1, i in zip(ids, xs, ys, zs, ra, dec, redshift, sconf, fs, range(len(ra))):
        for id2, x2, y2, z2, ra2, dec2, red2, sc2, f2 in zip(ids[i + 1:], xs[i + 1:], ys[i + 1:], zs[i + 1:], ra[i + 1:], dec[i + 1:], redshift[i + 1:], sconf[i + 1:], fs[i + 1]):
            print 'id1!', id1, 'id2!', id2
            theta, com_dist1, com_dist2, pi_Mpc, pi_v, proj_dist = distance(red1, ra1, dec1, red2, ra2, dec2)
            pairs.append(
                [int(id1), int(id2), theta, com_dist1, com_dist2, pi_Mpc, pi_v, proj_dist, ra1, dec1, red1, x1, y1, z1,
                 ra2, dec2, red2, x2, y2, z2, sc1, sc2])

    np.savetxt("%s/%s/cats/lae_pairs.dat" %(folder, fitcat), pairs,
               header='id1 id2 theta com_dist1 com_dist2 pi_Mpc pi_v proj_dist ra1 dec1 redshift1 x1 y1 z1 f1 ra2 dec2 redshift2 x2 y2 z2 f2 sconf1 sconf2',
               fmt='%.4e')
