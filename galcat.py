#!/usr/bin/env python
from sys import argv

import astropy.io.fits as fits  #
import numpy as np
from astropy import wcs
from pyfits import getdata

from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

fitcat = rdarg(argv, key='fitcat', type=str, default='EAGLE')  # 'UDF
snap = rdarg(argv, key='snap', type=int, default=12)  # 'UDF

def l2pix(l):
	"""Convert wavelength to pixels"""
	l0 = 4750.
	pixsize = 1.25  # 1 pixel = 1.25 Angstrom
	return (l - l0) / pixsize + 1


def lya_ob(z):
	"""Convert lya to observed wavelength"""
	lya_rest = 1215.67
	return lya_rest * (1 + z)


folder = '../../'

if fitcat == 'mosaic10':
	cat = '%s/%s/cats/laes.fits' % (folder, 'UDF')
elif fitcat == 'EAGLE':
	cat = '%s/EAGLE/cats/gals_snap%d.fits' % (folder, snap)
else:
	cat = '%s/%s/cats/laes.fits' % (folder, fitcat)
data = getdata(cat, 1)

if fitcat == 'HDFS':
	conf = data['confidence']
	good = np.where(conf >= 4)[0]
	sconf = data['Sofi_confidence'][good]
	ids = data['id'][good]
	ra = data['ra'][good]
	dec = data['dec'][good]
	m300 = data['m300'][good]
	x = data['x'][good]
	y = data['y'][good]
	z = data['z'][good]
	redshift = data['redshift'][good]
	l = data['LYALPHA_LUM'][good]
	l[np.isnan(l)] = 0

if fitcat == 'UDF' or fitcat == 'mosaic10':
	conf = data['CONFID']
	good = np.where(conf >= 1)[0]
	sconf = data['Sofi_confidence'][good]
	ids = data['ID'][good]
	ra = data['RA'][good]
	dec = data['DEC'][good]
	redshift = data['Z_MUSE'][good]
	wavelength = data['LYALPHA_LBDA_OBS'][good]
	if fitcat == 'mosaic10':
		filename = '/net/galaxy-data/export/galaxydata/gallegos/mosaic.z1300.fits'
		hdulist = fits.open(filename)
		ttt, header_data_cube = getdata(filename, 0, header=True)
		# Removing COMMENT key to avoid problems reading non-ascii characters!!!!!!!!!!!!
		cards = header_data_cube.cards
		bad = ['COMMENT' == b[0] for b in cards]
		for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
		w = wcs.WCS(header_data_cube, hdulist)
		x, y, z = w.all_world2pix(ra, dec, wavelength, 1)
		z = l2pix(wavelength)
	else:
		x = data['x'][good]
		y = data['y'][good]
		z = data['z'][good]
	l = data['LYALPHA_LUM'][good]
	# F606W = data['F606W'][good]
	l[np.isnan(l)] = 0

if fitcat == 'mosaic':
	conf = data['Soficonfid']
	good = np.where(conf > 0)[0]
	sconf = conf  # = data['Sofi_confidence'][good]
	ids = data['ID'][good]
	ra = data['RA'][good]
	dec = data['DEC'][good]
	redshift = data['Z_MUSE'][good]
	wavelength = data['LYALPHA_LBDA_OBS'][good]
	filename = '/net/galaxy-data/export/galaxydata/gallegos/mosaic.z1300.fits'
	hdulist = fits.open(filename)
	ttt, header_data_cube = getdata(filename, 0, header=True)
	# Removing COMMENT key to avoid problems reading non-ascii characters!!!!!!!!!!!!
	cards = header_data_cube.cards
	bad = ['COMMENT' == b[0] for b in cards]
	for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
	w = wcs.WCS(header_data_cube, hdulist)
	x, y, z = w.all_world2pix(ra, dec, wavelength, 1)
	z = l2pix(wavelength)
	l = data['LYALPHA_LUM'][good]
	F606W = data['HST_F606W'][good]
	l[np.isnan(l)] = 0

if fitcat == 'EAGLE':
	import params
	ids = data['ID']
	lcube = params.lcube
	coml = params.coml  # cMpc
	com2pix = params.com2pix
	comx = data['x']
	comy = data['y']
	comz = data['z']
	x = comx * com2pix
	y = comy * com2pix
	z = comz * com2pix
	vx = data['vx']
	vy = data['vy']
	vz = data['vz']
	# rmin = 3.007
	# rmax = 3.0325
	# redshift with respect to x coordinate!!!!
	# redshift = (rmax-rmin)*x/float(lcube)+rmin
	# l = z+1
	# kpc2arcsec = 0.12750223128904756 # 7.843  pkpc per arcsec
	# pix2kpc= 6.103515625 #~6 kpc per pixel
	pix2kpc = coml * 1000 / float(lcube)
	asec2kpc = params.asec2kpcs[snap]
	sred = params.redshifts[snap]
	nz = params.zlens[snap]
	pix2deg = pix2kpc / asec2kpc / (1 + sred) / 3600.  # initial units were comoving and so 1/(1+z) should be used
	pix2asec = pix2kpc / asec2kpc / (1 + sred)  # initial units were comoving and so 1/(1+z) should be used

if fitcat == 'EAGLE':
	unames = ('com_dist angx angy angz shearx sheary shearz thetax thetay thetaz xt yt zt').split()
	dnames = ('id x y z vx vy vz U G R I Z Y J H K gass_mass DM_mass stellar_mass xc yc zc').split()

else:
	unames = ('theta proj_pdist proj_cdist pi_Mpc pi_v redshift com_dist angle shear').split()
	dnames = ('id lum com_dist ra dec redshift x y z sconf').split()

n5ds = []
for i in ids:
	i0 = ids==i
	x0, y0, z0 = x[i0][0], y[i0][0], z[i0][0]

	if fitcat=='EAGLE':
		nd = (x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2
		nd = np.sort(nd)
		n5ds.append(np.sqrt(nd[5]*pix2kpc))
	else:
		nd = (x - x0) ** 2 + (y - y0) ** 2
		nd = np.sort(nd)
		n5ds.append(np.sqrt(nd[5])*.2)


try:  data.columns.del_col('n5d')
except: print 'n5d not there'


a = [data['ID']]
fmt = '%d'
h='ID'
for n in data.names[1:]:
	a.append(data[n])
	if n=='SubGroupNumber': fmt += ' %d'
	else: fmt += ' %.3f '
	h += ' ' + n
a.append(n5ds)
fmt += ' %.3f '
h += ' n5d'
np.savetxt('%s/%s/cats/n5d_snap%d.dat' % (folder, fitcat, snap), np.array(a).T, header=h, fmt=fmt)

if fitcat == 'EAGLE': outname = '%s/%s/cats/gals_snap%d.fits' % (folder, fitcat, snap)
else: outname = '%s/%s/cats/laes_n5d.fits' % (folder, fitcat)

#writeto(outname, data, clobber=True)
