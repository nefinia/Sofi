#!/usr/bin/env python
from math import pi, atan2
from sys import argv

import astropy.io.fits as fits  #
import numpy as np
from astropy import wcs
from pyfits import getdata, Column, ColDefs, BinTableHDU

from tools_sofi import distance, rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
snap = rdarg(argv, key='snap', type=int, default=12)  # 'UDF
unique = rdarg(argv, key='unique', type=bool, default=True)  # 'UDF


def rdarg(argv, key, type=None, default=None, listtype=int):
	# print argv, key, type, default

	if len(argv) > 1:
		opt = np.where([a == '-%s' % key for a in argv])[0] + 1
		if len(opt) > 0:
			name = argv[int(opt)]
			if type is list:
				name = name.split(',')
				if listtype == int:
					name = [int(i) for i in name]
				elif listtype == float:
					name = [float(i) for i in name]
			elif type is bool:
				print key, name
				name = eval(str(name))
			elif type is int:
				print key, name
				name = int(name)
			elif type is float:
				name = float(name)
			elif type is str:
				name = str(name)
			return name

	if default is not None:
		return default


def l2pix(l):
	"""Convert wavelength to pixels"""
	l0 = 4750.
	pixsize = 1.25  # 1 pixel = 1.25 Angstrom
	return (l - l0) / pixsize + 1


def lya_ob(z):
	"""Convert lya to observed wavelength"""
	lya_rest = 1215.67
	return lya_rest * (1 + z)


centerfix = False
theia = False

if theia:
	folder = '/scratch/gallegos/MUSE/'
else:
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

ngal = len(ids)
_data = {}
sdat = {}
nw = np.zeros(ngal)
# nw2 = np.zeros(ngal)
nps = np.zeros(ngal)
# np2 = np.zeros(ngal)
n5ds = np.zeros(ngal)

if fitcat == 'EAGLE':
	unames = ('com_dist angx angy angz shearx sheary shearz thetax thetay thetaz xt yt zt').split()
	dnames = ('id x y z vx vy vz U G R I Z Y J H K gass_mass DM_mass stellar_mass xc yc zc').split()

else:
	unames = ('theta proj_pdist proj_cdist pi_Mpc pi_v redshift com_dist angle shear').split()
	dnames = ('id lum com_dist ra dec redshift x y z sconf').split()
	_data['npairs1'] = []
	_data['npairs2'] = []

ufmt = ['D'] * len(unames)
dfmt = np.append(['J'], ['D'] * (len(dnames) - 1))

for name in unames: _data[name] = []
for name in dnames:
	if fitcat == 'EAGLE':
		if name == 'x' or name == 'y' or name == 'z':
			sdat[name] = data[name] * com2pix
		elif name == 'xc' or name == 'yc' or name == 'zc':
			sdat[name] = data[name[:-1]]
		else:
			sdat[name] = data[name]
	else:
		sdat[name] = data[name]
	_data['%s1' % name] = []
	_data['%s2' % name] = []

_data['nw1'] = []
_data['nw2'] = []

_data['n5d1'] = []
_data['n5d2'] = []

nn = range(ngal)
dens = 1. / float(ngal)

for i in nn:
	nd = []
	if unique:
		_nn = nn[i:]
	else:
		_nn = nn
	for j in _nn:
		if i != j:
			print 'id1!', ids[i], 'id2!', ids[j]
			if fitcat == 'EAGLE':
				xt = (x[j] - x[i] + lcube / 2) % lcube
				yt = (y[j] - y[i] + lcube / 2) % lcube
				zt = (z[j] - z[i] + lcube / 2) % lcube
				angx = atan2(zt, yt)
				angy = atan2(xt, zt)
				angz = atan2(yt, xt)
				shearx = (xt - lcube / 2) * nz / float(lcube)
				sheary = (yt - lcube / 2) * nz / float(lcube)
				shearz = (zt - lcube / 2) * nz / float(lcube)

				ra2 = lcube * pix2deg / 2
				dec2 = ra2
				ra1 = yt * pix2deg
				dec1 = zt * pix2deg
				thetax = 3600 * np.sqrt(
					((ra1 - ra2) * np.cos((dec1 + dec2) / 2. * pi / 180.)) ** 2 + (dec1 - dec2) ** 2)
				ra1 = zt * pix2deg
				dec1 = xt * pix2deg
				thetay = 3600 * np.sqrt(
					((ra1 - ra2) * np.cos((dec1 + dec2) / 2. * pi / 180.)) ** 2 + (dec1 - dec2) ** 2)
				ra1 = xt * pix2deg
				dec1 = yt * pix2deg
				thetaz = 3600 * np.sqrt(
					((ra1 - ra2) * np.cos((dec1 + dec2) / 2. * pi / 180.)) ** 2 + (dec1 - dec2) ** 2)
				com_dist = np.sqrt((xt - lcube / 2) ** 2 + (yt - lcube / 2) ** 2 + (zt - lcube / 2) ** 2) / com2pix

				udat = [com_dist, angx, angy, angz, shearx, sheary, shearz, thetax, thetay, thetaz, xt, yt, zt]
				nd.append(com_dist)

			else:
				ang = atan2(y[j] - y[i], x[j] - x[i])
				shear = z[i] - z[j]
				theta, com_dist1, com_dist2, pi_Mpc, pi_v, proj_dist = \
					distance(redshift[i], ra[i], dec[i], redshift[j], ra[j], dec[j])
				z1 = max(redshift[i], redshift[j])
				z2 = min(redshift[i], redshift[j])
				z12 = abs((1 + z1) / (1 + z2) - 1)
				proj_cdist = proj_dist * (1. + z2)
				com_dist = np.sqrt(proj_cdist * proj_cdist + pi_Mpc * pi_Mpc)
				udat = [theta, proj_dist, proj_cdist, pi_Mpc, pi_v, z12, com_dist, ang, shear]
				nd.append(com_dist)
			for name, d in zip(unames, udat): _data[name].append(d)
			for name in dnames:
				_data['%s1' % name].append(sdat[name][i])
				_data['%s2' % name].append(sdat[name][j])

			if 20 > com_dist > 0.5:
				if fitcat != 'EAGLE':
					if not unique: nps[i] += 1
				# np2[j] += 1
				if not unique: nw[i] += dens / com_dist
			# nw2[j] += dens/com_dist
	if not unique:
		nd = np.sort(nd)
		n5ds[i] = nd[4]
c = []

for name, df in zip(dnames, dfmt):
	print name
	n1 = '%s1' % name
	n2 = '%s2' % name
	c.append(Column(name=n1, format=df, array=_data[n1]))
	c.append(Column(name=n2, format=df, array=_data[n2]))
for name, uf in zip(unames, ufmt):
	c.append(Column(name=name, format=uf, array=_data[name]))

if not unique:
	for i in nn:
		for j in nn:
			if i != j:
				if fitcat != 'EAGLE':
					_data['npairs1'].append(nps[i])
					_data['npairs2'].append(nps[j])
				_data['nw1'].append(nw[i])
				_data['nw2'].append(nw[j])
				_data['n5d1'].append(n5ds[i])
				_data['n5d2'].append(n5ds[j])

	if fitcat != 'EAGLE':
		c.append(Column(name='npairs1', format='J', array=_data['npairs1']))
		c.append(Column(name='npairs2', format='J', array=_data['npairs2']))
	c.append(Column(name='nw1', format='D', array=_data['nw1']))
	c.append(Column(name='nw2', format='D', array=_data['nw2']))
	c.append(Column(name='d5th1', format='D', array=_data['n5d1']))
	c.append(Column(name='d5th2', format='D', array=_data['n5d2']))

coldefs = ColDefs(c)

tbhdu = BinTableHDU.from_columns(coldefs)
if fitcat == 'EAGLE':
	outname = '%s/%s/cats/lae_pairs_snap%d.fits' % (folder, fitcat, snap)
else:
	outname = '%s/%s/cats/lae_pairs%s.fits' % (folder, fitcat, '_unique'*unique)

tbhdu.writeto(outname, clobber=True)
