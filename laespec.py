#!/usr/bin/env python
import os
from sys import argv

import astropy.io.fits as fits
import numpy as np
from astropy import wcs
from pyfits import getdata, PrimaryHDU

from tools_sofi import rdarg


def l2pix(l):
	"""Convert wavelength to pixels"""
	l0 = 4750.
	pixsize = 1.25  # 1 pixel = 1.25 Angstrom
	return (l - l0) / pixsize


def lya_ob(z):
	"""Convert lya redshift to observed wavelength"""
	lya_rest = 1215.67
	return lya_rest * (1 + z)


simplestack = rdarg(argv, 'stack', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
hdu = PrimaryHDU()
pairs = []
makeim = rdarg(argv, 'makeim', bool, False)
makefit = rdarg(argv, 'makefit', bool, True)
mask = rdarg(argv, 'mask', bool, False)
gmask = rdarg(argv, 'gmask', bool, False)
model = rdarg(argv, 'model', str, 'HM01')
offset = rdarg(argv, 'offset', bool, False)
centerfix = rdarg(argv, 'centerfix', bool, False)
cubex = rdarg(argv, 'cubex', bool, False)
corr = rdarg(argv, 'corr', bool, False)
scorr = '.corr' * corr
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
folder = rdarg(argv, key='folder', type=str, default='/net/galaxy-data/export/galaxydata/gallegos/')
foldercat = rdarg(argv, key='folder', type=str, default='../../')
extraname = rdarg(argv, key='extraname', type=str, default='')
csub = rdarg(argv, 'csub', bool, True)
single = rdarg(argv, 'single', bool, False)
id0 = rdarg(argv, 'id', int, 1)
galrad = rdarg(argv, 'rad', int, 3)
smooth = rdarg(argv, 'smooth', bool, True)
parallel = rdarg(argv, 'parallel', bool, True)
random = rdarg(argv, 'random', bool, False)
nrand = rdarg(argv, 'nrand', int, 1)
line = rdarg(argv, 'line', str, 'lya')
type = rdarg(argv, 'type', str, 'LLS')
white = rdarg(argv, 'white', bool, False)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)

coord = rdarg(argv, 'coord', str, 'x')
snap = rdarg(argv, 'snap', int, 12)

verbose = rdarg(argv, 'verbose', int, 1)
ext = '.fits'
scsub = '.csub' * csub
if verbose < 2:
	vb = ' > /dev/null 2>&1'
else:
	vb = ''

dovar = True
xw = rdarg(argv, 'xw', int, 200)  # for EAGLE it was last time 500 (15.10.18)
yw = rdarg(argv, 'yw', int, 200)  # for EAGLE it was last time 500 (15.10.18)
zw = rdarg(argv, 'zw', int, 50)
vmin = rdarg(argv, 'vmin', float, -.5)
vmax = rdarg(argv, 'vmax', float, 1)
std = rdarg(argv, 'std', float, None)
binsize = 1
if xw%2 == 1: xw -= 1
if yw%2 == 1: yw -= 1
if zw%2 == 1: zw -= 1
xmin = -xw / 2
xmax = xw / 2
ymin = -yw / 2
ymax = yw / 2
zmin = -zw / 2  # min(z)
zmax = zw / 2  # max(z)
lenx = xw / binsize + 1
leny = yw / binsize + 1
lenz = zw + 1

hname = None
cubesmooth = False
if cubesmooth:
	ssmooth = ' -boxsm 3 '
	smooth = False
else:
	ssmooth = ''

if fitcat == 'EAGLE':
	csub = False
	smooth = False
	layers = {'8': 45, '9': 39, '10': 34, '11': 29, '12': 25, '13': 20, '14': 18, '15': 16}
	cubename = 'snap%s_%s_%s_512_4096_%d.%s.fits' % (snap, coord, model, layers[str(snap)], type)
	if mask: maskname = 'snap%s_%s.mask.fits' % (snap, coord)
	if gmask:
		gmaskname = 'snap%s_%s.galmask_%darcsec.fits' % (snap, coord, galrad)
	nconf = 1

if fitcat == 'HDFS':
	#cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s.fits' % '.csub' * csub
	cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s%s.fits' % (scsub, scorr)
	if mask: maskname = 'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id.fits'
	if gmask: gmaskname = 'HDFS.galmask_%darcsec.fits' % galrad
	if dovar:
		varname = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
	if mask: maskname = 'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id.fits'
	nconf = 1

if fitcat == 'UDF':
	#cubename = 'UDF.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'UDF10.z1300%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-10%s.fits' % '.csub' * csub
	if gmask: gmaskname = 'UDF.galmask_%darcsec.fits' % galrad
	cubename = 'DATACUBE_UDF-10%s%s.fits' % (scsub, scorr)
	print cubename
	hname = '%s/%s/%s' % (folder, fitcat, 'DATACUBE_UDF-10.fits')
	if mask: maskname = 'DATACUBE_UDF-10.IM.Objects_Id.fits'
	nconf = 1

if fitcat == 'mosaic' or fitcat == 'mosaic10':
	offset = False
	#cubename = '/net/galaxy-data/export/galaxydata/gallegos/mosaic/DATACUBE_UDF-MOSAIC.z1300%s.fits' % '.csub' * csub
	#cubename = 'mosaic.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-MOSAIC.z1300%s.corr.fits' % scsub
	cubename = 'DATACUBE_UDF-MOSAIC%s%s.fits' % (scsub, scorr)
	hname = '%s/%s/%s' % (folder, fitcat, 'DATACUBE_UDF-MOSAIC.fits')
	if mask: maskname = 'DATACUBE_UDF-MOSAIC.IM.Objects_Id.fits'
	if gmask: gmaskname = 'mosaic.galmask_%darcsec.fits' % galrad
	nconf = 1
	vmin = -3
	vmax = 10
	if std is None: std = 100

filename = '%s/%s/%s' % (folder, fitcat, cubename)

data_cube = getdata(filename, 0)
zlim, ylim, xlim = data_cube.shape
if mask:
	if fitcat == 'EAGLE': mcube = getdata('%s/%s/%s' % (folder, fitcat, maskname), 0)
	else: mcube = getdata('%s/%s/%s' % (folder, fitcat, maskname), 0)[0]
if gmask:
	gcube = getdata('%s/%s/%s' % (folder, fitcat, gmaskname), 0)
if fitcat == 'HDFS' or fitcat == 'mosaic' or fitcat == 'mosaic10':
	cut = 12
else:
	cut = 0

xpix = []
ypix = []
zpix = []
xr = np.arange(xw + 1)
yr = np.arange(yw + 1)
zr = np.arange(zw + 1)

if fitcat == 'mosaic10':
	cat = '%s/%s/cats/laes.fits' % (foldercat, 'UDF')
elif fitcat == 'EAGLE':
	cat = '%sEAGLE/cats/gals_snap%d.fits' % (foldercat, snap)
else:
	cat = '%s/%s/cats/laes.fits' % (foldercat, fitcat)
if centerfix:
	catout = '%s/%s/cats/laesfixed.txt' % (foldercat, fitcat)
	fout = open(catout, 'w')
	fout.write(
		'#id ra dec x y z flux_lya ra_CubEx dec_CubEx lambda_CubEx x_CubEx y_CubEx z_CubEx flux_lya_CubEx redshift_CubEx lya_lum_CubEx diff com_dist\n')
data = getdata(cat, 1)

if fitcat == 'HDFS':
	conf = data['confidence']
if fitcat == 'UDF' or fitcat == 'mosaic':
	conf = data['CONFID']


if fitcat == 'EAGLE':
	csub = False
	lcom = 25.  # comoving length of the simulation
	lcube = 4096
	conv = lcube / lcom
	ids = data['id']
	wavelength = np.zeros(len(ids)) + 1
	flya = np.zeros(len(ids)) + 1
	ra = np.zeros(len(ids)) + 1
	dec = np.zeros(len(ids)) + 1
	cs = ['x', 'y', 'z']
	cs.remove(coord)
	xs = np.round(data[cs[0]] * xlim / lcom).astype(int)
	ys = np.round(data[cs[1]] * ylim / lcom).astype(int)
	zs = (data[coord] * zlim / lcom).astype(int)

	com2pix = 163.84  # lcube/coml
	kpc2pix = lcube / float(lcom * 1e3)
	off = np.zeros(len(ids)).astype(int)

else:
	good = np.where(conf >= nconf)[0]

if fitcat == 'HDFS':
	m300 = data['m300'][good]
	redshift = data['z'][good]
	wavelength = lya_ob(redshift)
	flya = data['LYALPHA_FLUX'][good]
	ids = data['id'][good]
	ra = data['ra'][good]
	dec = data['dec'][good]
	if offset: off = data['offset'][good]
	else: off = ids - ids

if fitcat == 'UDF' or fitcat == 'mosaic':
	ids = data['ID'][good]
	ra = data['RA'][good]
	dec = data['DEC'][good]
	redshift = data['Z_MUSE'][good]
	if fitcat == 'mosaic':
		flya = data['LYALPHA_FLUX_SUM'][good]
	else:
		flya = data['LYALPHA_FLUX'][good]
	wavelength = data['LYALPHA_LBDA_OBS'][good]
	if offset: off = data['offset'][good]
	else: off = ids - ids

if hname is None: hname = filename

ttt, header_data_cube = getdata(hname, 0, header=True)  # .replace('.csub', '')

# Removing COMMENT key to avoid problems reading non-ascii characters
cards = header_data_cube.cards
bad = ['COMMENT' == b[0] for b in cards]
for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
hdulist = fits.open(filename)

if fitcat != 'EAGLE':
	w = wcs.WCS(header_data_cube, hdulist)
	xs, ys, zs = np.round(w.all_world2pix(ra, dec, [1]*len(ra), 1)).astype(int)
	zs = np.round(l2pix(wavelength)).astype(int)
	# xs, ys, zs = [data['x'][good], data['y'][good], data['z'][good]]
	if line == 'ovi':
		wav = (1 + redshift) * 1035
		zs = l2pix(wav)
	if line == 'halpha':
		wav = (1 + redshift) * 6562.8
		zs = l2pix(wav)

if fitcat == 'HDFS':
	xs = np.round(data['x'][good]).astype(int)
	ys = np.round(data['y'][good]).astype(int)
	zs = np.round(data['z'][good]).astype(int)

hdu = PrimaryHDU()

if random:
	ndata = len(xs)
	rnum = ndata*nrand
	_ids = ids
	_zs = zs
	_off = off
	rnums = np.zeros(ndata)
	xs = np.random.randint(cut, xlim-cut, size=rnum)
	ys = np.random.randint(cut, ylim-cut, size=rnum)
	for i in range(nrand-1):
		ids = np.concatenate((ids, _ids), 0)
		off = np.concatenate((off, _off), 0)
		#zs = np.concatenate((zs, _zs), 0)
		rnums = np.concatenate((rnums, np.zeros(ndata)+i+1), 0)
	zs = np.random.randint(zw, zlim-zw, size=rnum)


if fitcat == 'EAGLE':  cool = ids > 0
else: cool = (zs < (zlim - 5)) & (zs > 5)

if parallel:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	ncool = np.sum(cool)
	print "Parallel!, rank", rank, 'ncores', size
	r0, r1 = [ncool*rank/size, ncool*(rank+1)/size]
	print 'range', r0, r1
	cool = np.where(cool)[0][r0: r1]

if single:
	cool = ids == id0
#	_range = [0]
#if not single or random:
#	print 'aaa'
	#_range = range(len(xs[cool]))

extraname += '.rand'*random + ('.%s' % line) * (line != 'lya') + scsub + scorr
print extraname

for j in range(len(xs[cool])):

	x, y, z, i = [xs[cool][j], ys[cool][j], zs[cool][j]+off[cool][j], ids[cool][j]]

	# print 'Catalogue values', x, y, z, data_cube[
	#    int(z - 1), int(y - 1), int(x - 1)]  # , data_cube[z1 - d:z1 + d, y1 - d:y1 + d, x1 - d:x1 + d]

	if fitcat == 'EAGLE':
		name = "%s/%s/gals/snap%d_%s/%d%s.%s" % (folder, fitcat, snap, coord, i, extraname, type)
		if mask: mname = "%s/%s/gals/snap%d_%s/%d%s.mask" % (folder, fitcat, snap, coord, i, extraname)
		if gmask: gname = "%s/%s/gals/snap%d_%s/%d.gmask%d" % (folder, fitcat, snap, coord, i, galrad)
	else:
		name = "%s/%s/LAEs/%d%s" % (folder, fitcat, i, extraname)
		if mask: mname = "%s/%s/LAEs/%d%s.mask" % (folder, fitcat, i, extraname)
		if gmask: gname = "%s/%s/LAEs/%d.gmask%d" % (folder, fitcat, i, galrad)

	if random and (nrand > 1):
		nr = rnums[cool][j]
		name += '%d' % nr
		if mask: mname += '.%d' % nr
		if gmask: gname += '.%d' % nr

	isf = os.path.isfile(name+ext)
	if overwrite or not isf:
		print "--------------------"
		print 'Overwriting' * overwrite * isf, 'New' * (not isf), fitcat, i
		if makefit: flux = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if mask: flux2 = np.zeros([yw + 1, xw + 1]) + float('NaN')
		if gmask: flux3 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if fitcat == 'EAGLE':
			if mask:
				roll = np.roll(mcube, (ylim/2-y, xlim/2-x), (0, 1))
				flux2 = roll[ylim/2+ymin: ylim/2+ymax+1, xlim/2+xmin: xlim/2+xmax+1]
			if gmask:
				roll = np.roll(gcube, (zlim/2-z, ylim/2-y, xlim/2-x), (0, 1, 2))
				_cube = np.copy(roll[zlim/2+zmin: zlim/2+zmax+1, ylim/2+ymin: ylim/2+ymax+1, xlim/2+xmin: xlim/2+xmax+1])
				gal = _cube == i
				print 'gmask gal', np.sum(gal)
				_cube[gal] = 0
				flux3 = _cube
			if makefit:
				roll = np.roll(data_cube, (zlim/2-z, ylim/2-y, xlim/2-x), (0, 1, 2))
				flux = roll[zlim/2+zmin: zlim/2+zmax+1, ylim/2+ymin: ylim/2+ymax+1, xlim/2+xmin: xlim/2+xmax+1]
		else:
			_xmin = max(cut, x+xmin)
			dxmin = _xmin-x-xmin
			_xmax = min(xlim-cut-1, x+xmax)
			dxmax = _xmax-x-xmax

			_ymin = max(cut, y+ymin)
			dymin = _ymin-y-ymin
			_ymax = min(ylim-cut-1, y+ymax)
			dymax = _ymax-y-ymax

			if mask:
				_cube = mcube[_ymin: _ymax+1, _xmin: _xmax+1]
				flux2[dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube

			if makefit or gmask:
				_zmin = max(0, z+zmin)
				dzmin = _zmin-z-zmin
				_zmax = min(zlim-1, z+zmax)
				dzmax = _zmax-z-zmax
			if makefit:
				_cube = data_cube[_zmin: _zmax+1, _ymin: _ymax+1, _xmin: _xmax+1]
				flux[dzmin: zw+dzmax+1, dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube
			if gmask:
				_cube = np.copy(gcube[_zmin: _zmax + 1, _ymin: _ymax + 1, _xmin: _xmax + 1])
				gal = _cube == i
				print 'gmask gal', np.sum(gal)
				_cube[gal] = 0
				flux3[dzmin: zw+dzmax+1, dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube

		if makefit:
			hdu.data = flux
			hdu.writeto(name+ext, clobber=True)
		if makeim:
			hdu.data = np.nansum(flux, 0)
			hdu.writeto(name+'.IM'+ext, clobber=True)
		if mask:
			hdu.data = flux2
			hdu.writeto(mname+ext, clobber=True)
		if gmask:
			hdu.data = flux3
			hdu.writeto(gname+ext, clobber=True)

