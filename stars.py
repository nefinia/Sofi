#!/usr/bin/env python
import astropy.io.fits as fits
import os, params
import matplotlib.pyplot as plt
from astropy import wcs
from pyfits import getdata, PrimaryHDU
from sys import argv
import scipy.signal as sp
import scipy.ndimage as ndimage

import numpy as np

from tools_sofi import rdarg


def l2pix(l):
	"""Convert wavelength to pixels"""
	l0 = 4750.
	pixsize = 1.25  # 1 pixel = 1.25 Angstrom
	return (l - l0) / pixsize

def pix2l(p,l0=4750.,pixsize=.125):
	l0+p*pixsize

def l2z(l):
	lya_rest = 1215.67
	return l/lya_rest-1

def lya_ob(z):
	"""Convert lya redshift to observed wavelength"""
	lya_rest = 1215.67
	return lya_rest * (1 + z)

caseA = rdarg(argv, 'caseA', bool, True)
scase = '_CaseA' * caseA
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
extraname = rdarg(argv, key='extranamecat', type=str, default='')
csub = rdarg(argv, 'csub', bool, True)
single = rdarg(argv, 'single', bool, False)
id0 = rdarg(argv, 'id', int, 1)
galrad = rdarg(argv, 'rad', int, 6)
smooth = rdarg(argv, 'smooth', int, 0)

if smooth:
	ssmooth = '.smooth%d' % smooth
	import scipy.ndimage as ndimage
spec = rdarg(argv, 'spec', bool, False)
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
scsub = ''#''.csub' * csub
if verbose < 2:
	vb = ' > /dev/null 2>&1'
else:
	vb = ''

dovar = True
xw = rdarg(argv, 'xw', int, 200)  # for EAGLE it was last time 500 (15.10.18)
yw = rdarg(argv, 'yw', int, 200)  # for EAGLE it was last time 500 (15.10.18)
zw = rdarg(argv, 'zw', int, 100)
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

if fitcat == 'HDFS':
	#cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s.fits' % '.csub' * csub
	cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s%s.fits' % (scsub, scorr)
	if mask: maskname = 'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id.fits'
	if gmask: gmaskname = 'HDFS.galmask_%darcsec.fits' % galrad
	if dovar:
		varname = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
	if mask: maskname = 'DATACUBE-HDFS-1.35-PROPVAR.Objects_Id.fits'
	nconf = 1

if fitcat == 'UDF':
	#cubename = 'UDF.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'UDF10.z1300%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-10%s.fits' % '.csub' * csub
	if gmask: gmaskname = 'UDF.galmask_%darcsec.fits' % galrad
	cubename = 'DATACUBE_UDF-10%s%s.fits' % (scsub, scorr)
	print cubename
	hname = '%s/%s/%s' % (folder, fitcat, 'DATACUBE_UDF-10.fits')
	if mask: maskname = 'DATACUBE_UDF-10.Objects_Id.fits'
	nconf = 1

if fitcat == 'MXDF':
	#cubename = 'UDF.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'UDF10.z1300%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-10%s.fits' % '.csub' * csub
	if gmask: gmaskname = 'MXDF.galmask_%darcsec.fits' % galrad
	cubename = 'MXDF.bkgcorr.csub.fits'# % (scsub, scorr)
	hname = '%s/%s/MXDF_DATACUBE_FINAL_PROPVAR.csub.fits' % (folder, fitcat)
	print cubename
	if mask: maskname = 'MXDF_DATACUBE_FINAL_PROPVAR.Objects_Id.fits'
	nconf = 1

if fitcat == 'mosaic' or fitcat == 'mosaic10':
	offset = False
	#cubename = '/net/galaxy-data/export/galaxydata/gallegos/mosaic/DATACUBE_UDF-MOSAIC.z1300%s.fits' % '.csub' * csub
	#cubename = 'mosaic.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-MOSAIC.z1300%s.corr.fits' % scsub
	cubename = 'DATACUBE_UDF-MOSAIC%s%s.fits' % (scsub, scorr)
	hname = '%s/%s/%s' % (folder, fitcat, 'DATACUBE_UDF-MOSAIC.fits')
	if mask: maskname = 'DATACUBE_UDF-MOSAIC.Objects_Id.fits'
	if gmask: gmaskname = 'mosaic.galmask_%darcsec.fits' % galrad
	nconf = 1
	vmin = -3
	vmax = 10
	if std is None: std = 100

filename = '%s/%s/%s' % (folder, fitcat, cubename)

data_cube = getdata(filename, 0)
zlim, ylim, xlim = data_cube.shape
if mask:
	mcube = getdata('%s/%s/%s' % (folder, fitcat, maskname), 0)
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

cat = '%s/%s/cats/stars.fits' % (foldercat, fitcat)
data = getdata(cat, 1)

if fitcat == 'HDFS':
	conf = data['confidence']
if fitcat == 'UDF' or fitcat == 'mosaic' or fitcat == 'MXDF':
	conf = data['CONFID']

if fitcat == 'HDFS':
	m300 = data['m300']
	redshift = data['z']
	wavelength = lya_ob(redshift)
	flya = data['LYALPHA_FLUX']
	ids = data['id']
	ra = data['ra']
	dec = data['dec']
	if offset: off = data['offset']
	else: off = ids - ids

if fitcat == 'UDF' or fitcat == 'mosaic' or fitcat == 'MXDF':
	ids = data['ID']
	ra = data['RA']
	dec = data['DEC']
	redshift = data['Z_MUSE']
	if fitcat == 'mosaic' or fitcat == 'MXDF':
		flya = data['LYALPHA_FLUX_SUM']
	else:
		flya = data['LYALPHA_FLUX']
	wavelength = data['LYALPHA_LBDA_OBS']
	if offset: off = data['offset']
	else: off = ids - ids

if hname is None: hname = filename

ttt, header_data_cube = getdata(hname, 0, header=True)  # .replace('.csub', '')

# Removing COMMENT key to avoid problems reading non-ascii characters
cards = header_data_cube.cards
bad = ['COMMENT' == b[0] for b in cards]
for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
hdulist = fits.open(filename)

w = wcs.WCS(header_data_cube, hdulist)
xs, ys, zs = np.round(w.all_world2pix(ra, dec, [1]*len(ra), 1)).astype(int)

hdu = PrimaryHDU()

if random:
	ndata = len(xs)
	rnum = ndata*nrand
	_ids = ids
	_off = off
	_xs = xs
	_ys = ys
	rnums = np.zeros(ndata)
	#np.random.randint(cut, xlim-cut, size=rnum)
	#np.random.randint(cut, ylim-cut, size=rnum)
	for i in range(nrand-1):
		xs = np.concatenate((xs, _xs), 0)
		ys = np.concatenate((ys, _ys), 0)
		ids = np.concatenate((ids, _ids), 0)
		off = np.concatenate((off, _off), 0)
		rnums = np.concatenate((rnums, np.zeros(ndata)+i+1), 0)
	zs = np.random.randint(zw, zlim-zw, size=rnum)


cool = (xs>0) & (xs<xlim) & (ys>0) & (ys<ylim)

#for iii in ids[~cool]: os.system('rm %s/%s/stars/%d*.fits' % (folder, fitcat, iii))

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
fstack = []

for j in range(len(xs[cool])):

	x, y, z, i = [xs[cool][j], ys[cool][j], zs[cool][j]+off[cool][j], ids[cool][j]]

	# print 'Catalogue values', x, y, z, data_cube[
	#    int(z - 1), int(y - 1), int(x - 1)]  # , data_cube[z1 - d:z1 + d, y1 - d:y1 + d, x1 - d:x1 + d]
	name = "%s/%s/stars/%d%s" % (folder, fitcat, i, extraname)
	if mask: mname = "%s/%s/stars/%d.mask" % (folder, fitcat, i)
	if gmask: gname = "%s/%s/stars/%d.gmask%d" % (folder, fitcat, i, galrad)
	if spec: sname = "%s/%s/stars/spec/%d%s" % (folder, fitcat, i, extraname)

	if random and (nrand > 1):
		nr = rnums[cool][j]
		name += '%d' % nr
		if mask: mname += '.%d' % nr
		if gmask: gname += '.%d' % nr
		if spec: sname += '.%d' % nr

	isf = os.path.isfile(name+ext)
	if isf: flux = getdata(name + ext)

	if overwrite or not isf:
		print "--------------------"
		print 'Overwriting' * overwrite * isf, 'New' * (not isf), fitcat, i
		if makefit: flux = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if mask: flux2 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if gmask: flux3 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')

		_xmin = max(cut, x+xmin)
		dxmin = _xmin-x-xmin
		_xmax = min(xlim-cut-1, x+xmax)
		dxmax = _xmax-x-xmax

		_ymin = max(cut, y+ymin)
		dymin = _ymin-y-ymin
		_ymax = min(ylim-cut-1, y+ymax)
		dymax = _ymax-y-ymax

		_zmin = max(0, z + zmin)
		dzmin = _zmin - z - zmin
		_zmax = min(zlim - 1, z + zmax)
		dzmax = _zmax - z - zmax

		if makefit:
			_cube = data_cube[_zmin: _zmax+1, _ymin: _ymax+1, _xmin: _xmax+1]
			flux[dzmin: zw+dzmax+1, dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube

		if mask:
			_cube = np.copy(mcube[_zmin: _zmax+1, _ymin: _ymax+1, _xmin: _xmax+1])
			flux2[dzmin: zw+dzmax+1, dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube
			gid = flux2[zw/2, yw/2, xw/2]
			flux2[flux2 == gid] = 0

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
				hdu.data = np.nansum(flux[zw/2-3: zw/2+4], 0)
				hdu.writeto(name+'.IM'+ext, clobber=True)
		if mask:
			hdu.data = flux2
			hdu.writeto(mname+ext, clobber=True)
		if gmask:
			hdu.data = flux3
			hdu.writeto(gname+ext, clobber=True)
		if smooth:
			flux[flux2 > 0] = 0
			#flux[flux3 > 0] = 0
			flux[np.isnan(flux)] = 0
			flux4 = ndimage.gaussian_filter(flux, sigma=(0, smooth, smooth), order=0)
			hdu.data = flux4
			hdu.writeto(name+ssmooth+ext, clobber=True)
			if makeim:
				hdu.data = np.sum(flux4[zw/2-3: zw/2+4], 0)
				hdu.writeto(name+ssmooth+'.IM' + ext, clobber=True)

	if isf and makeim:
		zl0 = 3
		w0 = 30
		zl, yl, xl = flux.shape
		_fs = np.nanmean(flux[zl/2-zl0/2:zl/2+zl0/2+1, yl/2-w0/2:yl/2+w0/2+1, xl/2-w0/2:xl/2+w0/2+1], 0)
		_fs = ndimage.gaussian_filter(_fs, sigma=1, order=0)
		plt.figure()
		plt.title('%d%s'%(i,extraname))
		plt.imshow(_fs)
		plt.plot([w0/2], [w0/2], marker='x', color='red', markersize=10)
		plt.savefig(name.replace('/stars/', '/stars/images/')+'.png')
		plt.close()

	if spec: isf2 = os.path.isfile(sname + '.dat')

	if spec and isf:#(not isf or overwrite):
		zl0 = 10
		w0 = 5
		zl, yl, xl = flux.shape
		_fs = np.nanmean(flux[zl/2-zl0/2:zl/2+zl0/2+1, yl/2-w0/2:yl/2+w0/2+1, xl/2-w0/2:xl/2+w0/2+1], (1, 2))
		v = np.arange(zl0)-zl0/2
		np.savetxt(sname + '.dat', _fs)
		plt.figure(figsize=(12, 6))
		plt.xticks(v, v)
		plt.axvline(x=0, color='red', linestyle='--')
		plt.grid()
		plt.plot(v, sp.savgol_filter(_fs, 3, 1))
		plt.xlabel('offset')
		plt.ylabel('flux [1e-20 erg/s/cm^2/Angstrom]')
		plt.savefig(sname + '.1asec2.png')
		plt.close()
		w0 = 10
		_fs = np.nanmean(flux[zl/2-zl0/2:zl/2+zl0/2+1, yl/2-w0/2:yl/2+w0/2+1, xl/2-w0/2:xl/2+w0/2+1], (1, 2))
		v = np.arange(zl0)-zl0/2
		np.savetxt(sname + '.dat', _fs)
		plt.figure(figsize=(12, 6))
		plt.xticks(v, v)
		plt.axvline(x=0, color='red', linestyle='--')
		plt.grid()
		plt.plot(v, sp.savgol_filter(_fs, 3, 1))
		plt.xlabel('offset')
		plt.ylabel('flux [1e-20 erg/s/cm^2/Angstrom]')
		plt.savefig(sname + '.4asec2.png')
		plt.close()

	if simplestack:
		fstack.append(flux)


if simplestack:
	fstack = np.nanmedian(fstack, 0)
	hdu = PrimaryHDU()
	hdu.data = fstack
	hdu.writeto("%s/%s/stars/stack%s.fits" % (folder, fitcat, extraname), clobber=True)

