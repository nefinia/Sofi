#!/usr/bin/env python
import os
from sys import argv

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
import scipy.signal as sp
from astropy import wcs
from pyfits import getdata, PrimaryHDU

import params
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
mask2d = rdarg(argv, 'mask2d', bool, False)
gmask = rdarg(argv, 'gmask', bool, False)
model = rdarg(argv, 'model', str, 'HM12')
offset = rdarg(argv, 'offset', bool, False)
centerfix = rdarg(argv, 'centerfix', bool, False)
cubex = rdarg(argv, 'cubex', bool, False)
corr = rdarg(argv, 'corr', bool, False)
scorr = '.corr' * corr
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
folder = rdarg(argv, key='folder', type=str, default='/net/eos/scratch/gallegos/')#'/net/galaxy-data/export/galaxydata/gallegos/')#''
foldercat = rdarg(argv, key='foldercat', type=str, default='../../')
folderout = rdarg(argv, key='folderout', type=str, default='../../')
extraname = rdarg(argv, key='extraname', type=str, default='')
csub = rdarg(argv, 'csub', bool, True)
scsub = '.csub' * csub
single = rdarg(argv, 'single', bool, False)
id0 = rdarg(argv, 'id', int, 1)
galrad = rdarg(argv, 'rad', int, 3)
zw0 = rdarg(argv, 'zw0', int, 5)
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

if verbose < 2:
	vb = ' > /dev/null 2>&1'
else:
	vb = ''

dovar = True
xw = rdarg(argv, 'xw', int, 200)  # for EAGLE it was last time 500 (15.10.18)
yw = rdarg(argv, 'yw', int, 200)  # for EAGLE it was last time 500 (15.10.18)
zw = rdarg(argv, 'zw', int, 200)
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

if fitcat == 'EAGLE':
	csub = False
	#layers = {'8': 45, '9': 39, '10': 34, '11': 29, '12': 25, '13': 20, '14': 18, '15': 16}
	#_sgamma = '%.2E_%.2E_%.2E' % tuple(np.array(params.gamma_bkg[model][str(snap)]))

	#cubename = 'snap%s_%s_%s_512_4096_8_%d_%s_%.2E%s.%s.fits' % (snap, coord, model, params.zlens[str(snap)],
	#														_sgamma, params.nhSS[model][str(snap)], scase, type)
	cubename = 'snap%s_%s_%s_512_4096_12_%d_6.73E-03_CaseA.%s.fits' % (snap, coord, model, params.zlens[snap], type)

	if mask: maskname = 'snap%s_%s.mask.fits' % (snap, coord)
	if gmask:
		gmaskname = 'snap%s_%s.galmask_%darcsec_zw%d.fits' % (snap, coord, galrad, zw0)
	nconf = 1

if fitcat == 'HDFS':
	#cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s.fits' % '.csub' * csub
	cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s%s.fits' % (scsub, scorr)
	if mask2d: mask2dname = 'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id.fits'
	if gmask: gmaskname = 'HDFS.galmask_%darcsec_zw%d.fits' % (galrad, zw)
	if dovar:
		varname = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
	if mask: maskname = 'DATACUBE-HDFS-1.35-PROPVAR.Objects_Id.fits'
	nconf = 1

if fitcat == 'UDF':
	#cubename = 'UDF.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'UDF10.z1300%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-10%s.fits' % '.csub' * csub
	if gmask: gmaskname = 'UDF.galmask_%darcsec_zw%d.fits' % (galrad, zw0)
	cubename = 'DATACUBE_UDF-10%s%s.fits' % (scsub, scorr)
	print cubename
	hname = '%s/%s/%s' % ('/net/galaxy-data/export/galaxydata/gallegos/', fitcat, 'DATACUBE_UDF-10.fits')
	if mask: maskname = 'DATACUBE_UDF-10.Objects_Id.fits'
	if mask2d: mask2dname = 'DATACUBE_UDF-10.IM.Objects_Id.fits'
	nconf = 1

if fitcat == 'MXDF':
	#cubename = 'UDF.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'UDF10.z1300%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-10%s.fits' % '.csub' * csub
	if gmask: gmaskname = 'MXDF.galmask_%darcsec_zw%d.fits' % (galrad, zw0)
	cubename = 'DATACUBE_MXDF_ZAP_COR%s%s.fits' % (scsub, scorr)
	hname = '%s/%s/DATACUBE_MXDF_ZAP_COR.csub.fits' % (folder, fitcat)
	print cubename
	if mask: maskname = 'DATACUBE_MXDF_ZAP_COR.mask.fits'
	if mask2d: mask2dname = 'DATACUBE_MXDF_ZAP_COR.IM.Objects_Id.fits'
	nconf = 1

if fitcat == 'mosaic' or fitcat == 'mosaic10':
	offset = False
	#cubename = '/net/galaxy-data/export/galaxydata/gallegos/mosaic/DATACUBE_UDF-MOSAIC.z1300%s.fits' % '.csub' * csub
	#cubename = 'mosaic.bkgcorr%s.fits' % '.csub' * csub
	#cubename = 'DATACUBE_UDF-MOSAIC.z1300%s.corr.fits' % scsub
	cubename = 'DATACUBE_UDF-MOSAIC%s%s.fits' % (scsub, scorr)
	hname = '%s/%s/%s' % ('/net/galaxy-data/export/galaxydata/gallegos/', fitcat, 'DATACUBE_UDF-MOSAIC.fits')
	if mask: maskname = 'DATACUBE_UDF-MOSAIC.Objects_Id.fits'
	if mask2d: mask2dname = 'DATACUBE_UDF-MOSAIC.Objects_Id.fits'
	if gmask: gmaskname = 'mosaic.galmask_%darcsec_zw%d.fits' % (galrad, zw0)
	nconf = 1
	vmin = -3
	vmax = 10
	if std is None: std = 100

filename = '%s/%s/%s' % (folder, fitcat, cubename)

data_cube = getdata(filename, 0)
zlim, ylim, xlim = data_cube.shape
if mask:
	mcube = getdata('%s/%s/%s' % (folder, fitcat, maskname), 0)
if mask2d:
	m2dcube = getdata('%s/%s/%s' % (folder, fitcat, mask2dname), 0)[0]
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
if fitcat == 'UDF' or fitcat == 'mosaic' or fitcat == 'MXDF':
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

if fitcat == 'UDF' or fitcat == 'mosaic' or fitcat == 'MXDF':
	ids = data['ID'][good]
	ra = data['RA'][good]
	dec = data['DEC'][good]
	redshift = data['Z_MUSE'][good]
	
	if offset: off = data['offset'][good]
	else: off = ids - ids


if fitcat != 'EAGLE':

	savecoord = True
	#newcoord = True
	#if not newcoord:
	try:
		print 'Reading x,y,x coordinates!'
		xs = np.round(data['x'][good]).astype(int)
		ys = np.round(data['y'][good]).astype(int)
		zs = np.round(data['z'][good]).astype(int)
		# & (xs<xlim) & (ys>0) & (ys<ylim) (zs>0) & (zs<zlim)
	
	#else:
	except:
		print 'No x,y,x coordinates!'
		if hname is None: hname = filename
		ttt, header_data_cube = getdata(hname, 0, header=True)  # .replace('.csub', '')
		# Removing COMMENT key to avoid problems reading non-ascii characters
		cards = header_data_cube.cards
		bad = ['COMMENT' == b[0] for b in cards]
		for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
		hdulist = fits.open(filename)
		wavelength = 1215.67 * (1 + redshift)
		xs, ys, zs = np.zeros(len(wavelength)), np.zeros(len(wavelength)), np.zeros(len(wavelength))
		w = wcs.WCS(header_data_cube, hdulist)
		_xs, _ys, _zs = np.round(w.all_world2pix(ra, dec, [1] * len(ra), 1)).astype(int)
		_zs = np.round(l2pix(wavelength)).astype(int)
		np.savetxt('coords%s.dat' % fitcat, np.mat([ids, _xs, _ys, _zs]).T, header='id x y z',
		                         fmt='%d %d %d %d')
		cool = (xs > 0)
		xs[~cool] = _xs[~cool]
		ys[~cool] = _ys[~cool]
		zs[~cool] = np.round(_zs[~cool]).astype(int)
	
	xoff, yoff, zoff = 0, 0, 0
	
	if fitcat == 'MXDF': xoff, yoff, zoff = 0, 0, 38#30, -43, 83
	if fitcat == 'UDF' or fitcat == 'HDFS': xoff, yoff, zoff = 0, 0, 0
	if fitcat == 'mosaic': xoff, yoff, zoff = 0, 0, -1

	
	xs += xoff
	ys += yoff
	zs += zoff
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
else: cool = (xs>0) & (xs<xlim) & (ys>0) & (ys<ylim) & (zs < (zlim - 5)) & (zs > 5)

print 'Total LAEs', len(xs), 'cool', np.sum(cool)

#for iii in ids[~cool]: os.system('rm %s/%s/LAEs/%d*.fits' % (folder, fitcat, iii))

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

extraname += '.rand'*random + ('.%s' % line) * (line != 'lya') + scsub + scorr
print extraname
fstack = []

sblim=False
if sblim:
	if mask:
		c=np.copy(data_cube)
		c[mcube>0] = np.nan
		std = []
		for j in range(36):
			cool = c[j:(j + 1) * 100, :, :]
			for i in range(1000): std.append(np.nansum(np.random.choice(cool, 25)))
			print 'z', j * 100 - 50, '1sigma SB', np.nanstd(std)


for j in range(len(xs[cool])):
	x, y, z, i = np.array([xs[cool][j], ys[cool][j], zs[cool][j]+off[cool][j], ids[cool][j]]).astype(int)
	print "--------------------"
	print i, x, y, z
	# print 'Catalogue values', x, y, z, data_cube[
	#    int(z - 1), int(y - 1), int(x - 1)]  # , data_cube[z1 - d:z1 + d, y1 - d:y1 + d, x1 - d:x1 + d]

	if fitcat == 'EAGLE':
		name = "%s/%s/gals/snap%d_%s/%d%s.%s" % (folderout, fitcat, snap, coord, i, extraname, type)
		if mask: mname = "%s/%s/gals/snap%d_%s/%d%s.mask" % (folderout, fitcat, snap, coord, i, extraname)
		if gmask: gname = "%s/%s/gals/snap%d_%s/%d.gmask%d" % (folderout, fitcat, snap, coord, i, galrad)
		if spec: sname = "%s/%s/gals/snap%d_%s/spec/%d%s.%s" % (folderout, fitcat, snap, coord, i, extraname, type)
	else:
		name = "%s/%s/LAEs/%d%s" % (folderout, fitcat, i, extraname)
		if mask: mname = "%s/%s/LAEs/%d.mask" % (folderout, fitcat, i)
		if mask2d: m2dname = "%s/%s/LAEs/%d.mask2d" % (folderout, fitcat, i)
		if gmask: gname = "%s/%s/LAEs/%d.gmask%d" % (folderout, fitcat, i, galrad)
		if spec: sname = "%s/%s/LAEs/spec/%d%s" % (folderout, fitcat, i, extraname)

	if random and (nrand > 1):
		nr = rnums[cool][j]
		name += '%d' % nr
		if mask: mname += '.%d' % nr
		if gmask: gname += '.%d' % nr
		if spec: sname += '.%d' % nr

	isf = os.path.isfile(name + ext)
	

	if overwrite or not isf:
		print 'Overwriting' * overwrite * isf, 'New' * (not isf), fitcat, i
		if makefit: flux = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if mask: flux2 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if gmask: flux3 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if mask2d: flux2d = np.zeros([yw + 1, xw + 1]) + float('NaN')
		if fitcat == 'EAGLE':
			if mask:
				roll = np.roll(mcube, (zlim/2-z, ylim/2-y, xlim/2-x), (0, 1))
				flux2 = roll[zlim/2+zmin: zlim/2+zmax+1, ylim/2+ymin: ylim/2+ymax+1, xlim/2+xmin: xlim/2+xmax+1]
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
			if mask2d:
				_cube = np.copy(m2dcube[_ymin: _ymax+1, _xmin: _xmax+1])
				flux2d[dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube
				gid = flux2d[yw/2, xw/2]
				flux2d[flux2d == gid] = 0
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
		if mask2d:
			hdu.data = flux2d
			hdu.writeto(m2dname+ext, clobber=True)
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

	else: flux = getdata(name + ext)
	
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
		plt.savefig(name.replace('/LAEs/', '/LAEs/images/')+'.png')
		plt.close()

	if spec and isf:#(not isf or overwrite):
		zl0 = 10
		w0 = 5
		zl, yl, xl = flux.shape
		_fs = np.nanmean(flux[zl/2-zl0/2:zl/2+zl0/2+1, yl/2-w0/2:yl/2+w0/2+1, xl/2-w0/2:xl/2+w0/2+1], (1, 2))
		v = np.arange(zl0+1)-zl0/2
		np.savetxt(sname + '.1asec2.dat', _fs)
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
		v = np.arange(zl0+1)-zl0/2
		np.savetxt(sname + '.4asec2.dat', _fs)
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
	hdu.writeto("%s/%s/LAEs/stack%s.fits" % (folderout, fitcat, extraname), clobber=True)

