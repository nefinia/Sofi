#!/usr/bin/env python
__author__ = 'nefinia'

import h5py
import os
import time
from pyfits import getdata, PrimaryHDU
from sys import argv

import numpy as np


def rdarg(argv, key, type=None, default=None, listtype=int):
    # print argv, key, type, default

    if len(argv) > 1:
        opt = np.where([a == '-%s' % key for a in argv])[0] + 1
        if len(opt) > 0:
            name = argv[int(opt)]
            if type is list:
                name = name.split(',')
                if listtype==int: name = [int(i) for i in name]
                elif listtype==float: name = [float(i) for i in name]
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

def b2s(s):
	if s:
		return '.true.'
	else:
		return '.false.'

coord = rdarg(argv, 'coord', str, 'x')
csub = rdarg(argv, 'csub', bool, True)
cubetype = rdarg(argv, 'cubetype', str, 'SB')
daysmin = rdarg(argv, 'days', float, 10000)  # days = 1e-10#7
dmin = rdarg(argv, 'dmin', float, .5)
dmax = rdarg(argv, 'dmax', float, 20)
dims = rdarg(argv, 'dims', int, 3)
docat = rdarg(argv, 'docat', bool, True)
dorandom = rdarg(argv, 'random', bool, False)
dovar = rdarg(argv, 'dovar', bool, False)
extraname = rdarg(argv, 'extraname', str, '')
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
folder = rdarg(argv, 'folder', str, '/net/galaxy-data/export/galaxydata/gallegos/')  # '/scratch/gallegos/MUSE/'
foldercat = rdarg(argv, 'foldercat', str, '../../')  # '/scratch/gallegos/MUSE/'
folderout = rdarg(argv, 'folderout', str, '/net/galaxy-data/export/galaxydata/gallegos/')
galpos = rdarg(argv, 'galpos', bool, False)
i1 = rdarg(argv, 'id1', int, -1)
i2 = rdarg(argv, 'id2', int, -1)
idmin = rdarg(argv, 'idmin', int, 0)
idmax = rdarg(argv, 'idmax', int, 9e9)
imov = rdarg(argv, 'imov', bool, False)
laecentered = rdarg(argv, 'laecentered', bool, True)
hdf5 = rdarg(argv, 'hdf5', bool, True)
mask = rdarg(argv, 'mask', bool, False)
norm = rdarg(argv, 'norm', bool, False)
docat = rdarg(argv, 'docat', bool, False)
docalc = rdarg(argv, 'docalc', bool, True)
noshear = rdarg(argv, 'noshear', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
outtype = rdarg(argv, 'outtype', str, 'h5')
parallel = rdarg(argv, 'parallel', bool, True)
periodic = rdarg(argv, 'periodic', str, '0,0,0')
nmin = rdarg(argv, 'nmin', int, 0)
nmax = rdarg(argv, 'nmax', int, 9e9)
randdir = rdarg(argv, 'rdir', bool, False)
redmin = rdarg(argv, 'rmin', list, [0])
redmax = rdarg(argv, 'rmax', list, [99])
simon = rdarg(argv, 'cubical', bool, True)
single = rdarg(argv, 'single', bool, False)
sizemin = rdarg(argv, 'size', int, 1)
snap = rdarg(argv, 'snap', int, 12)
thetamin = rdarg(argv, 'thetamin', float, 6)
r = rdarg(argv, 'rad', int, 50)
xw = rdarg(argv, 'xw', int, 80)
xmore = rdarg(argv, 'xmore', int, 0)#100 for eagle
yw = rdarg(argv, 'yw', int, 80)
zw = rdarg(argv, 'zw', int, 6)
zw0 = rdarg(argv, 'zw0', int, 2)
xmin = rdarg(argv, 'xmin', int, 0)
ymin = rdarg(argv, 'ymin', int, 0)
xmax = rdarg(argv, 'xmax', int, 9e9)
ymax = rdarg(argv, 'ymax', int, 9e9)

if not extraname: extraname = ''
if noshear: extraname += '.noshear'
if not laecentered: extraname += '.pair'

smask = '.mask' * mask
extra0 = ('snap%d_%s' % (snap, coord)) * (fitcat == 'EAGLE')
if galpos: extraname = '.gals'+extraname

now = time.time()
byte2Mb = 9.5367431640625e-07

def doim(outname, px, py, pz, data, xwn, ywn, zwn, xmore=0, binsize=2, sb=True):
	f = np.zeros([2 * zwn + 1, 2 * ywn / binsize + 1, (2 * xwn + 1 + xmore) / binsize + 1])
	_x = ((px + xwn + .5) / binsize).astype(int)
	_y = ((py + ywn + .5) / binsize).astype(int)
	_z = (pz + zwn + .5).astype(int)
	for x, y, z, d in zip(_x, _y, _z, data): f[z, y, x] += d
	# Correction for SB
	if sb: f /= binsize ** 2
	hdu = PrimaryHDU()
	hdu.data = f
	hdu.writeto(outname, clobber=True)
	hdu.data = np.nansum(f[zw - zw0:zw + zw0 + 1, :, :], 0)
	hdu.writeto(outname.replace('.fits', '.IM.fits'), clobber=True)


def coordconv(id1, id2, x2, y2, z2, xwn=None, ywn=None, zwn=None, xmore=0, norm=False):
	def conv(u, d, xl, yl, zl):
		invd = 1. / np.linalg.norm(d[1:])
		px = np.zeros((zl, yl, xl)) + np.dot(u[1:], d[1:]) * invd
		py = np.zeros((zl, yl, xl)) + (u[1] * d[2] - u[2] * d[1]) * invd
		pz = u[0] - d[0] * px * invd
		return px, py, pz

	for i1 in id1:

		name = '%s/%s/gals/%s/%d%s%s.fits' % (folder, fitcat, extra0, i1, extraname, smask)
		if not os.path.isfile(name):
			print 'File %s does not exist' % name
		else:
			cube = getdata(name, 0)
			zl, yl, xl = cube.shape
			xw, yw, zw = [xl / 2, yl / 2, zl / 2]
			if xwn is None: xwn = xw
			if ywn is None: ywn = yw
			if zwn is None: zwn = zw
			C = np.ogrid[0:zl, 0:yl, 0:xl]
			Cl = np.array([zw, yw, xw])
			for x, y, z, i2 in zip(x2, y2, z2, id2):
				print 'Pair %d-%d' % (i1, i2)
				extra1 = ('SB_snap%d_%s' % (snap, coord)) * (fitcat == 'EAGLE')
				fitsname = '%s/%s/pairs/%s/%d-%d%s%s.fits' % (folder, fitcat, extra1, i1, i2, extraname, smask)
				if not os.path.isfile(fitsname) or imov:
					if docalc:
						Cn = np.array([z+zw, y+yw, x+xw])
						d = Cn - Cl
						u = C - Cl
						px, py, pz = conv(u, d, xl, yl, zl)
						cool = np.where((px <= xwn + xmore) & (px >= -xwn) & (abs(py) <= ywn) & (abs(pz) <= zwn))
	
						output = '%s/%s/cats/%s/%d-%d.%s' % (folder, fitcat, extra1, i1, i2, outtype)
						px = px[cool]
						py = py[cool]
						pz = pz[cool]
						sb = cube[cool]
						isfile = os.path.isfile(output)
						if (not isfile or overwrite) and docat:
							if not isfile: print 'Output %s file does not exist.' % outtype
							if overwrite: 'Overwriting output %s file.' % outtype
							print output
							if outtype == 'txt':
								out = open(output, 'w')
								for xx, yy, zz, ss in zip(px, py, pz, sb):
									out.write('%f %f %f %f\n' % (xx, yy, zz, ss))
								out.close()
							if outtype == 'h5':
								with h5py.File(output, 'w') as h:
									h.create_dataset('px', data=px)
									h.create_dataset('py', data=py)
									h.create_dataset('pz', data=pz)
									h.create_dataset('SB', data=sb)
						else:
							if isfile: print 'Output file does exist and we are not overwriting it.'
							if not docat: print 'No catalog saved.'
					else:
						data = h5py.File(output, 'r')
						px = data['px']
						py = data['py']
						pz = data['pz']
						sb = data['SB']
	
					fitsname = '%s/%s/pairs/%s/%d-%d%s%s.fits' % (folder, fitcat, extra1, i1, i2, extraname, smask)
					if not os.path.isfile(fitsname) or imov:
						print 'Making fits!', 2*xwn+1+xmore, 2*ywn+1, 2*zwn+1
						doim(fitsname, px, py, pz, sb, xwn, ywn, zwn, xmore)
				else:
					print 'Image aready exists.'


if fitcat == 'HDFS':
	cat = '%s/%s/cats/laes.fits' % (foldercat, fitcat)
	data = getdata(cat, 1)
	ids = data['ID']
	zs = data['redshift']
	ra = data['RaHMS']
	dec = data['DecDMS']
	sconf = data['Sofi_confidence']

if fitcat == 'UDF':
	cat = '%s/%s/cats/laes.fits' % (foldercat, fitcat)
	data = getdata(cat, 1)
	ids = data['ID']
	zs = data['Z_MUSE']
	sconf = data['Sofi_confidence']

if fitcat == 'mosaic':
	cat = '%s/%s/cats/laes.fits' % (foldercat, fitcat)
	data = getdata(cat, 1)
	ids = data['ID']
	zs = data['Z_MUSE']
	sconf = np.zeros(len(ids)) + 2  # data['Sofi_confidence']

if fitcat == 'mosaic10':
	cat = '%s/%s/cats/laes.fits' % (foldercat, 'UDF')
	data = getdata(cat, 1)
	ids = data['ID']
	zs = data['Z_MUSE']
	sconf = data['Sofi_confidence']

if fitcat == 'EAGLE':
	periodic = '0,0,0'  # ''1,1,1'#
	dims = 3
	nz = {'10': 34, '11': 29, '12': 25}
	lcube = 4096
	dovar = False
	csub = False
	redshift = 3.5  # 2+(15-snap)/3.
	cat = '%s/EAGLE/cats/gals_snap%d.fits' % (foldercat, snap)
	data = getdata(cat, 1)
	ids = data['ID']
	ngal = len(ids)
	zs = np.zeros(ngal) + redshift
	sconf = np.zeros(ngal) + 2
	cubename = '%s_snap%d_%s%s.fits' % (cubetype, snap, coord, mask * '.mask')  # 'LLS.fits'#'Fluxcube.fits'#
	dovar = False
	csub = False
	kpc2arcsec = 0.12750223128904756  # 7.843  pkpc per arcsec
	pix2kpc = 6.103515625  # ~6 kpc per pixel
	pix2arcsec = 0.1937778538085323  # 0.7782118608950657/(1+z)
	conv2SB = 1.6512175984082942  # 1/(pix2arcsec)**2

if fitcat == 'mosaic10':
	cat2 = '%s/%s/cats/lae_pairs_UDF10.fits' % (foldercat, 'mosaic')
	fitcat = 'mosaic'
elif fitcat == 'EAGLE':
	cat2 = '%s/%s/cats/lae_pairs_snap%d.fits' % (foldercat, fitcat, snap)
else:
	cat2 = '%s/%s/cats/lae_pairs.fits' % (foldercat, fitcat)
data = getdata(cat2, 1)

if fitcat == 'HDFS':
	cubename = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
	filevar = '%s/%s/%s' % (folder, fitcat, cubename)

if fitcat == 'UDF':
	# cubename = 'DATACUBEFINALuser_20141021T055145_72f74684.fits'
	cubename = 'UDF10.z1300.fits'  # 'UDF/DATACUBE-UDF10.fits'
	filevar = '%s/%s/%s' % (folder, fitcat, 'UDF10.z1300.fits')

if fitcat == 'mosaic' or fitcat == 'mosaic10':
	cubename = 'DATACUBE_UDF-MOSAIC.z1300.fits'  #####

filename = '%s/%s/%s' % (folder, fitcat, cubename.replace('.fits', '.csub' * csub + '.fits'))

id1 = data['id1'].astype(int)
id2 = data['id2'].astype(int)
x1 = np.round(data['x1']).astype(int)
y1 = np.round(data['y1']).astype(int)
z1 = np.round(data['z1']).astype(int)
x2 = data['x2']
y2 = data['y2']
z2 = data['z2']
if fitcat == 'EAGLE':
	#check this!!!
	xl, yl, zl = [501, 501, 21]
	# xc, yc, zc = [251, 251, 7]
	dist = data['com_dist']
	theta = data['theta%s' % coord]
	r1 = id1 - id1 + redshift
	r2 = r1
	xt = data['xt']
	yt = data['yt']
	zt = data['zt']
	shearx = data['shearx']
	sheary = data['sheary']
	shearz = data['shearz']
else:
	xl, yl, zl = [201, 201, 13]
	# xc, yc, zc = [101, 101, 7]
	r1 = data['redshift1']
	r2 = data['redshift2']
	dist = data['pi_Mpc']
	theta = data['theta']

idpairs = (dist <= dmax) & (dist >= dmin) & (theta > 6) \
		  & (r1 < redmax) & (r2 < redmax) & (r1 > redmin) & (r2 > redmin) \
		  & (id1 > idmin) & (id1 < idmax) & (x1 > xmin) & (y1 > ymin) & (x2 > xmin) & (y2 > ymin) \
		  & (x1 < xmax) & (y1 < ymax) & (x2 < xmax) & (y2 < ymax)  # & (r1 < 4) & (r2 < 4)

if i1 != -1: idpairs &= (id1 == i1)
if i2 != -1: idpairs &= (id2 == i2)

if fitcat != 'EAGLE':
	sconf1 = data['sconf1']
	sconf2 = data['sconf2']
	idpairs &= (sconf1 > 0) & (sconf2 > 0)

npairs = np.sum(idpairs)
print 'Npairs', npairs
_n0 = max(0, nmin)
_n1 = min(npairs, nmax)
ids1 = id1[idpairs][_n0:_n1]


if simon:
	sstype = 'cubical'
else:
	sstype = 'cylindrical'
if dorandom: sstype = 'random_' + sstype
if norm:
	ssnorm = 'normalized'
else:
	ssnorm = 'not-normalized'

uids1 = np.unique(ids1)
nuids = len(uids1)
print 'Unique ids', nuids

if parallel:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	print "Parallel!, rank", rank, 'ncores', size
	r0, r1 = [nuids*rank/size, nuids*(rank+1)/size]
	print 'range', r0, r1
	_ids = uids1[r0: r1]
else:
	_ids = uids1

ngen = 0
_cool = []
for i in _ids:
	cool = idpairs & (id1 == i)
	_cool.append(cool)
	ngen += np.sum(cool)

print 'Generating %d subcubes ' % ngen

for i, cool in zip(_ids, _cool):
	if fitcat == 'EAGLE':
		print 'Transforming coords centered on id', i
		print 'Number of neighbors', np.sum(cool)

		if coord == 'x': c1, c2, c3 = [yt[cool] - lcube / 2,  zt[cool] - lcube / 2, shearx[cool]]
		if coord == 'y': c1, c2, c3 = [xt[cool] - lcube / 2, zt[cool] - lcube / 2, sheary[cool]]
		if coord == 'z': c1, c2, c3 = [xt[cool] - lcube / 2, yt[cool] - lcube / 2, shearz[cool]]

	else:
		c1, c2, c3 = [x2[cool] - x1[cool], y2[cool] - y1[cool], z2[cool] - z1[cool]]

	coordconv([i], id2[cool], c1, c2, c3, xw, yw, zw, xmore)
