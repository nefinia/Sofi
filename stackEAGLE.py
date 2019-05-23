__author__ = 'gallegos'

import os
from os.path import isdir
from pyfits import getdata, PrimaryHDU
from sys import argv

import numpy as np

from tools_sofi import rdarg, pdfim


def sclipping(fits, nsigma, dim=None, mask=None, iter=3):
	# print 'Asign nan to values > nsigma in a fits array'
	for i in range(iter):
		if dim is None:
			stds = np.nanstd(fits[:, mask])
		else:
			stds = np.nanstd(fits[:, mask], dim)
		high_sigma = np.abs(fits) > nsigma * stds
		fits[high_sigma] = np.nan
	return fits, stds

doanalysis = rdarg(argv, 'analysis', bool, False)
binsize = rdarg(argv, 'binsize', int, 2)
contours = rdarg(argv, 'contours', bool, True)
dosclip = rdarg(argv, 'dosclip', bool, True)
extraname = rdarg(argv, 'extraname', str, '')
fitcat = rdarg(argv, key='fitcat', type=str, default='EAGLE')
folder = rdarg(argv, 'folder', str, '/net/galaxy-data/export/galaxydata/gallegos/')  # '/scratch/gallegos/MUSE/'
foldercat = rdarg(argv, 'foldercat', str, '../../')  # '/scratch/gallegos/MUSE/'
folderout = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
fshear = rdarg(argv, 'flipshear', bool, False)
highq = rdarg(argv, 'highq', bool, False)
imtype = rdarg(argv, 'imtype', str, 'mean')  # 'mean'
jin = rdarg(argv, 'jin', int, 0)
makeim = rdarg(argv, 'makeim', bool, True)
mask = rdarg(argv, 'mask', bool, False)
makepdf = rdarg(argv, 'makepdf', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
parallel = rdarg(argv, 'parallel', bool, True)
ncores = rdarg(argv, 'ncores', int, 2)
prename = rdarg(argv, 'prename', str, '')
propvar = rdarg(argv, 'pvar', bool, False)
statvar = rdarg(argv, 'svar', bool, False)
reject = rdarg(argv, 'reject', bool, False)
rotate = rdarg(argv, 'rotate', bool, False)
sb = rdarg(argv, 'sb', bool, False)
scalelims = rdarg(argv, 'scalelims', str, '-0.01 0.01')
sclip = rdarg(argv, 'sclip', int, 3)
scliptype = rdarg(argv, 'scliptype', int, 1)
simon = rdarg(argv, 'cubical', bool, True)
smooth = rdarg(argv, 'smooth', bool, False)
std = rdarg(argv, 'std', float, 2)  # 0.0707064#

vmin = rdarg(argv, 'vmin', float, 0)
vmax = rdarg(argv, 'vmax', float, 0.1)

xmin1 = rdarg(argv, 'xmin1', float, 0)
xmax1 = rdarg(argv, 'xmax1', float, 4096)
ymin1 = rdarg(argv, 'ymin1', float, 0)
ymax1 = rdarg(argv, 'ymax1', float, 4096)
xmin2 = rdarg(argv, 'xmin2', float, 0)
xmax2 = rdarg(argv, 'xmax2', float, 4096)
ymin2 = rdarg(argv, 'ymin2', float, 0)
ymax2 = rdarg(argv, 'ymax2', float, 4096)
yw0 = rdarg(argv, 'yw0', int, 2)
zw0 = rdarg(argv, 'zw0', int, 2)

if imtype == 'median': extraname = '.' + imtype
if not prename: prename = ''
if not extraname: extraname = ''
smask = '.mask' * mask
# d distance
# p proj distance
# r redshift
# los line of sight distance
# n number of neighbours
# vel velocity
dmin = rdarg(argv, 'dmin', float, .5)
dmax = rdarg(argv, 'dmax', float, 20)
d5min = rdarg(argv, 'd5min', float, .5)
d5max = rdarg(argv, 'd5max', float, 20)
pmin = rdarg(argv, 'pmin', float, 16)
pmax = rdarg(argv, 'pmax', float, 2000)
rmin = rdarg(argv, 'rmin', float, 2.9)
rmax = rdarg(argv, 'rmax', float, 4.0)
lmin = rdarg(argv, 'umin', float, 0)
lmax = rdarg(argv, 'umax', float, 99)
losmin = rdarg(argv, 'losmin', float, 0.5)
losmax = rdarg(argv, 'losmax', float, 20)
nmin = rdarg(argv, 'nmin', int, 0)
nmax = rdarg(argv, 'nmax', int, 1000)
nwmin = rdarg(argv, 'nwmin', float, 0)
nwmax = rdarg(argv, 'nwmax', float, 1)
velmin = rdarg(argv, 'velmin', float, 0)
velmax = rdarg(argv, 'velmax', float, 10000)
snap = rdarg(argv, 'snap', int, 10)
coord = rdarg(argv, 'coord', str, 'x')

cat2 = '%s/%s/cats/lae_pairs_snap%d.fits' % (foldercat, fitcat, snap)
data = getdata(cat2, 1)
ids1 = data['id1']
ids2 = data['id2']
laedata = getdata('%s/%s/cats/gals_snap%d.fits' % (foldercat, fitcat, snap), 1)
idlae = laedata['ID']
flae = idlae - idlae + 1
pdists = data['x2'] - data['x1']
zs1 = data['z1']
zs2 = data['z2']
dists = data['com_dist']  # data['pi_Mpc']
theta = data['theta%s' % coord]
fitcats = [fitcat] * len(data)
nw1 = data['nw1']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids1])
nw2 = data['nw2']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids2])
d5th1 = data['n5d1']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids1])
d5th2 = data['n5d2']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids2])
u1 = data['U1']
abu1 = np.abs(u1)
u2 = data['U2']
xs1 = data['x1']
xs2 = data['x2']
ys1 = data['y1']
ys2 = data['y2']

h = 'theta '
hdu = PrimaryHDU()
_close = np.where((dists <= dmax) & (dists > dmin) & (theta <= pmax)
				 & (theta > pmin) & (d5th1 > d5min) & (d5th1 < d5max) & (abu1 > lmin) & (abu1 <= lmax)
				 & (xs1 > xmin1) & (ys1 > ymin1) & (xs1 <= xmax1) & (ys1 <= ymax2)
				 & (xs2 > xmin2) & (ys2 > ymin2) & (xs2 <= xmax2) & (ys2 <= ymax2))[0]
nclose = len(_close)

if parallel:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	print "Parallel!, rank", rank, 'ncores', size
	r0, r1 = [nclose*rank/size, nclose*(rank+1)/size]
	print 'range', r0, r1
	close = _close[r0: r1]
else:
	close = _close

print 'Number of close neighbours', len(close)

if len(close) > 1:
	id1 = ids1[close]
	id2 = ids2[close]
	z1 = zs1[close]
	z2 = zs2[close]
	fcat = np.array(fitcats)[close]
	dist = dists[close]
	theta = theta[close]
	lst = []


	for i in range(len(id1)):
		pair = '%s/%s/pairs/SB_snap%d_%s/%d-%d%s%s.fits' % (folder, fcat[i], snap, coord, id1[i], id2[i], prename, smask)
		if os.path.isfile(pair): lst.append(pair)

	nstack = len(lst)
	print '\n Number of existing subcubes to be stacked', nstack

	foldername = '%s/%s/stacks/snap%d_%s_d%d-%d_pd%d-%d_d5th%.1f-%.1f_u%.1f-%.1f%s%s%s' % (
		folderout, fitcat, snap, coord, dmin, dmax,
		pmin, pmax, d5min, d5max, lmin, lmax, prename, extraname, smask)

	print "Output files in", foldername
	if not isdir(foldername): os.system('mkdir %s' % foldername)

	lstname = '%s/stack.lst' % (foldername)
	flst = open(lstname, 'w')
	for l in lst: flst.write('%s\n' % l)
	flst.close()
	title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$,' % (
		nstack, dmin, dmax, pmin, pmax) + '\n' + extraname

	###################
	## Stacking part ##
	###################
	print 'Stacking part!'
	stackname = '%s/stack.fits' % foldername
	nstack = len(lst)

	if not os.path.isfile(stackname) or overwrite:
		fits = np.array([getdata(l) for l in lst])

		if imtype == 'mean':
			f = np.nanmean(fits, 0)
			std = np.nanstd(fits, 0)
		if imtype == 'median': f = np.nanmedian(fits, 0)
		zl, yl, xl = f.shape

	else:
		f = getdata(stackname)
		zl, yl, xl = f.shape

	if parallel:
		fcomb = None
		stdcomb = None
		if rank == 0:
			fcomb = np.empty((size, zl, yl, xl), dtype='d')
			stdcomb = np.empty((size, zl, yl, xl), dtype='d')
		comm.Gather(f, fcomb, root=0)
		comm.Gather(std, stdcomb, root=0)
		if rank == 0:
			fcomb = np.nanmean(fcomb, 0)
			stdcomb = np.nanmean(stdcomb, 0)
			hdu.data = fcomb
			hdu.writeto(stackname, clobber=True)
			hdu.data = stdcomb
			hdu.writeto(stackname.replace('.fits', '.STD.fits'), clobber=True)
			hdu.data = np.nanmean(fcomb[zl / 2-zw0:zl/2+zw0+1, :, :], 0)
			hdu.writeto(stackname.replace('.fits', '.IM.fits'), clobber=True)
		if 0:
			if rank != 0:
				comm.send(f.shape, dest=0)
				comm.Send(f, dest=0)
			else:
				print 'Combining all %d cores stacks' % size
				fs = [f]
				for i in range(size-1):
					print 'Receiving file from %d' % rank
					shape = comm.recv(source=i+1)
					print 'Shape of data to receive:', shape
					_f = np.empty(shape, dtype='d')  # allocate space to receive the array
					comm.Recv(_f, source=0)
					fs.append(_f)
				if imtype == 'mean':
					f = np.nanmean(fs, 0)
					std = np.nanstd(fs, 0)
					hdu.data = std
					hdu.writeto(stackname.replace('.fits', '.STD.fits'), clobber=True)
				if imtype == 'median': f = np.nanmedian(fs, 0)
				zl, yl, xl = f.shape
				hdu.data = f
				hdu.writeto(stackname, clobber=True)
				hdu.data = np.nanmean(f[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
			hdu.writeto(stackname.replace('.fits', '.IM.fits'), clobber=True)
			
	print 'stacks shape', xl, yl, zl

	if 0:
		astroim(stackname, smooth=smooth, saveim=makeim,
				show=False, highq=highq,
				cbfrac=.08, pad=.006, dfig=(8, 10), contours=True,
				title='', vmin=vmin, vmax=vmax, gray=True,
				xmin=-xl / 2, xmax=xl / 2, ymin=-yl / 2, ymax=yl / 2,
				zpw=1, sb=True)

	if makepdf:
		imout = lstname.replace('.lst', '.pdf')
		if not os.path.isfile(imout) or overwrite:
			lstIM = [l for l in lst]
			pdfim(lstIM, fcats=fcat, imout=imout, contours=False)

