#!/usr/bin/env python
__author__ = 'gallegos'
import glob
#import h5py
# import eagleSqlTools
import matplotlib.pyplot as plt
import os
import scipy.interpolate as inter
from math import sqrt
#from matplotlib.cm import get_cmap
#from matplotlib.colors import LinearSegmentedColormap
from pyfits import getdata, PrimaryHDU
from sys import argv
import params
import numpy as np

from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,
coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
snaps = rdarg(argv, 'snap', list, [8, 9, 10, 11, 12, 13, 14], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
halfplot = rdarg(argv, 'halfplot', bool, False)
extraname = rdarg(argv, 'extraname', str, '')#'HM12')#
caseA = rdarg(argv, 'caseA', bool, True)#'HM12')#
scase = '_CaseA' * caseA
galcov = rdarg(argv, 'galcov', bool, False)
laemask = rdarg(argv, 'laemask', bool, False)
histover = rdarg(argv, 'histover', bool, False)
LLScubes = rdarg(argv, 'LLScubes', bool, False)
sql = rdarg(argv, 'sql', bool, False)
sphericalLLS = rdarg(argv, 'sphericalLLS', bool, False)
pairLLS = rdarg(argv, 'pairLLS', bool, False)
circlehist = rdarg(argv, 'circlehist', bool, False)
h2d = rdarg(argv, 'h2d', bool, False)
kde = rdarg(argv, 'kde', bool, False)
nhiprof = rdarg(argv, 'nhi', bool, False)
minres = rdarg(argv, 'minres', int, 512)
maxres = rdarg(argv, 'maxres', int, 4096)
model = rdarg(argv, 'model', str, 'HM01')#'HM12')#
npref = rdarg(argv, 'npref', int, 8)
rad = rdarg(argv, 'rad', int, 6)  # arcsec
radprof = rdarg(argv, 'radprof', bool, False)
sbhist = rdarg(argv, 'sbhist', bool, False)
_ssthr = rdarg(argv, 'ssthr', float, 6.73e-3)#None) or 1e10
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
type = rdarg(argv, 'type', str, 'NHI') #NHtot, f_NHI, NHII
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)
do_not = 0


cdict1 = {'red': ((0.0, 0.0, 0.0),
				  (0.5, 0.0, 0.1),
				  (1.0, 1.0, 1.0)),

		  'green': ((0.0, 0.0, 0.0),
					(1.0, 0.0, 0.0)),

		  'blue': ((0.0, 0.0, 1.0),
				   (0.5, 0.1, 0.0),
				   (1.0, 0.0, 0.0))
		  }
colors = [(0, 0, 1), (1, 1, 0), (1, 0, 0)]
#cm = LinearSegmentedColormap.from_list('sofi', colors, N=50)


gamma_bkg = params.gamma_bkg[model]
heat_bkg = params.heat_bkg[model]
snames = params.snames
redshifts = params.redshifts
asec2kpcs = params.asec2kpcs
ssthr = params.nhSS[model]
if _ssthr is not None:
	ssthr = {}
	for snap in snaps:
		ssthr[str(snap)] = _ssthr
dz = params.dz
sbpeaks = params.sbpeaks
zlens = params.zlens
dndz = params.dndz
fLLScorr = params.fLLScorr

nhi_fit = params.nhi_fit[model]

lcube = maxres # 4096
coml = 25  # cMpc
com2pix = 163.84  # lcube/coml
kpc2pix = lcube / float(coml * 1e3)
rads = np.array([0, 2, 4, 8, 12, 20, 30, 50, 100, 200])

lognhi = [10, 11, 12, 13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21, 22, 23, 30]
nhi = np.power(10, lognhi)
sb = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96, 0.98, 1, 1, 1, 1, 1, 1]

nhi2sb = inter.interp1d(nhi, sb)
sb2nhi = inter.interp1d(sb, nhi)


cubefolder = '/net/galaxy-data/export/galaxy/shared/MUSE/gallegos/EAGLE/'
#cubefolder = '../../EAGLE/'

do_inter = False

b = [0, 	0, 		0, 			0, 			0]
brun = rdarg(argv, 'run', list, b, int)
chombo, init_evol, col_dens, do_layers, do_dndx = brun
gal_mask, do_nhfrac, do_lls = 0, 0, 0

lmax = int(np.log(maxres/minres)/np.log(2)) #3 if res 4096, 2 for 2048, 1 for 1024, 0 for 512

hdu = PrimaryHDU()
y0 = -1.1549019599857431  # np.log10(.07)
y1 = -0.1549019599857431  # np.log10(.7)
y2 = 0  # np.log10(1)
x0 = 16.5
x1 = 18
x2 = 20

a1 = 0.666666  # (y1-y0)/(x1-x0)
a2 = .0775  # (y2-y1)/(x2-x1)
b1 = -12.154988  # y1-a1*x1
b2 = -1.55  # y2-a2*x2
NHI = {}
data = {}
fLLS = {}
minNHI = 17.3

if do_inter:
	plt.close()
	nhi2sb2 = inter.InterpolatedUnivariateSpline(nhi, sb)
	plt.plot(lognhi, sb, 'bo', label='Data')
	x = np.arange(13, 21, .01)
	plt.plot(x, nhi2sb2(10 ** x), 'k--', label='Spline fit')
	plt.plot(x, nhi2sb(10 ** x), 'k--', label='1d fit', color='red')
	plt.semilogy()
	plt.legend()
	plt.xlim(15, 21.5)
	plt.hlines(y=1, xmin=15, xmax=21.5)
	plt.ylim(.01, 2.1)
	plt.savefig('interp.png')
	plt.close()

	y = np.arange(.0001, 1, .001)
	plt.plot(sb, nhi, 'bo', label='Data')
	plt.plot(y, sb2nhi(y), 'k--', label='1d fit')
	plt.semilogy()
	plt.legend()
	plt.savefig('interp2.png')
	plt.close()

_cub = []
if do_dndx: mats = {}
for snap in snaps:
	s = str(snap)
	saeed = '/net/galaxy-data/export/galaxydata/saeed/EAGLE/RefL0025N0752/'

	print 'snapshot', snap
	red = redshifts[s]
	# sbpeak = sbpeaks[snap]
	asec2kpc = asec2kpcs[s]
	nl = zlens[s]#*3
	print nl, 'layers'
	asec2pix = asec2kpc * (1 + red) * kpc2pix
	print 'asec2pix', asec2pix
	sb2flux = asec2pix ** -2.
	#cat = getdata('../../EAGLE/cats/gals_snap%d.fits' % snap, 1)
	#size = cat['Mmean200']#in pkpc
	#ngals = len(size)
	#rad = np.max(([1]*ngals, size/asec2kpc), 0)*asec2pix
	_rad = rad*asec2pix
	fg = 2.
	# check lower or upper case E in the scientific notation!!!
	_sgamma = '%.2E_%.2E_%.2E' % tuple(np.array(gamma_bkg[s])*fg)
	_sheat = '%.2E_%.2E_%.2E' % tuple(np.array(heat_bkg[s])*fg)
	sgamma = '%.2E %.2E %.2E' % tuple(np.array(gamma_bkg[s])*fg)
	sheat = '%.2E %.2E %.2E' % tuple(np.array(heat_bkg[s])*fg)
	sname = snames[s]
	fssthr = 2.
	nHIssthr = ssthr[s] * fssthr

	fname = '%s/snapshot_%s/' % (saeed, sname)
	if chombo:
		print '\n\n\nChombo'
		chname = fname+'chombo/snap_%s_%d_%d_%d.chombo.hdf5' % \
				(sname, minres, maxres, npref)
		if not os.path.isfile(chname) or overwrite:
			if os.path.isfile(chname):
				print 'removing chombo file'
				os.system('rm %s' % chname)
			glob.os.chdir(fname + 'chombo')
			os.system('./chombo.sh %d %d %d' % (minres, maxres, npref))
			if do_not:
				notification.notify(
					title='Done!',
					message='Particle to Chombo',
					app_name='Nefinia',
					app_icon='~/notifsofi.png'
				)
		glob.os.chdir(scodes)
	if init_evol:
		print '\n\n\nInit evol'
		glob.os.chdir(fname+'init_evol')
		finit = 'SO.snap_%s_%d_%d_%d_%s_%s_%.2E.ts0000' % \
				(sname, minres, maxres, npref, _sgamma, _sheat, nHIssthr)
		if not os.path.isfile(finit) or overwrite:
			if os.path.isfile(finit):
				print 'removing init evol file'
				os.system('rm %s' % finit)
			srun = './init_evol.sh %d %d %d 25. %.3f "%s" "%s" %.2E %s' % \
				   (minres, maxres, npref, red, sgamma, sheat, nHIssthr, 'true'*caseA)
			print srun
			os.system(srun)
			if do_not:
				notification.notify(
					title='Done!',
					message='InitEvol',
					app_name='Nefinia',
					app_icon='~/notifsofi.png'
				)
		glob.os.chdir(scodes)
	if col_dens:
		print '\n\n\nColumn density!'
		glob.os.chdir(fname + 'column_density')
		finit = '../init_evol/SO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000' % \
				(sname, minres, maxres, npref, _sgamma, _sheat, nHIssthr, scase)
		if os.path.isfile(finit):
			nc = 39  # number of cores
			_s = saeed + '/snapshot_%s/column_density/layers/' % sname
			_fl = '%sSO.snap_%s_%d_%d_%d_%s_%s_%.2E.ts0000_var_%s_proj_%d_lmax_%d_l_' % \
							 (_s, sname, minres, maxres, npref, _sgamma, _sheat, nHIssthr, type, 1, lmax)
			isfiles = True
			for i in range(nl):
				_file = '%s%d_%d.fits' % (_fl, i+1, nl)
				_isfile = os.path.isfile(_file)
				if overwrite and _isfile: os.system('rm %s' % _file)
				isfiles &= _isfile
			if not isfiles or overwrite:
				srun = './column_density_layers.sh %s %s %d %d %d' % (finit, type, lmax, nl, nc)
				print srun
				os.system(srun)
				if do_not:
					notification.notify(
						title='Done!',
						message='Column density',
						app_name='Nefinia',
						app_icon='~/notifsofi.png'
					)
		else: print '%s not there?' % finit
		glob.os.chdir(scodes)

	sname = snames[s]
	_s = saeed + '/snapshot_%s/column_density/layers/' % sname

	for coord, proj in zip(['x'], [1]):#zip(coordnames, [1, 2, 3]):
		if gal_mask:
			zw = 3
			outname = '%s/snap%s_%s.galmask_%darcsec.fits' % (cubefolder, snap, coord, rad)
			if not os.path.isfile(outname) or overwrite:
				mask = np.zeros((nl, lcube, lcube))
				y, x = np.ogrid[0: lcube, 0: lcube]
				c = ['x', 'y', 'z']
				c.remove(coord)
				xc, yc, zc = [cat[c[0]] * com2pix, cat[c[1]] * com2pix, (cat[coord] * nl / coml).astype(int)]
				ids = cat['ID']
				cool = cat['U'] < 0
				for i in np.arange(len(ids))[cool]:
					print '%d %d %d %d' % (i, ids[i], xc[i], yc[i])
					gal = ((x - xc[i]) ** 2 + (y - yc[i]) ** 2 < _rad ** 2)
					mask[zc[i]-zw: zc[i]+zw+1, gal] = ids[i]
				hdu.data = mask
				hdu.writeto(outname, clobber=True)

		if do_layers:
			print 'Combinig layers!'
			outname = '%s/snap%d_%s_%s_%d_%d_%d_%d_%s_%.2E%s.%s%s.fits' % \
					  (cubefolder, snap, coord, model, minres, maxres, npref, nl, _sgamma, nHIssthr, scase, type, extraname)
			if not os.path.isfile(outname) or overwrite:
				cubes = []
				for i in range(nl):
						flayer = '%sSO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000_var_%s_proj_%d_lmax_%d_l_%d_%d.fits' % \
								 (_s, sname, minres, maxres, npref, _sgamma, _sheat, nHIssthr, scase, type, proj, lmax, i+1,
								  nl)
						cubes.append(getdata(flayer))

				cubes = np.array(cubes)
				hdu.data = cubes
				hdu.writeto(outname, clobber=True)
				print 'done %s cube' % type

			else:
				cubes = getdata(outname)

			if do_not:
				notification.notify(
					title='Done!',
					message='Layers',
					app_name='Nefinia',
					app_icon='~/notifsofi.png'
				)

			if type == 'density' or type == 'uniform-density':
				print 'Mean density (contrast)', np.nanmean(cubes)
		if do_nhfrac:
			cubes = getdata('%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E.%s.fits' %
						(cubefolder, snap, coord, model, minres, maxres, npref, nl, nHIssthr, type))
			logcube = np.log10(cubes)
			nhis = np.arange(15, 22, .1)
			cool = (logcube > nhis[0]) & (logcube < nhis[-1])
			hist, bins = np.histogram(logcube[cool], nhis)
			sumh = float(np.sum(hist))
			cumh = [np.sum(hist[i:])/sumh for i in range(len(hist))]
			_dndz = np.array(cumh)*nl/dz[s]

			plt.figure(figsize=(7, 7))
			plt.scatter(bins[:-1], np.log10(hist)-bins[:-1])
			plt.xlabel(type)
			plt.ylabel('log(f_%s)' % type)
			#plt.ylim([-4.2, 0.2])
			plt.savefig('../../EAGLE/%sfrac_%s_%s_%d_%d_%d_%.2E.jpg' %
						(type, coord, model, minres, maxres, nl, nHIssthr))
			plt.close()

		if do_lls:
			zs, ls, a = 3, 1.46, 1.7
			_ls, _a = .11, .22
			def l(z, zs=3, ls=1.46, a=1.70):
				return ls * ((1 + z) / (1 + zs)) ** a
			fnhiname = '%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.NHI.fits' % (cubefolder, snap, coord, model, minres, maxres, npref, nl, nHIssthr, scase)

			fllsname = '%s/snap%d_%s_%s_%d_%d_%d%s.LLS.fits' % (cubefolder, snap, coord, model, minres, maxres, nl, scase)
			fllsnamecorr = '%s/snap%d_%s_%s_%d_%d_%d%s.LLScorr.fits' % (cubefolder, snap, coord, model, minres, maxres, nl, scase)
			if not os.path.isfile(fllsname) or overwrite:
				cubes = getdata(fnhiname)
				_cubes = np.copy(cubes)
				llsNHI = 17.5
				lls0 = cubes > 10 ** llsNHI
				#lls = np.zeros(cubes.shape)
				#lls[cool] = 1
				hdu.data = lls0.astype(float)#lls
				hdu.writeto(fllsname, clobber=True)
				pre = 10
				lz = l(red, zs, ls, a)

				if 1:#not os.path.isfile(fllsnamecorr) or overwrite:
					if 1:
						dndz0 = np.mean(lls0) * nl / dz[s]
						dndz = dndz0
						diff = dndz / lz
						#nhi = cubes*lz/dndz
						#dndz2 = np.mean(nhi>10**llsNHI)*nl/dz[s]
						#print 'new dndz', dndz2


						while abs(diff-1)>1e-3:
							_cubes /= diff
							lls = _cubes > 10 ** llsNHI
							dndz = np.mean(lls) * nl / dz[s]
							diff = dndz / lz
							print 'diff', diff, 'dndz', dndz
					else:
						lls = cubes > 10 ** nhi_fit[s]

					hdu.data = _cubes
					hdu.writeto(fllsname.replace('.LLS.', '.NHIcorr.'), clobber=True)
					hdu.data = lls.astype(float)
					hdu.writeto(fllsnamecorr, clobber=True)
			else:
				lls0 = getdata(fllsname)

			#error propagation
			error = np.sqrt(((1 + red) / (1 + zs)) ** a * _ls**2 * ls * a * ((1 + red) / (1 + zs)) ** (a-1) * _a**2)
			lz = l(red, zs, ls, a)

			print 'dndz sim %.2f simcorr %.2f obs %.2f +-%.2f corr %.2f +-%.2f' % (dndz0, dndz, lz, error, lz / dndz, error / dndz)

		if do_dndx:
			print 'dndx!'
			_type = type#'NHI'
			cname = 'snap%d_%s_%s_%d_%d_%d_%d_%s_%.2E%s.%s' %\
					   (snap, coord, model, minres, maxres, npref, nl, _sgamma, nHIssthr, scase, _type)
			if _type=='NHIcorr':
				cname = 'snap%d_%s_%s_%d_%d_%d%s.%s' % \
						(snap, coord, model, minres, maxres, nl, scase, _type)
			cubename = '%s/%s.fits' % (cubefolder, cname)
			zw = 0#nl
			print cubename
			v = 2.99792458e5 * dz[s] * zw / float(nl) / (1 + red)
			fdat = '../../UVB/dndzdx_%s_snap%s_zw%d_%s_ssthr%.0e_%d_%d_%d%s.%s.dat' % (model, snap, zw, _sgamma, nHIssthr, minres, maxres, npref, scase, _type)
			if not os.path.isfile(fdat) or overwrite:
				cubes = getdata(cubename)
				if type == 'density':
					cubes *= params.meanbardens*coml*params.Mpc2cm/float(nl)/(1.+red)
				zl, yl, xl = cubes.shape

				zn = zl / (zw+1)
				if zw == 0: _c = cubes
				else:
					_c = []
					for i in range(zn): _c.append(np.sum(cubes[i:i+zw+1, :, :], 0))
					_c = np.array(_c)
				lc = np.log10(_c)
				dxdz = (1.+red)**2/sqrt(params.omega_m*(1+red)**3+params.omega_l)
				_lls = np.arange(15, 22.6, .25)
				dndx, dndz, cddf = [[], [], []]

				for i in range(len(_lls)-1):
					cool = (lc>=_lls[i])&(lc<_lls[i+1])
					dNHI = float(10.**_lls[i+1]-10**_lls[i])
					#lpos = np.nanmean(_c[cool])#(10.**_lls[i+1]+10**_lls[i])/2.
					#if np.isnan(lpos): lpos = 10**_lls[i]
					scool = np.mean(cool)*zn
					_dndz = scool/dz[s]
					_dndx = _dndz/dxdz
					_cddf = _dndx/dNHI
					print 'NHI %.2f %.2f' %(_lls[i], _lls[i+1]), 'dndz', _dndz, 'cddf', _cddf, 'dndx', _dndx
					dndx.append(_dndx)
					dndz.append(_dndz)
					cddf.append(_cddf)

				if do_not:
					notification.notify(
						title='Done!',
						message='dndx',
						app_name='Nefinia',
						app_icon='~/notifsofi.png'
					)
				lnhis = _lls[:-1]
				ldndx = np.log10(dndx)
				dndz = np.array(dndz)
				cddf = np.log10(cddf)
				#_dndz = [np.sum(dndzs[i:]) for i in range(len(dndzs))]
				sb = nhi2sb(10**lnhis)
				mat = np.array([lnhis, dndz, dndx, cddf, sb])
				np.savetxt(fdat, mat.T, header='logNHI dndz dndx cddf SBfit', fmt='%.3f')
			else: mat = np.loadtxt(fdat).T

