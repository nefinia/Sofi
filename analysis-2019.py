#!/usr/bin/env python
__author__ = 'gallegos'
import glob
# import h5py
# import eagleSqlTools
# import matplotlib.pyplot as plt
import os
from math import sqrt
from sys import argv

import numpy as np
import scipy.interpolate as inter
# from matplotlib.cm import get_cmap
# from matplotlib.colors import LinearSegmentedColormap
from pyfits import getdata, PrimaryHDU

import params
from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
snaps = rdarg(argv, 'snap', list, [8, 9, 10, 11, 12, 13, 14], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
caseA = rdarg(argv, 'caseA', bool, True)
scase = '_CaseA' * caseA
halfplot = rdarg(argv, 'halfplot', bool, False)
extraname = rdarg(argv, 'extraname', str, '')#'HM12')#
galcov = rdarg(argv, 'galcov', bool, False)
laemask = rdarg(argv, 'laemask', bool, False)
histover = rdarg(argv, 'histover', bool, False)
LLScubes = rdarg(argv, 'LLS', bool, False)
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
npref = rdarg(argv, 'npref', int, 12)
rad = rdarg(argv, 'rad', int, 3)  # arcsec
radprof = rdarg(argv, 'radprof', bool, False)
sbhist = rdarg(argv, 'sbhist', bool, False)
_ssthr = rdarg(argv, 'ssthr', float, None)#6.73e-3 or 1e10
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
type = rdarg(argv, 'type', str, 'NHI') #NHtot, f_NHI, NHII
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)
zw = rdarg(argv, 'zw', int, 5)  # arcsec
do_not = 0
if do_not: from plyer import notification


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


if temperature:
	res, lmax = 512, 0
	f0 = '/net/galaxy-data/export/galaxydata/saeed/EAGLE/RefL0025N0752/'
	colors = {'10': 'red', '11': 'green', '12': 'blue'}
	for snap in snaps:
		data = getdata('../../EAGLE/cats/gals_snap%d.fits' % snap, 1)
		ids = data['ID']
		if 0:
			xc, yc, zc, u, ids = data['x'], data['y'], data['z'], data['U'], data['ID']
			ngal = len(xc)
			cube = getdata('%s/snapshot_%s/fits_3d/snap_%s_%s_%s_8_var_uniform-temperature_lmax_%d.fits' %
						   (f0, snames[str(snap)], snames[str(snap)], res, res, lmax))
			zl, yl, xl = cube.shape
			flux2sb = zl / float(coml)
			z, y, x = np.ogrid[1: zl + 1, 1: yl + 1, 1: xl + 1]
			T = []  # 20, 40 and 80 kpc
			for i in range(ngal):
				print 'id', ids[i], '%d/%d' % (i, ngal)
				cool = ((.05 * flux2sb) < (z - zc[i] * flux2sb) ** 2 + (y - yc[i] * flux2sb) ** 2 + (x - xc[i] * flux2sb) ** 2) & \
					   ((z - zc[i] * flux2sb) ** 2 + (y - yc[i] * flux2sb) ** 2 + (x - xc[i] * flux2sb) ** 2 < (.08 * flux2sb))
				_t = np.nanmean(cube[cool])
				print 'Temperature between 50 and 80 kpc %f, U %f, number of pixels used %d' % (_t, u[i], np.sum(cool))
				T.append(_t)
			np.savetxt('../../EAGLE/cats/temp_snap%d' % snap, np.array([ids, T]).T, header='ID Temperature')
		if 1:
			i1, i2, fconns, frands = np.loadtxt('../../EAGLE/analysis/avsb_snap%d_x.dat' % snap).T
			fconn = []
			for i in ids:
				cool = i1 == i
				fconn.append(np.nanmean(frands[cool]))

			np.savetxt('../../EAGLE/cats/fconn_snap%d' % snap, np.array([ids, fconn]).T, header='ID fconn')

		else:
			T = data['Temperature']
			U = data['U']
			plt.scatter(U, T, label='z %.2f' % redshifts[str(snap)], color=colors[str(snap)])
	if 0:
		plt.semilogy()
		plt.xlabel(r'U')
		plt.ylabel(r'T[K]')
		plt.legend(scatterpoints=1)
		plt.show()

if LLScubes:
	cubefolder = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/'
	#cubefolder = '../../EAGLE/'
	#cubefolder = '/scratch/gallegos/EAGLE/'
	
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
		fg = 1.
		# check lower or upper case E in the scientific notation!!!
		_sgamma = '%.2E_%.2E_%.2E' % tuple(np.array(gamma_bkg[s])*fg)
		_sheat = '%.2E_%.2E_%.2E' % tuple(np.array(heat_bkg[s])*fg)
		sgamma = '%.2E %.2E %.2E' % tuple(np.array(gamma_bkg[s])*fg)
		sheat = '%.2E %.2E %.2E' % tuple(np.array(heat_bkg[s])*fg)
		sname = snames[s]

		fname = '%s/snapshot_%s/' % (saeed, sname)
		if chombo:
			print '\n\n\nChombo'
			chname = fname+'chombo/snap_%s_%d_%d_%d.chombo.h5' % \
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
			finit = 'SO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000' % \
					(sname, minres, maxres, npref, _sgamma, _sheat, ssthr[s], scase)
			if not os.path.isfile(finit) or overwrite:
				if os.path.isfile(finit):
					print 'removing init evol file'
					os.system('rm %s' % finit)
				srun = './init_evol.sh %d %d %d 25. %.3f "%s" "%s" %.2E' % \
					   (minres, maxres, npref, red, sgamma, sheat, ssthr[s])
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
					(sname, minres, maxres, npref, _sgamma, _sheat, ssthr[s], scase)
			if os.path.isfile(finit):
				nc = 39  # number of cores
				_s = saeed + '/snapshot_%s/column_density/layers/' % sname
				_fl = '%sSO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000_var_%s_proj_%d_lmax_%d_l_' % \
								 (_s, sname, minres, maxres, npref, _sgamma, _sheat, ssthr[s], scase, type, 1, lmax)
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
			else: print '%s 1 um' % finit
			glob.os.chdir(scodes)

		sname = snames[s]
		_s = saeed + '/snapshot_%s/column_density/layers/' % sname

		for coord, proj in zip(['x'], [1]):#zip(coordnames, [1, 2, 3]):
			if gal_mask:
				outname = '%s/snap%s_%s.galmask_%darcsec_zw%d.fits' % (cubefolder, snap, coord, rad, zw)
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
				print 'combinig layers!'
				outname = '%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.%s.fits' % (cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[s], scase, type)
				if not os.path.isfile(outname) or overwrite:
					cubes = []
					for i in range(nl):
							flayer = '%sSO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000_var_%s_proj_%d_lmax_%d_l_%d_%d.fits' % \
									 (_s, sname, minres, maxres, npref, _sgamma, _sheat, ssthr[s], scase, type, proj, lmax, i+1, nl)
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
			if do_nhfrac:
				cubes = getdata('%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E.%s.fits' %
							(cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[s], type))
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
							(type, coord, model, minres, maxres, nl, ssthr[s]))
				plt.close()

			if do_lls:
				zs, ls, a = 3, 1.46, 1.7
				_ls, _a = .11, .22
				def l(z, zs=3, ls=1.46, a=1.70):
					return ls * ((1 + z) / (1 + zs)) ** a
				fnhiname = '%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.NHI.fits' % (cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[s], scase)

				fllsname = '%s/snap%d_%s_%s_%d_%d_%d.LLS.fits' % (cubefolder, snap, coord, model, minres, maxres, nl)
				fllsnamecorr = '%s/snap%d_%s_%s_%d_%d_%d.LLScorr.fits' % (cubefolder, snap, coord, model, minres, maxres, nl)
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

				#_z,_y,_x = lls0.shape
				#_lls0 = np.zeros((nl/5, _y, _x))
				#for i in range(nl/5):
			#		_lls0[i] = np.sum(lls0[i: i+6], 0)

				#error propagation
				error = np.sqrt(((1 + red) / (1 + zs)) ** a * _ls**2 * ls * a * ((1 + red) / (1 + zs)) ** (a-1) * _a**2)
				lz = l(red, zs, ls, a)
				if 0:
					lzs = []
					lzs.append(l(red, zs, ls, a))
					lzs.append(l(red, zs, ls+_ls, a))
					lzs.append(l(red, zs, ls-_ls, a))
					lzs.append(l(red, zs, ls, a+_a))
					lzs.append(l(red, zs, ls, a-_a))
					lzs.append(l(red, zs, ls+_ls, a+_a))
					lzs.append(l(red, zs, ls+_ls, a-_a))
					lzs.append(l(red, zs, ls-_ls, a+_a))
					lzs.append(l(red, zs, ls-_ls, a-_a))
					lzmax, lzmin = np.amax(lzs), np.amin(lzs)
				print 'dndz sim %.2f simcorr %.2f obs %.2f +-%.2f corr %.2f +-%.2f' % (dndz0, dndz, lz, error, lz / dndz, error / dndz)

			if do_dndx:
				print 'dndx!'
				_type = type#'NHI'
				cname = 'snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.%s' %\
						   (snap, coord, model, minres, maxres, npref, nl, ssthr[s], scase, _type)
				if _type=='NHIcorr':
					cname = 'snap%d_%s_%s_%d_%d_%d.%s' % \
							(snap, coord, model, minres, maxres, nl, _type)
				cubename = '%s/%s.fits' % (cubefolder, cname)
				zw = 1#nl/5
				print cubename
				v = 2.99792458e5 * dz[s] * zw / float(nl) / (1 + red)
				fdat = '../../UVB/dndzdx_%s_snap%s_zw%d_ssthr%.0e_%d_%d_%d.%s.dat' % (model, snap, zw, ssthr[s], minres, maxres, npref, _type)
				if not os.path.isfile(fdat) or overwrite:
					cubes = getdata(cubename)
					if type == 'density':
						cubes *= (1+red)**2*coml*params.Mpc2cm/nl*params.meanbardens
					zl, yl, xl = cubes.shape

					if 0:
						lc = np.sum(cubes, 0)
						lc = np.log10(lc)
						dzdx = (1.+red)**2/sqrt((1-.714)*(1+red)**3+.714)
						_lls = np.arange(12, 23, .5)
						nhis, dndxs, dndzs = [[], [], []]

						for i in range(len(_lls)-1):
							cool = np.sum((lc>=_lls[i])&(lc<_lls[i+1]))
							dNHI = float(10.**_lls[i+1]-10**_lls[i])
							lpos = (10.**_lls[i+1]+10**_lls[i])/2.
							_dndz = cool/dz[s]/float(maxres)**2

							dndx = _dndz/dzdx/dNHI
							print 'NHI', _lls[i], _lls[i+1], 'dn/dz', _dndz, 'dn/dx', dndx
							nhis.append(lpos)
							dndxs.append(dndx)
							dndzs.append(_dndz)

						lnhis = np.log10(nhis)
						ldndx = np.log10(dndxs)

						_dndz = [np.sum(dndzs[i:]) for i in range(len(dndzs))]
						sb = nhi2sb(nhis)
						mat = np.array([lnhis, dndzs, _dndz, ldndx, sb])
						np.savetxt(fdat.replace('.dat', '.allz.dat'), mat.T, header='logNHI cddf  dndz logdndx SBfit')

					else:
						zn = zl / zw
						if zw == 1: _c = cubes
						
						else:
							_c = []
							for i in range(zn): _c.append(np.sum(cubes[i:i+zw, :, :], 0))
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

if sql:
	conn = eagleSqlTools.connect('sgallego', 'dpTT852J')

	for snap in snaps:
		sql = "SELECT \
				gal.CentreOfMass_x as x,\
				gal.CentreOfMass_y as y,\
				gal.CentreOfMass_z as z,\
				gal.Velocity_x as vx,\
				gal.Velocity_y as vy,\
				gal.Velocity_z as vz,\
				gal.Redshift as redshift,\
				gal.MassType_Star as stellar_mass,\
				gal.StarFormationRate as SFR,\
				gal.MassType_DM as DM_mass,\
				gal.MassType_Gas gass_mass,\
				sizes.R_halfmass30 as size,\
				gal.SubGroupNumber as SubGroupNumber,\
				mag.u_nodust as U,\
				mag.g_nodust as G,\
				mag.r_nodust as R,\
				mag.i_nodust as I,\
				mag.z_nodust as Z,\
				mag.Y_nodust as Y,\
				mag.J_nodust as J,\
				mag.H_nodust as H,\
				mag.K_nodust as K,\
				flux.Johnson_V as f_V,\
				flux.Johnson_R as f_R,\
				fof.Group_R_Mean200 as Rmean200,\
				fof.Group_M_Mean200 as Mmean200\
			  FROM\
				RefL0025N0752_SubHalo as gal,\
				RefL0025N0752_Magnitudes as mag,\
				RefL0025N0752_Aperture as ape,\
				RefL0025N0752_Sizes as sizes,\
				RefL0025N0752_FOF as fof,\
				RefL0025N0752_DustyFluxes as flux\
			  WHERE\
				gal.SnapNum = %d and\
				ape.Mass_Star > 1.0e8 and\
				ape.ApertureSize = 30 and\
				gal.GalaxyID = mag.GalaxyID and\
				gal.GalaxyID = ape.GalaxyID and\
				gal.GalaxyID = SIZES.GalaxyID and\
				gal.GalaxyID = flux.GalaxyID and\
				gal.GroupID = fof.GroupID" % snap

		data = eagleSqlTools.execute_query(conn, sql)
		np.savetxt('../../EAGLE/cats/gals_%d.dat' % snap, data,
				   header='x y z vx vy vz redshift stellar_mass SFR DM_mass gass_mass size SubGroupNumber '
						  'U G R I Z Y J H K f_V f_R Rmean200 Mmean200')

if sbhist:

	def sb_roi(cat, xmin, xmax, ymin, ymax, zmin, zmax):
		data = h5py.File(cat, 'r')
		pz = np.array(data['pz'])
		zgood = (pz >= zmin) & (pz <= zmax)
		f = np.array(data['SB'])[zgood]
		px = np.array(data['px'])[zgood]
		py = np.array(data['py'])[zgood]
		right = (px > xmin) & (px < xmax) & (py > ymin) & (py < ymax)
		left = (px < -xmin) & (px > -xmax) & (py > ymin) & (py < ymax)
		down = (px > ymin) & (px < ymax) & (py < -xmin) & (py > -xmax)
		up = (px > ymin) & (px < ymax) & (py > xmin) & (py < xmax)
		rand = left | down | up
		fright = np.nanmean(f[right])
		frand = np.nanmean(f[rand])
		return fright, frand


	folder = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/cats/'
	props = [['u1', 'd'], ['nw1', 'd']]
	xmin
	for snap in snaps:
		for coord in coordnames:
			catfolder = folder + 'SB_snap%d_%s/' % (snap, coord)
			glob.os.chdir(catfolder)
			cats = glob.glob('*.h5')
			cats.sort()
			id1, id2 = np.array([np.array(c[:-3].split('-')).astype(int) for c in cats]).T

if h2d:

	if 0:
		freexmax = False
		if freexmax:
			extraname = '.freexmax'  # ''_3-7asec'
		else:
			extraname = ''  # ''.d<6'
		from itertools import product

		zr = []  # np.arange(-2,3)

		# trying just u1 and d
		pnameshist = [['u1',
					   'd']]  # [['nw1', 'd'], ['nw1', 'd'], ['nw1', 'u1']]#, ['u2', 'd'], ['nw1', 'd']]  # [pnames[0], pnames[2]]
		# pnameshist = [['u1', 'u2'], ['u1', 'nw1']]#, ['nw1', 'u2']]  # [pnames[0], pnames[2]]

		thetamin = 20
		h_all = []
		n_all = []
		hsb_all = []
		superall = {}
		ssss = 0

		pn0, pn1 = pnameshist[0]
		txt0 = '../../EAGLE/analysis/histograms/hist_XXX_snapYY_%s-%s%s.dat' % (pn0, pn1, extraname)

		for snap in snaps:
			houts = {}
			txt = txt0.replace('YY', str(snap))
			houts['fright'] = txt.replace('XXX', 'fright')
			houts['frand'] = txt.replace('XXX', 'frand')
			houts['nsubs'] = txt.replace('XXX', 'nsubs')

			hbool = True
			for ho in houts: hbool &= os.path.isfile(houts[ho])

			if hbool:  # not hbool:
				for snap in snaps:
					sbpeak = sbpeaks[str(snap)]
					red = redshifts[str(snap)]
					asec2kpc = asec2kpcs[str(snap)]

					sb_dimming = np.power(1 + red, -4)
					all = {}
					thetamin = 16
					# height 2 arcsec
					# width 6 to 12 arcsec -> same as my paper!
					xasecmin = 6
					xasecmax = 12
					yasecmin = -1
					yasecmax = 1

					asec2pix = asec2kpc * (1 + red) * kpc2pix
					sb2flux = asec2pix ** -2.
					# limits in pixels
					xmin = xasecmin * asec2pix
					if not freexmax: xmax = xasecmax * asec2pix
					ymin = yasecmin * asec2pix
					ymax = yasecmax * asec2pix
					zmax = 2
					zmin = -zmax
					rad_asec = np.array([rads[i + 1] for i in range(len(rads) - 1)]) / asec2pix

					fsize = 20
					keys = ['right', 'left', 'top', 'down']

					# Use observable quantities if possible
					# dz = x1-x2
					pairs = getdata('../../EAGLE/cats/lae_pairs_snap%d.fits' % snap, 1)

					pnames = ['u1', 'u2', 'd', 'm1', 'm2', 'nw1']

					params = {'u1': pairs['U1'], 'u2': pairs['U2'], 'd': pairs['com_dist'],
							  'm1': np.log10(pairs['stellar_mass1']), 'm2': np.log10(pairs['stellar_mass2']),
							  'nw1': pairs['nw1']}

					umin = min(params['u1'])
					umax = max(params['u1'])

					urange = [-21.405001, -19.551001, -19.080999, -18.549, -18.066, -17.507999]
					# urange = [-21.405001, -20.125999, -19.551001, -19.306999, -19.080999, -18.822001, -18.549, -18.291,
					#          -18.066, -17.811001, -17.507999]
					pranges = {'u1': urange, 'u2': urange,
							   'd': np.arange(0, 21, 2), 'nw1': np.arange(0.08, .17, .02),
							   'm1': np.arange(10.5, 12.1, .17), 'm2': np.arange(10.5, 12.1, .17)}
					plabels = {'u1': u'$\mathrm{U_1}$', 'u2': u'$\mathrm{U_2}$', 'd': u'$\mathrm{d\,[cMpc]}$',
							   'm1': u'$\mathrm{M_{stars,\,1}}$', 'm2': u'$\mathrm{M_{stars,\,2}}$',
							   'nw1': u'$\mathrm{Number\,of\,neighbors}$'}

					props1 = [[[], []]] * 18
					props2 = [[]] * 4
					pnames = ['id', 'x', 'y', 'z', 'nw', 'U', 'G', 'R', 'I', 'Z', 'Y', 'J', 'H', 'K',
							  'gas_mass', 'DM_mass', 'size', 'stellar_mass']

					for i in range(len(props1)):
						props1[i] = [pairs['%s1' % pnames[i]], pairs['%s2' % pnames[i]]]

					id, x, y, z, nw, u, g, r, i, z, y, j, h, k, mgas, mdm, size, mstars = props1

					for coord in coordnames:

						pairsname = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/cats/SB_snap%d_%s/' \
									% (snap, coord)
						outname2 = '../../EAGLE/analysis/avsb_snap%d_%s%s.dat' % (snap, coord, extraname)

						pnames2 = ['theta%s' % coord, 'com_dist', 'ang%s' % coord, 'shear%s' % coord]

						for i in range(len(props2)): props2[i] = pairs[pnames2[i]]
						theta, dist, angle, shear = props2

						# props = [props1, props2]

						cool = []
						for i1, i2 in zip(id[0], id[1]):
							cool.append(os.path.isfile(pairsname + '%d-%d.h5' % (i1, i2)))
						cool &= (theta > thetamin) & (dist < 20) & (dist > .5)  # (dist<6)
						ncool = np.sum(cool)
						cool = np.where(cool)[0]

						if freexmax:
							xmax = (theta[cool] - 6) * asec2pix
						else:
							xmax = np.zeros(ncool) + xasecmax * asec2pix

						if 1:  # ncool>0:
							condition = (not os.path.isfile(outname2)) or overwrite
							if condition:

								fout = open(outname2, 'w')
								fout.write('#id1 id2 f_right f_rand\n')
								foutz = {}
								for _z in zr:
									foutz[str(_z)] = open(outname2.replace('.dat', '.%d.dat' % _z), 'w')
									foutz[str(_z)].write('#id1 id2 f_right f_rand\n')

								for i in range(ncool):
									print 'Pair %d-%d, snap%d %s. %d of %d' % (
										id[0][cool[i]], id[1][cool[i]], snap, coord, i, ncool)
									cat = '%s%d-%d.h5' % (pairsname, id[0][cool[i]], id[1][cool[i]])
									data = h5py.File(cat, 'r')
									pz = np.array(data['pz'])
									zgood = (pz >= zmin) & (pz <= zmax)
									f = np.array(data['SB'])[zgood]
									px = np.array(data['px'])[zgood]
									py = np.array(data['py'])[zgood]
									right = (px > xmin) & (px < xmax[i]) & (py > ymin) & (py < ymax)
									left = (px < -xmin) & (px > -xmax[i]) & (py > ymin) & (py < ymax)
									down = (px > ymin) & (px < ymax) & (py < -xmin) & (py > -xmax[i])
									up = (px > ymin) & (px < ymax) & (py > xmin) & (py < xmax[i])
									rand = left | down | up
									f_ = np.nanmean(f[right])
									# f_ = np.nanmedian(f[right])
									# f_n = np.nanmean(right)
									f_rand = np.nanmean(f[rand])
									# f_rand = np.nanmedian(f[rand])
									# f_rand_n = np.nanmean(rand)
									fout.write('%d %d %f %f\n' % (id[0][cool[i]], id[1][cool[i]], f_, f_rand))
									for _z in zr:
										# remember f[zgood].....
										rightz = right & (pz == _z)
										randz = rand & (pz == _z)
										f_ = np.nanmean(f[rightz])
										# f_n = np.nanmean(rightz)
										f_rand = np.nanmean(f[randz])
										# f_rand_n = np.nanmean(randz)
										foutz[str(_z)].write(
											'%d %d %f %f\n' % (id[0][cool[i]], id[1][cool[i]], f_, f_rand))

									print 'sb: right %1.3f rand %1.3f' % \
										  (f_, f_rand)
								fout.close()
								for _z in zr: foutz[str(_z)].close()

							else:
								print "Output file(s) already exist(s)", outname2

							_f = np.loadtxt(outname2)
							_fz = {}
							for _z in zr: _fz[str(_z)] = np.loadtxt(outname2.replace('.dat', '.%d.dat' % _z))

							for p in pnameshist:
								pn0, pn1 = p

								if pn0 != pn1:
									r0, r1 = [pranges[pn0], pranges[pn1]]
									n0, n1 = [len(r0) - 1, len(r1) - 1]
									id1 = list(_f[:, 0])
									id2 = list(_f[:, 1])
									cool2 = []
									nid0 = len(id[0])
									nid1 = len(id1)
									for _j in range(len(id[0])):
										end = False
										iii = 0
										while not end:
											if iii >= len(id1):
												end = True
											else:
												if (id1[iii] == id[0][_j]) and (id2[iii] == id[1][_j]):
													id1.remove(id1[iii])
													id2.remove(id2[iii])
													cool2.append(_j)
													end = True
											iii += 1
									_p0, _p1 = [params[pn0][cool2], params[pn1][cool2]]
									f_right = _f[:, 2]
									# f_right_n = _f[:, 3]
									f_rand = _f[:, 3]
									# f_rand_n = _f[:, 5]
									print 'Combining %s with %s' % (pn0, pn1)

									outs = ['../../EAGLE/analysis/histograms/hist_fright_snap%d_%s_%s-%s%s.dat' % (
										snap, coord, pn0, pn1, extraname),
											'../../EAGLE/analysis/histograms/hist_frand_snap%d_%s_%s-%s%s.dat' % (
												snap, coord, pn0, pn1, extraname),
											# '../../EAGLE/analysis/histograms/hist_n_snap%d_%s_%s-%s%s.dat' % (snap, coord, pn0, pn1, extraname),
											'../../EAGLE/analysis/histograms/hist_nsubs_snap%d_%s_%s-%s%s.dat' % (
												snap, coord, pn0, pn1, extraname)]

									isfiles = True
									for out in outs: isfiles &= os.path.isfile(out)

									if isfiles and not overwrite and not histover:
										print 'h2d files exist'
										all[coord, pn0, pn1] = [np.loadtxt(out) for out in outs]
									else:
										if histover and isfiles: 'Overwriting h2d files!'
										if not isfiles: print 'No h2d files!'
										# n = np.zeros((n0, n1))
										nsubs = np.zeros((n0, n1))
										hfright = np.zeros((n0, n1))
										hfrand = np.zeros((n0, n1))

										for i in range(n0):
											for j in range(n1):
												coolh = (_p0 > r0[i]) & (_p0 <= r0[i + 1]) & (_p1 > r1[j]) & (
													_p1 <= r1[j + 1])
												# n[i, j] = np.nansum(f_right_n[coolh])
												nsubs[i, j] = np.nansum(coolh)  # /n[i, j]
												hfright[i, j] = np.nanmean(f_right[coolh])
												hfrand[i, j] = np.nanmean(f_rand[coolh])  # / np.nansum(f_rand_n[coolh])

										all[coord, pn0, pn1] = [hfright, hfrand, nsubs]

										for out, _h in zip(outs, all[coord, pn0, pn1]):
											print 'Saving file', out
											np.savetxt(out, _h)

									if 0:
										def doplot(h, name, labelname, vmin=None, vmax=None):
											# from matplotlib.colors import LogNorm
											if pn0 == 'd': plt.figure(figsize=(7, 5))
											if pn1 == 'd': plt.figure(figsize=(6, 7))
											plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin,
													   vmax=vmax)
											dx = (r0[1] - r0[0])
											dy = (r1[1] - r1[0])
											plt.xticks(2 * np.arange(n0 / 2) + .5,
													   np.arange(r0[0] + dx, r0[-1], 2 * dx))
											plt.yticks(2 * np.arange(n1 / 2) + .5,
													   np.arange(r1[0] + dy, r1[-1], 2 * dy))
											plt.xlabel(plabels[pn0], fontsize=fsize)
											plt.ylabel(plabels[pn1], fontsize=fsize)
											cbar = plt.colorbar()
											cbar.set_label(labelname, size=fsize)
											plt.savefig(
												'../../EAGLE/analysis//plots/coords/hist_%s_snap%d_%s_%s-%s%s.png'
												% (name, snap, coord, pn0, pn1, extraname))
											plt.close()


										# Number of neighbors per bin
										doplot(n, 'N', u'$\mathrm{n}$')  # , vmin=0, vmax=50)

										# sb
										doplot(hfright, 'SB', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=.4, vmax=.7)

										# sb rand
										doplot(hfrand, 'SB_rand', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=.4, vmax=.7)

										# sb - sb rand
										doplot((hfright - hfrand) * 100 / hfrand, 'SBminusSBrand',
											   u'$\mathrm{SB\,\,increase\,[percentage]}$', vmin=0, vmax=20)

										# sb wrt sbpeak
										doplot(hfright / sbpeak, 'SBvsUVB',
											   u'$\mathrm{SB/SB_{UVB}}$', vmin=0, vmax=.5)

						else:
							print 'No subcubes for %s coordinate' % coord

					for p in pnameshist:  # list(product(pnames, pnames)):
						superall['%d' % snap] = all
						ssss += 1
						print 'ssss', ssss


						def doplot(h, name, labelname, vmin=None, vmax=None, title=''):
							# from matplotlib.colors import LogNorm
							if pn0 == 'd': plt.figure(figsize=(7, 5))
							if pn1 == 'd': plt.figure(figsize=(5.5, 6.21))
							plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin, vmax=vmax,
									   cmap=cm)
							dx = (r0[1] - r0[0])
							dy = (r1[1] - r1[0])
							plt.title(title)
							plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
							plt.yticks([-.5, 1.5, 3.5, 5.5, 7.5, 9.5], [0, 4, 8, 12, 16, 20])
							# plt.xticks(2 * np.arange(n0 / 2) + .5, np.arange(r0[0] + dx, r0[-1], 2 * dx))
							# plt.yticks(2 * np.arange(n1 / 2) + .5, np.arange(r1[0] + dy, r1[-1], 2 * dy))
							plt.xlabel(plabels[pn0], fontsize=fsize)
							plt.ylabel(plabels[pn1], fontsize=fsize)
							cbar = plt.colorbar()
							cbar.set_label(labelname, size=fsize)
							plt.savefig(
								'../../EAGLE/analysis/plots/hist_%s_snap%d_%s-%s%s.png'
								% (name, snap, pn0, pn1, extraname))
							plt.close()


						pn0, pn1 = p
						print 'Combining %s with %s for all coords and snap %d' % (pn0, pn1, snap)

						r0, r1 = [pranges[pn0], pranges[pn1]]
						n0, n1 = [len(r0), len(r1)]
						outs = ['../../EAGLE/analysis/histograms/hist_fright_snap%d_%s-%s%s.dat' % (
							snap, pn0, pn1, extraname),
								'../../EAGLE/analysis/histograms/hist_frand_snap%d_%s-%s%s.dat' % (
									snap, pn0, pn1, extraname),
								# '../../EAGLE/analysis/histograms/hist_n_snap%d_%s-%s%s.dat' % (snap, pn0, pn1, extraname),
								'../../EAGLE/analysis/histograms/hist_nsubs_snap%d_%s-%s%s.dat' % (
									snap, pn0, pn1, extraname)]
						isfiles = True
						for out in outs: isfiles &= os.path.isfile(out)

						if isfiles and not overwrite and not histover:
							print 'h2d files exist'
							_all = [np.loadtxt(out) for out in outs]
						else:
							print 'No h2d files!'
							_all = []
							for coord in coordnames: _all.append(all[coord, p[0], p[1]])
							_all = np.nanmean(_all, 0)
							for out, _h in zip(outs, _all):
								print 'Saving file', out
								np.savetxt(out, _h)
							hfright, hfrand, nsubs = _all

							# from matplotlib.colors import LogNorm

							# Number of pixels per bin
							# doplot(np.log10(n), 'N', u'$\mathrm{log(pixels)}$')#, vmin=0, vmax=40)
							vmin, vmax = [.04, .13]  # [None, None]#
							doplot(hfright, 'SBvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=vmin, vmax=vmax,
								   title=u'snapshot %d, z=%.3f' % (snap, red))

							# sb rand
							doplot(hfrand, 'SBrandvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=vmin, vmax=vmax,
								   title=u'snapshot %d, z=%.3f' % (snap, red))

						if 0:
							# Number of neighbors per bin
							doplot(nsubs, 'Nsub', u'$\mathrm{Number\,\,of\,\,orientations}$')  # , vmin=0, vmax=40)
							# sb
							doplot(hfright * sbpeak, 'SB', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=0, vmax=.14)

							doplot(hfright, 'SBvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=0, vmax=.14)

							# sb rand
							doplot(hfrand * sbpeak, 'SB_rand', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=0, vmax=.14)

							# sb - sb rand
							doplot((hfright - hfrand) * 100 / hfrand, 'SBminusSBrand',
								   u'$\mathrm{SB\,\,increase\,[percentage]}$', vmin=-30, vmax=80)

							# SNR
							doplot(hfright * np.sqrt(nsubs), 'SNR', u'$\mathrm{Relative\,\,SNR}$', vmin=0.2, vmax=1.5)

			else:

				houts2 = ['../../EAGLE/analysis/histograms/hist_fright_%s-%s%s.dat' % (pn0, pn1, extraname),
						  '../../EAGLE/analysis/histograms/hist_frand_%s-%s%s.dat' % (pn0, pn1, extraname),
						  '../../EAGLE/analysis/histograms/hist_nsubs_%s-%s%s.dat' % (pn0, pn1, extraname)]
				hbool = True
				for ho in houts2: hbool &= os.path.isfile(ho)

				hright, hrand, nsubs = [[], [], []]

				if hbool:
					print 'aa'
					hists = [hright, hrand, nsubs]
					for ho, h in zip(houts2, hists): h = np.loadtxt(ho)
				else:
					for ho in houts:
						h = np.loadtxt(ho)


	def doplot(h, name, labelname, p0='', p1='', r0=None, r1=None, vmin=None, vmax=None, title='',
			   xlabel='HST_F606W', ylabel='d [cMpc]', cticks=None):

		# from matplotlib.colors import LogNorm
		# if pn0 == 'd': plt.figure(figsize=(7, 5))
		plt.figure(figsize=(5, 6.2))  # if pn1 == 'd':
		plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin, vmax=vmax,
				   cmap=cm)
		if p0 == 'u1':
			plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
		else:
			plt.xticks([0, 1.5, 3, 4.5], [26, 28, 30, 32])
		# plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
		plt.yticks([-.5, 1.5, 3.5, 5.5, 7.5, 9.5], [0, 4, 8, 12, 16, 20])

		if r0 is not None:
			dx = (r0[1] - r0[0])
			plt.xticks(2 * np.arange(n0 / 2) + .5, np.arange(r0[0] + dx, r0[-1], 2 * dx))
		if r1 is not None:
			dy = (r1[1] - r1[0])
			plt.yticks(2 * np.arange(n1 / 2) + .5, np.arange(r1[0] + dy, r1[-1], 2 * dy))
		plt.xlabel(xlabel, fontsize=fsize)
		plt.ylabel(ylabel, fontsize=fsize)
		plt.title(title)
		cbar = plt.colorbar()
		cbar.set_label(labelname, size=fsize)
		if cticks is not None: cbar.set_ticks(cticks)
		plt.savefig(
			'../../EAGLE/analysis/plots/hist_%s_%s-%s%s.png'
			% (name, p0, p1, extraname))
		plt.close()


	# recheck all of this part!!

	if 0:
		print 'Final analysis'
		uvb = None
		huvb = None
		uvbrand = None
		iii = 0
		for k in superall.keys():
			sa = superall[k]
			for kk in sa.keys():
				iii += 1
				_all = sa[kk]
				if huvb is None:
					# n = _all[2]
					nsub = _all[2]
					hfright = _all[0] * sbpeaks[k]  # *nsub
					hfrand = _all[1] * sbpeaks[k]  # *nsub
					huvb = _all[0]  # *nsub
					uvb = np.nanmean(_all[0])
					uvbrand = np.nanmean(_all[1])
					nuvb = _all[2]
					print iii, k, kk, 'oriented %.3f random %.3f' % (uvb, uvbrand)

				else:
					# n += _all[2]
					_nsub = _all[2]
					nsub += _nsub
					_hright = _all[0]  # *_all[2]
					hfright += _hright * sbpeaks[k]
					_hfrand = _all[1]  # *_all[2]
					hfrand += _hfrand * sbpeaks[k]
					huvb += _all[0]
					_uvb = np.nanmean(_hright)
					_uvbrand = np.nanmean(_hfrand)
					uvb += _uvb
					uvbrand += _uvbrand
					nuvb += _nsub
					print iii, k, kk, 'oriented %.4f random %.4f' % (_uvb, _uvbrand)
		ntot = 9.  # np.sum(n)
		# uvb = uvb/ntot/sbpeaks['11']
		# uvbrand = uvbrand/ntot/sbpeaks['11']

		hfright /= ntot
		hfrand /= ntot
		huvb /= ntot

	else:
		urange = [-21.405001, -19.551001, -19.080999, -18.549, -18.066, -17.507999]
		drange = np.arange(0, 21, 2)
		n0 = len(urange)
		n1 = len(drange)
		fsize = 14
		extraname = ''
		hright = {}
		hrand = {}
		huvb = {}
		for snap in snaps:
			ss = str(snap)
			vmin, vmax = [.04, .16]
			_f = getdata('../../EAGLE/cats/lae_pairs_snap%d.fits' % snap, 1)
			fright = np.concatenate([_f['f_right_x'], _f['f_right_y'], _f['f_right_z']])
			fright[np.isnan(fright)] = 0
			frand = np.concatenate([_f['f_rand_x'], _f['f_rand_y'], _f['f_rand_z']])
			frand[np.isnan(frand)] = 0
			u1 = np.concatenate([_f['U1'], _f['U1'], _f['U1']])
			dist = np.concatenate([_f['com_dist'], _f['com_dist'], _f['com_dist']])
			n = np.histogram2d(u1, dist, [urange, drange], weights=fright > 0)[0]
			h = np.histogram2d(u1, dist, [urange, drange], weights=fright)[0]
			huvb[ss] = h / n
			hright[ss] = huvb[ss] * sbpeaks[ss]
			n = np.histogram2d(u1, dist, [urange, drange], weights=frand > 0)[0]
			h = np.histogram2d(u1, dist, [urange, drange], weights=frand)[0]
			hrand[ss] = h * sbpeaks[ss] / n
			doplot(huvb[ss], 'SBvsUVB_snap%d' % snap, u'$\mathrm{SB/SB_{UVB}}$', 'u1', 'd',  # urange, drange,
				   xlabel=u'$\mathrm{U}}$', ylabel=u'$\mathrm{d\,[cMpc]}}$',
				   title='snapshot %d, z=%.3f' % (snap, redshifts[ss]), vmin=vmin, vmax=vmax)
			doplot(hright[ss], 'SB_snap%d' % snap, u'$\mathrm{SB\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', 'u1', 'd',
				   # urange, drange,
				   xlabel=u'$\mathrm{U}}$', ylabel=u'$\mathrm{d\,[cMpc]}}$',
				   title='snapshot %d, z=%.3f' % (snap, redshifts[ss]), vmin=.04, vmax=.26,
				   cticks=[.05, .1, .15, .2, .25])
		hfright = np.nanmean([hright['10'], hright['11'], hright['12']], 0)
		hfrand = np.nanmean([hrand['10'], hrand['11'], hrand['12']], 0)
		huvb = np.nanmean([huvb['10'], huvb['11'], huvb['12']], 0)

	# print 'total SB/SB_UVB: oriented %.4f random %.4f' % (uvb, uvbrand)

	noise1 = 1.  # assuming 1e-20 e.g. UDF -> but 0.44 with the combination of 390 subcubes
	noise2 = noise1 * np.sqrt(3)  # assuming 1e-20 e.g. UDF

	udf = getdata('../../UDF/cats/lae_pairs.fits', 1)
	mosaic = getdata('../../mosaic/cats/lae_pairs-2018.fits', 1)

	d1 = udf['com_dist']
	t1 = udf['theta']
	f1 = udf['HST_F606W']
	d2 = mosaic['com_dist']
	f2 = mosaic['HST_F606W']
	t2 = mosaic['theta']

	g1 = (d1 > .5) & (d1 < 20) & (t1 > 16) & (f1 > 0)
	g2 = (d2 > .5) & (d2 < 20) & (t2 > 16) & (f2 > 0)

	dbins = np.arange(0, 21, 2)
	# fbins = [25.5921, 27.3024, 27.9615, 28.5004, 28.9428, 29.343,
	#         29.724,  30.0098,  30.5226,  31.3483,  32.046]
	fbins = [25.5921, 27.9615, 28.9428, 29.724, 30.5226, 32.046]

	p0 = 'u1'
	p1 = 'd'
	# n0 = len(pranges[p0])
	# n1 = len(pranges[p1])
	r0 = fbins
	r1 = dbins

	h1 = np.histogram2d(f1, d1, [fbins, dbins])[0]
	h2 = np.histogram2d(f2, d2, [fbins, dbins])[0]
	doplot(h1, 'Nsub_UDF', u'$\mathrm{Number\,\,of\,\,orientations}$', 'HST_F606W', p1, r0, r1)  # , vmin=.2, vmax=1.5)
	doplot(h2, 'Nsub_mosaic', u'$\mathrm{Number\,\,of\,\,orientations}$', 'HST_F606W', p1, r0,
		   r1)  # , vmin=.2, vmax=1.5)#, vmin=.2, vmax=.5)

	SNR1 = hfright * np.sqrt(h1) / noise1
	SNRrand1 = hfrand * np.sqrt(h1) / noise1
	SNR2 = hfright * np.sqrt(h2) / noise2
	SNRrand2 = hfrand * np.sqrt(h2) / noise2
	vmin, vmax = [None, None]  # [.3, 1.6]
	doplot(huvb, 'SBvsUVBtot', u'$\mathrm{SB/SB_{UVB}}$', 'u1', p1, vmin=.06, vmax=.14,
		   xlabel=r'$\mathrm{U}$')  # r0, r1,
	doplot(SNR1, 'SNR_UDF', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)
	doplot(SNRrand1, 'SNRrand_UDF', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin,
		   vmax=vmax)  # , vmin=.2, vmax=.5)
	doplot((SNR1 - SNRrand1) * 100 / SNRrand1, 'OrientedvsRandom_UDF', u'$\mathrm{SB\,\,increase\,[percentage]}$',
		   'HST_F606W', p1, r0, r1, vmin=0, vmax=40)  # , vmin=.2, vmax=.5)
	doplot(SNR2, 'SNR_mosaic', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)  # , vmin=.2, vmax=.5)
	doplot(SNRrand2, 'SNRrand_mosaic', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin,
		   vmax=vmax)  # , vmin=.2, vmax=.5)
	doplot((SNR2 - SNRrand2) * 100 / SNRrand2, 'OrientedvsRandom_mosaic', u'$\mathrm{SB\,\,increase\,[percentage]}$',
		   'HST_F606W', p1, r0, r1, vmin=0, vmax=40)  # , vmin=.2, vmax=.5)

	nh1 = np.nansum(h1)
	fconn1 = np.nansum((huvb * h1).reshape(-1)) / nh1
	SNR_tot1 = np.nansum((SNR1 * h1).reshape(-1)) / nh1 / noise1
	SNRrand_tot1 = np.nansum((SNRrand1 * h1).reshape(-1)) / nh1 / noise1

	nh2 = np.nansum(h2)
	fconn2 = np.nansum((huvb * h2).reshape(-1)) / nh2
	SNR_tot2 = np.nansum((SNR2 * h2).reshape(-1)) / nh2 / noise2
	SNRrand_tot2 = np.nansum((SNRrand2 * h2).reshape(-1)) / nh2 / noise2

	print 'Total SNR UDF %.3f rand %.3f' % (SNR_tot1, SNRrand_tot1)
	print 'f_conn UDF %.3f' % (fconn1)
	print 'Total SNR MOSAIC %.3f rand %.3f' % (SNR_tot2, SNRrand_tot2)
	print 'f_conn MOSAIC %.3f' % (fconn2)

if nhiprof:
	s='11'
	H0=69.3
	Ol = 0.714
	Om = 1-0.714
	red = redshifts[s]
	f1 = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/'
	f2 = '/net/abnoba/scratch2/gallegos/Research/MUSE/EAGLE/'
	dat = getdata(f2+'cats/gals_snap11.fits', 1)
	m = dat['DM_mass']
	v = {}
	offset = {}
	for c in coordnames:
		v[c] = dat['v%s' % c]
		offset[c] = v[c]*(1+red)/H0/np.sqrt(Om*(1+red)**3+Ol)

	print 'min %.3e, max %.3e, mean %.3e, median %.3e' % (np.amin(m), np.amax(m), np.nanmean(m), np.nanmax(m))
	id = dat['ID']
	f = []
	for c in coordnames:
		for i in id:
			print i, c
			f.append(getdata(f1+'gals/snap11_%s/%d.NHI.fits' % (c, i)))
	mmin = 3e9
	mmax = 5e11
	lm = np.log10(m)
	lm = np.concatenate([lm, lm, lm], 0)
	r = np.arange(9.5, 12.5, .5)
	fm = []
	f = np.array(f)
	for i in range(len(r) - 1):
		print r[i], r[i + 1]
		fm.append(np.nanmean(f[(lm > r[i]) & (lm < r[i + 1]), :, :, :], 0))
	fm = np.array(fm)

	if 1:
		hdu = PrimaryHDU()
		for i in range(len(r) - 1):
			hdu.data = fm[i]
			hdu.writeto('radprof_logm%.1f-%.1f.fits', clobber=True)
		for i in range(len(r) - 1):
			hdu.data = fm[i]
			hdu.writeto('radprof_logm%.1f-%.1f.fits' % (r[i], r[i + 1]), clobber=True)

	n, zl, yl, xl = fm.shape
	y, x = np.ogrid[0: yl, 0: xl]
	rs = np.array([0, 1, 2, 4, 6, 8, 10, 15, 20, 30, 40])
	nr = len(rs)
	asec2pix = asec2kpcs[s] * (1 + redshifts[s]) * kpc2pix
	rp = rs/asec2pix
	prof = []
	for i in range(nr-1):
		cool = ((x-xl/2)**2+(y-yl/2)**2<rp[i+1]) & ((x-xl/2)**2+(y-yl/2)**2>rp[i])
		prof.append(np.mean(fm[:, zl/2, cool]))
	prof = np.array(prof)
	plt.figure()
	x = (rp[1:]+rp[:-1])/2.
	for i in range(len(r) - 1):
		plt.plot(x, prof[:, i])
	plt.show()

if radprof:
	feagle = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/'
	rads = np.array([[0, 2], [2, 4], [4, 6], [6, 8], [8, 12], [12, 16],
			[16, 20], [20, 24], [24, 30], [30, 40], [40, 60], [60, 80], [80, 100], [100, 150], [150, 200]])
	feagles = []
	l = 601
	hdu = PrimaryHDU()
	yo, xo = np.ogrid[0:lcube, 0:lcube]
	for s in snaps:
		asec2kpc, red = [asec2kpcs[str(s)], redshifts[str(s)]]
		data = getdata('../../EAGLE/cats/gals_snap%d.fits' % s, 1)
		asec2pix = asec2kpc*(1+red)*kpc2pix
		radpix = rads*asec2pix
		rcubes = []
		for c in coordnames:
			cube = getdata('%s/snap%d_%s.NHI.fits' % (feagle, s, c))
			zl, yl, xl = cube.shape
			ids = data['ID']
			cs = ['x', 'y', 'z']
			cs.remove(c)
			xs = np.round(data[cs[0]] * xl / coml).astype(int)
			ys = np.round(data[cs[1]] * yl / coml).astype(int)
			zs = np.round(data[c] * zl / coml).astype(int)

			for i, x, y, z in zip(ids, xs, ys, zs):
				print i, x, y, z
				radname = '%s/gals/snap%d_%s/%d.RADPROF.fits' % (feagle, s, c, i)
				if not os.path.isfile(radname) or overwrite:
					pos = (xo - x) ** 2 + (yo - y) ** 2
					frad = []
					for r in radpix:
						cool = (pos > (r[0])**2) & \
							   (pos < (r[1])**2)
						fcool = np.nanmean(cube[:, cool], 1)
						frad.append(np.roll(fcool, zl/2-z))
					rcube = np.array(frad)
					hdu.data = rcube
					hdu.writeto(radname, clobber=True)
				else: rcube = getdata(radname)
				rcubes.append(rcube)
		hdu.data = np.nanmean(rcubes, 0)
		hdu.writeto('%s/gals/snap%d.RADPROF.fits' % (feagle, s), clobber=True)

if galcov:
	zw0 = 2
	folder = '/net/abnoba/scratch2/gallegos/Research/MUSE/'

	cov = {}
	rads = [[0.1, 2], [2, 4], [4, 8], [8, 12], [12, 16], [16, 20], [20, 30]]#, [30, 60], [60, 100]]
	nr = len(rads)
	for ff in ['HDFS']:#, 'mosaic']:
		pairs = getdata('../../%s/cats/lae_pairs.fits' % ff, 1)
		theta = pairs['theta']
		cool = (abs(pairs['shear']) <= zw0)# & (pairs['redshift1']<=3) & (pairs['redshift2']<=3)
		theta = theta[cool]
		cov[ff] = np.zeros(nr)
		for i in range(nr):
			area = np.pi*(rads[i][1]**2-rads[i][0]**2)
			inside = (theta > rads[i][0]) & (theta < rads[i][1])
			_cov = np.sum(inside)/area
			print 'catalog', ff, 'rads', rads[i], 'fgal %.3f area %.3e arcsec^2' % (_cov, area)
			cov[ff][i] = _cov

	ce = {}
	for coord in coordnames:
		for snap in snaps:
			k = '%d_%s' % (snap, coord)
			ce[k] = np.zeros(nr)
			cs = ['x', 'y', 'z']
			cs.remove(coord)
			pairs = getdata('../../EAGLE/cats/lae_pairs_snap11.fits', 1)
			theta = pairs['theta%s' % coord]
			cool = abs(pairs['shear%s' % coord]) <= zw0
			for i in range(nr):
				area = np.pi*(rads[i][1]**2-rads[i][0]**2)
				inside = (theta > rads[i][0]) & (theta < rads[i][1])
				_cov = np.sum(inside)/area
				print k, 'r', rads[i], 'fgal %.3f area %.3e' % (_cov, area)
				ce[k][i] = _cov
	cov['EAGLE'] = np.nanmean([ce[_ce] for _ce in ce.iterkeys()], 0)
	plt.figure()
	rm = np.mean(rads, 1)
	import itertools
	marker = itertools.cycle(('*', 'X', '^', 'h', 'd'))

	for k in cov.iterkeys():
		plt.plot(rm, cov[k])
		plt.scatter(rm, cov[k], label=k, marker=marker.next(), s=30)
	plt.legend()
	plt.semilogy()
	plt.xlabel(r'distance [arcsec]')
	plt.ylabel(r'$\mathrm{f_{gal}\,[arcsec^{-2}\,\Delta v^{-1}}]$')
	plt.savefig('../../Figures/Gal_covfrac.jpg')

if cubecorr:
	r = 20
	zw0 = 3
	folder = '/net/galaxy-data/export/galaxydata/gallegos/'
	fitcat = ['mosaic', 'UDF', 'HDFS']
	names = ['DATACUBE_UDF-MOSAIC', 'DATACUBE_UDF-10', 'DATACUBE-HDFS-1.35-PROPVAR']
	mnames = ['DATACUBE_UDF-MOSAIC.IM.Objects_Id', 'DATACUBE_UDF-10.IM.Objects_Id',
			  'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id']
	gnames = ['mosaic.galmask_%darcsec' % r, 'UDF.galmask_%darcsec' % r,
			  'HDFS.galmask_%darcsec' % r]
	ext = '.fits'
	do_csub = True
	cut = 20
	dim = (1, 2)
	for ft, name, mname, gname in [zip(fitcat, names, mnames, gnames)[0]]:
		print ft
		spec = open('%sspec.dat' % ft, 'w')
		spec.write('#layer mean_1l mean_%dl\n' % (zw0*2+1))
		fin = '%s%s/%s' % (folder, ft, name)
		if do_csub: fin += '.csub'
		fout = '%s.corr' % fin
		fmask = '%s%s/%s' % (folder, ft, mname)
		fgmask = '%s%s/%s' % (folder, ft, gname)

		print 'get original %scube' % ('csub ' * do_csub)
		cube = getdata(fin+ext)
		cube[cube == -999] = np.nan
		z, y, x = cube.shape
		_f = np.copy(cube)

		print 'mask continuum objects'
		mask = getdata(fmask+ext)
		bad = mask > 0
		_f[:, bad[0]] = np.nan

		print 'mask 3d objects'
		mask = getdata(fmask.replace('.IM', '')+ext)
		bad = mask > 0
		_f[bad] = np.nan

		print 'mask LAE halos'
		gmask = getdata(fgmask + ext)
		#gmask2 = np.nansum(gmask, 0)
		halos = gmask > 0
		_f[halos] = np.nan

		print 'find high bkg layers'

		smooth = False
		print 'Calculate mean bkg for non high layers, not smoothed!!'
		m = np.nanmean(_f[:, cut: y - cut, cut: x - cut], dim)
		for zz in range(z):
			if smooth:
				zmin = max(0, zz-zw0)
				zmax = min(z, zz+zw0)
				m5 = np.nanmean(_f[zmin: zmax+1, cut: y - cut, cut: x - cut])
			else:
				m5 = np.nanmean(_f[zz, cut: y - cut, cut: x - cut])
			print '%d %f %f' % (zz, m[zz], m5)
			spec.write('%d %f %f\n' % (zz, m[zz], m5))
			cube[zz, :, :] -= m5
		mm = np.nanmean(m)
		print 'mean sb value', mm, 'sum of all layers with their non masked pixels', mm*z
		std = np.nanstd(m)
		high = abs(m-mm) > 4*std
		print np.sum(high), 'high bkg layers of', z, 'std', std
		cube[high, :, :] = np.nan

		hdu = PrimaryHDU()
		hdu.data = cube
		hdu.writeto(fout+ext, clobber=True)

if laemask:
	import astropy.io.fits as fits
	from astropy import wcs

	def red2pix(z):
		"""Convert redshift to pixel"""
		lya_rest = 1215.67
		l0 = 4750.
		pixsize = 1.25  # 1 pixel = 1.25 Angstrom
		return (lya_rest * (1 + z) - l0) / pixsize + 1
	fin = '../../'
	fout = '../../'#'/net/galaxy-data/export/galaxydata/gallegos/'
	ext = '.fits'
	fitcat = ['mosaic', 'UDF', 'HDFS', 'MXDF']
	names = ['DATACUBE_UDF-MOSAIC', 'DATACUBE_UDF-10', 'DATACUBE-HDFS-1.35-PROPVAR',
			 'DATACUBE_MXDF_ZAP_COR']
	hnames = ['DATACUBE_UDF-MOSAIC', 'UDF10.z1300', 'DATACUBE-HDFS-1.35-PROPVAR',
			  'DATACUBE_MXDF_ZAP_COR']
	mnames = ['DATACUBE_UDF-MOSAIC.IM.Objects_Id', 'UDF10.z1300.IM.Objects_Id',
			  'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id',
			  'DATACUBE_MXDF_ZAP_COR.IM.Objects_Id']
	gnames = ['%s.galmask_%darcsec' % (fc, rad) for fc in fitcat]
	for i in range(len(fitcat)):
		ft = fitcat[i]
		name = names[i]
		cat = getdata(fin + '%s/cats/laes.fits' % ft, 1)
		cubename = '/net/galaxy-data/export/galaxydata/gallegos/%s/%s%s' % (ft, name, ext)
		_cube = getdata(cubename)
		zl, yl, xl = _cube.shape
		hname = '/net/galaxy-data/export/galaxydata/gallegos/%s/%s%s' % (ft, hnames[i], ext)
		ttt, header_data_cube = getdata(hname, 0, header=True)
		ra = cat['RA']
		dec = cat['DEC']
		try: redshift = cat['z_muse']
		except: redshift = cat['redshift']
		ids = cat['ID']
		# Removing COMMENT key to avoid problems reading non-ascii characters
		cards = header_data_cube.cards
		bad = ['COMMENT' == b[0] for b in cards]
		for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
		hdulist = fits.open(cubename)
		w = wcs.WCS(header_data_cube, hdulist)
		x, y, z = np.round(w.all_world2pix(ra, dec, [1] * len(ra), 1)).astype(int)
		z = np.round(red2pix(redshift)).astype(int)
		zw0 = 5
		_y, _x = np.ogrid[0:yl, 0:xl]
		cube = np.zeros((zl, yl, xl))

		rpix = rad*5
		for i in range(len(z)):
			gal = (x[i]-_x)**2+(y[i]-_y)**2 < rpix**2
			print i, np.sum(gal)
			zmin = max(0, z[i] - zw0)
			zmax = min(zl, z[i] + zw0 + 1)
			cube[zmin: zmax, gal] += ids[i]
		hdu = PrimaryHDU()
		hdu.data = cube
		hdu.writeto(fout+'%s/%s.galmask_%darcsec.fits' % (ft, ft, rad), clobber=True)

if snr:
	def sclipping(fits, nsigma, dim=(0), mask=None):
		# print 'Asign nan to values > nsigma in a fits array'
		for i in range(len(fits)):
			fits[i, :, np.where(mask)[i]] = np.nan
		stds = np.nanstd(fits, dim)
		high_sigma = np.abs(fits) > nsigma * stds
		fits[high_sigma] = np.nan
		return fits

	folder = '/net/abnoba/scratch2/gallegos/Research/MUSE/'
	flaes = '/net/galaxy-data/export/galaxydata/gallegos/'
	fstack = folder + '/all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100/'
	scat = ''
	for i in fitcat: scat += '_' + i
	ncores = 10
	nrand = 10
	zw0 = rdarg(argv, 'zw', int, 2)
	ptype = rdarg(argv, 'ptype', str, '')
	do_laes = 0
	do_eagle_gals, do_eagle, do_eplot = 0, 0, 0
	vcorr = 1
	do_muse, do_rand, do_mplot = 0, 1, 0
	do_spectra = 0
	do_snr = 0
	do_mask = mask
	do_sclip = 0
	fix_red = 0
	zoff = 0
	spectype = [1, 2, 3, 4]#[2]
	sspec = '.spec'
	for _sp in spectype: sspec += '%d' % _sp
	if ptype=='laes':
		do_laes = 1
		do_rand = 0
		do_muse = 1
		do_spectra = 1
		zw0 = 0
		fix_red = 1
	if ptype=='spectra':
		do_muse = 1
		do_spectra = 1
		do_rand = 0
		zw0 = 0

	if ptype=='muse':
		do_muse = 1
		do_mplot = 1
		zw0 = 1
	if ptype=='snr':
		do_eagle = 1
		do_muse = 1
		do_snr = 1
		zw0 = 0
		zoff = 0
	if ptype=='eagle':
		do_eagle =1
		do_eplot = 1
		zw0 = 2
	zw = 2*zw0+1
	SB = 1.27  # SB at z=3.5
	type = 'LLS'  # 'NHI'
	nsigma = 3

	# 1 spaxel is (0.2 arcsec)^2
	# flux2sb comes from 1.25 Angstrom wavelength width times 1/(0.2)**2 from flux2sbersion to the correct aperture
	flux2sb = 31.25

	if do_eagle_gals:
		fall = []
		fstds = []
		hdu = PrimaryHDU()
		for s in snaps:
			cat = getdata('../../EAGLE/cats/gals_snap%d.fits' % s, 1)
			pairs = getdata('../../EAGLE/cats/lae_pairs_snap%d.fits' % s, 1)
			ids = cat['ID']
			us = cat['U']
			mstar = cat['stellar_mass']
			x, y, z = [cat['x'], cat['y'], cat['z']]
			n5d = []
			#fit = np.zeros([lcube, lcube, lcube])
			for i, _x, _y, _z in zip(ids, x, y, z):
				d = np.sqrt(((x - _x + 12.5) % 25 - 12.5) ** 2 + ((y - _y + 12.5) % 25 - 12.5) ** 2 + (
						(z - _z + 12.5) % 25 - 12.5) ** 2)
				dsort = np.argsort(d)
				d.sort()
				n5d.append(d[5])
				cool = pairs['id1'] == i
				id2, x2, y2, z2, cd2 = np.array([pairs['id2'][cool], pairs['x2'][cool], pairs['y2'][cool],
												 pairs['z2'][cool], pairs['com_dist'][cool]])
				#cls = cd2 < 2
				#for cl in cls:
			#		fit[x2[cl], y2[cl], z2[cl]] += id2[cl]

			# n5d = cat['n5d']
			med_u = np.nanmedian(us)
			med_n5d = np.nanmedian(n5d)
			med_mstar = np.nanmedian(mstar)

			feagles = []
			feagle_med = {'n5d_high': [], 'n5d_low': [],  'u_high': [], 'u_low': []}
			try:
				glob.os.makedirs('../../EAGLE/simplestacks/')
			except:
				pass
			if 1:

				for c in coordnames:
					for i, u, m in zip(ids, us, n5d):
						ffits = flaes + '/EAGLE/gals/snap%d_%s/%d.%s.fits' % (s, c, i, type)
						print ffits
						if os.path.isfile(ffits):
							fit = getdata(ffits)
							feagles.append(fit)
							if u < med_u: feagle_med['u_low'].append(fit)
							if u > med_u: feagle_med['u_high'].append(fit)
							if m < med_n5d: feagle_med['n5d_low'].append(fit)
							if m > med_n5d: feagle_med['n5d_high'].append(fit)

			if 1:
				feagle = np.nanmean(feagles, 0)
				fstd = np.nanstd(feagles, 0)
				fall.append(feagle)
				fstds.append(fstd)

				hdu.data = feagle
				hdu.writeto('../../EAGLE/simplestacks/snap%d.%s.fits' % (s, type), clobber=True)
				hdu.data = fstd
				hdu.writeto('../../EAGLE/simplestacks/snap%d.STD.%s.fits' % (s, type), clobber=True)

				for k in feagle_med.keys():
					feagle_med[k] = np.nanmean(feagle_med[k], 0)
					hdu.data = feagle_med[k]
					hdu.writeto('../../EAGLE/simplestacks/snap%d_%s.%s.fits' % (s, k, type), clobber=True)
		fcomb = np.nanmean(fall, 0)
		hdu.data = fcomb
		hdu.writeto('../../EAGLE/simplestacks/snaps10-11-12.%s.fits' % type, clobber=True)

	fdat = 'muse-vs-eagle%s.dat' % scat
	_fin = None
	if do_eagle:
		do_emask = mask
		cat = getdata('../../EAGLE/cats/gals_snap11.fits', 1)
		f2 = '/net/abnoba/scratch2/gallegos/Research/MUSE/EAGLE/'
		dat = getdata(f2 + 'cats/gals_snap11.fits', 1)
		m = dat['DM_mass']
		s = '11'
		H0 = 69.3
		Ol = 0.714
		Om = 1 - 0.714
		red = redshifts[s]
		offset = dat['vx'] * (1 + red) / H0 / np.sqrt(Om * (1 + red) ** 3 + Ol)
		offpix = np.round(offset*zlens['11']/coml).astype(int)

		pnames = ['U', 'SFR', 'Temperature', 'n5d', 'Mmean200']
		prop = {}
		pmed = {}
		for p in pnames:
			prop[p] = cat[p]
			pmed[p] = np.nanmedian(prop[p])
		cename = '../../EAGLE/simplestacks/snap11_x%s%s.fits' % ('.masked'*do_emask, '.vcorr'*vcorr)
		if not os.path.isfile(cename) or overwrite:
			stack = []
			mstack = []
			for i, off in zip(cat['ID'], offpix):
				print i, off
				_s = getdata(flaes + '/EAGLE/gals/snap11_x/%d.%s.fits' % (i, type))
				if vcorr: _s = np.roll(_s, off, 0)
				if mask:
					_m = getdata(flaes + '/EAGLE/gals/snap11_x/%d.mask.fits' % i)
					gal = _m > 0
					_s[:, gal] = np.nan
					mstack.append(_m)
				stack.append(_s)
			stack = np.array(stack)
			_stack = np.nanmean(stack, 0)
			hdu = PrimaryHDU()
			hdu.data = _stack
			hdu.writeto(cename, clobber=True)

			if mask:
				mstack = np.array(mstack)
				_mstack = np.nansum(mstack, 0)
				_mstack2 = np.nanmean(mstack>0, 0)
				hdu.data = _mstack
				hdu.writeto(cename.replace('.fits', '.gals.fits'), clobber=True)
				hdu.data = _mstack2
				hdu.writeto(cename.replace('.fits', '.galcov.fits'), clobber=True)
			for p in pnames:
				high = prop[p] > pmed[p]
				shigh = np.nanmean(stack[high], 0)
				slow = np.nanmean(stack[~high], 0)
				hdu.data = shigh
				hdu.writeto(cename.replace('.fits', '.high-%s.fits' % p), clobber=True)
				hdu.data = slow
				hdu.writeto(cename.replace('.fits', '.low-%s.fits' % p), clobber=True)
			#del stack, mstack

		feagle = getdata(cename)
		fhigh, flow = {},{}
		for p in pnames:
			fhigh[p] = getdata(cename.replace('.fits', '.high-%s.fits' % p))
			flow[p] = getdata(cename.replace('.fits', '.low-%s.fits' % p))
		zl, yl, xl = feagle.shape
		zle0 = zl
		if ptype!='snr':
			_fite = []
			zrange_e = [zl/2, zl/2+2, zl/2+4, zl/2+8]#np.arange(zl/2-(zl/zw0*zw-zw0, zl/2+(zl/zw0)*zw, zw)
			#zrange = zrange[(zrange > 0) & (zrange + zw < zl)]
			for zz in zrange_e:
				_fite.append(np.nansum(feagle[zz-zw0: zz+zw0+1], 0))
			feagle = np.array(_fite)
			del _fite
		zle, yle, xle = feagle.shape
		asec2pix_eagle = .1836  # average btw snapshots
	asec2pix_muse = .2  # muse
	fall = []
	fouts = []
	n5d = []
	try:
		glob.os.makedirs('../../all/simplestacks/')
	except:
		pass

	wmax = 201/2
	fall = []
	mall = []
	fall2 = []
	frand = {}
	mrand = {}

	if do_muse:

		if overwrite:
			for nr in range(nrand): frand['%d' % nr], mrand['%d' % nr] = [[], []]
			for _fitcat in fitcat:
				cat = getdata('../../%s/cats/laes.fits' % _fitcat, 1)
				if _fitcat == 'UDF':
					fwhm = cat['LYALPHA_FWHM_OBS']

				ids = cat['ID']
				zs = cat['redshift']
				#ls = np.round((1215.67*(1+zs)-4750.)/1.25).astype(int)
				ls = 1215.67 * (1 + zs)
				sconf = cat['sconf']
				spec = cat['spec']

				scool = False
				for _sp in spectype: scool |= spec == _sp
				cool = (sconf >= 2) & (zs < 4) & scool
				fitin = []
				corrupt = 0
				for i, l in zip(ids[cool], ls):
					ffits = flaes + '/%s/LAEs/%d.fits' % (_fitcat, i)
					if os.path.isfile(ffits):
						if do_rand:
							for nr in range(nrand):
								try:
									fr = getdata(ffits.replace('.fits', '.rand%d.fits' %nr))
									mr = getdata(ffits.replace('.fits', '.rand.mask%d.fits' %nr))
									zr, yr, xr = fr.shape
									if ptype=='snr':
										_fr = fr[:, yr/2-wmax: yr/2+wmax+1, xr/2-wmax: xr/2+wmax+1]
										frand['%d' % nr].append(_fr)
									else:
										_fr = fr[zr/2-zw0: zr/2+zw0+1, yr/2-wmax: yr/2+wmax+1, xr/2-wmax: xr/2+wmax+1]
										frand['%d' % nr].append(np.nansum(np.nanmean(_fr, (1,2))))
									_mr = mr[yr/2-wmax: yr/2+wmax+1, xr/2-wmax: xr/2+wmax+1]
									mrand['%d' % nr].append(_mr)
								except:
									#frand['%d' % nr].append(frand['0'])
									#mrand['%d' % nr].append(mrand['0'])
									print 'Missing rand?', nr
									pass
						try:
							fit = getdata(ffits)
							#fall2.append(fit)
							zlm, ylm, xlm = fit.shape
							zlm0 = zlm
							off = 0

							_fit = []
							if zw0 == 0:
								_fit = fit[:, ylm/2-wmax: ylm/2+wmax+1, xlm/2-wmax: xlm/2+wmax+1]
								zrange = np.arange(zlm)
							else:
								#zrange = [zlm/2, zlm/2+zw, zlm/2-zw, zlm/2+2*zw, zlm/2-2*zw, zlm/2+3*zw, zlm/2-3*zw]#np.arange(zlm/2-(zlm/zw0)*zw-zw0, zlm/2+(zlm/zw0)*zw, zw)
								zrange = [zlm/2+zoff, zlm/2+3, zlm/2-3, zlm/2+7, zlm/2+10, zlm/2+15, zlm/2+20]#np.arange(zlm/2-(zlm/zw0)*zw-zw0, zlm/2+(zlm/zw0)*zw, zw)
								if do_spectra:
									zrange = np.arange(zlm / 2 - (zlm / zw0) * zw - zw0, zlm / 2 + (zlm / zw0) * zw, zw0)
									zrange = zrange[(zrange-zw0>0) & (zrange+zw0+1<zlm)]
								for zz in zrange:
									fcool = fit[zz-zw0: zz+zw0+1, ylm/2-wmax: ylm/2+wmax+1, xlm/2-wmax: xlm/2+wmax+1]
									_fit.append(fcool)
							_fit = np.array(_fit)
							fall.append(_fit)
							if do_mask:
								mit = getdata(ffits.replace('.fits', '.mask.fits'))
								_mit = mit[ylm/2-wmax: ylm/2+wmax+1, xlm/2-wmax: xlm/2+wmax+1]
								mall.append(_mit)

							if do_laes:
								ym, xm = np.ogrid[0: ylm, 0: xlm]
								rgal = 1/asec2pix_muse
								rcgm = np.array([2, 5])/asec2pix_muse
								cm = (xm - xlm / 2.) ** 2 + (ym - ylm / 2.) ** 2
								gal = (cm < rgal ** 2)
								cgm = (cm >= rcgm[0] ** 2) & (cm < rcgm[1] ** 2)
								if mask: _fit[:, np.array(_mit) > 0] = np.nan
								fgal = np.nanmedian(_fit[:, gal], 1)
								fcgm = np.nanmedian(_fit[:, cgm], 1)
								fnan = np.where(np.isnan(fgal))[0]
								nnan = len(fnan)

								fgal[fnan] = 0
								fcgm[fnan] = 0
								galstd = np.nanstd(fgal)
								cgmstd = np.nanstd(fcgm)
								galmax = np.nanmax(fgal[zlm/2-6: zlm/2+6])
								cgmmax = np.nanmax(fcgm[zlm/2-6: zlm/2+6])
								xmax = np.where(fgal == galmax)[0]-zlm/2
								xmax2 = np.where(fcgm == cgmmax)[0]-zlm/2
								print i, 'fgal %.3f std %.3f' % (galmax, galstd), (xmax+1)*68.5
								print 'fcgm %.3f std %.3f' % (cgmmax, cgmstd), (xmax2+1)*68.5
								fig, ax = plt.subplots(figsize=(10, 5))
								xz = (zrange-zlm0/2)
								xlim = 22
								_xz = xz*68.5
								_xlim = xlim*68.5
								plt.plot(_xz, [0]*len(xz), color='gray')
								plt.xlim(-_xlim, _xlim)
								plt.grid()

								plt.xlabel('v [km/s]')
								plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg/s/cm^2/arcsec^2]}}$')
								ax2 = ax.twiny()
								plt.xlabel(r'wavelength [$\mathrm{\AA}$]')
								_xz = xz*1.25+l
								_xlim = xlim*1.25

								plt.plot(_xz, fcgm * flux2sb, label='CGM')
								plt.scatter(_xz, fcgm*flux2sb)
								#plt.scatter(xmax, galmax*flux2sb*1.1, label=r'Ly$\alpha$ peak', marker='v')
								for ff in fnan:
									plt.axvline(_xz[ff], color='red', linewidth=10, alpha=.2)
								plt.xlim(l-_xlim, l+_xlim)
								plt.legend()
								plt.savefig('../../%s/plots/%d_CGMspec.jpg' % (_fitcat, i))
								plt.plot(_xz, fgal * flux2sb, label='Galaxy')
								plt.scatter(_xz, fgal * flux2sb)
								plt.savefig('../../%s/plots/%d_spec.jpg' % (_fitcat, i))
								plt.close()
						except:
							corrupt += 1
							print 'Problem!?', ffits
							pass
			if ptype != 'laes':
				fall = np.array(fall)
				# median or mean with sclipping?
				# fmuse = np.nanmedian(fall, 0)
				if do_mask:
					bad = np.array(mall) > 0
					for i in range(len(fall)):
						fall[i, :, bad[i]] = np.nan
				if do_sclip:
					stds = np.nanstd(fall, 0)
					high_sigma = np.abs(fall) > nsigma * stds
					fall[high_sigma] = np.nan
				# return fits
				# fall = sclipping(fall, nsig, mask=bad)
				fmuse = np.nanmean(fall, 0)

				# fall2 = np.array(fall2)
				# fall2 = sclipping(fall2, 5, 0)[0]
				# fmuse2 = np.nanmean(fall2, 0)
				hdu = PrimaryHDU()
				hdu.data = fmuse
				hdu.writeto('../../all/simplestacks/stack%s%s%s%s.fits' % (scat, '.mask' * do_mask, sspec,
																		   extraname), clobber=True)

				if do_rand:
					for nr in range(nrand):
						frand['%d' % nr] = np.array(frand['%d' % nr])
						if do_mask:
							for i in range(len(frand['%d' % nr])):
								bad = np.array(mrand['%d' % nr][i] > 0)
								if ptype == 'snr':
									frand['%d' % nr][i, :, bad] = np.nan
								else:
									frand['%d' % nr][i, bad] = np.nan
						if do_sclip:
							stds = np.nanstd(frand['%d' % nr], 0)
							high_sigma = np.abs(frand['%d' % nr]) > nsigma * stds
							frand['%d' % nr][high_sigma] = np.nan
						# frand['%d' % nr] = sclipping(np.array(frand['%d' % nr]), nsig, mask=bad)
						fr_mean = np.nanmean(frand['%d' % nr], 0)
						if ptype == 'snr':
							zlr, ylr, xlr = fr_mean.shape
						else:
							ylr, xlr = fr_mean.shape


		else:
			ffff= '/net/abnoba/scratch2/gallegos/Research/MUSE/all/simplestacks/'
			fmuse = getdata(ffff+'stack_UDF.gmask.sclip4.fits')
			fr = []
			for nr in range(nrand):
				_fr = getdata(ffff+'randstack_UDF.gmask.sclip.%d.fits' % nr)
				frand['%d' % nr] = _fr
				fr.append(_fr)
			fr_mean = np.nanmean(fr, 0)

		zlm, ylm, xlm = fmuse.shape
		zlm0 = zlm
		if do_rand: zlr, ylr, xlr = fr_mean.shape

		if ptype == 'snr':
			fout = open('UVB.dat', 'w')
			fout.write('#r0 r1 zw zoff SB SB_std SB_rand fLLS Gamma Gamma_std\n')
			rpix_eagle = np.array(rads) / asec2pix_eagle
			y, x = np.ogrid[0: yle, 0: xle]
			ce = (x - xle / 2.) ** 2 + (y - yle / 2.) ** 2

			ym, xm = np.ogrid[0: ylm, 0: xlm]
			cm = (xm - xlm / 2.) ** 2 + (ym - ylm / 2.) ** 2
			if do_rand:
				yr, xr = np.ogrid[0: ylr, 0: xlr]
				cr = (xr - xlr / 2.) ** 2 + (yr - ylr / 2.) ** 2

			r0 = np.arange(8, 15, 1)
			r1 = np.arange(12, 28, 1)
			zoffs = [0, -1, -2, -3, 1]
			zws = [0, 1, 2]
			g = {}
			gstd = {}
			sbstd = {}
			sb = {}
			sbrand = {}
			_nm = fall.shape[0]
			for i in r0:
				for j in r1:
					if i!=j:
						rm0 = i / asec2pix_muse
						rm1 = j / asec2pix_muse
						re0 = i / asec2pix_eagle
						re1 = j / asec2pix_eagle
						print 'Between %.1f and %.1f arcsec' % (i, j)
						inside_muse = (cm >= rm0 ** 2) & (cm < rm1 ** 2)
						nin = np.sum(inside_muse)
						for _zoff in zoffs:
							for _zw in zws:
								print 'zw %d zoff %d' % (_zw*2+1, _zoff)
								_fmuse = fmuse[zlm/2-_zw+_zoff: zlm/2+_zw+_zoff+1, inside_muse]
								fin = np.nanmean(np.nansum(_fmuse, 0))
								inside_rand = (cr >= rm0 ** 2) & (cr < rm1 ** 2)
								fr = fr_mean[zlr/2-_zw: zlr/2+_zw+1, inside_rand]
								fr = np.nanmean(np.nansum(fr, 0))
								fstd = []
								for nr in range(nrand):
									#_n = frand['%d' % nr].shape[0]
									#corr = np.sqrt(_n)/np.sqrt(_nm)
									fr2 = frand['%d' % nr][:, zlr/2-_zw: zlr/2+_zw+1, inside_rand]
									fstd.append(np.nanmean(np.nansum(fr2, 1)))
								fstd = np.nanstd(fstd)

								inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
								fe = np.nanmean(np.nansum(feagle[zle/2-_zw: zle/2+_zw+1, inside_eagle],0))
								fy = fin*flux2sb/fe
								sb2gamma = .519 # for z=3.5 expected SB is 1.31e-20 for a gamma .684e-12
								gamma = fy*sb2gamma
								gamma_std = fstd * flux2sb / fe * sb2gamma
								fout.write('%.3f %.3f %d %d %.3f %.3f %.3f %.3f %.3f %.3f\n' %
										   (i, j, _zw, _zoff, fin*flux2sb, fstd*flux2sb, fr*flux2sb, fe, gamma, gamma_std))
								if 0:#fin > 2*fstd+fr:
									sb['%.1f-%.1f_%d_%d' % (i, j, _zw, _zoff)] = fin*flux2sb
									sbstd['%.1f-%.1f_%d_%d' % (i, j, _zw, _zoff)] = fstd*flux2sb
									sbrand['%.1f-%.1f_%d_%d' % (i, j, _zw, _zoff)] = fr*flux2sb
								print 'SB %.3f, random %.3f. 2sigma noise %.3f. fLLS %.3f. SB EAGLE HM12 %.3f. SBUVB %.3f. Gamma %.3f std %.3f' % \
										  (fin*flux2sb, fr*flux2sb, 2*fstd*flux2sb, fe, fe*SB, fy, gamma, gamma_std)
			fout.close()
			if 0:
				_g = np.array(g.values())
				_gstd = np.array(gstd.values())
				_k = np.array(g.keys())
				gmin = np.amin(_g)
				gstdmin = np.amin(_gstd)
				min1 = _g==gmin
				min2 = _gstd == gstdmin
				print 'Min value gamma',_k[min1], _g[min1], _gstd[min1]
				print 'Min value gamma std', _k[min2], _g[min2], _gstd[min2]



			if 1:
				# Selected values
				_r0, _r1, _zoff, _zw = [8, 19, -1, 2]
				rm0 = _r0 / asec2pix_muse
				rm1 = _r1 / asec2pix_muse
				re0 = _r0 / asec2pix_eagle
				re1 = _r1 / asec2pix_eagle
				inside_rand = (cr >= rm0 ** 2) & (cr < rm1 ** 2)
				inside_muse = (cm >= rm0 ** 2) & (cm < rm1 ** 2)
				inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
				_fmuse = fmuse[zlm / 2 - _zw + _zoff: zlm / 2 + _zw + _zoff + 1, inside_muse]
				#zrange = [zlm/2, zlm/2+zoff, zlm/2-zoff, zlm/2+20]
				zrange = np.arange(-10, 11)
				_f = []
				_fr = []
				_std = []
				_fe = []
				for zz in zrange:
					fcool = fit[zlm/2+zz-_zw: zlm/2+zz+_zw+1, inside_muse]
					_f.append(np.nanmean(np.nansum(fcool, 0)))
					inside_rand = (cr >= rm0 ** 2) & (cr < rm1 ** 2)
					fr = fr_mean[zlr/2+zz-_zw: zlr/2+zz+_zw+1, inside_rand]
					fr = np.nanmean(fr)
					fstd = []
					for nr in range(nrand):
						_n = frand['%d' % nr].shape[0]
						corr = np.sqrt(_n) / np.sqrt(_nm)
						fstd.append(corr * np.nanmean(frand['%d' % nr][:, inside_rand]))
					_std.append(np.nanstd(fstd))

					inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
					_fe.append(np.nanmean(feagle[0, inside_eagle]))

		else:
			#rads = [[0, .2], [.2, .4], [.4, 1], [1, 2], [2, 4], [4, 6], [6, 8], [8, 12], [12, 16], [16, 20], [20, 24], [24, 30], [40, 50]]#,	[50, 60], [60, 70], [70, 80], [80, 90], [90, 100]]
			rads = [[0, 2], [2, 4], [4, 6], [6, 10], [10, 16], [16, 24], [24, 30],#, [40, 50], [50, 60]]#,
					[60, 70], [70, 80], [80, 90], [90, 100]]
			if do_snr: rads = [[6, 12], [6, 14], [6, 16], [6, 20], [6, 30], [8, 14], [10, 20], [8, 20]]
			#rads = [[0, 1], [4, 12], [8, 12], [8, 16], [10, 20], [14, 20], [18, 30]]
			if do_spectra: rads = [[0, 2], [2, 4], [4, 8], [8, 12], [12, 20]]

			if do_muse: rpix_muse = np.array(rads) / asec2pix_muse
			nrad = len(rads)

			if do_eagle:
				rpix_eagle = np.array(rads) / asec2pix_eagle
				y, x = np.ogrid[0: yle, 0: xle]
				ce = (x - xle / 2.) ** 2 + (y - yle / 2.) ** 2
				feh = {}
				fel = {}
				for p in pnames:
					feh[p] = [[]]*nrad
					fel[p] = [[]]*nrad

			ym, xm = np.ogrid[0: ylm, 0: xlm]
			cm = (xm - xlm / 2.) ** 2 + (ym - ylm / 2.) ** 2
			if do_rand:
				yr, xr = np.ogrid[0: ylr, 0: xlr]
				cr = (xr - xlr / 2.) ** 2 + (yr - ylr / 2.) ** 2
			fconn = []
			fstds = []
			cmap = get_cmap('tab10')
			color = 0

			_fin = []
			_fe = []
			_fstd = []
			_fe_std = []
			_fy = []
			_fr = []

			for k in range(nrad):
				rm = rpix_muse[k]
				if do_eagle: re = rpix_eagle[k]
				r = rads[k]
				print 'Between %.1f and %.1f arcsec' % (r[0], r[1])

				inside_muse = (cm >= rm[0] ** 2) & (cm < rm[1] ** 2)
				nin = np.sum(inside_muse)
				# print 'Number of pixels', nin
				_fmuse = fmuse[:, inside_muse]
				fin = np.nanmean(_fmuse, 1)
				if do_rand:
					inside_rand = (cr >= rm[0] ** 2) & (cr < rm[1] ** 2)
					fr = fr_mean[:, inside_rand]
					fr = np.nansum(fr)
					fstd = []
					for nr in range(nrand): fstd.append(np.nansum(np.nanmean(frand['%d' % nr][:, inside_rand], 1)))
					fstd = np.nanstd(fstd)
				fe, fe_std, fy = None, None, None

				if do_eagle:
					inside_eagle = (ce >= re[0] ** 2) & (ce < re[1] ** 2)
					tfe = feagle[:, inside_eagle]
					fe = np.nanmean(tfe, 1)
					# fjack = []
					# for i in range(len(feagles)): fjack.append(np.nanmean(np.delete(feagles, i, 0)[:, inside_eagle]))
					# del feagles
					fe_std = np.nanstd(tfe, 1)
					for p in pnames:
						feh[p][k] = np.nanmean(fhigh[p][:, inside_eagle], 1)
						fel[p][k] = np.nanmean(flow[p][:, inside_eagle], 1)

				if do_eagle and do_muse and do_rand:
					fy = (fin[0]-fr)*flux2sb/fe[0]
					gamma_std = fstd*flux2sb/fe[0]*1.9/4.27
					print 'SB %.3f, random %.3f. 2sigma noise %.3f. fLLS %.3f. SB EAGLE HM12 %.3f. SBUVB %.3f. Gamma %.3f std %.3f' % \
						  (fy*fe[0], fr*flux2sb, 2*fstd*flux2sb, fe[0], fe[0]*SB, fy, 1.9*fy/4.27, gamma_std)
					if color == 0:
						lbs = ['Measured', 'Expected', 'Upper limit']
					else:
						lbs = [None, None, None]
				_fin.append(fin)
				_fe.append(fe)
				_fe_std.append(fe_std)
				_fy.append(fy)
				if do_rand:
					_fstd.append(fstd)
					_fr.append(fr)
			_fin = np.array(_fin)
			_fe = np.array(_fe)
			_fe_std = np.array(_fe_std)
			_fy = np.array(_fy)
			if do_rand:
				_fr = np.array(_fr)
				_fstd = np.array(_fstd)

			if do_spectra:
				nrad = len(rads)
				nr = range(nrad)

				zmin, zmax = 0, zlm + 1#zlm/2-10, zlm/2+11
				v = 2.99792458e5 * (np.arange(zlm)*zlm0/zlm-zlm0/2) * 1.25 / 1215.67 / 4.5
				for i in nr:
					fig, ax = plt.subplots(figsize=(10, 4))
					r = rads[i]

					f = _fin[i][zmin: zmax]# - _fr[i]
					if do_rand: fstd = _fstd[i]
					ax.plot(v[zmin: zmax], f*flux2sb)
					#ax.plot(v[zmin: zmax], [2*fstd*flux2sb]*len(f), label=r'2$\sigma$ noise level')
					ax.plot(v[zmin: zmax], [0]*len(f), color='gray')
					#ax.scatter([0], [1], color='red')
					plt.xlabel('v[km/s]')
					plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg/s/cm^2/arcsec^2]}}$')
					plt.grid()

					#ax2 = ax.twinx()
					#ax2.plot(v[zmin: zmax], f/fstd, alpha=0)
					#plt.ylabel(r'SNR')

					plt.savefig('../../Figures/SB_spectra%s_r%d-%d_zw%d%s%s%s.png' % (scat, r[0], r[1], zw, '_mask'*do_mask, sspec, extraname))
					# plt.show()
					plt.close()

				fig, ax = plt.subplots(figsize=(10, 4))
				ax.plot(v[zmin: zmax], [0] * zlm, color='gray')
				plt.ylim(-zw*.35, zw*1.6)
				plt.xlabel('v[km/s]')
				plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg/s/cm^2/arcsec^2]}}$')
				for i in nr[1:]:
					r = rads[i]
					f = _fin[i][zmin: zmax]# - _fr[i]
					if do_rand: fstd = _fstd[i]
					ax.plot(v[zmin: zmax], f*flux2sb, label='%d < r < %d arcsec' % (r[0], r[1]))
				plt.legend()
				plt.grid()
				plt.savefig('../../Figures/SB_spectra%s_zw%d%s%s%s.png' % (scat, zw, '_mask'*do_mask, sspec, extraname))
				plt.rcParams["figure.figsize"] = (8, 4)
				plt.xlim(-1200, 1200)
				plt.savefig('../../Figures/SB_spectra%s_zw%d%s%s%s_2.png' % (scat, zw, '_mask' * do_mask, sspec, extraname))
				# plt.show()
				plt.close()
			import itertools

			marker = itertools.cycle(('*', 'X', '^', 'h', 'd', 'o', 'v'))
			scale = 1  # e-20
			offset = 0
			rm = np.mean(rads, 1)
			if do_mplot:
				fig, ax = plt.subplots(figsize=(7, 4))

				for i in range(len(zrange))[:3]:
					_y = (_fin[:, i]) * flux2sb * scale + offset
					vmin = 2.99792458e5 * (zrange[i]-zlm0/2-zw0) * 1.25 / 1215.67 / 4.5
					vmax = 2.99792458e5 * (zrange[i]-zlm0/2+zw0) * 1.25 / 1215.67 / 4.5
					ax.plot(rm, _y, c=cmap(i))
					ax.scatter(rm, _y, c=cmap(i), marker=marker.next(), s=30, label=r'$%.f<v<%.f$ km/s' % (vmin, vmax))

				_y = (2*_fstd+_fr)*flux2sb
				_ymin = list(_fstd - 10)
				ax.plot(rm, _y, c='black')
				ax.scatter(rm, _y, c='black', marker='.', s=15, label=r'2$\sigma$ noise level')
				plt.fill_between(rm, _ymin, list(_y), facecolor='black', alpha=0.3, lw=0, edgecolor='none')

				plt.xlim(2, 50)
				#plt.ylim(-.2, 2)
				plt.xlabel('distance [arcsec]')
				plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg/s/cm^2/arcsec^2]}}$')
				plt.legend()
				plt.savefig('../../Figures/SB_MUSE_layers%s_zw%d%s%s.pdf' % (scat, zw, '_mask'*do_mask, extraname), ext='pdf', pdi=200)
				plt.close()

			if do_eplot:
				fig, ax = plt.subplots(figsize=(7, 4))
				#a = [list(f) for f in _fe]
				#_fe = np.array(a).T
				#_fe = np.nanmean([_fe, _fe[::-1]], 0)
				_y = _fe*scale
				#a = [list(f) for f in _fe_std]
				#_fe_std = np.array(a).T
				#_fe_std = np.nanmean([_fe_std, _fe_std[::-1]], 0)
				_ystdmin = (_fe[:, 0]-2*_fe_std[:, 0])*scale
				_ystdmax = (_fe[:, 0]+2*_fe_std[:, 0])*scale

				for i in range(len(zrange_e)):
					v0 = 69.3*np.sqrt((1+redshifts['11'])**3*(1-0.714)+0.714)*coml/float(1+redshifts['11'])/zlens['11']
					v = v0*(zrange_e[i]-zle0/2)
					#if zrange_e[i] == zle0/2:
					#	plt.fill_between(rm, _ystdmin[i], _ystdmax[i], facecolor=cmap(i), alpha=0.3, lw=0, edgecolor='none')

					ax.scatter(rm, _y[:, i], c=cmap(i), marker=marker.next(), s=30, label=r'$\Delta$v=%.0f km/s' % v)
					ax.plot(rm, _y[:, i], c=cmap(i))
				plt.semilogy()
				plt.ylim(.3e-2, 1)
				#plt.yticks()
				plt.xlim(1, 95)
				#plt.xticks([.1, .3, 1, 3, 10, 30, 80])
				plt.xlabel('distance [arcsec]')
				plt.legend()
				plt.semilogy()
				plt.ylabel(r'$\mathrm{f_{LLS}}$')
				yti = np.array([1e-2, 1e-1, 1])
				plt.yticks(yti, [.01, .1, 1])

				#ax2 = ax.twinx()
				#plt.semilogy()
				#plt.ylim(.3e-2, 1)
				#plt.ylabel('dndz')
				#plt.yticks(yti, [2, 20, 200])

				# leg = ax.get_legend()
				# for ls in leg.legendHandles: ls.set_color('black')
				plt.savefig('../../Figures/SB_EAGLE_layers_zw%d%s.pdf' % (zw, '_masked'*do_emask), ext='pdf', pdi=200)
				#plt.show()
				plt.close()


				for p in pnames:
					fig, ax = plt.subplots(figsize=(7, 4))
					a = [list(f) for f in feh[p]]
					_feh = np.array(a).T
					a = [list(f) for f in fel[p]]
					_fel = np.array(a).T

					ax.scatter(rm, _feh[zle0/2], c='red', marker=marker.next(), s=30, label=p+' high')
					ax.plot(rm, _feh[zle0/2], c='red')
					ax.scatter(rm, _fel[zle0/2], c='blue', marker=marker.next(), s=30, label=p+' low')
					ax.plot(rm, _fel[zle0/2], c='blue')
					plt.semilogy()
					plt.ylim(.2e-2, 1)
					# plt.yticks()
					plt.xlim(1, 95)
					plt.xlabel('distance [arcsec]')
					plt.legend()
					plt.semilogy()
					plt.ylabel(r'$\mathrm{f_{LLS}}$')
					yti = np.array([1e-2, 1e-1, 1])
					plt.yticks(yti, [.01, .1, 1])

					# leg = ax.get_legend()
					# for ls in leg.legendHandles: ls.set_color('black')
					plt.savefig('../../Figures/SB_EAGLE_layers_zw%d.%s.pdf' % (zw, p), ext='pdf', pdi=200)
					# plt.show()
					plt.close()

if superstack:

	def sclipping(fits, nsigma, dim=None, mask=None):
		# print 'Asign nan to values > nsigma in a fits array'
		if dim is None:
			stds = np.nanstd(fits[:, mask])
		else:
			stds = np.nanstd(fits[:, mask], dim)
		high_sigma = np.abs(fits) > nsigma * stds
		fits[high_sigma] = np.nan
		return fits, stds


	ds = [[0, 2], [1, 3], [2, 4], [3, 6], [5, 8], [7, 10], [9, 15], [14, 20]]
	spmin = {}
	spmax = {}
	sbmin = {}
	sbmax = {}
	srand = {}
	stdmin = {}
	stdmax = {}
	zw0 = 5
	folder = '../../EAGLE/stacks/'
	f2 = '../..all/stacks/'
	f3 = '../../EAGLE/cats/'
	hdu = PrimaryHDU()
	imtype = 'mean'  # 'median'#
	mask = False
	smask = '.mask' * mask
	stype = ('.' + imtype) * (imtype == 'median')
	for dd in ds:
		dmin = dd[0]
		dmax = dd[1]
		all = []
		ncores = dmax
		for snap in snaps:
			red = redshifts[str(snap)]
			asec2kpc = asec2kpcs[str(snap)]
			thetamin = 16
			# height 2 arcsec, width 6 to 12 arcsec -> same as my paper!
			xasecmin = 6
			xasecmax = 12
			yasecmin = -1
			yasecmax = 1
			bin = 2.  # binning after creating subcubes
			asec2pix = asec2kpc * (1 + red) * kpc2pix / bin
			sb2flux = asec2pix ** -2.
			# limits in pixels
			xmin = int(xasecmin * asec2pix)
			xmax = int(xasecmax * asec2pix)
			ymin = int(yasecmin * asec2pix)
			ymax = int(yasecmax * asec2pix)
			# limits in pixels
			# xmin =

			rad_asec = np.array([rads[i + 1] for i in range(len(rads) - 1)]) / asec2pix
			data = getdata(f3 + 'lae_pairs_snap%d.fits' % snap, 1)
			dist = data['com_dist']

			props = ['U1', 'n5d1']
			pr = props[1]
			doall = True
			prop = data[pr]
			cool = (dist > dmin) & (dist < dmax)
			pmedian = abs(np.median(prop[cool]))
			print '%s median' % pr, pmedian
			if pr == props[0]:
				smax = 'd%s-%d_pd16-2000_d5th0.5-20.0_u0.0-%.1f' % (dmin, dmax, pmedian)
				smin = 'd%s-%d_pd16-2000_d5th0.5-20.0_u%.1f-99.0' % (dmin, dmax, pmedian)
				pmin = 'umin'
				pmax = 'umax'
				sprop = 'U'
			if pr == props[1]:
				smax = 'd%s-%d_pd16-2000_d5th0.5-%.1f_u0.0-99.0' % (dmin, dmax, pmedian)
				smin = 'd%s-%d_pd16-2000_d5th%.1f-20.0_u0.0-99.0' % (dmin, dmax, pmedian)
				pmin = 'd5min'
				pmax = 'd5max'
				sprop = 'd5th'

			stacks_pmin = []
			stacks_pmax = []
			std_pmin = []
			std_pmax = []
			for coord in coordnames:  # ['x']:#
				# os.system(
				#		'mpirun -n %d python stackEAGLE.py -snap %d -coord %s -overwrite True -imtype %s -mask %s'
				#		% (ncores, snap, coord, imtype, mask))
				print 'd %d-%d snap %d coord %s' % (dmin, dmax, snap, coord)
				sname = folder + 'snap%d_%s_%s%s%s/stack.fits' % (
					snap, coord, smax, stype, smask)
				if not os.path.isfile(sname) or overwrite:
					os.system(
						'mpirun -n %d python stackEAGLE.py -snap %d -coord %s -dmin %d -dmax %d -%s %f -overwrite True -imtype %s -mask %s'
						% (ncores, snap, coord, dmin, dmax, pmax, pmedian, imtype, mask))
				_stack = getdata(sname)
				zl, yl, xl = _stack.shape
				stacks_pmax.append(_stack[:, :, :yl])
				if doall:
					_std = getdata(sname.replace('.fits', '.STD.fits'))
					std_pmax.append(_std[:, :, :yl])
				sname = folder + 'snap%d_%s_%s%s%s/stack.fits' % (
					snap, coord, smin, stype, smask)
				if not os.path.isfile(sname) or overwrite:
					os.system(
						'mpirun -n %d python stackEAGLE.py -snap %d -coord %s -dmin %d -dmax %d -%s %f -overwrite True -imtype %s -mask %s' %
						(ncores, snap, coord, dmin, dmax, pmin, pmedian, imtype, mask))

				_stack = getdata(sname)
				zl, yl, xl = _stack.shape
				stacks_pmin.append(_stack[:, :, :yl])
				if doall:
					_std = getdata(sname.replace('.fits', '.STD.fits'))
					std_pmin.append(_std[:, :, :yl])
			stack_pmin = np.nanmean(stacks_pmin, 0)
			stack_pmax = np.nanmean(stacks_pmax, 0)
			if doall:
				std_pmin = np.nanmean(std_pmin, 0)
				std_pmax = np.nanmean(std_pmax, 0)
			zl, yl, xl = stack_pmin.shape
			hdu.data = stack_pmax
			fs = folder + 'snap%d_d%d-%d_pd16-2000_d5th0.5-4.0_u0.0-%.1f%s%s/' % (
				snap, dmin, dmax, pmedian, stype, smask)
			if not os.path.isdir(fs): glob.os.makedirs(fs)
			if not os.path.isfile(fs + 'stack.IM.fits'):
				hdu.writeto(fs + 'stack.fits', clobber=True)
				stackim = np.nanmean(stack_pmax[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
				hdu.data = stackim
				hdu.writeto(fs + 'stack.IM.fits', clobber=True)
			hdu.data = stack_pmin
			fs = folder + 'snap%d_d%d-%d_pd16-2000_d5th0.5-4.0_u%.1f-99.0%s%s/' % (
				snap, dmin, dmax, pmedian, stype, smask)
			if not os.path.isdir(fs): glob.os.makedirs(fs)
			if not os.path.isfile(fs + 'stack.IM.fits'):
				hdu.writeto(fs + 'stack.fits', clobber=True)
				stackim = np.nanmean(stack_pmin[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
				hdu.data = stackim
				hdu.writeto(fs + 'stack.IM.fits', clobber=True)

			x0, x1 = xl / 2 + xmin, xl / 2 + xmax + 1
			_x0, _x1 = xl / 2 - xmax, xl / 2 - xmin + 1
			y0, y1 = yl / 2 + ymin, yl / 2 + ymax + 1
			_y0, _y1 = yl / 2 - ymax, yl / 2 - ymin + 1
			z0, z1 = zl / 2 - zw0, zl / 2 + zw0 + 1
			spnhi = np.power(10, sb2nhi(stack_pmin[z0: z1, y0: y1, x0: x1]))
			spmin[('%d' % snap, '%d' % dmax)] = np.nanmean(nhi2sb(np.log10(np.nansum(spnhi, 0))))
			spnhi = np.power(10, sb2nhi(stack_pmax[z0: z1, y0: y1, x0: x1]))
			spmax[('%d' % snap, '%d' % dmax)] = np.nanmean(nhi2sb(np.log10(np.nansum(spnhi, 0))))
			if doall:
				# npix to properly propagate the std
				npix = 1. / np.sqrt(3 * (xmax - xmin) * (ymax - ymin) * zw0)

				stdmin[('%d' % snap, '%d' % dmax)] = np.nanmean(std_pmin[z0: z1, y0: y1, x0: x1]) * npix
				stdmax[('%d' % snap, '%d' % dmax)] = np.nanmean(std_pmax[z0: z1, y0: y1, x0: x1]) * npix

				# srand, random means average over other directions (up, left, bottom)
				sr = np.nanmean(stack_pmin[z0: z1, y0: y1, _x0: _x1])
				sr += np.nanmean(stack_pmin[z0: z1, _y0: _y1, x0: x1])
				sr += np.nanmean(stack_pmin[z0: z1, _y0: _y1, _x0: _x1])
				sr += np.nanmean(stack_pmax[z0: z1, y0: y1, _x0: _x1])
				sr += np.nanmean(stack_pmax[z0: z1, _y0: _y1, x0: x1])
				sr += np.nanmean(stack_pmax[z0: z1, _y0: _y1, _x0: _x1])
				srand[('%d' % snap, '%d' % dmax)] = sr * (2 * zw0 + 1) / 6.

			if 0:
				all.append(stack * sbpeaks[str(snap)])
				stack = np.nanmean(all, 0)
				zl, yl, xl = stack.shape
				hdu.data = stack
				fs = folder + 'd%d-%d_pd16-2000_nw0-1000%s%s/' % (d, d2, stype, smask)
				if not os.path.isdir(fs): glob.os.makedirs(fs)
				hdu.writeto(fs + '/stack.fits', clobber=True)
				stackim = np.nansum(stack[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
				hdu.data = stackim
				hdu.writeto(fs + 'stack.IM.fits', clobber=True)
		if 0:

			fn = f2 + 'd0-%d_los0-20_pd16-2000_z2.9-4.0_l0-2000_n1-1000_v0-10000/' % d
			nname = fn + 'stack.NPIX.IM.fits'
			if not os.path.isfile(nname): os.system('python stack.py -dmax %d -overwrite True' % d)
			npix = getdata(nname)
			yl, xl = npix.shape
			noise = 10.  # assuming 1e-19 cgs noise
			snr = stackim[:, :xl] * np.sqrt(npix) / noise / (
				2 * zw0 + 1.)  # for the image stack I am doing a sum therefore the noise does not decrease, so I added the factor 1/(2*zw0+1.)
			hdu.data = snr
			hdu.writeto(fs + 'stack.SNR.fits', clobber=True)

	types = ['o', '^', 's', 'D', 'p']
	colors = ['red', 'green', 'blue', 'orange', 'magenta']
	SB = [.73, 1.27, 2.49, 5.14, 6.97]
	fLLS = [0.01329666667, 0.008964333333, 0.005651, 0.004793716667, 0.003494351]
	fconn = [0.01999333333, 0.01287333333, 0.007092333333, 0.004155633333, 0.0030469]
	dpi = 200
	ext = 'pdf'
	zs = [redshifts['%d' % s] for s in snaps]
	plt.figure()
	plt.ylim([0, .021])
	plt.yticks([0.000, .005, .01, .015, .02])
	plt.xlabel(r'$\mathrm{redshift}$')
	# plt.xticks(zs)
	plt.scatter(zs, fLLS, color='red', label=r'$\mathrm{f_{LLS}}$')
	plt.plot(zs, fLLS, color='red')
	plt.scatter(zs, fconn, color='blue', label=r'$\mathrm{SB/SB_{UVB}}$')
	plt.plot(zs, fconn, color='blue')
	plt.legend()
	plt.savefig('../../EAGLE/plots/fLLS_fconn-vs-redshift.%s' % ext, dpi=dpi, format=ext)
	plt.close()
	sb2cMpc = []

	for i in range(len(snaps)):
		snap = snaps[i]
		t = types[i]
		c = colors[i]
		sb = SB[i]
		z = redshifts['%d' % snap]
		print snap
		fmin = np.array([spmin[('%d' % snap, '%d' % dd[1])] for dd in ds])
		fmax = np.array([spmax[('%d' % snap, '%d' % dd[1])] for dd in ds])
		f = (fmin + fmax) / 2.
		x = np.array([(dd[0] + dd[1]) / 2. for dd in ds])
		f0 = np.zeros(len(fmax))
		xlabel = r'neighbour distance [cMpc]'
		xlim = [1, 17]

		plt.figure(1)
		plt.ylim([-.25, .25])
		plt.xlim(xlim)
		plt.xlabel(xlabel)
		plt.ylabel(r'$\Delta\mathrm{f_{conn,\,%s}}$' % sprop, fontsize=13)
		plt.plot(x, (fmin - fmax) / f, color=c)
		plt.scatter(x, (fmin - fmax) / f, color=c, marker=t, label=r'$z= %.1f$' % z)  # , label='Snap %d' % snap)

		plt.figure(5)
		plt.xlim()
		plt.xlabel(r'$\mathrm{redshift}$')
		plt.ylabel(r'$\mathrm{SB\,[10^{-20}erg/s/cm^2/arcsec^2]}}$', fontsize=13)
		plt.scatter(z, f[0] * sb, color=c, marker=t)
		sb2cMpc.append(f[0] * sb)

		if doall:
			frand = np.array([srand[('%d' % snap, '%d' % dd[1])] for dd in ds])
			smin = np.array([stdmin[('%d' % snap, '%d' % dd[1])] for dd in ds])
			smax = np.array([stdmax[('%d' % snap, '%d' % dd[1])] for dd in ds])
			s = (smin + smax) / 2.
			plt.figure(2)
			# plt.ylim([0, .5])
			# plt.yticks([0, .03, .06, .09, .12, .15, .18, .21])
			plt.xlim(xlim)
			plt.xlabel(xlabel)
			plt.ylabel(r'$\mathrm{f_{conn}}$', fontsize=13)
			plt.plot(x, f, color=c)
			plt.scatter(x, f, color=c, marker=t, label=r'$z=%.1f$' % z)
			plt.fill_between(x, np.max([f0, f - s], 0), f + s, facecolor=c, alpha=0.3, lw=0, edgecolor='none')

			plt.figure(3, figsize=(7, 5))
			plt.ylim([-.07, .17])
			plt.xlim(xlim)
			plt.xlabel(xlabel)
			plt.ylabel(r'$\mathrm{(f_{conn}-f_{rand})/f_{rand}}$', fontsize=13)
			plt.plot(x, (f - frand) / frand, color=c)
			plt.scatter(x, (f - frand) / frand, color=c, marker=t, label=r'$z=%.1f$' % z)
			# plt.fill_between(x, np.max([f0, f-frand-s], 0), f-frand+s, facecolor=c, alpha=0.3, lw=0, edgecolor='none')

			plt.figure(4)
			# plt.ylim([0, .2])
			# plt.yticks([0, .03, .06, .09, .12, .15, .18, .21])
			plt.xlim(xlim)
			plt.xlabel(xlabel)
			plt.ylabel(r'$\mathrm{SB\,[10^{-20}erg/s/cm^2/arcsec^2]}}$', fontsize=13)
			plt.plot(x, f * sb, color=c)
			plt.scatter(x, f * sb, color=c, marker=t, label=r'$z=%.1f$' % z)
		# plt.fill_between(x, np.max([f0, (f-s)*sb], 0), (f+s)*sb, facecolor=c, alpha=0.3, lw=0, edgecolor='none')

	plt.figure(1)
	plt.legend()
	plt.savefig('../../EAGLE/plots/%s_diff-vs-d.%s' % (sprop, ext), dpi=dpi, format=ext)
	plt.close()
	plt.figure(5)
	plt.plot(zs, sb2cMpc, color='black', zorder=-1)
	plt.savefig('../../EAGLE/plots/SB-vs-z.%s' % ext, dpi=dpi, format=ext)
	plt.close()

	if doall:
		plt.figure(2)
		plt.legend()
		plt.savefig('../../EAGLE/plots/fconn-vs-d.%s' % ext, dpi=dpi, format=ext)
		plt.close()
		plt.figure(3)
		plt.legend()
		plt.savefig('../../EAGLE/plots/fconn-frand.%s' % ext, dpi=dpi, format=ext)
		plt.close()
		plt.figure(4)
		plt.legend()
		plt.savefig('../../EAGLE/plots/SB-vs-d.%s' % ext, dpi=dpi, format=ext)
		plt.close()

if sbprof:
	zw = 5
	yw = 2
	zoff = -1
	dw = 1
	c = (2 * zw + 1) * 7.8125
	asec2kpc = 7.47
	folder = '../../EAGLE/stacks/snap11_x_d0-20_pd16-2000_nw0.000-1.000/'
	folder2 = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100/'
	stack = getdata(folder + 'stack.fits')
	stack2 = getdata(folder2 + 'stack.fits')
	srand = getdata(folder2 + 'random_stack.fits')
	frs = []
	for i in range(200):
		hdur = getdata(folder2 + 'randoms/random_stack.%d.fits' % i)
		frs.append(hdur)
	zl, yl, xl = stack.shape
	zl2, yl2, xl2 = stack2.shape
	print xl, xl2, yl, yl2, zl, zl2
	xmin = yl / 2
	t = np.concatenate([np.arange(xmin, xmin + 12, 1), [xmin + 13, xmin + 18, xmin + 27, xmin + 41]]).astype('int')
	zi = zl / 2 - zw / 2
	zf = zl / 2 + zw / 2 + 1
	zi2 = zl2 / 2 - zw / 2 + zoff
	zf2 = zl2 / 2 + zw / 2 + 1 + zoff
	if 1:
		dir = ['right', 'top', 'left', 'bottom']
		yymin = yl / 2 - yw / 2
		yymax = yl / 2 + yw / 2 + 1
		sb = []

		for i in range(len(t) - 1):
			yis = [yymin, t[i] - 1 - dw, yymin, yl - t[i + 1] - dw]
			yfs = [yymax, t[i + 1] + dw, yymax, yl - t[i] + 1 + dw]
			xis = [t[i] - 1 - dw, yymin, yl - t[i + 1] - dw, yymin]
			xfs = [t[i + 1] + dw, yymax, yl - t[i] + 1 + dw, yymax]
			sb.append((t[i] + t[i + 1]) * .5)
			# for yi, yf, xi, xf, ds in zip(yis, yfs, xis, xfs, dir):
			yi, yf, xi, xf, ds = zip(yis, yfs, xis, xfs, dir)[0]
			if 1:
				sb.append(np.nanmean(stack[zi:zf, yi:yf, xi:xf]))
				sb.append(np.nanmean(stack2[zi2:zf2, yi:yf, xi:xf]) * c)
				sb.append(np.nanmean(srand[zi2:zf2, yi:yf, xi:xf]) * c)
				sb.append(np.nanstd([np.nanmean(ff[zi:zf, yi:yf, xi:xf]) * c for ff in frs]))
	lx = len(t) - 1
	# ld = len(dir)*3 + 1
	ld = 5
	sb = np.array(sb).reshape((lx, ld))
	np.savetxt('../../EAGLE/plots/fconn.txt', sb,
			   header='dpix sim muse rand rand_std')

	x1 = 2
	x2 = 9
	y1 = 0.0001
	y2 = 12
	hm = False
	hmSB = 1.14
	y1b = -4
	y2b = 6
	fsize = 30
	width = 4
	sigmas = 2

	plt.figure(figsize=(20, 12))
	plt.loglog()
	ax = plt.axes()
	cax = plt.gca()
	plt.xlim([x1, x2])
	xrange = np.arange(2, 9, 2)
	plt.xticks(xrange, xrange, fontsize=fsize)
	plt.xlabel(r"$\theta$ [arcsec]", fontsize=fsize)
	# plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize + 2)
	yticks = [0, 2, 4, 6, 8, 10]
	plt.yticks(yticks, yticks, fontsize=fsize)
	plt.minorticks_off()
	plt.twiny(ax=None)
	plt.loglog()
	plt.xlim([x1, x2])
	plt.xlabel(r'$r_p\,\rm{[pkpc]}$', fontsize=fsize)
	kpc = np.arange(15, 90, 15)
	plt.xticks(kpc / asec2kpc, kpc.astype('int'), fontsize=fsize)
	plt.ylim([y1, y2])
	plt.plot((2, x2), (0, 0), '--', label='', color="gray", linewidth=width + 1)
	if hm: plt.plot((x1, x2), (hmSB, hmSB), '--', label=r"LLS Fluorescence from HM12", color="green",
					linewidth=width + 1)

	x = (t[:-1] + (t[1:] - t[:-1]) / 2. - yl / 2) * .4  # 0.4 arcsec per pixel

	# plt.plot(x, sb[:, 2], label='MUSE oriented', lw=width, color='black')
	# plt.scatter(x, sb[:, 2], s=width * 10, color='black')
	sb = np.array(sb)
	sb[sb < 0] = 1e-10
	plt.plot(x, sb[:, 3], label='SB MUSE', lw=width, color='dodgerblue')
	plt.scatter(x, sb[:, 3], s=width * 10, color='dodgerblue')
	ymin = sb[:, 3] - sigmas * sb[:, 4]
	ymin[ymin < 0] = 1e-10
	plt.fill_between(x, ymin, sb[:, 3] + sigmas * sb[:, 4], facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')

	sbuvmean = 0.8
	sbuvup = 1.2
	sbuvdown = 0.6
	plt.plot(x, sb[:, 1] * sbuvmean, label='SB EAGLE+HM12', lw=width, color='red')
	plt.scatter(x, sb[:, 1] * sbuvmean, s=width * 10, color='red')
	plt.fill_between(x, sb[:, 1] * sbuvdown, sb[:, 1] * sbuvup, facecolor='red', alpha=0.3, lw=0, edgecolor='none')
	plt.minorticks_off()
	gamma = 3.16e13  # ionizing photon rate from the galaxy /1e40
	kpc2cm = 30.86  # kpc2cm /1e20
	zmean = 3.5
	sb_diff = sb[:, 3] - sb[:, 1] * sbuvmean
	E_lya = 1.64e-11
	asec2rad_sqd = 2.35e-11
	fesc = 4 * np.pi * (1 + zmean) ** 4 * (
		x * asec2kpc * kpc2cm) ** 2 * sb_diff * 1e-20 / E_lya / .6 / gamma / asec2rad_sqd
	gamma_std = 3e13
	obs_std = sb[:, 4]
	sim_std = 0.3
	sb_std = np.sqrt(np.power(obs_std, 2) + sim_std ** 2)
	tot_std = np.sqrt((gamma_std / gamma) ** 2 + (sb_std / sb_diff) ** 2) * fesc
	print 'fesccccc', fesc
	print 'sb std', sb_std
	print 'tot std', tot_std
	fesc[fesc < 0] = 1e-10
	plt.plot(x, fesc, label='Escape fraction', lw=width, color='black')
	ymin = fesc - tot_std
	ymin[ymin < 0] = 1e-10
	plt.fill_between(x, ymin, fesc + tot_std, facecolor='gray', alpha=0.3, lw=0, edgecolor='none')

	plt.legend(fontsize=fsize, loc='best')  # (3.5,2))

	plt.savefig('../../analysis/muse-vs-eagle.png')

if unique:
	# I'll use pixelssss
	yw = 5
	yl = 141
	n = 10000
	ngal = 303
	nsub = 2016
	nfrac = nsub / ngal
	tot = 0

	# in asec
	rads = [2, 4, 6, 10, 20]

	# for z 3
	fconn = [.17, .11, .105, .10, .11]
	# for z 3.5
	fconn = [.14, .10, .9, .85, .9]
	frand = [.89, .96, .99, 1.02, 1]
	noise = np.sqrt(3)
	SB = [.73, 1.27, 2.49, 5.14, 6.97]
	asec2pix = 5
	Nsub_udf = [46, 80, 140, 208, 352]
	Ngal_udf = [30, 48, 61, 65, 75]

	Nsub_mosaic = [375, 1413, 2571, 4120, 7091]
	Ngal_mosaic = [207, 292, 311, 324, 328]

	for d, f, fr, nsubm, ngalm, nsubu, ngalu in zip(rads, fconn, frand, Nsub_mosaic, Ngal_mosaic, Nsub_udf, Ngal_udf):
		tot = 0
		npix = int(2 * np.pi * d)
		nfrac = nsubm / ngalm
		print 'Mosaic d', d, 'nsub', nsubm, 'ngal', ngalm
		print npix, 'pixels at rad', d
		for j in range(n):
			xrand = np.random.randint(0, npix, nfrac)
			tot += len(np.unique(xrand))
		lun = tot / float(n * nfrac)
		print 'less unique', lun
		print 'Efective signal', SB[1] * f * np.sqrt(lun * nsubm) / np.sqrt(3) / 10.
		print 'Random prof signal', SB[1] * f * fr * np.sqrt(lun * nsubm * npix) / np.sqrt(3) / 10.

		tot = 0
		nfrac = nsubu / ngalu
		print 'UDF d', d, 'nsub', nsubu, 'ngal', ngalu
		print npix, 'pixels at rad', d
		for j in range(n):
			xrand = np.random.randint(0, npix, nfrac)
			tot += len(np.unique(xrand))
		lun = tot / float(n * nfrac)
		print 'less unique', lun
		print 'Efective signal', SB[1] * f * np.sqrt(lun * nsubu) / 10.
		print 'Random prof signal', SB[1] * f * fr * np.sqrt(lun * nsubu * npix) / 10.
