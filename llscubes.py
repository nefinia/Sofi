#!/usr/bin/env python
__author__ = 'gallegos'
import glob
import os
from math import sqrt
from sys import argv

# import eagleSqlTools
import numpy as np
import scipy.interpolate as inter
# from matplotlib.cm import get_cmap
# from matplotlib.colors import LinearSegmentedColormap
from pyfits import getdata, PrimaryHDU

import params
from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
snaps = rdarg(argv, 'snap', list, [6, 7, 8, 9, 10, 11, 12, 13, 14], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
caseA = rdarg(argv, 'caseA', bool, True)
scase = '_CaseA' * caseA
extraname = rdarg(argv, 'extraname', str, '')#'HM12')#
galcov = rdarg(argv, 'galcov', bool, False)
laemask = rdarg(argv, 'laemask', bool, False)
nhiprof = rdarg(argv, 'nhi', bool, False)
minres = rdarg(argv, 'minres', int, 512)
maxres = rdarg(argv, 'maxres', int, 4096)
model = rdarg(argv, 'model', str, 'HM12')#'HM01')#
npref = rdarg(argv, 'npref', int, 12)
rad = rdarg(argv, 'rad', int, 3)  # arcsec
zw0 = rdarg(argv, 'zw', int, 5)  # +- zw0 layers for galmask
_ssthr = rdarg(argv, 'ssthr', float, 6.73e-3)# 1e10 6.73e-3 or None
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
type = rdarg(argv, 'type', str, 'NHI') #NHtot, f_NHI, NHII
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)

gamma_bkg = params.gamma_bkg[model]
heat_bkg = params.heat_bkg[model]
snames = params.snames
redshifts = params.redshifts
asec2kpcs = params.asec2kpcs
ssthr = params.nhSS[model]
if _ssthr is not None:
	ssthr = {}
	for snap in snaps:
		ssthr[snap] = _ssthr

dz = params.dz
sbpeaks = params.sbpeaks
zlens = params.zlens
dndz = params.dndz
fLLScorr = params.fLLScorr


lcube = maxres # 4096
coml = 25  # cMpc
com2pix = 163.84  # lcube/coml
kpc2pix = lcube / float(coml * 1e3)

nhi_fit = params.nhi_fit[model]

coord = 'x'
proj = 1
h = params.h
coml = params.coml
lh = coml * h
cuts = 4
size = lh / float(cuts)
ncores = 10

rads = np.array([0, 2, 4, 8, 12, 20, 30, 50, 100, 200])

lognhi = [10, 11, 12, 13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21, 22, 23, 30]
nhi = np.power(10, lognhi)
sb = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96, 0.98, 1, 1, 1, 1, 1, 1]

nhi2sb = inter.interp1d(nhi, sb)
sb2nhi = inter.interp1d(sb, nhi)


#cubefolder = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/'
cubefolder = '/net/eos/scratch/gallegos/EAGLE/'

do_inter = False
b = [0, 0, 0, 0, 0]
brun = rdarg(argv, 'run', list, b, int)
chombo, init_evol, col_dens, do_layers, do_dndx = brun
gal_mask, do_nhfrac, do_lls = 1, 0, 0
lmax = int(np.log(maxres / minres) / np.log(2))  # 3 if res 4096, 2 for 2048, 1 for 1024, 0 for 512

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

_cub = []
if do_dndx: mats = {}
for snap in snaps:
	s = str(snap)
	saeed = '/net/eos/scratch/gallegos/EAGLE/'#'/net/galaxy-data/export/galaxydata/saeed/EAGLE/RecalL0025N0752/'#
	
	print 'snapshot', snap
	red = redshifts[snap]
	# sbpeak = sbpeaks[snap]
	asec2kpc = asec2kpcs[snap]
	nl = zlens[snap]  # *3
	print nl, 'layers'
	asec2pix = asec2kpc * (1 + red) * kpc2pix
	print 'asec2pix', asec2pix
	sb2flux = asec2pix ** -2.
	# cat = getdata('../../EAGLE/cats/gals_snap%d.fits' % snap, 1)
	# size = cat['Mmean200']#in pkpc
	# ngals = len(size)
	# rad = np.max(([1]*ngals, size/asec2kpc), 0)*asec2pix
	_rad = rad * asec2pix
	fg = 1.
	# check lower or upper case E in the scientific notation!!!
	_sgamma = '%.2E_%.2E_%.2E' % tuple(np.array(gamma_bkg[snap]) * fg)
	_sheat = '%.2E_%.2E_%.2E' % tuple(np.array(heat_bkg[snap]) * fg)
	sgamma = '%.2E %.2E %.2E' % tuple(np.array(gamma_bkg[snap]) * fg)
	sheat = '%.2E %.2E %.2E' % tuple(np.array(heat_bkg[snap]) * fg)
	sname = snames[snap]
	
	fname = '%s/snapshot_%s/' % (saeed, sname)
	if chombo:
		print '\n\n\nChombo'
		chname = fname + 'chombo/snap_%s_%d_%d_%d.chombo.hdf5' % \
		         (sname, minres, maxres, npref)
		if not os.path.isfile(chname) or overwrite:
			if os.path.isfile(chname):
				print 'removing chombo file'
				os.system('rm %s' % chname)
			glob.os.chdir(fname + 'chombo')
			os.system('./chombo.sh %d %d %d' % (minres, maxres, npref))
		glob.os.chdir(scodes)
	if init_evol:
		print '\n\n\nInit evol'
		glob.os.chdir(fname + 'init_evol')
		finit = 'SO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000' % \
		        (sname, minres, maxres, npref, _sgamma, _sheat, ssthr[snap], scase)
		if not os.path.isfile(finit) or overwrite:
			if os.path.isfile(finit):
				print 'removing init evol file'
				os.system('rm %s' % finit)
			srun = './init_evol.sh %d %d %d 25. %.3f "%s" "%s" %.2E' % \
			       (minres, maxres, npref, red, sgamma, sheat, ssthr[snap])
			print srun
			os.system(srun)
		glob.os.chdir(scodes)
	if col_dens:
		print '\n\n\nColumn density!'
		glob.os.chdir(fname + 'column_density')
		finit = '../init_evol/SO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000' % \
		        (sname, minres, maxres, npref, _sgamma, _sheat, ssthr[snap], scase)
		if os.path.isfile(finit):
			nc = 39  # number of cores
			_s = saeed + '/snapshot_%s/column_density/layers/' % sname
			_fl = '%sSO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000_var_%s_proj_%d_lmax_%d_l_' % \
			      (_s, sname, minres, maxres, npref, _sgamma, _sheat, ssthr[snap], scase, type, 1, lmax)
			isfiles = True
			for i in range(nl):
				_file = '%s%d_%d.fits' % (_fl, i + 1, nl)
				_isfile = os.path.isfile(_file)
				if overwrite and _isfile: os.system('rm %s' % _file)
				isfiles &= _isfile
			if not isfiles or overwrite:
				srun = './column_density_layers.sh %s %s %d %d %d' % (finit, type, lmax, nl, nc)
				print srun
				os.system(srun)
		else:
			print '%s 1 um' % finit
		glob.os.chdir(scodes)
	
	sname = snames[snap]
	_s = saeed + '/snapshot_%s/column_density/layers/' % sname
	
	for coord, proj in zip(['x'], [1]):  # zip(coordnames, [1, 2, 3]):
		if gal_mask:
			cat = getdata('../../EAGLE/cats/gals_snap%d.fits' % snap, 1)
			outname = '%s/snap%s_%s.galmask_%darcsec_zw%d.fits' % (cubefolder, snap, coord, rad, zw0)
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
					mask[zc[i] - zw0: zc[i] + zw0 + 1, gal] = ids[i]
				hdu.data = mask
				hdu.writeto(outname, clobber=True)
		
		if do_layers:
			print 'combinig layers!'
			outname = '%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.%s.fits' % (
			cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[snap], scase, type)
			if not os.path.isfile(outname) or overwrite:
				cubes = []
				for i in range(nl):
					flayer = '%sSO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000_var_%s_proj_%d_lmax_%d_l_%d_%d.fits' % \
					         (_s, sname, minres, maxres, npref, _sgamma, _sheat, ssthr[snap], scase, type, proj, lmax,
					          i + 1, nl)
					cubes.append(getdata(flayer))
				
				cubes = np.array(cubes)
				hdu.data = cubes
				hdu.writeto(outname, clobber=True)
				print 'done %s cube' % type
			
			else:
				cubes = getdata(outname)
			
		if do_nhfrac:
			cubes = getdata('%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E.%s.fits' %
			                (cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[snap], type))
			logcube = np.log10(cubes)
			nhis = np.arange(15, 22, .1)
			cool = (logcube > nhis[0]) & (logcube < nhis[-1])
			hist, bins = np.histogram(logcube[cool], nhis)
			sumh = float(np.sum(hist))
			cumh = [np.sum(hist[i:]) / sumh for i in range(len(hist))]
			_dndz = np.array(cumh) * nl / dz[snap]
			
			import matplotlib.pyplot as plt
			
			plt.figure(figsize=(7, 7))
			plt.scatter(bins[:-1], np.log10(hist) - bins[:-1])
			plt.xlabel(type)
			plt.ylabel('log(f_%s)' % type)
			# plt.ylim([-4.2, 0.2])
			plt.savefig('../../EAGLE/%sfrac_%s_%s_%d_%d_%d_%.2E.jpg' %
			            (type, coord, model, minres, maxres, nl, ssthr[snap]))
			plt.close()
		
		if do_lls:
			zs, ls, a = 3, 1.46, 1.7
			_ls, _a = .11, .22
			
			
			def l(z, zs=3, ls=1.46, a=1.70):
				return ls * ((1 + z) / (1 + zs)) ** a
			
			
			fnhiname = '%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.NHI.fits' % (
				cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[snap], scase)
			fllsname = '%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.LLS.fits' % (
				cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[snap], scase)
			fllsnamecorr = '%s/snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.LLScorr.fits' % (
				cubefolder, snap, coord, model, minres, maxres, npref, nl, ssthr[snap], scase)
			if not os.path.isfile(fllsname) or overwrite:
				cubes = getdata(fnhiname)
				_cubes = np.copy(cubes)
				llsNHI = 17.5
				lls0 = cubes > 10 ** llsNHI
				# lls = np.zeros(cubes.shape)
				# lls[cool] = 1
				hdu.data = lls0.astype(float)  # lls
				hdu.writeto(fllsname, clobber=True)
				pre = 10
				lz = l(red, zs, ls, a)
				dndz0 = np.mean(lls0) * nl / dz[snap]
				dndz = dndz0
				diff = dndz / lz
				# nhi = cubes*lz/dndz
				# dndz2 = np.mean(nhi>10**llsNHI)*nl/dz[snap]
				# print 'new dndz', dndz2
				
				while abs(diff - 1) > 1e-3:
					_cubes /= diff
					lls = _cubes > 10 ** llsNHI
					dndz = np.mean(lls) * nl / dz[snap]
					diff = dndz / lz
					print 'diff', diff, 'dndz', dndz
					
				hdu.data = _cubes
				hdu.writeto(fnhiname.replace('.NHI.', '.NHIcorr.'), clobber=True)
				hdu.data = lls.astype(float)
				hdu.writeto(fllsnamecorr, clobber=True)
			else:
				lls0 = getdata(fllsname)
			
			# error propagation
			error = np.sqrt(
				((1 + red) / (1 + zs)) ** a * _ls ** 2 * ls * a * ((1 + red) / (1 + zs)) ** (a - 1) * _a ** 2)
			lz = l(red, zs, ls, a)
			print 'dndz sim %.2f simcorr %.2f obs %.2f +-%.2f corr %.2f +-%.2f' % (
			dndz0, dndz, lz, error, lz / dndz, error / dndz)
		
		if do_dndx:
			print 'dndx!'
			_type = type  # 'NHI'
			cname = 'snap%d_%s_%s_%d_%d_%d_%d_%.2E%s.%s' % \
			        (snap, coord, model, minres, maxres, npref, nl, ssthr[snap], scase, _type)
			cubename = '%s/%s.fits' % (cubefolder, cname)
			zw = 1  # nl/5
			print cubename
			v = 2.99792458e5 * dz[snap] * zw / float(nl) / (1 + red)
			fdat = '../../UVB/dndzdx_%s_snap%s_zw%d_ssthr%.2e_%d_%d_%d.%s.dat' % (
			model, snap, zw, ssthr[snap], minres, maxres, npref, _type)
			if not os.path.isfile(fdat) or overwrite:
				cubes = getdata(cubename)
				if type == 'density':
					cubes *= (1 + red) ** 2 * coml * params.Mpc2cm / nl * params.meanbardens
				zl, yl, xl = cubes.shape
				
				zn = zl / zw
				if zw == 1:
					_c = cubes
				
				else:
					_c = []
					for i in range(zn): _c.append(np.sum(cubes[i:i + zw, :, :], 0))
					_c = np.array(_c)
				lc = np.log10(_c)
				dxdz = (1. + red) ** 2 / sqrt(params.omega_m * (1 + red) ** 3 + params.omega_l)
				_lls = np.arange(15, 22.6, .25)
				dndx, dndz, cddf = [[], [], []]
				
				for i in range(len(_lls) - 1):
					cool = (lc >= _lls[i]) & (lc < _lls[i + 1])
					dNHI = float(10. ** _lls[i + 1] - 10 ** _lls[i])
					# lpos = np.nanmean(_c[cool])#(10.**_lls[i+1]+10**_lls[i])/2.
					# if np.isnan(lpos): lpos = 10**_lls[i]
					scool = np.mean(cool) * zn
					_dndz = scool / dz[snap]
					_dndx = _dndz / dxdz
					_cddf = _dndx / dNHI
					print 'NHI %.2f %.2f' % (_lls[i], _lls[i + 1]), 'dndz (diff)', _dndz, 'cddf', _cddf, 'dndx', _dndx
					dndx.append(_dndx)
					dndz.append(_dndz)
					cddf.append(_cddf)
				
				lnhis = _lls[:-1]
				ldndx = np.log10(dndx)
				dndz = np.array(dndz)
				cddf = np.log10(cddf)
				# _dndz = [np.sum(dndzs[i:]) for i in range(len(dndzs))]
				sb = nhi2sb(10 ** lnhis)
				mat = np.array([lnhis, dndz, dndx, cddf, sb])
				np.savetxt(fdat, mat.T, header='logNHI dndz dndx cddf SBfit', fmt='%.3f')
			else:
				mat = np.loadtxt(fdat).T

