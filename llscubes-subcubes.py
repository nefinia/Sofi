#!/usr/bin/env python
__author__ = 'gallegos'
import glob
# import eagleSqlTools
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
LLScubes = rdarg(argv, 'LLScubes', bool, False)
sql = rdarg(argv, 'sql', bool, False)
sphericalLLS = rdarg(argv, 'sphericalLLS', bool, False)
pairLLS = rdarg(argv, 'pairLLS', bool, False)
circlehist = rdarg(argv, 'circlehist', bool, False)
h2d = rdarg(argv, 'h2d', bool, False)
kde = rdarg(argv, 'kde', bool, False)
nhiprof = rdarg(argv, 'nhi', bool, False)
res = rdarg(argv, 'res', int, 512)
model = rdarg(argv, 'model', str, 'HM01')#'HM12')#
npref = rdarg(argv, 'npref', int, 20)
rad = rdarg(argv, 'rad', int, 6)  # arcsec
radprof = rdarg(argv, 'radprof', bool, False)
sbhist = rdarg(argv, 'sbhist', bool, False)
_ssthr = rdarg(argv, 'ssthr', float, None)#6.73e-3)# 1e10 or None
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
#cubefolder = '../../EAGLE/'
#if you run it in eos
cubefolder = '/scratch/gallegos/EAGLE/'

do_inter = False
b = [0, 	0, 		0, 			0]
brun = rdarg(argv, 'run', list, b, int)
basic, combine, do_lls, do_dndx = brun

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
	#data = getdata('/net/abnoba/scratch2/gallegos/Research/MUSE/EAGLE/cats/gals_snap%d.fits' % snap, 1)
	#xs, ys, zs, ids = data['x'], data['y'], data['z'], data['ID']
	red = redshifts[s]
	print 'snapshot', snap, 'redshift', red
	# sbpeak = sbpeaks[snap]
	#size and nl starting from the galaxy coordinate -> total cube size size*2,size*2,nl*2
	#size = .0256*(1+red)*h*asec2kpcs[s] # this is to ensure I have 0.2 arcsec per pixel
	#size2 = size*2/h
	nl = zlens[s]/cuts
	#sizel = 3.5*1.028e-3*lh/dz[s] # in order to have 7 layers of 1.25 Angstrom

	fg = 1.
	# check lower or upper case E in the scientific notation!!!
	_sgamma = '%.2E_%.2E_%.2E' % tuple(np.array(gamma_bkg[s])*fg)
	_sheat = '%.2E_%.2E_%.2E' % tuple(np.array(heat_bkg[s])*fg)
	sgamma = '%.2E %.2E %.2E' % tuple(np.array(gamma_bkg[s])*fg)
	sheat = '%.2E %.2E %.2E' % tuple(np.array(heat_bkg[s])*fg)
	sname = snames[s]
	fname = '%s/snapshot_%s/' % (cubefolder, sname)
	
	cubeall = np.zeros((nl*cuts, res*cuts, res*cuts))
	exe = './P2C/build/P2C_omp'
	fin = '/scratch/gallegos/EAGLE/RecalL0025N0752/snapshot_%s/snap_%s.%%d.hdf5' % (sname, sname)
	
	fchombo = fname + 'chombo/snap_%s_%d_%d_%d.chombo.hdf5' % \
	          (sname, res, res, npref)
	finit_base = fname + '/init_evol/SO.snap_%s_%d_%d_%d_%s_%s_%.2E%s' % \
	        (sname, res, res, npref, _sgamma, _sheat, ssthr[s], scase)
	_s = cubefolder + '/snapshot_%s/column_density/layers/' % sname
	flayer_base = '%sSO.snap_%s_%d_%d_%d_%s_%s_%.2E%s.ts0000_var_%s_proj_%d_lmax_0_l_' % \
	              (_s, sname, res, res, npref, _sgamma, _sheat, ssthr[s], scase, type, proj)

	if basic:
		for i in range(cuts):
			for j in range(cuts):
				for k in range(cuts):
					print i, j, k
					isfiles = True
					for l in range(nl):
						_file = '%s%d_%d.%d%d%d.fits' % (flayer_base, l + 1, nl, i, j, k)
						_isfile = os.path.isfile(_file)
						if overwrite and _isfile: os.system('rm %s' % _file)
						isfiles &= _isfile
					if not isfiles or overwrite:
						finit = finit_base+'.ts0000'
						finit2 = finit_base+'.%d%d%d.ts0000' % (i, j, k)
						
						if not os.path.isfile(finit2) or overwrite:
							# if chombo:
							print '\n\n\nChombo'
							glob.os.chdir(fname + 'chombo')
							os.system('%s -base_grid %d -lmax 0 -npref %d -cut_le "%.2f,%.2f,%.2f" -cut_re "%.2f,%.2f,%.2f" ' \
						     '-inp %s -out %s -file_fmt eaglefmt | tee -a %s' \
						     % (exe, res, npref, i*size, j*size, k*size, (i+1)*size, (j+1)*size, (k+1)*size, fin, fchombo,
						        '%d%d%d.log' % (i, j, k)))
			
							#if init_evol:
							print '\n\n\nInit evol'
							glob.os.chdir(fname + 'init_evol')
							if os.path.isfile(finit2):
								print 'removing init evol file'
								os.system('rm %s' % finit2)
							srun = './init_evol.sh %d %d %d %.3f %.3f "%s" "%s" %.2E' % \
								   (res, res, npref, size, red, sgamma, sheat, ssthr[s])
							print srun
							os.system(srun)
							if os.path.isfile(finit):
								os.system('rm %s' % fchombo)
								os.system('mv %s %s' % (finit, finit2))
					if isfiles:
						#if col_dens:
						print '\n\n\nColumn density!'
						glob.os.chdir(fname + 'column_density')
						nc = 20  # number of cores
						srun = './column_density_layers.sh %s %s 0 %d %d' % (finit2, type, nl, nc)
						print srun
						os.system(srun)
						os.system('echo done')
			
					else: print 'subcube exists'
	
	if combine:
		glob.os.chdir(scodes)
		print 'combinig layers!'
	
		for i in range(cuts):
			for j in range(cuts):
				for k in range(cuts):
					print i, j, k
					outname = '%s/snap%d_%s_%s_%d_%d_%d_%.2E%s.%s.%d%d%d.fits' % \
					          (cubefolder, snap, coord, model, res, res, npref, ssthr[s], scase, type, i, j, k)
					#finit = finit_base + '.ts0000'
					#finit2 = '%s.%d%d%d' % (finit_base, i, j, k)
					#if os.path.isfile(finit2): os.system('rm %s'% finit)
					if not os.path.isfile(outname) or overwrite:
						cubes = []
						for l in range(nl):
							flayer = flayer_base.replace('ts0000', 'ts0000.%d%d%d') + '%d_%d.%d%d%d.fits' % \
							         (l + 1, nl, i, j, k)
							cubes.append(getdata(flayer))
						
						cubes = np.array(cubes)
						hdu.data = cubes
						hdu.writeto(outname, clobber=True)
						print 'done %s cube' % type
						#if os.path.isfile(outname):
						#	for l in range(nl):
						#		flayer = '%s%d_%d.%d%d%d.fits' % (flayer_base, l + 1, nl, i, j, k)
								#if os.path.isfile(flayer): os.system('rm %s' % flayer)
					else: cubes = getdata(outname)
					cubeall[i * nl:(i + 1) * nl, k * res:(k + 1) * res, j * res:(j + 1) * res] = cubes
		
		foutall = '%s/snap%d_%s_%s_%d_%d_%d_%.2E%s.%s.fits' % \
				  (cubefolder, snap, coord, model, res*cuts, res*cuts, npref, ssthr[s], scase, type)
		hdu.data = cubeall
		hdu.writeto(foutall, clobber=True)

	if do_lls:
		_nl = nl*cuts
		_res = res*cuts
		zs, ls, a = 3, 1.46, 1.7
		_ls, _a = .11, .22
		def l(z, zs=3, ls=1.46, a=1.70):
			return ls * ((1 + z) / (1 + zs)) ** a
		fnhiname = '%s/snap%d_%s_%s_%d_%d_%d_%.2E%s.NHI.fits' % (cubefolder, snap, coord, model, _res, _res, npref,
		                                                         ssthr[s], scase)

		fllsname = '%s/snap%d_%s_%s_%d_%d_%d.LLS.fits' % (cubefolder, snap, coord, model, _res, _res, _nl)
		fllsnamecorr = '%s/snap%d_%s_%s_%d_%d_%d.LLScorr.fits' % (cubefolder, snap, coord, model, _res, _res, _nl)
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
		glob.os.chdir(scodes)
		print 'dndx!'
		_res = res*cuts
		_nl = nl*cuts
		type = 'NHI'
		cname = 'snap%d_%s_%s_%d_%d_%d_%.2E%s.%s' % \
		        (snap, coord, model, _res, _res, npref, ssthr[s], scase, type)
		#cname = 'snap%d_%s_%s_%d_%d_%d.%s' % \
		#        (snap, coord, model, _res, _res, _nl, type)
		cubename = '%s/%s.fits' % (cubefolder, cname)
		zw = 1  # nl/5
		print cubename
		#v = 2.99792458e5 * dz[s] * zw / float(_nl) / (1 + red)
		fdat = '../../UVB/dndzdx_%s_snap%s_zw%d_ssthr%.0e_%d_%d_%d.%s.dat' % (
		model, snap, zw, ssthr[s], _res, _res, npref, type)
		if not os.path.isfile(fdat) or overwrite:
			cubes = getdata(cubename)
			if type == 'density':
				cubes *= (1 + red) ** 2 * coml * params.Mpc2cm / _nl * params.meanbardens
			zl, yl, xl = cubes.shape
			
			zn = zl / zw
			if zw == 1: _c = cubes
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
				scool = np.mean(cool) * zn
				_dndz = scool / dz[s]
				_dndx = _dndz / dxdz
				_cddf = _dndx / dNHI
				print 'NHI %.2f %.2f' % (_lls[i], _lls[i + 1]), 'dndz', _dndz, 'cddf', _cddf, 'dndx', _dndx
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