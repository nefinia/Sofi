#!/usr/bin/env python
__author__ = 'gallegos'
import glob
import h5py
# import eagleSqlTools
import matplotlib.pyplot as plt
import os
import scipy.interpolate as inter
from math import sqrt
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap
from pyfits import getdata, PrimaryHDU
from sys import argv
import params

import numpy as np

from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
snaps = rdarg(argv, 'snap', list, [6, 7, 8, 9, 10, 11, 12, 13, 14], int)
ns = len(snaps)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
caseA = rdarg(argv, 'caseA', bool, True)#'HM12')#
scase = '_CaseA' * caseA
halfplot = rdarg(argv, 'halfplot', bool, False)
extraname = rdarg(argv, 'extraname', str, '')  # 'HM12')#
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
npref = rdarg(argv, 'npref', int, 12)
sbhist = rdarg(argv, 'sbhist', bool, False)
_ssthr = rdarg(argv, 'ssthr', float, 6.73E-3)#6.73E-3 or 1e10 or None
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
types = rdarg(argv, 'type', list, ['NHI'], str)  # NHtot, f_NHI, NHII
type = rdarg(argv, 'type', str, 'NHI')  # NHtot, f_NHI, NHII
radprof = rdarg(argv, 'radprof', bool, False)
rad = rdarg(argv, 'rad', int, 8)
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)
model = rdarg(argv, 'model', str, 'HM12')  # 'HM01' or 'HM12')#
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)
zw = rdarg(argv, 'zw', int, 1)

lognhi = [12, 13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21,
		  22, 23, 30]
nhi = np.power(10, lognhi)
sb = [0.0000, 0.00001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96,
	  0.98, 1, 1, 1, 1, 1, 1]

nhi2sb = inter.interp1d(nhi, sb)
lognhi2sb = inter.interp1d(lognhi, sb)
sb2nhi = inter.interp1d(sb, nhi)

do_eagle = 1

logNHI = [12.875, 13.125, 13.375, 13.625, 13.875, 14.125, 14.375, 14.625, 14.875, 15.25, 15.75, 16.25, 16.75, 17.25,
		  17.75]
f = [-11.077, -11.437, -11.795, -12.163, -12.507, -12.992, -13.459, -13.781, -14.29, -14.933, -15.97, -16.674, -17.219,
	 -17.896, -19.174]

logNHI = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
f = [-10.1, -11, -13, -15, -17, -17.5, -19.5, -20.5, -21.5, -23, -24.5]

logNHI = [12, 12.875, 13.125, 13.375, 13.625, 13.875, 14.125, 14.375, 14.625, 14.875, 15.25, 15.75, 16.25, 16.75, 17.25,
		  17.75, 18, 19, 20, 21, 22]
f = [-10.1, -11.077, -11.437, -11.795, -12.163, -12.507, -12.992, -13.459, -13.781, -14.29, -14.933, -15.97, -16.674,
	 -17.219,
	 -17.896, -19.174, -19.5, -20.5, -21.5, -23, -24.5]

nhi2f = inter.interp1d(logNHI, f)
nhirange = np.arange(15, 21, .1)
f2 = nhi2f(nhirange)

sb = lognhi2sb(nhirange)

nhifrac = np.power(10, f2) * np.power(10, nhirange)
nhisum = np.sum(nhifrac)
sbc = nhifrac * sb
sbmax = np.amax(sbc)

# low = nhirange < 17.5
# sblow = np.sum(sbc[low])
# sbhigh = np.sum(sbc[~low])
# print 'low', sblow, 'high', sbhigh

sbcum = np.array([np.sum(sbc[i:]) for i in range(len(sbc))])
sbsum = np.sum(sbc)
# sbcum = [np.sum(_sb[i:])/sbtot for i in range(len(cddf))]

# from scipy import integrate
# integrate.quad(lognhi2sb, 12, 21, args=(alpha, mstar))
colors = ['red', 'orange', 'green', 'pink', 'blue', 'purple']
zcolors = {14: 'black', 13: 'darkgrey', 12: 'darkorchid', 11: 'dodgerblue', 10: 'darkturquoise', 9: 'limegreen', 8: 'khaki', 7: 'orange', 6: 'red'}

fsize = 14
plt.figure()
plt.plot(nhirange, sbc, label='SB contribution')
plt.plot(nhirange, sb, label='SB conversion')
plt.plot(nhirange, nhifrac / np.amax(nhifrac), label='NHI fraction')
# plt.plot(nhirange, f2, label='f conversion')
plt.semilogy()
plt.legend()
# plt.xlim(14, 21)
# plt.ylim(-27, -9)
# plt.xticks([16, 18, 20, 22])
# plt.yticks([-25, -23, -21, -19, -17, -15, -13, -11, -9])
plt.xlabel(r'log$_{10}\mathrm{N_{HI}[cm^{-2}]}$]')
# plt.ylabel(r'SB $\rm{[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}$')
plt.grid()
plt.savefig('../../Figures/sbcont_obs.jpg')
plt.close()

if do_eagle:
	gamma_bkg = params.gamma_bkg[model]

	heat_bkg = params.heat_bkg[model]
	snames = params.snames
	redshifts = params.redshifts
	asec2kpcs = params.asec2kpcs

	dz = params.dz
	sbpeaks = params.sbpeaks
	zlens = params.zlens
	dndz = params.dndz
	fLLScorr = params.fLLScorr

	nhi_fit = params.nhi_fit[model]

	lcube = params.lcube
	coml = params.coml
	com2pix = params.com2pix
	kpc2pix = params.kpc2pix

	mats = {}

	fdats = 'dndzdx_HM01_snap8_zw9_v458.dat', 'dndzdx_HM01_snap9_zw7_v398.dat', 'dndzdx_HM01_snap10_zw6_v371.dat', \
			'dndzdx_HM01_snap11_zw5_v342.dat', 'dndzdx_HM01_snap12_zw5_v380.dat', 'dndzdx_HM01_snap13_zw4_v358.dat', \
			'dndzdx_HM01_snap14_zw3_v261.dat'

	for snap in snaps:
		if _ssthr is not None: ssthr = _ssthr
		else: ssthr = params.nhSS[model][snap]
		fg = 1
		sgamma = '%.2E_%.2E_%.2E' % tuple(np.array(gamma_bkg[snap]) * fg)
		fdat = '../../UVB/dndzdx_%s_snap%d_zw%d_ssthr%.2e_%d_%d_%d.%s.dat' %\
			   (model, snap, zw, ssthr, minres, maxres, npref, type)
		mat = np.loadtxt(fdat).T
		mats[snap] = mat

	fig, ax = plt.subplots(1, figsize=(7, 8))
	fsize = 16
	i = 0
	Mpc2cm = 3.086e24
	for snap in snaps:
		s = str(snap)
		lnhis, dndz, dndx, cddf, sb = mats[snap]
		llspos = lnhis >= 17.5
		dlapos = lnhis >= 20
		# _dndx = np.sum(10 ** cddf[llspos] * (10 ** (lnhis[llspos] + .1) - 10 ** lnhis[llspos]))
		# _dndz = _dndx * (1 + redshifts[snap]) ** 2 / np.sqrt(params.omega_m * (1 + redshifts[snap]) ** 3 - params.omega_l)
		_dndz = np.nansum(dndz[llspos])
		_dndx = np.nansum(dndx[llspos])
		dndzcum = np.array([np.sum(dndz[k:]) for k in range(len(dndz))])
		llsfit = inter.interp1d(dndzcum, lnhis)
		#lcool = llsfit(params.dndz[snap])

		print 'snap', snap, 'z', redshifts[snap], 'dndx', _dndx, 'dndz', _dndz#, 'NHI match', lcool
		# cddf = np.log10(dndz / (10 ** (lnhis + .1) - 10 ** lnhis))
		# print '%d dndz LLS %.3f DLA %.3f' % (snap, np.sum(dndz[llspos]), np.sum(dndz[dlapos])),\
		#		'dndx LLS %.3f DLA %.3f' %(np.sum(dndx[llspos]), np.sum(dndx[dlapos]))

		#ln2 = np.arange(13, 23, .125)
		#N2n = 10 ** ln2 * zlens[snap] / 25 / (1 + redshifts[snap]) / Mpc2cm
		#gamma_eff = gamma_bkg[snap][0] * \
	#				(0.98 * (1. + (N2n / ssthr[snap]) ** 1.64) ** (-2.28) + \
#					 0.02 * (1. + N2n / ssthr[snap]) ** (-0.84))

		# print 'NHI\t\tdensities\tgamma eff'
		# for j in range(len(ln2)):
		#	print '%.3f\t\t%.3e\t\t%.3e' % (ln2[j], N2n[j], gamma_eff[j])

		# cddf2 = np.log10(dndz/np.sqrt((1+redshifts[snap])/.3)/(10**(lnhis+.25)-10**(lnhis)))
		plt.plot(lnhis, cddf, label=r'$z=%.2f$' % redshifts[snap], color=zcolors[snap], linestyle='--',
				 dashes=(3, len(snaps) - i))
		# plt.plot(lnhis, cddf2, label=r'$z=%.2f$' % redshifts[snap], color=zcolors[len(snaps)-i-1], linestyle='--', dashes=(3, len(snaps)-i))
		# lnhis, cddf, dndz, ldndx, sb = mats2[snap]
		# plt.plot(lnhis, ldndx, label=r'$z=%.2f$  all z' % redshifts[snap], color=zcolors[len(snaps)-i-1], linestyle='--', dashes=(3, len(snaps)-i))
		i += 1

	lnhi = np.array([19., 19.5, 20.3, 20.65, 21., 21.35, 21.7])

	f = np.array([-20.5, -21.2, -21.9, -22.5, -23., -23.8])
	fmin = np.array([-20.4, -21.1, -21.9, -22.4, -22.9, -23.6])
	fmax = np.array([-20.8, -21.4, -22, -22.6, -23.2, -24.1])
	x = (lnhi[1:] + lnhi[:-1]) / 2.
	xerr = np.array([(lnhi[i] - lnhi[i + 1]) / 2 for i in range(len(x))])
	yerr = fmax - fmin
	plt.errorbar(x, f, xerr=xerr, yerr=yerr, label=r'Peroux+05 $1.78<z<3.5$', fmt="none", fontsize=fsize, color='orange')

	f = np.array([-20.2, -21.1, -21.6, -22.5, -23., -23.8])
	fmin = np.array([-20.4, -21.1, -21.9, -22.4, -22.9, -23.6])
	fmax = np.array([-20.8, -21.4, -22, -22.6, -23.2, -24.1])

	yerr = fmax - fmin
	plt.errorbar(x, f, xerr=xerr, yerr=yerr, label=r'Peroux+05 $3.5<z<5.0$', fmt="none", fontsize=fsize, color='brown')

	x = np.arange(13.2, 17.3, .5)
	y = 10.322 - 1.65 * x
	plt.errorbar(x, y, xerr=.25, yerr=.017 * x, label=r'Rudie+13 $2.0<z<2.8$', fmt="none", fontsize=fsize, color='green')


	xmin = np.array([20.00, 20.10, 20.20, 20.30, 20.40, 20.50, 20.60, 20.70, 20.80, 20.90, 21.00, 21.10, 21.20, 21.30, 21.40,
	        21.50, 21.60, 21.70, 21.80, 21.90, 22.00, 22.20])
	xmax = np.array([20.10, 20.20, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21, 21.2, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8,
	        21.9, 22, 22.2, 22.4])
	x = (xmin + xmax)/2.
	xerr = (xmax - xmin) / 2.
	y = [-21.44, -21.47, -21.59, -21.68, -21.82, -21.98, -22.14, -22.32, -22.51, -22.67, -22.91, -23.11, -23.28, -23.58,
	     -23.81, -24.01, -24.20, -24.62, -24.85, -25.60, -26.05, -26.25]
	yerr = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.05, 0.06, 0.07, 0.08, 0.08, 0.12,
	        0.18, 0.53, 0.53, 0.53]
	plt.errorbar(x, y, xerr=xerr, yerr=yerr, label=r'Noterdaeme+12 $\langle z=2.5\rangle$', fmt="none", fontsize=fsize,
	             color='blue')



	if 0:
		x = np.array(
			[12.11, 12.21, 12.38, 12.63, 12.84, 13.03, 13.28, 13.49, 13.72, 13.93, 14.18, 14.43, 14.7, 14.91, 15.05,
			 15.26,
			 15.47, 15.74, 15.97, 16.2, 16.41, 16.54, 16.76, 16.95, 17.08, 17.25, 17.47, 17.79, 18.08, 18.35, 18.6,
			 18.92,
			 19.21, 19.5, 19.79, 20.08, 20.34, 20.56, 20.79, 21.04, 21.23, 21.46, 21.53, 21.67, 21.73, 21.82])
		y = np.array(
			[-10.45, -10.32, -10.36, -10.53, -10.79, -11.04, -11.34, -11.68, -12.11, -12.41, -12.75, -13.09, -13.47,
			 -13.9,
			 -14.24, -14.67, -15.05, -15.43, -15.82, -16.12, -16.41, -16.84, -17.18, -17.4, -17.57, -17.78, -18.08,
			 -18.46,
			 -18.89, -19.23, -19.61, -19.82, -20.25, -20.54, -20.93, -21.31, -21.74, -22.08, -22.46, -22.93, -23.4,
			 -24.09,
			 -24.77, -25.33, -26.01, -26.57])
		plt.plot(x, y, label=r'Nastasha z=3.5')
		dndxs = [np.sum(10 ** y[i:-1] * (10 ** x[i + 1:] - 10 ** x[i:-1])) for i in range(len(x) - 1)]
		fit = inter.interp1d(x[:-1], dndxs)
		dndx = fit(17.5)
		dndz = dndx * 4.5 ** 2 / np.sqrt(params.omega_m * 4.5 ** 3 - params.omega_l)
		print 'nastasha z=3.5 dndx', dndx, 'dndz', dndz

		x = np.array(
			[15.76, 15.92, 16.1, 16.3, 16.48, 16.67, 16.85, 17.03, 17.2, 17.44, 17.57, 17.74, 17.95, 18.15, 18.39,
			 18.57,
			 18.78, 19.01, 19.23, 19.46, 19.72, 19.92, 20.13, 20.29, 20.44, 20.64, 20.87, 20.99, 21.16, 21.3, 21.47,
			 21.6,
			 21.68, 21.82, 21.9, 22, 22.12, 22.18, 22.3])
		y = np.array(
			[-15.01, -15.33, -15.69, -16.05, -16.37, -16.74, -17.08, -17.34, -17.58, -17.92, -18.18, -18.33, -18.65,
			 -18.91, -19.19, -19.45, -19.69, -19.95, -20.22, -20.52, -20.9, -21.18, -21.42, -21.65, -21.91, -22.17,
			 -22.49,
			 -22.77, -23.03, -23.31, -23.72, -24.1, -24.42, -24.84, -25.14, -25.48, -25.81, -26.07, -26.43])
		import scipy.signal as sp

		plt.plot(x, y, label=r'Ali z=5')
		dndxs = [np.sum(10 ** y[i:-1] * (10 ** x[i + 1:] - 10 ** x[i:-1])) for i in range(len(x) - 1)]
		fit = inter.interp1d(x[:-1], dndxs)
		dndx = fit(17.5)
		dndz = dndx * 6 ** 2 / np.sqrt(params.omega_m * 6 ** 3 - params.omega_l)
		print 'Ali z=5 dndx', dndx, 'dndz', dndz

		x = np.array([15.21, 15.39, 15.61, 15.84, 16.07, 16.31, 16.56, 16.82, 17.04, 17.33, 17.58, 17.88, 18.11, 18.42,
					  18.68, 19, 19.27, 19.56, 19.86, 20.13, 20.22, 20.51, 20.74, 20.96, 21.19, 21.32, 21.48, 21.57,
					  21.73, 21.8, 21.9, 22.09, 22.17])
		y = np.array(
			[-14.4, -14.72, -15.18, -15.62, -16.04, -16.42, -16.88, -17.26, -17.66, -18.06, -18.42, -18.8, -19.12,
			 -19.48,
			 -19.82, -20.2, -20.54, -20.9, -21.26, -21.64, -21.74, -22.14, -22.46, -22.82, -23.16, -23.38, -23.74,
			 -24.08,
			 -24.54, -24.84, -25.26, -25.66, -25.9
			 ])
		plt.plot(x, y, label=r'Ali z=4')
		dndxs = [np.sum(10 ** y[i:-1] * (10 ** x[i + 1:] - 10 ** x[i:-1])) for i in range(len(x) - 1)]
		fit = inter.interp1d(x[:-1], dndxs)
		dndx = fit(17.5)
		dndz = dndx * 5 ** 2 / np.sqrt(params.omega_m * 5 ** 3 - params.omega_l)
		print 'Ali z=4 dndx', dndx, 'dndz', dndz

		x = np.array(
			[15.21, 15.29, 15.45, 15.56, 15.74, 15.89, 16.09, 16.26, 16.52, 16.68, 16.86, 17.07, 17.25, 17.44, 17.62,
			 17.81, 18.02, 18.3, 18.51, 18.83, 19.01, 19.24, 19.52, 19.73, 20.03, 20.23, 20.52, 20.8, 21.03, 21.25,
			 21.41, 21.53, 21.61, 21.68, 21.86, 21.98, 22.09, 22.19])

		y = np.array(
			[-14.9, -15.02, -15.3, -15.52, -15.78, -16.1, -16.42, -16.72, -17.06, -17.3, -17.56, -17.86, -18.16, -18.42,
			 -18.7, -18.92, -19.22, -19.54, -19.88, -20.18, -20.44, -20.7, -20.94, -21.26, -21.52, -21.82, -22.1,
			 -22.48, -22.84, -23.2, -23.4, -23.76, -24.1, -24.46, -25, -25.34, -25.68, -25.94])
		plt.plot(x, y, label=r'Ali z=3')
		dndxs = [np.sum(10 ** y[i:-1] * (10 ** x[i + 1:] - 10 ** x[i:-1])) for i in range(len(x) - 1)]
		fit = inter.interp1d(x[:-1], dndxs)
		dndx = fit(17.5)
		dndz = dndx * 4 ** 2 / np.sqrt(params.omega_m * 4 ** 3 - params.omega_l)
		print 'Ali z=3 dndx', dndx, 'dndz', dndz

		x = np.array(
			[15.27, 15.48, 15.73, 16.03, 16.43, 16.65, 16.95, 17.21, 17.46, 17.61, 17.81, 17.99, 18.18, 18.33, 18.57,
			 18.77,
			 19.04, 19.25, 19.43, 19.62, 19.83, 19.99, 20.07, 20.21, 20.41, 20.64, 20.93, 21.07, 21.26, 21.42, 21.55,
			 21.65,
			 21.75, 21.92, 22, 22.12, 22.22, 22.35])
		y = np.array(
			[-15.27, -15.71, -16.07, -16.55, -17.14, -17.54, -17.88, -18.3, -18.65, -18.89, -19.09, -19.27, -19.59,
			 -19.79,
			 -20.03, -20.28, -20.58, -20.78, -20.94, -21.18, -21.38, -21.54, -21.58, -21.77, -21.99, -22.27, -22.59,
			 -22.89,
			 -23.29, -23.68, -24.04, -24.36, -24.72, -25.14, -25.62, -25.89, -26.27, -26.59])
		plt.plot(x, y, label=r'Ali z=2')
		dndxs = [np.sum(10 ** y[i:-1] * (10 ** x[i + 1:] - 10 ** x[i:-1])) for i in range(len(x) - 1)]
		fit = inter.interp1d(x[:-1], dndxs)
		dndx = fit(17.5)
		dndz = dndx * 3 ** 2 / np.sqrt(params.omega_m * 3 ** 3 - params.omega_l)
		print 'Ali z=2 dndx', dndx, 'dndz', dndz

		x = np.array(
			[15.23, 15.41, 15.61, 15.9, 16.19, 16.54, 16.76, 16.96, 17.25, 17.5, 17.89, 18.23, 18.59, 18.9, 19.17, 19.4,
			 19.65, 19.85, 20.22, 20.46, 20.67, 20.9, 20.98, 21.11, 21.27, 21.47, 21.64, 21.84, 22.05, 22.19, 22.33])
		y = np.array(
			[-15.5, -15.78, -16.08, -16.54, -16.98, -17.5, -17.8, -18.1, -18.5, -18.86, -19.34, -19.76, -20.16, -20.48,
			 -20.74, -21, -21.26, -21.46, -21.84, -22.1, -22.38, -22.72, -23.02, -23.4, -23.72, -24.18, -24.52, -24.98,
			 -25.66, -26.06, -26.56])

		plt.plot(x, y, label=r'Ali z=1')
		dndxs = [np.sum(10 ** y[i:-1] * (10 ** x[i + 1:] - 10 ** x[i:-1])) for i in range(len(x) - 1)]
		fit = inter.interp1d(x[:-1], dndxs)
		dndx = fit(17.5)
		dndz = dndx * 2 ** 2 / np.sqrt(params.omega_m * 2 ** 3 - params.omega_l)
		print 'Ali z=1 dndx', dndx, 'dndz', dndz

		# NHI_SS = 4e17*2.49e-18/sigmaHI
		# for snap in snaps
		#

		# y = 7.562-1.425*x
		# plt.errorbar(x, y, xerr=.25, yerr=.033*x, label=r'Rudie+13 $2.0<z<2.8$ CGM', fmt="none")

	plt.xlim(15, 22.6)
	plt.ylim(-26.5, -14.5)
	plt.xticks([16, 17, 18, 19, 20, 21, 22], fontsize=fsize)
	plt.yticks([ -26, -24, -22, -20, -18, -16], fontsize=fsize)
	plt.xlabel(r'log$_{10}\mathrm{(N_{HI})[cm^{-2}]}$', fontsize=fsize)
	plt.ylabel(r'log$_{10}[\mathrm{f(N_{HI}, X)}]$', fontsize=fsize)
	h, l = ax.get_legend_handles_labels()
	plt.legend(h, l, fontsize=fsize)
	leg1 = ax.legend(h[:ns], l[:ns], loc=[.68, 0.55], fontsize=fsize)
	#leg1.set_title('EAGLE', prop={'size': fsize})
	leg2 = ax.legend(h[ns:], l[ns:], loc=[.02, .02], fontsize=fsize)
	ax.add_artist(leg1)
	ssnap = ''
	for snap in snaps:
		if _ssthr is None: sssthr = '-%.0E' % params.nhSS[model][snap]
		else: sssthr = ''
		ssnap += '_%d%s' % (snap, sssthr)
	if _ssthr is not None: sssthr += '_ssthr%.2E' % _ssthr

	plt.tight_layout()

	plt.savefig('../../Figures/cddf%s_%d_%d_%d_%s%s.%s%s.pdf' % (ssnap, minres, maxres, npref, sssthr, scase, type, extraname))
	plt.close()
