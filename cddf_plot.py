#!/usr/bin/env python
__author__ = 'gallegos'
from sys import argv

# import eagleSqlTools
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter

import params
from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
snaps = rdarg(argv, 'snap', list, [8, 9, 10, 11, 12, 13, 14], int)

scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
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
npref = rdarg(argv, 'npref', int, 8)
sbhist = rdarg(argv, 'sbhist', bool, False)
ssthr = rdarg(argv, 'ssthr', float, 1e10)  # 6.73E-3
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
types = rdarg(argv, 'type', list, ['NHI'], str)  # NHtot, f_NHI, NHII
radprof = rdarg(argv, 'radprof', bool, False)
rad = rdarg(argv, 'rad', int, 8)
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)
model = rdarg(argv, 'model', str, 'HM01')  # 'HM12')#
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)

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
zcolors = ['darkorchid', 'dodgerblue', 'darkturquoise', 'limegreen', 'khaki', 'orange', 'red']

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


gamma_bkg = params.gamma_bkg[model]
heat_bkg = params.heat_bkg[model]
snames = params.snames
redshifts = params.redshifts
asec2kpcs = params.asec2kpcs
ssthr = params.nhSS[model]
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
mats2 = {}

zw = 1

fdats = ['dndzdx_HM01_snap8_zw1_ssthr5e-02_512_4096_8.NHI.dat',
		 'dndzdx_HM01_snap14_zw1_ssthr1e-03_512_4096_8.NHI.dat']
labels = ['z=5 gammax10', 'z=2.2 gammax10']
ns = len(fdats)
fig, ax = plt.subplots(1, figsize=(8, 6))
i = 0
for snap, _fdat, l in zip(['11','11', '11', '11'], fdats, labels):
	s = str(snap)
	fdat = '../../UVB/%s' % _fdat
	mat = np.loadtxt(fdat).T
	mats[s] = mat

	fsize = 14

	Mpc2cm = 3.086e24
	s = str(snap)

	lnhis, dndz, dndx, cddf, sb = mats[s]
	llspos = lnhis >= 17.5
	dlapos = lnhis >= 20
	# _dndx = np.sum(10 ** cddf[llspos] * (10 ** (lnhis[llspos] + .1) - 10 ** lnhis[llspos]))
	# _dndz = _dndx * (1 + redshifts[s]) ** 2 / np.sqrt(params.omega_m * (1 + redshifts[s]) ** 3 - params.omega_l)
	_dndz = np.nansum(dndz[llspos])
	_dndx = np.nansum(dndx[llspos])
	dndzcum = np.array([np.sum(dndz[k:]) for k in range(len(dndz))])
	llsfit = inter.interp1d(dndzcum, lnhis)
	lcool = llsfit(params.dndz[s])

	print 'snap', s, l, 'dndx', _dndx, 'dndz', _dndz, 'NHI match', lcool

	plt.plot(lnhis, cddf, label=l, color=zcolors[len(fdats) - i - 1], linestyle='--',
			 dashes=(3, len(fdats) - i))

	i += 1

if 1:
	if 0:
		lnhi = np.array([19., 19.5, 20.3, 20.65, 21., 21.35, 21.7])

		f = np.array([-20.5, -21.2, -21.9, -22.5, -23., -23.8])
		fmin = np.array([-20.4, -21.1, -21.9, -22.4, -22.9, -23.6])
		fmax = np.array([-20.8, -21.4, -22, -22.6, -23.2, -24.1])
		x = (lnhi[1:] + lnhi[:-1]) / 2.
		xerr = np.array([(lnhi[i] - lnhi[i + 1]) / 2 for i in range(len(x))])
		yerr = fmax - fmin
		plt.errorbar(x, f, xerr=xerr, yerr=yerr, label=r'Peroux+05 $1.78<z<3.5$', fmt="none")

		f = np.array([-20.2, -21.1, -21.6, -22.5, -23., -23.8])
		fmin = np.array([-20.4, -21.1, -21.9, -22.4, -22.9, -23.6])
		fmax = np.array([-20.8, -21.4, -22, -22.6, -23.2, -24.1])

		yerr = fmax - fmin
		plt.errorbar(x, f, xerr=xerr, yerr=yerr, label=r'Peroux+05 $3.5<z<5.0$', fmt="none")

		x = np.arange(13.2, 17.3, .5)
		y = 10.322 - 1.65 * x
		plt.errorbar(x, y, xerr=.25, yerr=.017 * x, label=r'Rudie+13 $2.0<z<2.8$', fmt="none")

	x = np.array(
		[12.11, 12.21, 12.38, 12.63, 12.84, 13.03, 13.28, 13.49, 13.72, 13.93, 14.18, 14.43, 14.7, 14.91, 15.05, 15.26,
		 15.47, 15.74, 15.97, 16.2, 16.41, 16.54, 16.76, 16.95, 17.08, 17.25, 17.47, 17.79, 18.08, 18.35, 18.6, 18.92,
		 19.21, 19.5, 19.79, 20.08, 20.34, 20.56, 20.79, 21.04, 21.23, 21.46, 21.53, 21.67, 21.73, 21.82])
	y = np.array(
		[-10.45, -10.32, -10.36, -10.53, -10.79, -11.04, -11.34, -11.68, -12.11, -12.41, -12.75, -13.09, -13.47, -13.9,
		 -14.24, -14.67, -15.05, -15.43, -15.82, -16.12, -16.41, -16.84, -17.18, -17.4, -17.57, -17.78, -18.08, -18.46,
		 -18.89, -19.23, -19.61, -19.82, -20.25, -20.54, -20.93, -21.31, -21.74, -22.08, -22.46, -22.93, -23.4, -24.09,
		 -24.77, -25.33, -26.01, -26.57])

	if 0:#nastasha
		plt.plot(x, y, label=r'Nastasha z=3.5 with SS')
		dndx = np.sum(10 ** y[26:-1] * (10 ** x[27:] - 10 ** x[26:-1]))
		dndz = dndx * 4.5 ** 2 / np.sqrt(params.omega_m * 4.5 ** 3 - params.omega_l)
		print 'nastasha z=3.5 dndx SS', dndx, 'dndz', dndz

		x = np.array(
			[12.01, 12.34, 12.74, 13.03, 13.32, 13.56, 13.83, 14.1, 14.39, 14.71, 15.11, 15.48, 15.96, 16.45, 16.64, 16.89,
			 17.16, 17.37, 17.58, 17.81, 18.04, 18.25, 18.52, 18.81, 19.11, 19.42, 19.69, 19.92, 20.2, 20.47, 20.76, 20.99,
			 21.22, 21.37, 21.64, 21.87, 22.04, 22.21])

		y = np.array(
			[-10.53, -10.41, -10.66, -11, -11.39, -11.73, -12.16, -12.59, -12.97, -13.53, -14.34, -14.98, -15.79, -16.56,
			 -16.95, -17.2, -17.72, -18.02, -18.32, -18.74, -19, -19.38, -19.68, -20.15, -20.58, -21.01, -21.52, -21.91,
			 -22.25, -22.63, -23.15, -23.7, -24.34, -25.03, -25.58, -26.18, -26.52, -27])

		plt.plot(x, y, label=r'Nastasha z=3.5 no SS')
		dndx = np.sum(10 ** y[26:-1] * (10 ** x[27:] - 10 ** x[26:-1]))
		dndz = dndx * 4.5 ** 2 / np.sqrt(params.omega_m * 4.5 ** 3 - params.omega_l)
		print 'nastasha z=3.5 dndx no SS', dndx, 'dndz', dndz

	if 1:#Ali
		x = np.array(
			[15.76, 15.92, 16.1, 16.3, 16.48, 16.67, 16.85, 17.03, 17.2, 17.44, 17.57, 17.74, 17.95, 18.15, 18.39, 18.57,
			 18.78, 19.01, 19.23, 19.46, 19.72, 19.92, 20.13, 20.29, 20.44, 20.64, 20.87, 20.99, 21.16, 21.3, 21.47, 21.6,
			 21.68, 21.82, 21.9, 22, 22.12, 22.18, 22.3])
		y = np.array(
			[-15.01, -15.33, -15.69, -16.05, -16.37, -16.74, -17.08, -17.34, -17.58, -17.92, -18.18, -18.33, -18.65,
			 -18.91, -19.19, -19.45, -19.69, -19.95, -20.22, -20.52, -20.9, -21.18, -21.42, -21.65, -21.91, -22.17, -22.49,
			 -22.77, -23.03, -23.31, -23.72, -24.1, -24.42, -24.84, -25.14, -25.48, -25.81, -26.07, -26.43])
		plt.plot(x, y, label=r'Ali z=5')
		dndx = np.sum(10 ** y[9:-1] * (10 ** x[10:] - 10 ** x[9:-1]))
		dndz = dndx * 6 ** 2 / np.sqrt(params.omega_m * 6 ** 3 - params.omega_l)
		print 'Ali z=5 dndx', dndx, 'dndz', dndz

		x = np.array(
			[15.27, 15.48, 15.73, 16.03, 16.43, 16.65, 16.95, 17.21, 17.46, 17.61, 17.81, 17.99, 18.18, 18.33, 18.57, 18.77,
			 19.04, 19.25, 19.43, 19.62, 19.83, 19.99, 20.07, 20.21, 20.41, 20.64, 20.93, 21.07, 21.26, 21.42, 21.55, 21.65,
			 21.75, 21.92, 22, 22.12, 22.22, 22.35])
		y = np.array(
			[-15.27, -15.71, -16.07, -16.55, -17.14, -17.54, -17.88, -18.3, -18.65, -18.89, -19.09, -19.27, -19.59, -19.79,
			 -20.03, -20.28, -20.58, -20.78, -20.94, -21.18, -21.38, -21.54, -21.58, -21.77, -21.99, -22.27, -22.59, -22.89,
			 -23.29, -23.68, -24.04, -24.36, -24.72, -25.14, -25.62, -25.89, -26.27, -26.59])
		plt.plot(x, y, label=r'Ali z=2')
		dndx = np.sum(10 ** y[8:-1] * (10 ** x[9:] - 10 ** x[8:-1]))
		dndz = dndx * 3 ** 2 / np.sqrt(params.omega_m * 3 ** 3 - params.omega_l)
		print 'Ali z=2 dndx', dndx, 'dndz', dndz




plt.xlim(15.5, 23.)
plt.ylim(-26.5, -14.5)
plt.xticks([16, 17, 18, 19, 20, 21, 22])
plt.yticks([-26, -24, -22, -20, -18, -16])
plt.xlabel(r'log$_{10}\mathrm{(N_{HI})[cm^{-2}]}$')
plt.ylabel(r'log$_{10}[\mathrm{f(N_{HI}, X)}]$')
h, l = ax.get_legend_handles_labels()
plt.legend(h, l)
leg1 = ax.legend(h[:ns], l[:ns], loc='upper right', title='EAGLE')
leg2 = ax.legend(h[ns:], l[ns:], loc='lower left')
ax.add_artist(leg1)

plt.savefig('../../Figures/cddf_%d_%d_%d_noSS.pdf' % (minres, maxres, npref))
plt.close()
