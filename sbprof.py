#!/usr/bin/python

__author__ = 'gallegos'
from sys import argv

# import eagleSqlTools
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
from matplotlib.cm import get_cmap

import params
from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['UDF', 'mosaic'], str)#HDFS', '
snaps = rdarg(argv, 'snap', list, [8, 9, 10, 11, 12], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
stack = rdarg(argv, 'stack', bool, False)
corr = rdarg(argv, 'corr', bool, False)
scorr = '.corr' * corr
csub = rdarg(argv, 'csub', bool, True)
scsub = '.csub'*csub
halfplot = rdarg(argv, 'halfplot', bool, False)
extraname = rdarg(argv, 'extraname', str, '')#'HM12')#
prename = rdarg(argv, 'prename', str, '')#'HM12')#
galcov = rdarg(argv, 'galcov', bool, False)
histover = rdarg(argv, 'histover', bool, False)
LLScubes = rdarg(argv, 'LLScubes', bool, False)
sql = rdarg(argv, 'sql', bool, False)
sphericalLLS = rdarg(argv, 'sphericalLLS', bool, False)
pairLLS = rdarg(argv, 'pairLLS', bool, False)
circlehist = rdarg(argv, 'circlehist', bool, False)
h2d = rdarg(argv, 'h2d', bool, False)
kde = rdarg(argv, 'kde', bool, False)
nhiprof = rdarg(argv, 'nhi', bool, False)
nconf = rdarg(argv, 'conf', int, 2)
minres = rdarg(argv, 'minres', int, 512)
maxres = rdarg(argv, 'maxres', int, 4096)
npref = rdarg(argv, 'npref', int, 8)
galrad = rdarg(argv, 'rad', int, 6)
sbhist = rdarg(argv, 'sbhist', bool, False)
ssthr = rdarg(argv, 'ssthr', float, 1e10)#6.73E-3
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
type = rdarg(argv, 'type', str, 'LLS') #NHtot, f_NHI, NHII
types = rdarg(argv, 'type', list, ['NHI'], str) #NHtot, f_NHI, NHII
radprof = rdarg(argv, 'radprof', bool, False)
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, True)
mask2d = rdarg(argv, 'mask2d', bool, True)
gmask = rdarg(argv, 'gmask', bool, False)
model = rdarg(argv, 'model', str, 'HM01')#'HM12')#
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)
parallel = rdarg(argv, 'parallel', bool, False)

gamma_bkg = {}
heat_bkg = {}
sgamma = {}
sheat = {}
gamma_bkg['HM12'] = params.gamma_bkg['HM12']
heat_bkg['HM12'] = params.heat_bkg['HM12']
gamma_bkg['HM01'] = params.gamma_bkg['HM01']
heat_bkg['HM01'] = params.heat_bkg['HM01']
snames = params.snames
redshifts = params.redshifts
asec2kpcs = params.asec2kpcs
dz = params.dz
sbpeaks = params.sbpeaks
zlens = params.zlens
dndz = params.dndz
nhi_fit = {}
nhi_fit['HM01'] = params.nhi_fit['HM01']
nhi_fit['HM12'] = params.nhi_fit['HM12']
fLLScorr = params.fLLScorr


lcube = maxres # 4096
coml = 25  # cMpc
com2pix = 163.84  # lcube/coml
kpc2pix = lcube / float(coml * 1e3)
rads = np.array([0, 2, 4, 8, 12, 20, 30, 50, 100, 200])

lognhi = [13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21]
nhi = np.power(10, lognhi)
sb = [0.0001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96, 0.98, 1, 1, 1]

nhi2sb = inter.interp1d(nhi, sb)
sb2nhi = inter.interp1d(sb, nhi)

def sclipping(fits, nsigma, dim=(0), mask=None):
	# print 'Asign nan to values > nsigma in a fits array'
	for i in range(len(fits)):
		fits[i, :, np.where(mask)[i]] = np.nan
	stds = np.nanstd(fits, dim)
	high_sigma = np.abs(fits) > nsigma * stds
	fits[high_sigma] = np.nan
	return fits

folder = '../../'
flaes = '../../'

scat = ''
for i in fitcat: scat += '_' + i
ncores = 10
nrand = 10
zw0 = rdarg(argv, 'zw', int, 2)
zoff = rdarg(argv, 'zoff', int, 0)
zw = 2*zw0+1
do_eagle = 1
ov_eagle = False
vcorr = 1
do_mask = mask
do_mask2d = mask2d
do_gmask = gmask
smask = ('.%smask' % (('g%d' % galrad)*do_gmask))*do_mask
do_sclip = 1
nsigma = 3
sclip = ('.sclip%d' % nsigma)*do_sclip
do_rand = 0
do_sbcorr = 0
fix_red = 0
SB = 1.27  # SB at z=3.5

fontsize = 14
figsize = (8, 5)
vlim = 2500

cmap = get_cmap('gist_rainbow')

#colors = [(0, 0, 1), (1, 1, 0), (1, 0, 0)]
#cmap = LinearSegmentedColormap.from_list('sofi', colors, N=50)
# 1 spaxel is (0.2 arcsec)^2
# flux2sb comes from 1.25 Angstrom wavelength width times 1/(0.2)**2 from flux2sbersion to the correct aperture
flux2sb = 31.25



zreal = [3.70]
zprob = [[3.65, 3.75]]
sncat = ['_18_45']

zreal = [3.26]
zprob = [[2.9, 3.6]]
sncat = ['_46_214']

zreal = [3.70]
zprob = [[3.65, 3.75]]
sncat = ['_18_45']

zreal = [3.35]
zprob = [[2.9, 3.9]]
sncat = ['_59_263']


zprob = [[2.9, 3.4]]
zreal = [3.15]
sncat = ['_33_140']



zprob = [[2.9, 3.4], [3.4, 4.5], [4.5, 5.5]]
zreal = [3.15, 3.89, 4.91]
sncat = ['_33_140', '_46_242', '_43_162']



zprob = [[3.4, 4.5]]
zreal = [3.89]
sncat = ['_46_242']

extraname += prename

zlm, ylm, xlm = 201, 201, 201
zrange = np.arange(zlm)
w = ylm / 2
zz, yy, xx = np.ogrid[0:zlm, 0:ylm, 0:xlm]
cm = ((yy - w) ** 2 + (xx - w) ** 2)[0]
mat = {}
zrealp = {}
ncats = {}
fsize = 14
colors2 = 'm', 'darkorange', 'g', 'blue', 'black', 'purple', 'red', 'orange', 'gray', 'gold', 'purple', 'brown', 'pink'
nn = 5
nnn = zlm/nn

col = ['blue', 'green', 'red']



for i in range(len(zprob)):
	fig, ax = plt.subplots(figsize=figsize)
	_zp = zprob[i]
	zp = (_zp[0]+_zp[1])/2.
	zrealp['%.2f' % zp] = zp
	_zreal = zreal[i]


	wred = _zp[1]-_zp[0]
	print 'Probing %.1f < z < %.1f' % (_zp[0], _zp[1])
	fdat = 'muse-vs-eagle%s_z%.2f.dat' % (scat, zp)
	_fin = None

	asec2pix_muse = 5  # muse .2 arcsec per pixel
	
	odat = '../../UVB/SNR/UVB_zr%.2f%s%s%s%s%s_z%.1f-%.1f_ngal%s%s.dat' \
	       % (_zreal, scat, scsub, scorr, smask, sclip, _zp[0], _zp[1], sncat[i], extraname)


	mat['%.2f' % zp] = np.loadtxt(odat)
	r, r0, r1, zwi, z, sb, std, sbmean, fc, flls, flls_uplim, gammac, gamma_uplimc, gmean, snr, sbr, fe, gamma, gamma_uplim, sbhm12, v = mat['%.2f' % zp].T
	#vmin, vmax = 2.99792458e5*zmin*1.25/1215.67/(1+_zreal), 2.99792458e5*zmax*1.25/1215.67/(1+_zreal)
	vmin, vmax = -vlim, vlim
	v0 = v[(zwi==0) & (z==zoff-zw0)][0]
	v1 = v[(zwi==0) & (z==zoff+zw0)][0]
	c = (zwi==zw0) & (z==zoff)

	ignore = [8, 10], [10, 12], [12, 16], [6, 10], [6, 12], [12, 18], [10, 16], [6, 20]

	for ig in ignore: c &= ~((r0 == ig[0]) & (r1 == ig[1]))

	sbm = (sbmean[c]*zwi[c])[0]
	print 'Velocity window %.2f +%.2f dv %.2f km/s' % (v0, v1, v1-v0)
	print 'Base SB for this window', sbm
	x, y, xerr, yerr = r[c], sb[c]+sbm, [r[c]-r0[c], r1[c]-r[c]], [std[c], std[c]]
	y[y<=0] = 1e-3
	ax.scatter(x, y, label=r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}}$', color='green')
	ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='none', color='green')
	y = fc[c]
	ax.scatter(x, y, label=r'$\mathrm{f_{LLS}}$', color='blue')
	gm = (gmean[c]*zwi[c])[0]
	y, yerr = gamma[c]+gm, [gamma_uplim[c]/2., gamma_uplim[c]/2.]
	y[y <= 0] = 1e-3
	ax.scatter(x, y, label=r'$\mathrm{\Gamma_{HI}\,[10^{-12}\,s^{-1}}]}$', color='gold', zorder=0)
	ax.errorbar(x, y, yerr=yerr, fmt='none', color='gold')
	plt.xlabel(r'distance [arcsec]', fontsize=fsize)
	
	w=1./gamma_uplim[c][5:]
	wsum=np.sum(w)
	_gmean = np.nansum(gamma[c][5:]*w)/wsum
	gstd = np.std(gamma[c][5:])
	for r0_, r1_, rr, g, gup in zip(r0[c], r1[c], r[c], gamma[c], gamma_uplim[c]/2.):
		print 'r%d-%d %.2f G %.2f +- %.2f' % (r0_, r1_, rr, g+gm, gup)
	print 'Base G %.2e, Mean weighted G %.2f +- %.2f' % (gm, _gmean+gm, gstd)
	
	#plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}}$', fontsize=fsize)
	plt.xlim([0, 20])
	plt.ylim([5e-2, 6e1])
	plt.xticks(fontsize=fsize-2)
	plt.yticks(fontsize=fsize-2)
	plt.title('z=%.1f' % _zreal, fontsize=fsize)

	plt.grid()
	plt.legend(prop={"size": fsize})
	#plt.tight_layout()
	plt.semilogy()
	# ax2 = ax.twinx()
	# ax2.plot(v[zmin: zmax], f/fstd, alpha=0)
	# plt.ylabel(r'SNR')
	plt.savefig('../../Figures/SB_profile%s_redshift%d_zw%d_zoff%d%s%s.pdf' % (scat, _zreal * 10, zw, zoff, smask, extraname))
	plt.close()


sb620 = []
std620 = []
#colors = 'blue', 'green', 'gold', 'red', 'purple', 'orange', 'black', 'brown', 'pink', 'yellow'

if 0:#for rr, ylim in zip(rsel2, ylims):
	fig, ax = plt.subplots(figsize=figsize)
	_zw = 0
	i = 0
	sred = ''
	for zp in zprob:
		zpm = (zp[0]+zp[1])/2.
		_zp = zrealp['%.2f' % zpm]
		ncat = ncats['%.2f' % zpm]
		sred += '_%d' % (_zp*10)
		r, r0, r1, zwi, z, sb, std, fc, flls, flls_uplim, gammac, gamma_uplimc, snr, sbr, fe, gamma, gamma_uplim, sbhm12, v = mat['%.2f' % zpm].T
		c = (zwi == _zw)&(r0 == rr[0])&(r1 == rr[1])&(np.abs(v) <= vlim)
		ax.plot(v[c], sp.savgol_filter(sb[c]/1.25, 5, 1), label=r'$z=%.1f$' % _zp, color=colors[i])
		cout = (zwi == _zw) & (r0 == rr[0]) & (r1 == rr[1]) & (np.abs(v) > 1000)
		noise = np.std(sb[cout]/1.25)
		print noise
		ax.errorbar(-vlim*(.9-i*.1), ylim[1]*.5, yerr=noise, color=colors[i], capsize=3)
		#ax.plot(v[c], sb[c], label=r'$z=%.1f$' % _zp, color=colors[i])
		i += 1
		# ax.plot(v[zmin: zmax], [2*fstd*flux2sb]*len(f), label=r'2$\sigma$ noise level')
		#ax.plot(z[c]*1.25, [0] * len(r[c]), color='gray')
		# ax.scatter([0], [1], color='red')
		cool = (z==-2) & (zwi==2) & (r0==rr[0]) & (r1==rr[1])
		print 'z', zp, '%.2f' % _zp, '%d<r<%d' % (rr[0], rr[1]), ncat
		print 'SB', sb[cool], 'std', std[cool], '\nG', gamma[cool], 'corr', gammac[cool], 'uplim', gamma_uplim[cool],\
			'corr', gamma_uplimc[cool], '\nfLLS_eagle', fe[cool], \
			'corr', fc[cool], 'fLLS_muse', flls[cool], 'uplim', flls_uplim[cool]
		
	plt.xlabel(r'v [$\mathrm{km\,s^{-1}}$]', fontsize=fsize)
	plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}\AA^{-1}]}}$', fontsize=fsize)
	plt.xlim([-vlim, vlim])
	plt.ylim(ylim)
	plt.xticks(fontsize=fsize-2)
	plt.yticks(fontsize=fsize-2)
	plt.legend(prop={"size":fsize})
	if rr[0] < .3: plt.title('r<%.1f"' % rr[1])
	elif rr[0] < 1: plt.title('%.1f"<r<%.1f"' % (rr[0], rr[1]))
	else: plt.title('%d"<r<%d"' % (rr[0], rr[1]), fontsize=fsize)
	plt.grid()
	plt.tight_layout()
	plt.savefig('../../Figures/SB_spectra%s%s_r%d-%d_zw%d%s%s.pdf' % (
		sred, scat, rr[0], rr[1], _zw, smask, extraname))
	# plt.show()
	plt.close()


