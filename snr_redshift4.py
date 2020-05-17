#!/usr/bin/python

__author__ = 'gallegos'
import glob
# import eagleSqlTools
import matplotlib.pyplot as plt
import os
import scipy.interpolate as inter
from matplotlib.cm import get_cmap
from pyfits import getdata, PrimaryHDU
from sys import argv

import numpy as np

from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
snaps = rdarg(argv, 'snap', list, [12], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
halfplot = rdarg(argv, 'halfplot', bool, False)
extraname = rdarg(argv, 'extraname', str, '')#'HM12')#
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
minres = rdarg(argv, 'minres', int, 512)
maxres = rdarg(argv, 'maxres', int, 4096)
npref = rdarg(argv, 'npref', int, 8)
sbhist = rdarg(argv, 'sbhist', bool, False)
ssthr = rdarg(argv, 'ssthr', float, 1e10)#6.73E-3
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
types = rdarg(argv, 'type', list, ['NHI'], str) #NHtot, f_NHI, NHII
radprof = rdarg(argv, 'radprof', bool, False)
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)
model = rdarg(argv, 'model', str, 'HM01')#'HM12')#
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)
parallel = rdarg(argv, 'parallel', bool, False)

snames = {'10': '010_z003p984', '11': '011_z003p528', '12': '012_z003p017', '13': '013_z002p478', '14': '014_z002p237',
		  '15': '015_z002p012'}
reds = {'10': 984, '11': 528, '12': 17}
redshifts = {'10': 3.984, '11': 3.528, '12': 3.017, '13': 2.478, '14': 2.237, '15': 2.012}
asec2kpcs = {'10': 7.842, '11': 7.449, '12': 7.108, '13': 8.241, '14': 8.396, '15': 8.516}

gamma_bkg = {}
heat_bkg = {}
sgamma = {}
sheat = {}
gamma_bkg['HM12'] = {'10': [5.71E-13, 3.28E-13, 2.31E-16], '11': [6.77E-13, 3.89E-13, 8.70E-16],
			 '12': [7.66E-13, 4.42E-13, 2.10E-15], '13': [9.50E-13, 5.55E-13, 9.06E-15],
			 '14': [9.64E-13, 5.67E-13, 1.13E-14]}

heat_bkg['HM12'] = {'10': [2.27E-12, 2.18E-12, 7.24E-15],
			'11': [2.68E-12, 2.62E-12, 2.33E-14],
			'12': [3.02E-12, 3.05E-12, 5.01E-14],
			'13': [3.75E-12, 4.22E-12, 1.78E-13],
			'14': [3.81E-12, 4.42E-12, 2.18E-13]}


gamma_bkg['HM01'] = {'10': [.9E-12, 3.28E-13, 2.31E-16],
					 '11': [1E-12, 3.89E-13, 8.70E-16],
					 '12': [1.5E-12, 4.42E-13, 2.10E-15],
					 '13': [1.7E-12, 5.55E-13, 9.06E-15],
					 '14': [1.9E-12, 5.67E-13, 1.13E-14]}

heat_bkg['HM01'] = {'10': [3E-12, 2.18E-12, 7.24E-15],
					'11': [4E-12, 2.62E-12, 2.33E-14],
					'12': [5E-12, 3.05E-12, 5.01E-14],
					'13': [5.5E-12, 4.22E-12, 1.78E-13],
					'14': [6E-12, 4.42E-12, 2.18E-13]}


dz = {'10': 0.0255, '11': .03, '12': .035}

sbpeaks = {'10': .739, '11': 1.293, '12': 2.579}  # for an aperture of 1 asec^2!!! I am converting to flux later on
zlens = {'10': 34, '11': 29, '12': 25, '13': 20, '14': 18, '15': 16}
dndz = {'10': 3., '11': 2.3, '12': 2., '13': 1.3, '14': 1., '15': .5} # check last ones, high extrapolation w/r to Prochaska+10
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

folder = '/net/abnoba/scratch2/gallegos/Research/MUSE/'
flaes = '/net/galaxy-data/export/galaxydata/gallegos/'

#folder = '../../'
#flaes = '../../'

scat = ''
for i in fitcat: scat += '_' + i
ncores = 10
nrand = 10
zw0 = rdarg(argv, 'zw', int, 2)
zw = 2*zw0+1
do_eagle = 1
ov_eagle = False
vcorr = 1
do_mask = 1
do_gmask = 1
smask = ('.%smask' % ('g'*do_gmask))*do_mask
do_sclip = 1
nsigma = 4
sclip = ('.sclip%d' % nsigma)*do_sclip
do_csub = False
scsub = '.csub'*do_csub
do_rand = 0
fix_red = 0
zoff = 0
SB = 1.27  # SB at z=3.5
type = 'LLS'  # 'NHI'

cmap = get_cmap('gist_rainbow')

#colors = [(0, 0, 1), (1, 1, 0), (1, 0, 0)]
#cmap = LinearSegmentedColormap.from_list('sofi', colors, N=50)
# 1 spaxel is (0.2 arcsec)^2
# flux2sb comes from 1.25 Angstrom wavelength width times 1/(0.2)**2 from flux2sbersion to the correct aperture
flux2sb = 31.25

fdat = 'muse-vs-eagle%s.dat' % scat
_fin = None

if do_eagle:

	do_emask = mask
	cat = getdata('%s/EAGLE/cats/gals_snap12.fits' % folder, 1)
	m = cat['DM_mass']
	s = '12'
	H0 = 69.3
	Ol = 0.714
	Om = 1 - 0.714
	red = redshifts[s]
	offset = cat['vx'] * (1 + red) / H0 / np.sqrt(Om * (1 + red) ** 3 + Ol)
	offpix = np.round(offset*zlens['12']/coml).astype(int)

	pnames = ['U', 'SFR', 'Mmean200']
	dolog = [0, 1, 1]
	prop = {}
	pmed = {}
	phist = {}
	bins = [[-22, -20, -19, -17, -16], [-3, -.8, -.4, -.2, 0, .5, 3],
	        [10.44616846, 10.78110581, 11.11604317, 11.45098053, 11.78591788]]
	for p, dl, nh in zip(pnames, dolog, bins):
		if p == 'U':
			p0 = cat[p]
			noU = np.isnan(p0)
			p1 = -2.5*np.log10(cat['SFR'])-18.6
			prop[p] = p0
			prop[p][noU] = p1[noU]
		else: prop[p] = cat[p]
		pmed[p] = np.nanmedian(prop[p])
		if dl: phist[p] = np.histogram(np.log10(prop[p]), nh)
		else: phist[p] = np.histogram(prop[p], nh)

	cename = '../../EAGLE/simplestacks/snap12_x%s%s.fits' % ('.masked'*do_emask, '.vcorr'*vcorr)
	if not os.path.isfile(cename) or overwrite and ov_eagle:
		stack = []
		mstack = []
		for i, off in zip(cat['ID'], offpix):
			print i, off
			_s = getdata(flaes + '/EAGLE/gals/snap12_x/%d.%s.fits' % (i, type))
			if vcorr: _s = np.roll(_s, off, 0)
			if mask:
				_m = getdata(flaes + '/EAGLE/gals/snap12_x/%d.mask.fits' % i)
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
			_mstack2 = np.nanmean(mstack > 0, 0)
			hdu.data = _mstack
			hdu.writeto(cename.replace('.fits', '.gals.fits'), clobber=True)
			hdu.data = _mstack2
			hdu.writeto(cename.replace('.fits', '.galcov.fits'), clobber=True)
			
		for p, dl, nh in zip(pnames, dolog, bins):
			if dl:
				_p = np.log10(prop[p])
				_pm = np.log10(pmed[p])
			else:
				_p = prop[p]
				_pm = pmed[p]
			high = _p > _pm
			shigh = np.nanmean(stack[high], 0)
			slow = np.nanmean(stack[~high], 0)
			hdu.data = shigh
			hdu.writeto(cename.replace('.fits', '.high-%s.fits' % p), clobber=True)
			hdu.data = slow
			hdu.writeto(cename.replace('.fits', '.low-%s.fits' % p), clobber=True)
			for i in range(len(nh)-1):
				cool = (_p > phist[p][1][i]) & (_p < phist[p][1][i+1])
				hdu.data = np.nanmean(stack[cool], 0)
				hdu.writeto(cename.replace('.fits', '.%s_%.2f_%.2f.fits' % (p, phist[p][1][i], phist[p][1][i+1])), clobber=True)

	rads = [[0, 2], [2, 4], [4, 6], [6, 10], [10, 16], [16, 24], [24, 30], [40, 50], [50, 60],
	        [60, 70], [70, 80], [80, 90], [90, 100]]
	nr = len(rads)
	asec2pix_eagle = .1836
	rpix_eagle = np.array(rads) / asec2pix_eagle
	
	stack = getdata(cename)
	zle, yle, xle = stack.shape
	yy, xx = np.ogrid[0: yle, 0: xle]
	ce = (xx - xle / 2.) ** 2 + (yy - yle / 2.) ** 2
	
	y = {}
	label = ['U', 'log SFR', r'log M$_{200}$']
	colors = ['red', 'orange', 'green', 'pink', 'blue', 'purple']
	xl = [.9, 100]
	yl = [5e-3, .3]
	
	for p, l, nh in zip(pnames, label, bins):
		plt.figure()
		ph = phist[p][1]
		_p = (ph[1:]+ph[:-1])/2.
		x = [(r[0]+r[1])/2. for r in rads]
		for i in range(len(nh)-1):
			y['%s_%d' % (p, _p[i])] = []
			st = getdata(cename.replace('.fits', '.%s_%.2f_%.2f.fits' % (p, ph[i], ph[i + 1])))
			if np.nansum(st) > 0:
				for r in rads:
					cool = (ce >= r[0] ** 2) & (ce < r[1] ** 2)
					_y = np.nanmean(st[zle/2-zw0:zle/2+zw0+1, cool])
					if np.isnan(_y): _y = None
					y['%s_%d' % (p, _p[i])].append(_y)
				if p == 'Mmean200': _l = r'$%.1f\times10^{%d}\,\rm{M_\odot}$' % (10 ** _p[i] / 10. ** int(_p[i]), int(_p[i]))
				if p == 'SFR': _l = r'$%.1f\times10^{%d}\,\rm{M_\odot\,yr^{-1}}}$' % (10 ** _p[i] / 10. ** int(_p[i]), int(_p[i]))
				if p == 'U': _l = r'$\rm{M_U}=%.1f$' % (_p[i])
				plt.plot(x, y['%s_%d' % (p, _p[i])], label=_l, color=colors[i])
				plt.scatter(x, y['%s_%d' % (p, _p[i])], color=colors[i])
		plt.loglog()
		plt.xlim(xl)
		plt.ylim(yl)
		plt.xlabel('distance [arcsec]')
		plt.ylabel(r'f$_{\rm{LLS}}$')
		plt.legend()
		plt.savefig('%s/Figures/LLSprof-%s.jpg' % (folder, p))
		plt.close()

		
	feagle = getdata(cename)
	fhigh, flow = {},{}
	for p in pnames:
		fhigh[p] = getdata(cename.replace('.fits', '.high-%s.fits' % p))
		flow[p] = getdata(cename.replace('.fits', '.low-%s.fits' % p))
	zl, yl, xl = feagle.shape
	zle0 = zl
	zle, yle, xle = feagle.shape
	asec2pix_eagle = .1836  # average btw snapshots

asec2pix_muse = .2  # muse
fall = []
fouts = []
n5d = []
try:
	glob.os.makedirs('%s/all/simplestacks/' % folder)
except:
	pass

wmax = 201/2
fall = []
mall = []
gall = []
_id = []
fall2 = []
frand = {}
mrand = {}
grand = {}

for nr in range(nrand): frand['%d' % nr], mrand['%d' % nr], grand['%d' % nr] = [[], [], []]

for _fitcat in fitcat:
	cat = getdata('%s/%s/cats/laes.fits' % (folder, _fitcat), 1)
	ids = cat['ID']
	zs = cat['redshift']
	ls = 1215.67 * (1 + zs)
	sconf = cat['sconf']
	#spec = cat['spec']
	cool = (sconf >= 2) & (zs < 4)
	corrupt = 0
	for i, l in zip(ids[cool], ls):
		ffits = flaes + '/%s/LAEs/%d.fits' % (_fitcat, i)
		if os.path.isfile(ffits):
			if do_rand:
				for nr in range(nrand):
					fr = getdata(ffits.replace('.fits', '.rand%d.fits' % nr))
					zr, yr, xr = fr.shape
					_fr = fr[:, yr/2-wmax: yr/2+wmax+1, xr/2-wmax: xr/2+wmax+1]
					frand['%d' % nr].append(_fr)
					if do_mask:
						mr = getdata(ffits.replace('.fits', '.rand.mask%d.fits' % nr))
						_mr = mr[yr / 2 - wmax: yr / 2 + wmax + 1, xr / 2 - wmax: xr / 2 + wmax + 1]
						mrand['%d' % nr].append(_mr)
					if do_gmask:
						gr = getdata(ffits.replace('.fits', '.rand.gmask%d.fits' % nr))
						_gr = gr[:, yr/2-wmax: yr/2+wmax+1, xr/2-wmax: xr/2+wmax+1]
						grand['%d' % nr].append(_gr)

			fit = getdata(ffits)
			zlm, ylm, xlm = fit.shape
			zlm0 = zlm
			_fit = fit[:, ylm/2-wmax: ylm/2+wmax+1, xlm/2-wmax: xlm/2+wmax+1]
			zrange = np.arange(zlm)

			_fit = np.array(_fit)
			fall.append(_fit)
			if do_mask:
				mit = getdata(ffits.replace('.fits', '.mask.fits'))
				_mit = mit[ylm/2-wmax: ylm/2+wmax+1, xlm/2-wmax: xlm/2+wmax+1]
				mall.append(_mit)
			if do_gmask:
				git = getdata(ffits.replace('.fits', '.gmask.fits'))
				_git = git[:, ylm/2-wmax: ylm/2+wmax+1, xlm/2-wmax: xlm/2+wmax+1]
				_git = np.array(_git)
				gall.append(_git)
				_id.append(i)

fall = np.array(fall)

if do_mask:
	bad = np.array(mall) > 0
	if do_gmask: gall = np.array(gall)
	for i in range(len(fall)):
		if do_gmask:
			bad2 = (gall[i] > 0) & (gall[i] != _id[i])
			fall[i, bad2] = np.nan
		fall[i, :, bad[i]] = np.nan

if do_sclip:
	stds = np.nanstd(fall, 0)
	high_sigma = np.abs(fall) > nsigma * stds
	fall[high_sigma] = np.nan
fmuse = np.nanmean(fall, 0)
stackname = '%s/all/simplestacks/stack%s%s%s%s%s.fits' % (folder, scsub, scat, smask, sclip, extraname)
hdu = PrimaryHDU()
if not os.path.isfile(stackname) or overwrite:
	hdu.data = fmuse
	hdu.writeto(stackname, clobber=True)

zlm, ylm, xlm = fmuse.shape

if do_rand:
	for nr in range(nrand):
		frand['%d' % nr] = np.array(frand['%d' % nr])
		if do_mask:
			for i in range(len(frand['%d' % nr])):
				bad = np.array(mrand['%d' % nr][i] > 0)
				frand['%d' % nr][i, :, bad] = np.nan
				if do_gmask:
					bad2 = (grand['%d' % nr][i] > 0)
					frand['%d' % nr][i, bad2] = np.nan
		if do_sclip:
			stds = np.nanstd(frand['%d' % nr], 0)
			high_sigma = np.abs(frand['%d' % nr]) > nsigma * stds
			frand['%d' % nr][high_sigma] = np.nan
		fr_mean = np.nanmean(frand['%d' % nr], 0)
		zlr, ylr, xlr = fr_mean.shape

		randstackname = '../../all/simplestacks/randstack%s%s%s%s.%d.fits' % (scat, smask, sclip, extraname, nr)
		if not os.path.isfile(randstackname) or overwrite:
			hdu.data = fr_mean
			hdu.writeto(randstackname, clobber=True)

if do_eagle:
	rpix_eagle = np.array(rads) / asec2pix_eagle
	y, x = np.ogrid[0: yle, 0: xle]
	ce = (x - xle / 2.) ** 2 + (y - yle / 2.) ** 2

ym, xm = np.ogrid[0: ylm, 0: xlm]
cm = (xm - xlm / 2.) ** 2 + (ym - ylm / 2.) ** 2
if do_rand:
	yr, xr = np.ogrid[0: ylr, 0: xlr]
	cr = (xr - xlr / 2.) ** 2 + (yr - ylr / 2.) ** 2


rw = 3
extraname = ''#'.rw%d_zoff-5_5' % rw
fout = open('UVB%s%s%s%s.dat' % (scsub, smask, sclip, extraname), 'w')

if do_eagle: fout.write('#r0 r1 zw zoff SB SB_std SB_rand fLLS Gamma Gamma_std\n')
else: fout.write('#r0 r1 zw zoff SB SB_std SB_rand\n')

r0 = np.arange(0, 20, 1)
r1 = np.arange(2, 28, 1)
zoffs = np.arange(-50, 51)#[1, 0, -1, -2]
zws = [0, 1, 2, 3]
g = {}
gstd = {}
sbstd = {}
sb = {}
sbrand = {}
_nm = fall.shape[0]

for i in r0:
	for j in r1:
		if (j-i) >= rw:
			rm0 = i / asec2pix_muse
			rm1 = j / asec2pix_muse
			if do_eagle:
				re0 = i / asec2pix_eagle
				re1 = j / asec2pix_eagle
			print 'Between %.1f and %.1f arcsec' % (i, j)
			inside_muse = (cm >= rm0 ** 2) & (cm < rm1 ** 2)
			nin = np.sum(inside_muse)
			for _zw in zws:
				if do_eagle:
					inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
					fe = np.nanmean(np.nansum(feagle[zle/2-_zw: zle/2+_zw+1, inside_eagle], 0))
				if do_rand:
					inside_rand = (cr >= rm0 ** 2) & (cr < rm1 ** 2)
					fr = fr_mean[zlr/2-_zw: zlr/2+_zw+1, inside_rand]
					fr = np.nanmean(np.nansum(fr, 1))
					fstd = []
					for nr in range(nrand):
						fr2 = frand['%d' % nr][:, zlr / 2 - _zw: zlr / 2 + _zw + 1, inside_rand]
						fstd.append(np.nanmean(np.nansum(np.nanmean(fr2, 2), 1)))
					fstd = np.nanstd(fstd)
				else:
					fm = [fmuse[zz: zz+zw+1, inside_muse] for zz in [0, zw, 2*zw, zlm-3*zw, zlm-2*zw, zlm-zw]]
					_fr = [np.nansum(np.nanmean(_fm, 1)) for _fm in fm]
					fr = np.nanmean(_fr)
					fstd = np.nanstd(_fr)

				for _zoff in zoffs:
					print 'zw %d zoff %d' % (_zw*2+1, _zoff)
					zmin = max(0, zlm/2-_zw+_zoff)
					zmax = min(zlm, zlm/2+_zw+_zoff+1)
					_fmuse = fmuse[zmin: zmax, inside_muse]
					fin = np.nansum(np.nanmean(_fmuse, 1))
					if do_eagle:
						fy = (fin)*flux2sb/fe
						sb2gamma = .519 #expected SB 1.31e-20 for a gamma of .684e-12 at z=3.5
						gamma = fy*sb2gamma
						gamma_std = fstd * flux2sb / fe * sb2gamma
						fout.write('%.3f %.3f %d %d %.3f %.3f %.3f %.3f %.3f %.3f\n' %
								   (i, j, _zw, _zoff, fin*flux2sb, fstd*flux2sb, fr*flux2sb, fe, gamma, gamma_std))

						print 'SB %.3f, random %.3f. 2sigma noise %.3f. fLLS %.3f. SB EAGLE HM12 %.3f. SBUVB %.3f. Gamma %.3f std %.3f' % \
								  (fin*flux2sb, fr*flux2sb, 2*fstd*flux2sb, fe, fe*SB, fy, gamma, gamma_std)
					else:
						fout.write('%.3f %.3f %d %d %.3f %.3f %.3f\n' %
								   (i, j, _zw, _zoff, fin*flux2sb, fstd*flux2sb, fr*flux2sb))
						print 'SB %.3f, random %.3f. 2sigma noise %.3f.' % \
								  (fin*flux2sb, fr*flux2sb, 2*fstd*flux2sb)

fout.close()


if 0:
	# Selected values
	_r0, _r1, _zoff, _zw = [13, 19, -3, 2]#[11, 17, -2, 2]
	rm0 = _r0 / asec2pix_muse
	rm1 = _r1 / asec2pix_muse

	inside_muse = (cm >= rm0 ** 2) & (cm < rm1 ** 2)
	if do_eagle:
		re0 = _r0 / asec2pix_eagle
		re1 = _r1 / asec2pix_eagle
		inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
	zrange = np.arange(-50, 51)

	if do_rand:
		yr, xr = np.ogrid[0: ylr, 0: xlr]
		cr = (xr - xlr / 2.) ** 2 + (yr - ylr / 2.) ** 2
		inside_rand = (cr >= rm0 ** 2) & (cr < rm1 ** 2)
		fr = fr_mean[zlr / 2-_zw: zlr/2+_zw+1, inside_rand]
		fr = np.nanmean(fr)
		fstd = []
		for nr in range(nrand):
			fr2 = np.nansum(frand['%d' % nr][:, zlr/2-_zw: zlr/2+_zw+1, inside_rand], 1)
			fstd.append(np.nanmean(np.nanmean(fr2, 1)))
		_std = np.nanstd(fstd)
	else:
		fm = [fmuse[zz: zz + zw + 1, inside_muse] for zz in [0, zw, 2 * zw, zlm - 3 * zw, zlm - 2 * zw, zlm - zw]]
		_fr = [np.nansum(np.nanmean(_fm, 1)) for _fm in fm]
		fr = np.nanmean(_fr)
		fstd = np.nanstd(_fr)

	if do_eagle:
		inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
		fe = np.nanmean(np.nansum(feagle[zle/2-_zw: zle/2+_zw+1, inside_eagle], 0))

	f = np.nanmean(fmuse[:, inside_muse], 1)
	x = zrange*68.5
	sb = np.array(f-fr)*flux2sb
	sbstd = sb-sb + np.mean(fstd)*flux2sb
	sbrand = sb-sb + fr*flux2sb
	plt.figure(figsize=(10, 5))
	plt.xlabel('v [km/s]')
	plt.ylabel('SB')
	plt.plot(x, sb-sbrand,)
	plt.plot(x, 2*sbstd, label=r'$2\sigma$ level')
	plt.legend()
	plt.savefig('../../Figures/%d-%d_zw%d%s%s.jpg' % (_r0, _r1, _zw*2+1, smask, sclip))
	#plt.show()





if 0:
	# Selected values
	_r0, _r1, _zoff, _zw = [13, 19, -3, 2]#[11, 17, -2, 2]
	rm0 = _r0 / asec2pix_muse
	rm1 = _r1 / asec2pix_muse

	inside_muse = (cm >= rm0 ** 2) & (cm < rm1 ** 2)
	if do_eagle:
		re0 = _r0 / asec2pix_eagle
		re1 = _r1 / asec2pix_eagle
		inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
	zrange = np.arange(-50, 51)

	if do_rand:
		yr, xr = np.ogrid[0: ylr, 0: xlr]
		cr = (xr - xlr / 2.) ** 2 + (yr - ylr / 2.) ** 2
		inside_rand = (cr >= rm0 ** 2) & (cr < rm1 ** 2)
		fr = fr_mean[zlr / 2-_zw: zlr/2+_zw+1, inside_rand]
		fr = np.nanmean(fr)
		fstd = []
		for nr in range(nrand):
			fr2 = np.nansum(frand['%d' % nr][:, zlr/2-_zw: zlr/2+_zw+1, inside_rand], 1)
			fstd.append(np.nanmean(np.nanmean(fr2, 1)))
		_std = np.nanstd(fstd)
	else:
		fm = [fmuse[zz: zz + zw + 1, inside_muse] for zz in [0, zw, 2 * zw, zlm - 3 * zw, zlm - 2 * zw, zlm - zw]]
		_fr = [np.nanmean(np.nansum(_fm, 0)) for _fm in fm]
		fr = np.nanmean(_fr)
		fstd = np.nanstd(_fr)

	if do_eagle:
		inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
		fe = np.nanmean(np.nansum(feagle[zle/2-_zw: zle/2+_zw+1, inside_eagle], 0))

	f = np.nanmean(fmuse[:, inside_muse], 1)
	x = zrange*68.5
	sb = np.array(f-fr)*flux2sb
	sbstd = sb-sb + fstd*flux2sb
	sbrand = sb-sb + fr*flux2sb
	plt.figure(figsize=(10, 5))
	plt.xlabel('v [km/s]')
	plt.ylabel('SB')
	plt.plot(x, sb-sbrand,)
	plt.plot(x, 2*sbstd, label=r'$2\sigma$ level')
	plt.legend()
	plt.savefig('../../Figures/%d-%d_zw%d%s%s.jpg' % (_r0, _r1, _zw*2+1, smask, '.sclip'*do_sclip))
	#plt.show()