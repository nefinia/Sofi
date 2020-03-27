#!/usr/bin/python

__author__ = 'gallegos'
import glob
import os
from sys import argv

# import eagleSqlTools
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
from matplotlib.cm import get_cmap
from pyfits import getdata, PrimaryHDU

from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
snaps = rdarg(argv, 'snap', list, [8, 9, 10, 11, 12], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
stack = rdarg(argv, 'stack', bool, False)
corr = rdarg(argv, 'corr', bool, True)
scorr = '.corr' * corr
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
galrad = rdarg(argv, 'rad', int, 3)
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

gamma_bkg = {}
heat_bkg = {}
sgamma = {}
sheat = {}
gamma_bkg['HM12'] = {'7': [3.60E-13, 2.14E-13, 3.23E-18],
					 '8': [4.26E-13, 2.45E-13, 1.33E-17],
					 '9': [4.86E-13, 2.81E-13, 5.60E-17],
			 		 '10': [5.71E-13, 3.28E-13, 2.31E-16],
					 '11': [6.77E-13, 3.89E-13, 8.70E-16],
					 '12': [7.66E-13, 4.42E-13, 2.10E-15],
					 '13': [9.50E-13, 5.55E-13, 9.06E-15],
					 '14': [9.64E-13, 5.67E-13, 1.13E-14]}

heat_bkg['HM12'] = {'7': [3.60E-13, 2.14E-13, 3.23E-18],
					'8': [4.26E-13, 2.45E-13, 1.33E-17],
					'9': [4.86E-13, 2.81E-13, 5.60E-17],
					'10': [2.27E-12, 2.18E-12, 7.24E-15],
					'11': [2.68E-12, 2.62E-12, 2.33E-14],
					'12': [3.02E-12, 3.05E-12, 5.01E-14],
					'13': [3.75E-12, 4.22E-12, 1.78E-13],
					'14': [3.81E-12, 4.42E-12, 2.18E-13]}


gamma_bkg['HM01'] = {'7': [3.60E-13, 2.14E-13, 3.23E-18],
					 '8': [.6E-12, 2.45E-13, 1.33E-17],
					 '9': [.7E-12, 2.81E-13, 5.60E-17],
					 '10': [.9E-12, 3.28E-13, 2.31E-16],
					 '11': [1E-12, 3.89E-13, 8.70E-16],
					 '12': [1.5E-12, 4.42E-13, 2.10E-15],
					 '13': [1.7E-12, 5.55E-13, 9.06E-15],
					 '14': [1.9E-12, 5.67E-13, 1.13E-14]}

heat_bkg['HM01'] = {'7': [1.48E-12, 1.46E-12, 1.56E-16],
					'8': [2.8E-12, 2.45E-13, 1.33E-17],
					'9': [3E-12, 2.81E-13, 5.60E-17],
					'10': [3.5E-12, 2.18E-12, 7.24E-15],
					'11': [4E-12, 2.62E-12, 2.33E-14],
					'12': [5E-12, 3.05E-12, 5.01E-14],
					'13': [5.5E-12, 4.22E-12, 1.78E-13],
					'14': [6E-12, 4.42E-12, 2.18E-13]}

snames = {'7': '007_z005p487', '8': '008_z005p037', '9': '009_z004p485', '10': '010_z003p984', '11': '011_z003p528',
		  '12': '012_z003p017', '13': '013_z002p478', '14': '014_z002p237', '15': '015_z002p012'}
redshifts = {'7': 5.487, '8': 5.037, '9': 4.485, '10': 3.984, '11': 3.528, '12': 3.017, '13': 2.478, '14': 2.237, '15': 2.012}
asec2kpcs = {'7': 6.092, '8': 6.393, '9': 6.754, '10': 7.842, '11': 7.449, '12': 7.108, '13': 8.241, '14': 8.396, '15': 8.516}
dz = {'7': 0.0518, '8': 0.0462, '9': 0.0406,'10': 0.035, '11': .03, '12': .0255, '13': .0208, '14': .01695}
# sbpeaks = [2.396]
# sbpeaks = {'10': 1.011, '11': 1.484, '12':2.396} # for an aperture of 1 asec^2!!! I am converting to flux later on
#sbpeaks calculated with python plots.py -hm True
sbpeaks = {'7': .165, '8': .26, '9': .43,'10': .73, '11': 1.27, '12': 2.49, '13': 5.14, '14': 6.97, '15': 9.12}  #for an aperture of 1 asec^2!!! I am converting to flux later on
zlens = {'7': 50, '8': 45, '9': 39,'10': 34, '11': 29, '12': 25, '13': 20, '14': 18, '15': 16}
dndz = {'7': 6., '8': 5., '9': 4.,'10': 3., '11': 2.3, '12': 2., '13': 1.3, '14': 1., '15': .5} #check first and last ones, high extrapolation w/r to Prochaska+10

cenames = {'10': 'snap10_x.gmasked%d.vcorr.U_-18.23_-17.57' % galrad, '11': 'snap11_x.gmasked%d.vcorr.U_-18.33_-17.67' % galrad,
		   '12': 'snap12_x.gmasked%d.vcorr.U_-17.95_-16.36' % galrad, '9': 'snap9_x.gmasked%d.vcorr.U_-18.19_-17.02' % galrad,
		   '8': 'snap8_x.gmasked%d.vcorr.U_-18.36_-17.59' % galrad, '13': 'snap13_x.gmasked%d.vcorr.U_-18.39_-17.56' % galrad}

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
nsigma = 3
sclip = ('.sclip%d' % nsigma)*do_sclip
do_csub = 1
scsub = '.csub'*do_csub
do_rand = 0
do_sbcorr = 0
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


#zprob = [[2.7, 3.6], [3.6, 3.8], [3.8, 4.5], [3.8, 4.8], [4.5, 5.], [4.5, 5.5]]
zprob = [[2.7, 3.2], [3.2, 3.8], [3.8, 5]]
protos = [[3.467, 3.476], [3.687, 3.731], [4.496, 4.53]]
sred = 'z'
#extraname += '.protocluster'

zlm, ylm, xlm = 101, 201, 201
zrange = np.arange(zlm)
w = ylm / 2
zz, yy, xx = np.ogrid[0:zlm, 0:ylm, 0:xlm]
cm = ((yy - w) ** 2 + (xx - w) ** 2)[0]
mat = {}
zrealp = {}

for _zp in zprob:
	zp = (_zp[0]+_zp[1])/2.
	wred = _zp[1]-_zp[0]
	print 'Probing %.1f < z < %.1f' % (_zp[0], _zp[1])
	fdat = 'muse-vs-eagle%s_z%.2f.dat' % (scat, zp)
	_fin = None

	asec2pix_muse = 5  # muse .2 arcsec per pixel
	fouts = []
	n5d = []
	try:
		glob.os.makedirs('%s/all/simplestacks/' % folder)
	except:
		pass

	_id = []
	frand = {}
	fm = {}

	stackname = '%s/all/simplestacks/stack%s%s%s%s%s%s.z%.2f-%.2f.fits' % (folder, scat, scsub, smask, sclip, scorr, extraname, _zp[0], _zp[1])
	if do_rand: randstackname = '../../all/simplestacks/randstack%s%s%s%s%s.%d.fits' % (scat, smask, sclip, scorr, extraname, nr)

	zreal = {}
	ncat = {}
	hdu = PrimaryHDU()
	for _fitcat in fitcat:
		_stackname = '%s/all/simplestacks/stack_%s%s%s%s%s%s.z%.2f-%.2f.fits' % (
		folder, _fitcat, scsub, smask, sclip, scorr, extraname, _zp[0], _zp[1])

		if do_rand: _randstackname = '../../all/simplestacks/randstack%s%s%s%s%s.%d.fits' % (
		_fitcat, smask, sclip, scorr, extraname, nr)

		if do_gmask: _gstackname = '../../all/simplestacks/gstack_%s%s%s%s.z%.2f-%.2f.r%d.fits' % (
		_fitcat, sclip, scorr, extraname, _zp[0], _zp[1], galrad)

		fall = []
		gall = []
		print '%s catalog' % _fitcat,
		cat = getdata('%s/%s/cats/laes.fits' % (folder, _fitcat), 1)
		ids = cat['ID']
		zs = cat['redshift']
		ls = 1215.67 * (1 + zs)
		sconf = cat['sconf']
		cool = (sconf >= 2) & (zs < _zp[1]) & (zs > _zp[0])
		for proto in protos: cool &= ((zs < proto[0])|(zs > proto[1]))
		if fitcat == 'mosaic':
			udf10 = cat['IN_UDF10']
			cool &= ~udf10
		print 'cool gals', np.sum(cool)
		_zreal = np.nanmean(zs[cool])
		_ncat = np.sum(cool)
		print 'desired redshift', zp, 'real mean %.3f' % _zreal, 'ngal', _ncat
		zreal[_fitcat] = _zreal
		ncat[_fitcat] = _ncat
		if not os.path.isfile(_stackname) or stack:
			if do_rand:
				for nr in range(nrand): frand['%d' % nr] = []

			for i, l in zip(ids[cool], ls[cool]):
				ffits = flaes + '/%s/LAEs/%d%s%s.fits' % (_fitcat, i, scsub, scorr)
				if os.path.isfile(ffits):
					if do_rand:
						for nr in range(nrand):
							fr = getdata(ffits.replace('.fits', '.rand%d.fits' % nr))
							zr, yr, xr = fr.shape
							_fr = fr[:, yr/2-w: yr/2+w+1, xr/2-w: xr/2+w+1]
							if do_mask:
								mr = getdata(ffits.replace('.fits', '.rand.mask%d.fits' % nr))
								_mr = mr[yr / 2 - w: yr / 2 + w + 1, xr / 2 - w: xr / 2 + w + 1]
								bad = np.array(_mr) > 0
								_fr[:, bad[0]] = np.nan
							if do_gmask:
								gr = getdata(ffits.replace('.fits', '.rand.gmask%d.fits' % nr))
								_gr = gr[:, yr/2-w: yr/2+w+1, xr/2-w: xr/2+w+1]
								bad2 = (_gr > 0) & (_gr != i)
								_fr[bad2] = np.nan
							frand['%d' % nr].append(_fr)

					fit = getdata(ffits)
					_zl, _yl, _xl = fit.shape
					_fit = fit[:, _yl/2-w: _yl/2+w+1, _xl/2-w: _xl/2+w+1]

					if do_mask:
						mit = getdata(ffits.replace('.fits', '.mask.fits'))
						_yl, _xl = mit.shape
						_mit = mit[_yl/2-w: _yl/2+w+1, _xl/2-w: _xl/2+w+1]
						bad = np.array(_mit) > 0
						_fit[:, bad[0]] = np.nan

					if do_gmask:
						#out = (np.abs(zlm/2-zz) <= zw0) & (cm > (5*asec2pix_muse)**2)
						fgmask = ffits.replace('.fits', '.gmask%d.fits' % galrad)
						if os.path.isfile(fgmask):
							git = getdata(fgmask)
							_zl, _yl, _xl = git.shape
							_git = git[:, _yl/2-w: _yl/2+w+1, _xl/2-w: _xl/2+w+1]
							bad2 = (_git > 0)# & out
							_fit[bad2] = np.nan
							gall.append(bad2)

					_fit = np.array(_fit)

					fall.append(_fit)
				else:
					print 'No fits?', ffits
					ncat[_fitcat] -= 1
			fall = np.array(fall)

			if do_sclip:
				stds = np.nanstd(fall, 0)
				high_sigma = np.abs(fall) > nsigma * stds
				fall[high_sigma] = np.nan
			fm[_fitcat] = np.nanmean(fall, 0)


			hdu.data = fm[_fitcat]
			hdu.writeto(_stackname, clobber=True)
			if do_gmask:
				gall = np.nansum(gall, 0)
				hdu.data = gall
				hdu.writeto(_gstackname, clobber=True)
			if do_rand:
				if do_sclip:
					stds = np.nanstd(frand['%d' % nr], 0)
					high_sigma = np.abs(frand['%d' % nr]) > nsigma * stds
					frand['%d' % nr][high_sigma] = np.nan
				frand[_fitcat] = np.nanmean(frand['%d' % nr], 0)
				hdu.data = frand[_fitcat]
				hdu.writeto(randstackname, clobber=True)
				zlr, ylr, xlr = frand[_fitcat].shape
		else:
			fm[_fitcat] = getdata(_stackname)
			if do_rand:
				frand[_fitcat] = getdata(_randstackname)

	print 'Desired redshift', zp, 'mean redshift(s)', zreal, 'ngal', ncat

	_z, _n, _ss, sncat = 0., 0., {}, ''
	far = cm > (20 * asec2pix_muse) ** 2
	for _f_ in fitcat:
		if ncat[_f_] > 0:

			_ss[_f_] = np.nanstd(fm[_f_][:zlm-zw0, far])
			_z += zreal[_f_] / _ss[_f_]
			_n += ncat[_f_]
		sncat += '_%d' % ncat[_f_]
	_stdsum = np.sum([1./sss for sss in _ss.values()])
	_zreal = _z / _stdsum
	zrealp['%.2f' % zp] = _zreal
	sred += '_%.2f' % _zreal

	if not os.path.isfile(stackname) or stack:
		_fits = []
		for _f_ in fitcat:
			if ncat[_f_] > 0:
				if do_sbcorr:
					foff = np.nanmean(fm[_f_][:, far], 1)
					ff = []
					for i in zrange:
						ff.append(fm[_f_][i, :, :]-foff[i])
					ff = np.array(ff)
				else: ff = fm[_f_]
				_fits.append(ff/_ss[_f_])
		_f = np.nansum(_fits, 0)/_stdsum
		hdu.data = _f
		hdu.writeto(stackname, clobber=True)
		fmuse = _f
	else:
		fmuse = getdata(stackname)

	if do_eagle:
		rval = np.array(redshifts.values())
		rkey = np.array(redshifts.keys())
		if _zreal > 5: zhi, zlo = 5.037, 4.485
		else:
			zlo = np.amax(rval[rval < _zreal])
			zhi = np.amin(rval[rval > _zreal])
		m = (_zreal - zlo) / (zhi - zlo)
		print m
		s = [rkey[rval == _ze][0] for _ze in [zlo, zhi]]
		# cat = [getdata('%s/EAGLE/cats/gals_snap%s.fits' % (folder, _s), 1) for _s in s]
		cename = ['../../EAGLE/simplestacks/%s.fits' % cenames[_s] for _s in s]
		_f = [getdata(_ce) for _ce in cename]
		feagle = _f[0] * (1 - m) + _f[1] * m
		sbpeak = sbpeaks[s[0]] * (1 - m) + sbpeaks[s[1]] * m
		zle, yle, xle = feagle.shape
		zle0 = zle
		_a = [asec2kpcs[_s] * (1 + redshifts[_s]) * kpc2pix for _s in s]
		asec2pix_eagle = _a[0] * (1 - m) + _a[1] * m
		rpix_eagle = np.array(rads) * asec2pix_eagle
		y, x = np.ogrid[0: yle, 0: xle]
		ce = (x - xle / 2.) ** 2 + (y - yle / 2.) ** 2
		sb2gfac = 0.00127 # expected SB 1.31e-20 for a gamma of .684e-12 at z=3.5 (accounting for SB dimming)
		#_sb2g = [gamma_bkg['HM12'][_s][0] * 1e12 / sbpeaks[_s] for _s in s]
		_sb2g = [sb2gfac*(1+redshifts[_s])**4 for _s in s]
		sb2gamma = _sb2g[0] * (1 - m) + _sb2g[1] * m
		print 'SB to Gamma for z %.3f' % _zreal, sb2gamma

	if do_rand:
		yr, xr = np.ogrid[0: ylr, 0: xlr]
		cr = (xr - xlr / 2.) ** 2 + (yr - ylr / 2.) ** 2
	
	
	odat = '../../UVB/SNR/UVB_z%.1f-%.1f_zreal%.2f%s_ngal%s%s%s%s%s%s.dat' \
		   % (_zp[0], _zp[1], _zreal, scat, sncat, scsub, scorr, smask, sclip, extraname)
	mat['%.2f' % zp] = []
	if not os.path.isfile(odat) or overwrite:
		rw = 6
		extraname = ''#'.rw%d_zoff-5_5' % rw
		#fout = open(odat, 'w')
	
		#if do_eagle: fout.write('#r r0 r1 zw zoff SB SB_std SB_rand feagle fLLS SB_HM12 Gamma Gamma_uplim Npix\n')
		#else: fout.write('#r r0 r1 zw zoff SB SB_std SB_rand\n')
	
		#if _fitcat =='mosaic': rmin0, rmax0, rmin1, rmax1 = 0, 100, 10, 200
		#else: rmin0, rmax0, rmin1, rmax1 = 0, 15, 5, 20
	
		#rmin0, rmax0, rmin1, rmax1 = 0, 21, 5, 31
	
		#r0 = np.arange(rmin0, rmax0, rmin1-rmin0)
		#r1 = np.arange(rmin1, rmax1, rmin1-rmin0)
		rads = [[0, 0.5], [0.5, 1],[0, 1], [1, 2], [2, 4], [4, 6], [6, 8], [6, 10], [6, 11], [6, 12], [6, 13], [6, 14], [6, 15], [6, 18], [6, 20], [12, 20], [8, 12], [8, 14], [8, 15], [8, 20], [20, 25], [25, 30], [20, 30], [5, 20]]
		zoffs = np.arange(-50, 51)#
		zws = [1, 2, 3]
		g = {}
		gstd = {}
		sbstd = {}
		sb = {}
		sbrand = {}
	
		#for i in r0:
		#	for j in r1:
		for rad in rads:
			i, j = rad
			rm0 = i * asec2pix_muse
			rm1 = j * asec2pix_muse
			inside_muse = (cm >= rm0 ** 2) & (cm < rm1 ** 2)
			rr = np.sqrt(np.nanmean(cm[inside_muse]))/asec2pix_muse
			_rw = (j-i)/2.
			if do_eagle:
				re0 = i * asec2pix_eagle
				re1 = j * asec2pix_eagle
			print 'Between %d and %d arcsec. Real r %.1f' % (i, j, rr)
	
			if do_eagle: inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
			nin = np.sum(inside_muse)
			for _zw in zws:
				if do_eagle:
					zmin = max(0, zle/2-_zw)
					zmax = min(zle, zle/2+_zw+1)
					_feagle = feagle[zmin: zmax, inside_eagle]
					fe = np.nansum(np.nanmean(_feagle, 1))
	
				if do_rand:
					inside_rand = (cr >= rm0 ** 2) & (cr < rm1 ** 2)
					fr = frand[zlr/2-_zw: zlr/2+_zw+1, inside_rand]
					fr = np.nanmean(np.nansum(fr, 1))
					fstd = []
					for nr in range(nrand):
						fr2 = frand['%d' % nr][:, zlr / 2 - _zw: zlr / 2 + _zw + 1, inside_rand]
						fstd.append(np.nanmean(np.nansum(np.nanmean(fr2, 2), 1)))
					fstd = np.nanstd(fstd)
				else:
					_fr = []
					for zz in np.concatenate([zrange[_zw:zlm/2-3*_zw], zrange[zlm/2+3*_zw:-_zw]]):
						zmin = max(0, zlm/2-_zw+zz)
						zmax = min(zlm, zlm/2+_zw+zz+1)
						fm = fmuse[zz-_zw: zz+_zw+1, inside_muse]
						if do_sclip: fm[abs(fm)>nsigma*np.nanstd(fm)] = np.nan
						_fr.append(np.nansum(np.nanmean(fm, 1)))
					fr = np.nanmean(_fr)
					fstd = np.nanstd(_fr)
	
				noise = fstd*flux2sb
				if do_eagle:
					gamma_uplim = 2*noise*sb2gamma/fe
					sbhm12 = fe * sbpeak
	
				for _zoff in zoffs:
					#print 'zw %d zoff %d' % (_zw*2+1, _zoff)
					zmin = max(0, zlm/2-_zw+_zoff)
					zmax = min(zlm, zlm/2+_zw+_zoff+1)
					_fmuse = fmuse[zmin: zmax, inside_muse]#-np.nanmean(fmuse[zmin: zmax, far])
					if do_sclip: _fmuse[abs(_fmuse) > nsigma * np.nanstd(_fmuse)] = np.nan

					fin = np.nansum(np.nanmean(_fmuse, 1))-fr
					sb = fin*flux2sb
					_gamma = gamma_bkg['HM12'][s[0]][0]*(1-m)+gamma_bkg['HM12'][s[1]][0]*m
					flls = 2*noise*sb2gamma*1e-12/_gamma
	
					if do_eagle:
						nmuse = np.sum(inside_muse)
						fy = fin*flux2sb/fe
						gamma = fy*sb2gamma
						mat['%.2f' % zp].append([rr, i, j, _zw, _zoff, sb, noise, sb/noise, fr*flux2sb, fe, flls, sbhm12, gamma, gamma_uplim])
					else:
						mat['%.2f' % zp].append([(i+j)/2., i, j, _zw, _zoff, fin*flux2sb, noise, fr*flux2sb])

		mat['%.2f' % zp] = np.array(mat['%.2f' % zp])
		if do_eagle: header = 'r r0 r1 zw zoff SB SB_std SNR SB_rand feagle fLLS SB_HM12 Gamma Gamma_uplim'
		else: header = 'r r0 r1 zw zoff SB SB_std SB_rand'

		np.savetxt(odat, mat['%.2f' % zp], header=header)

	else:
		mat['%.2f' % zp] = np.loadtxt(odat)

	fig, ax = plt.subplots(figsize=(10, 5))

	r, r0, r1, zwi, z, sb, std, fe, fr, flls, sbhm12, gamma, gamma_uplim, nmuse = mat['%.2f' % zp].T
	zmin, zmax = -50, 50
	v = 2.99792458e5*z*1.25/1215.67/(1+_zreal)
	vmin, vmax = 2.99792458e5*zmin*1.25/1215.67/(1+_zreal), 2.99792458e5*zmax*1.25/1215.67/(1+_zreal)
	rsel = [2, 4], [4, 6], [6, 12], [12, 20]
	_zw_ = 1
	for rr in rsel:
		c = (zwi==_zw_)&(r0==rr[0])&(r1==rr[1])&(z>=zmin)&(z<=zmax)
		ax.plot(v[c], sb[c]/(2*_zw_+1), label=r'$%.1f<r<%.1f$' % (rr[0], rr[1]))
		# ax.plot(v[zmin: zmax], [2*stf]*len(r[c]), label=r'2$\sigma$ noise level')
		#ax.plot(z[c]*1.25, [0] * len(r[c]), color='gray')
		# ax.scatter([0], [1], color='red')
	plt.xlabel(r'v [km/s]')
	plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}}$')
	plt.xlim([vmin, vmax])
	plt.legend()
	plt.grid()

	# ax2 = ax.twinx()
	# ax2.plot(v[zmin: zmax], f/fstd, alpha=0)
	# plt.ylabel(r'SNR')

	plt.savefig('../../Figures/SB_spectra%s_redshift%.1f_zw%d%s%s%s.png' % (
	scat, _zreal, _zw_, scsub, smask, extraname))
	# plt.show()
	plt.close()

#rsel = [0, .5], [.5, 1], [0, 1], [1, 2], [2, 4], [4, 6], [6, 12], [8, 12], [10, 14]
#zmin, zmax = -20, 20
#zprob = [3.5, 4, 4.5]
#zprob = [3.72]
#if fmos: zprob = [3.5, 3.75, 4]
vlim = 3000
for rr in rsel:
	fig, ax = plt.subplots(figsize=(10, 5))
	_zw = 1
	for zp in zprob:
		zpm = (zp[0]+zp[1])/2.
		_zp = zrealp['%.2f' % zpm]
		r, r0, r1, zwi, z, sb, std, fe, fr, flls, sbhm12, gamma, gamma_uplim, nmuse = mat['%.2f' % zpm].T

		v = 2.99792458e5*z*1.25/1215.67/(1+_zp)
		c = (zwi == _zw)&(r0 == rr[0])&(r1 == rr[1])&(np.abs(v) <= vlim)
		ax.plot(v[c], sb[c]/(2*_zw+1.), label=r'$z=%.1f$' % _zp)
		# ax.plot(v[zmin: zmax], [2*fstd*flux2sb]*len(f), label=r'2$\sigma$ noise level')
		#ax.plot(z[c]*1.25, [0] * len(r[c]), color='gray')
		# ax.scatter([0], [1], color='red')

	plt.xlabel(r'v [km/s]')
	plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}}$')
	plt.xlim([-vlim, vlim])
	plt.legend()
	plt.grid()
	plt.savefig('../../Figures/SB_spectra%s%s_r%.1f-%.1f_zw%d%s%s%s.png' % (
		sred, scat, rr[0], rr[1], _zw, scsub, smask, extraname))
	# plt.show()
	plt.close()


if 0:
	fig, ax = plt.subplots(figsize=(10, 5))
	for zp in zprob:
		zpm = (zp[0] + zp[1]) / 2.
		_zp = zrealp['%.2f' % zpm]
		r, r0, r1, zwi, z, sb, std, fe, fr, flls, sbhm12, gamma, gamma_uplim, nmuse = mat['%.2f' % zpm].T
		v = 2.99792458e5 * z * 1.25 / 1215.67 / (1 + _zp)
		c = (zwi == _zw) & (r0 == rr[0]) & (r1 == rr[1]) & (np.abs(v) <= vlim)
		ax.plot(v[c], sb[c]*(1+zrealp['%.2f' % zpm])**4./(_zw+1.), label=r'$z=%.1f$' % _zp)
	# ax.plot(v[zmin: zmax], [2*fstd*flux2sb]*len(f), label=r'2$\sigma$ noise level')
	# ax.plot(z[c]*1.25, [0] * len(r[c]), color='gray')
	# ax.scatter([0], [1], color='red')

	plt.xlabel(r'v [km/s]')
	plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}}$')
	plt.xlim([-vlim, vlim])
	plt.legend()
	plt.grid()
	plt.savefig('../../Figures/SB_spectra_zcorr_%s%s_r%.1f-%.1f_zw%d%s%s%s.png' % (
		sred, scat, rr[0], rr[1], 0, scsub, smask, extraname))
	# plt.show()
	plt.close()