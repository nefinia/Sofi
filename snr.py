#!/usr/bin/python

__author__ = 'gallegos'
import glob
import os
from sys import argv

# import eagleSqlTools
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
import scipy.signal as sp
from matplotlib.cm import get_cmap
from pyfits import getdata, PrimaryHDU

import params
from tools_sofi import rdarg, UVBcalc  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

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
mask = rdarg(argv, 'mask', bool, False)
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

cenames = {'10': 'snap10_x.gmasked%d.vcorr.LLS' % (galrad),#'snap10_x.gmasked%d.vcorr.LLS.U_-18.79_-18.23' % (galrad),
           '11': 'snap11_x.gmasked%d.vcorr.LLS' % (galrad),#.U_-18.33_-17.67
		   '12': 'snap12_x.gmasked%d.vcorr.LLS' % (galrad),#.U_-17.95_-16.36
		   '9': 'snap9_x.gmasked%d.vcorr.LLS' % (galrad),#.U_-18.19_-17.02
		   '8': 'snap8_x.gmasked%d.vcorr.LLS' % (galrad),#.U_-18.36_-17.59
		   '13': 'snap13_x.gmasked%d.vcorr.LLS' % (galrad)}#.U_-18.39_-17.56

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
flaes = '/net/galaxy-data/export/galaxy/shared/MUSE/gallegos/'

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
do_mask = 0
do_gmask = 1
smask = ('.%smask' % (('g%d' % galrad)*do_gmask))*do_mask
do_sclip = 1
nsigma = 3
sclip = ('.sclip%d' % nsigma)*do_sclip
do_rand = 0
do_sbcorr = 0
fix_red = 0
zoff = 0
SB = 1.27  # SB at z=3.5

fontsize = 14
figsize = (8, 4.5)

cmap = get_cmap('gist_rainbow')

#colors = [(0, 0, 1), (1, 1, 0), (1, 0, 0)]
#cmap = LinearSegmentedColormap.from_list('sofi', colors, N=50)
# 1 spaxel is (0.2 arcsec)^2
# flux2sb comes from 1.25 Angstrom wavelength width times 1/(0.2)**2 from flux2sbersion to the correct aperture
flux2sb = 31.25


#zprob = [[2.7, 3.6], [3.6, 3.8], [3.8, 4.5], [3.8, 4.8], [4.5, 5.], [4.5, 5.5]]
zprob = [[2.8, 3.5], [3.2, 3.8], [3.2, 4.2], [3.5, 4.5], [3.8, 5], [4.5, 5.5]]#[2.8, 3.2]
#zprob = [[4.2, 4.5]]
protos = [[3.45, 3.5], [3.65, 3.75], [4.5, 4.6]]
#zprob = [[5.12, 5.14]]#, [4.76, 4.78], [3.467, 3.476], [3.687, 3.731], [4.496, 4.53]]
#zprob = [[3.467, 3.476], [3.687, 3.731], [4.496, 4.53]]
#protos = []
#zprob = [[3.6, 4.6], [3.8, 4.8]]
#zprob = [[2.8, 3.5], [3.5, 4.2], [4.2, 4.9], [4.9, 5.6]]#[2.8, 3.2]
zprob = [[2.8, 3.4], [3.4, 4.4], [4.4, 5.4], [5.4, 6.6]]#[2.8, 3.2]
#zprob = [[2.9, 3.4], [3.4, 4.4], [4.4, 5.4], [5.4, 6.6]]#[2.8, 3.2]
#zprob = [[4.4, 4.6]]#[2.8, 3.2]
#zprob = [[4.5, 6.6]]#[2.8, 3.2]
#zprob = [[3., 4.], [4., 5.], [5, 6.6]]#[2.8, 3.2]
#protos = []
#zprob = [[3.45, 3.5], [3.65, 3.75], [4.45, 4.55]]

sred = 'z'
#extraname += '.protocluster'

zlm, ylm, xlm = 101, 201, 201
zrange = np.arange(zlm)
w = ylm / 2
zz, yy, xx = np.ogrid[0:zlm, 0:ylm, 0:xlm]
cm = ((yy - w) ** 2 + (xx - w) ** 2)[0]
mat = {}
zrealp = {}


odats = ['UVB_zr3.16_UDF_mosaic.csub.corr.g6mask.sclip3_z2.8-3.4_ngal_40_104.dat',
		 'UVB_zr3.90_UDF_mosaic.csub.corr.g6mask.sclip3_z3.4-4.4_ngal_37_163.dat',
		 'UVB_zr4.88_UDF_mosaic.csub.corr.g6mask.sclip3_z4.4-5.4_ngal_34_156.dat',
		 'UVB_zr5.82_UDF_mosaic.csub.corr.g6mask.sclip3_z5.4-6.6_ngal_24_65.dat']
zreals = [3.16, 3.9, 4.88, 5.82]

fsize = 14
colors2 = 'm', 'darkorange','g'

for i in range(4):

	_zp = zprob[i]
	odat = '../../UVB/SNR/'+odats[i]
	_zreal = zreals[i]
	zp = (_zp[0]+_zp[1])/2.
	zrealp['%.2f' % zp] = _zreal


	wred = _zp[1]-_zp[0]
	print 'Probing %.1f < z < %.1f' % (_zp[0], _zp[1])
	fdat = 'muse-vs-eagle%s_z%.2f.dat' % (scat, zp)
	_fin = None

	asec2pix_muse = 5  # muse .2 arcsec per pixel


	if 1:
		fouts = []
		n5d = []
		try:
			glob.os.makedirs('%s/all/simplestacks/' % folder)
		except:
			pass

		_id = []
		frand = {}
		fm = {}
		gm = {}

		stackname = '%s/all/simplestacks/stack%s%s%s%s%s%s.z%.2f-%.2f.fits' % (folder, scat, scsub, smask, sclip, scorr, extraname, _zp[0], _zp[1])
		gstackname = '%s/all/simplestacks/gstack%s%s.z%.2f-%.2f.fits' % (folder, scat, extraname, _zp[0], _zp[1])
		if do_rand: randstackname = '../../all/simplestacks/randstack%s%s%s%s%s.%d.fits' % (scat, smask, sclip, scorr, extraname, nr)

		zreal = {}
		ncat = {}
		hdu = PrimaryHDU()
		for _fitcat in fitcat:
			_stackname = '%s/all/simplestacks/stack_%s%s%s%s%s%s.z%.2f-%.2f.fits' % (
			folder, _fitcat, scsub, smask, sclip, scorr, extraname, _zp[0], _zp[1])

			if do_rand: _randstackname = '../../all/simplestacks/randstack%s%s%s%s%s.%d.fits' % (
			_fitcat, smask, sclip, scorr, extraname, nr)

			if do_gmask: _gstackname = '../../all/simplestacks/gstack_%s%s.z%.2f-%.2f.r%d.fits' % (
			_fitcat, extraname, _zp[0], _zp[1], galrad)

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
									mr = getdata(ffits.replace('.fits', '.rand.mask.%d.fits' % nr))
									_mr = mr[yr / 2 - w: yr / 2 + w + 1, xr / 2 - w: xr / 2 + w + 1]
									bad = np.array(_mr) > 0
									_fr[:, bad[0]] = np.nan
								if do_gmask:
									gr = getdata(ffits.replace('.fits', '.rand.gmask%d.%d.fits' % (galrad, nr)))
									_gr = gr[:, yr/2-w: yr/2+w+1, xr/2-w: xr/2+w+1]
									bad2 = (_gr > 0) & (_gr != i)
									_fr[bad2] = np.nan
								frand['%d' % nr].append(_fr)

						fit = getdata(ffits)
						_zl, _yl, _xl = fit.shape
						_fit = fit[:, _yl/2-w: _yl/2+w+1, _xl/2-w: _xl/2+w+1]

						if do_mask:
							fmask = flaes + '/%s/LAEs/%d.mask.fits' % (_fitcat, i)
							mit = getdata(fmask)
							_zl, _yl, _xl = mit.shape
							#_mit = mit[:, _yl/2-w: _yl/2+w+1, _xl/2-w: _xl/2+w+1]
							#gid = mit[_zl/2, _yl/2, _xl/2]
							#mit[mit == gid] == 0
							bad = mit > 0
							_fit[bad] = np.nan

						if do_gmask:
							fgmask = flaes + '/%s/LAEs/%d.gmask%d.fits' % (_fitcat, i, galrad)
							#out = (np.abs(zlm/2-zz) <= zw0) & (cm > (5*asec2pix_muse)**2)
							if os.path.isfile(fgmask):
								git = getdata(fgmask)
								_zl, _yl, _xl = git.shape
								#_git = git[:, _yl/2-w: _yl/2+w+1, _xl/2-w: _xl/2+w+1]
								bad2 = (git > 0)# & out
								_fit[bad2] = np.nan
								gall.append(bad2)
							else: print 'No gmask?', fgmask

						bkgcorr = True
						if bkgcorr:
							out = ((xx-xlm/2)**2+(yy-ylm/2)**2 > 100**2) & (abs(zz-zlm/2)>10)
							fout = _fit[out]
							fstd = np.nanstd(fout)
							fout[abs(fout) > nsigma*fstd] = np.nan
							bkg = np.nanmean(fout)
							_fit -= bkg

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
					gm[_fitcat] = np.nanmean(gall, 0)
					hdu.data = gm[_fitcat]
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
		if len(fitcat)>1:
			for _f_ in fitcat:
				if ncat[_f_] > 0:
					_ss[_f_] = np.nanstd(fm[_f_][:zlm-zw0, far])
					_z += zreal[_f_] / _ss[_f_]
					_n += ncat[_f_]
				sncat += '_%d' % ncat[_f_]
			_stdsum = np.sum([1./sss for sss in _ss.values()])
			_zreal = _z / _stdsum
		else: sncat = '_%d' % ncat[fitcat[0]]
		zrealp['%.2f' % zp] = _zreal
		sred += '_%.2f' % _zreal

		if not os.path.isfile(stackname) or stack:
			_fits = []
			if do_gmask: _gits = []
			if len(fitcat) > 1:
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
						if do_gmask: _gits.append(gm[_f_] /_ss[_f_])
				_f = np.nansum(_fits, 0)/_stdsum
				if do_gmask: _g = np.nansum(_gits, 0) / _stdsum
			else:
				_f = fm[fitcat[0]]
				if do_gmask: _g = gm[fitcat[0]]
			hdu.data = _f
			hdu.writeto(stackname, clobber=True)
			fmuse = _f
			if do_gmask:
				hdu.data = _g
				hdu.writeto(gstackname, clobber=True)
				gmuse = _g
		else:
			fmuse = getdata(stackname)
			gmuse = getdata(gstackname)

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
			_f = [getdata('../../EAGLE/simplestacks/%s.fits' % cenames[_s]) for _s in s]
			_fcorr = [getdata('../../EAGLE/simplestacks/%s.fits' % cenames[_s].replace('.LLS', '.LLScorr')) for _s in s]
			#_fcorr = fLLScorr[s[0]] * (1 - m) + fLLScorr[s[1]] * m
			feagle = _f[0] * (1 - m) + _f[1] * m
			fcorr = _fcorr[0] * (1 - m) + _fcorr[1] * m
			sbpeak = sbpeaks[s[0]] * (1 - m) + sbpeaks[s[1]] * m
			_gamma12 = gamma_bkg['HM12'][s[0]][0] * (1 - m) + gamma_bkg['HM12'][s[1]][0] * m
			_gamma01 = gamma_bkg['HM01'][s[0]][0] * (1 - m) + gamma_bkg['HM01'][s[1]][0] * m
			zle, yle, xle = feagle.shape
			zle0 = zle
			_a = [asec2kpcs[_s] * (1 + redshifts[_s]) * kpc2pix for _s in s]
			asec2pix_eagle = _a[0] * (1 - m) + _a[1] * m
			rpix_eagle = np.array(rads) * asec2pix_eagle
			y, x = np.ogrid[0: yle, 0: xle]
			ce = (x-xle/2.)**2 + (y-yle/2.) ** 2
			#sb2gfac = .001268 #.677/1.27/4.528**4, for HM12
			RHM12, SBHM12, gamma_HM12 = UVBcalc(_zreal)
			sb2gamma = gamma_HM12*1e-8/SBHM12#sb2gfac*(1+_zreal)**4
			print 'SB to Gamma for z %.3f' % _zreal, sb2gamma

		if do_rand:
			yr, xr = np.ogrid[0: ylr, 0: xlr]
			cr = (xr - xlr / 2.) ** 2 + (yr - ylr / 2.) ** 2


		odat = '../../UVB/SNR/UVB_zr%.2f%s%s%s%s%s_z%.1f-%.1f_ngal%s%s.dat' \
			   % (_zreal, scat, scsub, scorr, smask, sclip, _zp[0], _zp[1], sncat, extraname)
		mat['%.2f' % zp] = []

		if not os.path.isfile(odat) or overwrite:

			rw = 6
			extraname = ''#'.rw%d_zoff-5_5' % rw

			#if _fitcat =='mosaic': rmin0, rmax0, rmin1, rmax1 = 0, 100, 10, 200
			#else: rmin0, rmax0, rmin1, rmax1 = 0, 15, 5, 20

			#rmin0, rmax0, rmin1, rmax1 = 0, 21, 5, 31

			#r0 = np.arange(rmin0, rmax0, rmin1-rmin0)
			#r1 = np.arange(rmin1, rmax1, rmin1-rmin0)
			#rads = [[0, .3], [0.3, 1], [0, 1], [1, 2], [2, 4], [4, 6], [6, 8], [8, 12], [12, 20], [6, 20], [20, 30], [8, 20], [8, 16]]
			rads = [[1, 2], [2, 4], [4, 6], [6, 8], [6, 20]]
			zoffs = np.arange(-50, 51)
			zws = [0, 2, 4, 6]
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
				print 'Between %.1f and %.1f arcsec. Real r %.1f' % (i, j, rr)

				if do_eagle: inside_eagle = (ce >= re0 ** 2) & (ce < re1 ** 2)
				nin = np.sum(inside_muse)
				for _zw in zws:
					if do_eagle:
						zmin = max(0, zle/2-_zw)
						zmax = min(zle, zle/2+_zw+1)
						_feagle = feagle[zmin: zmax, inside_eagle]
						_fcorr = fcorr[zmin: zmax, inside_eagle]
						fe = np.nansum(np.nanmean(_feagle, 1))
						fc = np.nansum(np.nanmean(_fcorr, 1))

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
						for _zz in np.concatenate([zrange[_zw:zlm/2-10*_zw], zrange[zlm/2+10*_zw:-_zw]]):
							zmin = max(0, zlm/2-_zw+_zz)
							zmax = min(zlm, zlm/2+_zw+_zz+1)
							fm = fmuse[_zz-_zw: _zz+_zw+1, inside_muse]
							if do_sclip: fm[abs(fm)>nsigma*np.nanstd(fm)] = np.nan
							_fr.append(np.nansum(np.nanmean(fm, 1)))
						fr = np.nanmean(_fr)
						fstd = np.nanstd(_fr)

					noise = fstd*flux2sb
					if do_eagle:
						gamma_uplim = 2*noise*sb2gamma/fe
						sbhm12 = fe*sbpeak
						gamma_uplimc = 2*noise*sb2gamma/fc
					for _zoff in zoffs:
						#print 'zw %d zoff %d' % (_zw*2+1, _zoff)
						zmin = max(0, zlm/2-_zw+_zoff)
						zmax = min(zlm, zlm/2+_zw+_zoff+1)
						v = 2.99792458e5 * (_zoff) * 1.25 / 1215.67 / (1 + _zreal)
						_fmuse = fmuse[zmin: zmax, inside_muse]#-np.nanmean(fmuse[zmin: zmax, far])
						if do_sclip: _fmuse[abs(_fmuse) > nsigma * np.nanstd(_fmuse)] = np.nan

						fin = np.nansum(np.nanmean(_fmuse, 1))-fr
						sb = fin*flux2sb
						snr = sb/noise
						sbr = fr*flux2sb

						flls = sb*sb2gamma*1e-12/_gamma12
						flls_uplim = 2*noise*sb2gamma*1e-12/_gamma12

						if do_eagle:
							nmuse = np.sum(inside_muse)
							fy = fin*flux2sb/fe
							fyc = fin*flux2sb/fc
							gamma = fy*sb2gamma
							gammac = fyc*sb2gamma
							mat['%.2f' % zp].append([rr, i, j, _zw, _zoff, sb, noise, fc, flls, flls_uplim, gammac, gamma_uplimc, snr, sbr, fe, gamma, gamma_uplim, sbhm12, v])
						else:
							mat['%.2f' % zp].append([(i+j)/2., i, j, _zw, _zoff, fin*flux2sb, noise, fr*flux2sb])

			mat['%.2f' % zp] = np.array(mat['%.2f' % zp])
			if do_eagle: header = 'r r0 r1 zw z SB SB_std fcorr fLLS fLLS_uplim Gamma Gamma_uplim SNR SB_rand fEAGLE GammaEAGLE Gamma_uplimEAGLE SB_HM12 v'
			else: header = 'r r0 r1 zw z SB SB_std SB_rand'

			np.savetxt(odat, mat['%.2f' % zp], header=header)

		else:
			mat['%.2f' % zp] = np.loadtxt(odat)

	else:
		mat['%.2f' % zp] = np.loadtxt(odat)

	fig, ax = plt.subplots(figsize=figsize)

	r, r0, r1, zwi, z, sb, std, fc, flls, flls_uplim, gammac, gamma_uplimc, snr, sbr, fe, gamma, gamma_uplim, sbhm12, v = mat['%.2f' % zp].T
	zmin, zmax = -50, 50
	#vmin, vmax = 2.99792458e5*zmin*1.25/1215.67/(1+_zreal), 2.99792458e5*zmax*1.25/1215.67/(1+_zreal)
	vmin, vmax = -2000, 2000
	rsel = [2, 4], [4, 6], [6, 20]#, [8, 20], [8, 14]
	_zw_ = 0
	i = 0
	for rr in rsel:
		c = (zwi==_zw_)&(r0==rr[0])&(r1==rr[1])&(z>=zmin)&(z<=zmax)
		ax.plot(v[c], sp.savgol_filter(sb[c], 3, 1), label=r"$%d''<r<%d''$" % (rr[0], rr[1]), color=colors2[i])
		cout = (zwi == _zw_) & (r0 == rr[0]) & (r1 == rr[1]) & (np.abs(v) > 1000)
		noise = np.std(sb[cout])
		print noise
		ax.errorbar(vmax * (.6 + i * .1), .3, yerr=noise, color=colors2[i], capsize=3)
		i += 1
		#ax.plot(v[c], sb[c], label=r"$%d''<r<%d''$" % (rr[0], rr[1]))
		# ax.plot(v[zmin: zmax], [2*stf]*len(r[c]), label=r'2$\sigma$ noise level')
		#ax.plot(z[c]*1.25, [0] * len(r[c]), color='gray')
		# ax.scatter([0], [1], color='red')
	plt.xlabel(r'v [$\mathrm{km\,s^{-1}}$]', fontsize=fsize)
	plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}}$', fontsize=fsize)
	plt.xlim([vmin, vmax])
	plt.ylim([-.3, .6])
	plt.xticks(fontsize=fsize-2)
	plt.yticks(fontsize=fsize-2)
	plt.legend(prop={"size":fsize})
	plt.title('z=%.1f' % _zreal, fontsize=fsize)
	plt.grid()

	# ax2 = ax.twinx()
	# ax2.plot(v[zmin: zmax], f/fstd, alpha=0)
	# plt.ylabel(r'SNR')

	plt.savefig('../../Figures/SB_spectra%s_redshift%d_zw%d.pdf' % (
	scat, _zreal*10, _zw_))
	# plt.show()
	plt.close()

#rsel = [0, .5], [.5, 1], [0, 1], [1, 2], [2, 4], [4, 6], [6, 12], [8, 12], [10, 14]
#zmin, zmax = -20, 20
#zprob = [3.5, 4, 4.5]
#zprob = [3.72]
#if fmos: zprob = [3.5, 3.75, 4]
vlim = 2000
#rsel = [0, .3], [0.3, 1], [1, 2], [2, 4], [4, 6], [6, 20]
rsel = [1, 2], [2, 4], [4, 6], [6, 20]
#ylims = [[-10, 60], [-5, 30], [-1, 5], [-.4, .6], [-.3, .4], [-.1, .15]]
ylims = [-1, 4], [-.4, .6], [-.3, .4], [-.1, .15]
sb620 = []
std620 = []
colors = 'blue', 'green', 'gold', 'red'
#colors = 'blue', 'green', 'red'#, 'gold'

for rr, ylim in zip(rsel, ylims):
	fig, ax = plt.subplots(figsize=figsize)
	_zw = 0
	i = 0
	sred = ''
	for zp in zprob:
		zpm = (zp[0]+zp[1])/2.
		_zp = zrealp['%.2f' % zpm]
		sred += '_%d' % (_zp*10)
		r, r0, r1, zwi, z, sb, std, fc, flls, flls_uplim, gammac, gamma_uplimc, snr, sbr, fe, gamma, gamma_uplim, sbhm12, v = mat['%.2f' % zpm].T
		c = (zwi == _zw)&(r0 == rr[0])&(r1 == rr[1])&(np.abs(v) <= vlim)
		ax.plot(v[c], sp.savgol_filter(sb[c], 3, 1), label=r'$z=%.1f$' % _zp, color=colors[i])
		cout = (zwi == _zw) & (r0 == rr[0]) & (r1 == rr[1]) & (np.abs(v) > 1000)
		noise = np.std(sb[cout])
		print noise
		ax.errorbar(-vlim*(.9-i*.1), ylim[1]*.5, yerr=noise, color=colors[i], capsize=3)
		#ax.plot(v[c], sb[c], label=r'$z=%.1f$' % _zp, color=colors[i])
		i += 1
		# ax.plot(v[zmin: zmax], [2*fstd*flux2sb]*len(f), label=r'2$\sigma$ noise level')
		#ax.plot(z[c]*1.25, [0] * len(r[c]), color='gray')
		# ax.scatter([0], [1], color='red')
		cool = (z==-2) & (zwi==2) & (r0==rr[0]) & (r1==rr[1])
		print 'z', zp, _zp, '%d<r<%d' % (rr[0], rr[1]), 'SB', sb[cool], 'std', std[cool], 'G', gammac[cool], 'uplim', gamma_uplimc[cool], 'fLLS', fc[cool]
	plt.xlabel(r'v [$\mathrm{km\,s^{-1}}$]', fontsize=fsize)
	plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}}$', fontsize=fsize)
	plt.xlim([-vlim, vlim])
	plt.ylim(ylim)
	plt.xticks(fontsize=fsize-2)
	plt.yticks(fontsize=fsize-2)
	plt.legend(prop={"size":fsize})
	if rr[0] < .3: plt.title('r<%.1f"' % rr[1])
	elif rr[0] < 1: plt.title('%.1f"<r<%.1f"' % (rr[0], rr[1]))
	else: plt.title('%d"<r<%d"' % (rr[0], rr[1]), fontsize=fsize)
	plt.grid()
	plt.savefig('../../Figures/SB_spectra%s%s_r%d-%d_zw%d.pdf' % (
		sred, scat, rr[0], rr[1], _zw))
	# plt.show()
	plt.close()

