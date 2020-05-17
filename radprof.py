#!/usr/bin/python

__author__ = 'gallegos'
import glob
# import eagleSqlTools
import matplotlib.pyplot as plt
import os
import scipy.interpolate as inter
from matplotlib.cm import get_cmap as cmap
from pyfits import getdata, PrimaryHDU
from sys import argv
import astropy.io.fits as fits
import params
import scipy.signal as sp
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
circlehist = rdarg(argv, 'circlehist', bool, False)
do_delaunay = rdarg(argv, 'delaunay', bool, False)
extraname = rdarg(argv, 'extraname', str, '')#'HM12')#
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
galcov = rdarg(argv, 'galcov', bool, True)
galrad = rdarg(argv, 'rad', int, 6)
h2d = rdarg(argv, 'h2d', bool, False)
halfplot = rdarg(argv, 'halfplot', bool, False)
histover = rdarg(argv, 'histover', bool, False)
kde = rdarg(argv, 'kde', bool, False)
LLScubes = rdarg(argv, 'LLScubes', bool, False)
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
nhiprof = rdarg(argv, 'nhi', bool, False)
minres = rdarg(argv, 'minres', int, 512)
maxres = rdarg(argv, 'maxres', int, 4096)
mask = rdarg(argv, 'mask', bool, False)
model = rdarg(argv, 'model', str, 'HM01')  # 'HM12')#
npref = rdarg(argv, 'npref', int, 12)
overwrite = rdarg(argv, 'overwrite', bool, False)
pairLLS = rdarg(argv, 'pairLLS', bool, False)
parallel = rdarg(argv, 'parallel', bool, False)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
snaps = rdarg(argv, 'snap', list, [6, 7, 8, 9, 10, 11, 12, 13, 14], int)
sql = rdarg(argv, 'sql', bool, False)
sphericalLLS = rdarg(argv, 'sphericalLLS', bool, False)
radprof = rdarg(argv, 'radprof', bool, False)
sbhist = rdarg(argv, 'sbhist', bool, False)
ssthr = rdarg(argv, 'ssthr', float, 1e10)#6.73E-3
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)
type = rdarg(argv, 'type', str, 'LLS') #NHtot, f_NHI, NHII
unique = rdarg(argv, 'unique', bool, False)


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

#folder = '/net/abnoba/scratch2/gallegos/Research/MUSE/'
#flaes = '/net/galaxy-data/export/galaxy/shared/MUSE/gallegos/'

folder = '../../'
flaes = '../../'

scat = ''
for i in fitcat: scat += '_' + i
ncores = 10
nrand = 10
xl = rdarg(argv, 'xl', int, 301)
yl = rdarg(argv, 'yl', int, 301)
zl = rdarg(argv, 'zl', int, 15)
zw0 = rdarg(argv, 'zw', int, 7)#3
zw = 2*zw0+1
sbexp = 0
do_eagle = 1
ov_eagle = True
do_muse = 0
ov_muse = False
do_p = 1
vcorr = 1
do_mask = 0
do_gmask = 0
smask = ('.%smask' % (('g%d' % galrad)*do_gmask))*do_mask
do_sclip = 1
nsigma = 3
sclip = ('.sclip%d' % nsigma)*do_sclip
do_csub = False
scsub = '.csub'*do_csub
do_rand = 0
fix_red = 0
zoff = 0
SB = 1.27  # SB at z=3.5

#colors = [(0, 0, 1), (1, 1, 0), (1, 0, 0)]
#cmap = LinearSegmentedColormap.from_list('sofi', colors, N=50)
# 1 spaxel is (0.2 arcsec)^2
# flux2sb comes from 1.25 Angstrom wavelength width times 1/(0.2)**2 from conversion to the correct aperture
flux2sb = 31.25
smodel = '_HM12' * (model == 'HM12')
fdat = 'muse-vs-eagle%s.dat' % scat
_fin = None
_flls, _fcum = {}, {}

xticks = [3, 6, 9, 12, 15, 18, 36, 54, 72]
xts = ['%d' % xt for xt in xticks]
zcolors = ['black', 'darkorchid', 'dodgerblue', 'darkturquoise', 'limegreen', 'khaki', 'orange', 'red', 'brown', 'darkred']
fsize = 14


rads = [[0, 1], [1, 2], [2, 4], [4, 6], [6, 10], [10, 16], [16, 24], [24, 30], [30, 50]]
x = [(r[0] + r[1]) / 2. for r in rads]
nr = len(rads)
asec2pix_eagle = asec2kpcs[11] * (1 + redshifts[11]) * kpc2pix
rpix_eagle = np.array(rads) * asec2pix_eagle
f = []



lognhi = [0, 12, 13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21,
		  22, 23, 30]
nhi = np.power(10, lognhi)
sb = [0, 0, 0.00001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96,
	  0.98, 1, 1, 1, 1, 1, 1]

nhi2sb = inter.interp1d(nhi, sb)
lognhi2sb = inter.interp1d(lognhi, sb)
sb2nhi = inter.interp1d(sb, nhi)

logNHI = [12, 12.875, 13.125, 13.375, 13.625, 13.875, 14.125, 14.375, 14.625, 14.875, 15.25, 15.75, 16.25, 16.75, 17.25,
		  17.75,18,19,20,21,22]
f = [-10.1, -11.077, -11.437, -11.795, -12.163, -12.507, -12.992, -13.459, -13.781, -14.29, -14.933, -15.97, -16.674, -17.219,
	 -17.896, -19.174,-19.5,-20.5,-21.5,-23,-24.5]

nhi2f = inter.interp1d(logNHI, f)
nhirange = np.arange(15, 21, .1)
f2 =nhi2f(nhirange)

sb = lognhi2sb(nhirange)

fdat = '../../fLLS/NHI_prof.fits'
fdat2 = '../../fLLS/NHI_frac.fits'

nhis = np.arange(14, 22, .25)#np.array([14.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5])
_nhi = nhis[:-1]
klls = 2
cat = fits.getdata('%s/EAGLE/cats/gals_snap11.fits' % (folder), 1)
nid = len(cat['ID'])
nr = len(rads)
nh = len(_nhi)

if sbexp:
	fmean = getdata('../../EAGLE/simplestacks/snap11_x.vcorr.NHIcorr.fits')
	zl, yl, xl = fmean.shape
	zz, yy, xx = np.ogrid[0:zl, 0:yl, 0:xl]
	ce = (xx - xl / 2.) ** 2 + (yy - yl / 2.) ** 2
	type = 'NHI'
	if not os.path.isfile(fdat) or overwrite:
		hdu = PrimaryHDU()
		fn = np.zeros((nid, nr, nh))
		nn = np.zeros((nid, nr, nh))
		for i in range(nid):
			print i, nid
			f = fits.getdata(flaes + '/EAGLE/gals/snap11_x/%d%s.%s.fits' % (cat['ID'][i], smodel, type))
			m = fits.getdata(flaes + '/EAGLE/gals/snap11_x/%d%s.gmask6.fits' % (cat['ID'][i], smodel))
			_zl, _xl, _yl = f.shape
			gals = (m > 0) & (m != cat['ID'][i])
			f[gals] = 0
			f = f[_zl/2-zw0: _zl/2+zw0+1, _yl/2-yl/2: _yl/2+yl/2+1, _xl/2-xl/2: _xl/2+xl/2+1]
			for j in range(nr):
				r = rpix_eagle[j]
				cool = (ce >= r[0] ** 2) & (ce < r[1] ** 2)
				scool = float(np.sum(cool))
				for k in range(nh):
					n0 = nhis[k]
					n1 = nhis[k+1]
					fsum = np.sum(f[:, cool[0]], 0)
					a = (fsum >= 10**n0) & (fsum < 10**n1)
					fn[i, j, k] = np.nansum(a)/scool
					#nn[i, j, k] = np.nanmean(fsum[a])
		hdu.data = fn
		hdu.writeto(fdat2, clobber=True)
		#hdu.data = nn
		#hdu.writeto(fdat, clobber=True)
	else:
		fn = getdata(fdat2)
		#nn = getdata(fdat)

	f = np.nanmean(fn, 0)
	#np.log10(np.nanmean(nn, 0))
	ds = .5
	plt.figure(figsize=(6, 6))

	cm = cmap('jet')

	sbs = f*lognhi2sb(_nhi)
	poslls = (_nhi >= 17.5)
	ysum = np.nansum(sbs, 1)
	flls = np.nansum(f[:, poslls], 1)
	sblls = np.nansum(sbs[:, poslls], 1)
	ycool = np.zeros(nr)
	ya = np.zeros(nr)

	norm = BoundaryNorm(_nhi, cm.N)

	colors = [cm(i*9) for i in range(nh)[::-1]]
	labels = [r'$\rm{log(N_{HI})}=%.1f$' % _n for _n in _nhi]
	_y =f*lognhi2sb(_nhi)
	_y[np.isnan(_y)] = 0
	y = np.vstack([_y[:, i] for i in range(nh)[::-1]])
	plt.stackplot(x, y, colors=colors)

	plt.imshow(np.array([[_nhi[0], _nhi[-1]]]), cmap=cm, aspect='auto', vmin=14, vmax=22, norm=norm)
	#plt.plot(x, ysum, label=r'Total SB')
	#plt.plot(x, ycool, label=r'SB LLS', linestyle=':')
	plt.plot(x, flls, label=r'$\rm{f_{LLS}}$', linestyle=':', fontsize=fsize)
	#plt.plot(x, sblls, label=r'$\rm{SB_{LLS}}$', linestyle=':')
	plt.title(r'$z=%.2f$' % redshifts['11'], fontsize=fsize)
	plt.colorbar(ticks=np.arange(14, 22), label=r'$\rm{log(N_{HI})}$', fontsize=fsize)
	plt.legend()
	plt.xlim([6, 40])
	plt.ylim([0, .3])
	plt.xlabel('distance [arcsec]', fontsize=fsize)
	plt.ylabel(r'$\rm{SB/SB_{max}}$', fontsize=fsize)
	plt.savefig('../../Figures/SB_expected_snap11.pdf')
	plt.close()


if do_eagle:
	yg = {}
	ycum = {}
	
	for snap in snaps:
		s = str(snap)
		asec2pix_eagle = asec2kpcs[snap] * (1 + redshifts[snap]) * kpc2pix
		zz, yy, xx = np.ogrid[0:zl, 0:yl, 0:xl]
		#out = (np.abs(zl / 2 - zz) <= zw0) & ((yy - yl / 2) ** 2 + (xx - xl / 2) ** 2 > (5 * asec2pix_eagle) ** 2)
		cat = fits.getdata('%s/EAGLE/cats/gals_snap%s.fits' % (folder, s), 1)
		_x, _y, _z = (cat['x']*lcube/coml).astype(int), (cat['y']*lcube/coml).astype(int), (cat['z']*lcube/coml).astype(int)
		#data_cube = getdata('%s/EAGLE/snap%d_x_%s_512_4096_%d.%s.fits' % (flaes, snap, model, zlens[str(snap)], type))
		#_zl, _yl, _xl = data_cube.shape
		#gcube = getdata('%s/EAGLE/snap%d_x.galmask_%darcsec.fits' % (flaes, snap, galrad))
		m = cat['DM_mass']
		H0 = 69.3
		Ol = 0.714
		Om = 1 - 0.714
		red = redshifts[snap]
		offset = cat['vx'] * (1 + red) / H0 / np.sqrt(Om * (1 + red) ** 3 + Ol)
		offpix = np.round(offset*zlens[snap]/coml).astype(int)

		if do_p:
			pnames = ['U', 'SFR', 'Mmean200', 'n5d']
			dolog = [0, 0, 1, 0]
			prop = {}
			pmed = {}
			phist = {}
			#old bins, just suitable for snap11
			bins = [[-22, -20, -19, -17, -16], [-3, -.8, -.4, -.2, 0, .5, 3],
					[10.44616846, 10.78110581, 11.11604317, 11.45098053, 11.78591788]]
			for p, dl in zip(pnames, dolog):
				if p == 'U':
					p0 = cat[p]
					noU = np.isnan(p0) | np.isinf(p0)
					p1 = -2.5*np.log10(cat['SFR'])-18.6
					prop[p] = p0
					prop[p][noU] = p1[noU]
				else: prop[p] = cat[p]
				if dl: _p = np.log10(prop[p])
				else: _p = prop[p]
				m0 = np.nanmedian(_p)
				pmin = np.nanmin(_p)
				pmax = np.nanmax(_p)
				phi = _p > m0
				plo = _p < m0
				mhi = np.nanmedian(_p[phi])
				mlo = np.nanmedian(_p[plo])
				pmed[p] = [pmin, mlo, m0, mhi, pmax]
				phist[p] = np.histogram(_p, pmed[p])

		cename = '../../EAGLE/simplestacks/snap%s_x%s%s.%s.fits' % (s, ('.gmasked%d' % galrad)*do_gmask, '.vcorr'*vcorr, type)
		gname = '../../EAGLE/simplestacks/snap%s_x.galcov.fits' % s
		docalc = not os.path.isfile(cename)
		hdu = PrimaryHDU()
		
		if do_p:
			for p, dl in zip(pnames, dolog):
				if dl:
					_p = np.log10(prop[p])
				else:
					_p = prop[p]
				_pm = pmed[p]
				nbins = len(_pm)
				for i in range(nbins - 1):
					docalc |= not os.path.isfile(cename.replace('.fits', '.%s_%.2f_%.2f.fits' % (p, _pm[i], _pm[i + 1])))
	
			if docalc or overwrite:
				stack = []
				mstack = []
				std = []
				for p, dl in zip(pnames, dolog):
					if dl: _p = np.log10(prop[p])
					else: _p = prop[p]
					_pm = pmed[p]
					nbins = len(_pm)
					for i in range(nbins-1):
						_fout = cename.replace('.fits', '.%s_%.2f_%.2f.fits' % (p, _pm[i], _pm[i + 1]))
						if not os.path.isfile(_fout) or overwrite:
							print 'bins', p, _pm[i], _pm[i+1]
							cool = (_p > _pm[i]) & (_p < _pm[i+1])
							#_st, _mt, _nt = np.zeros([zl, yl, xl]), np.zeros([zl, yl, xl]), np.zeros([zl, yl, xl])
							_st, _mt = [], []
							for j, off, x, y, z in zip(cat['ID'][cool], offpix[cool], _x[cool], _y[cool], _z[cool]):
								print j, off
	
								#roll = np.roll(data_cube, (_zl/2-z, _yl/2-y, _xl/2-x), (0, 1, 2))
								#_s = roll[_zl/2-zl/2: _zl/2+zl/2+1, _yl/2-yl/2: _yl/2+yl/2+1, _xl/2-xl/2: _xl/2+xl/2+1]
	
								_s = fits.getdata(flaes + '/EAGLE/gals/snap%s_x/%d%s.%s.fits' % (s, j, smodel, type))
								zz, yy, xx = _s.shape
								_s = _s[zz/2-zw0:zz/2+zw0+1]
								if do_gmask:
									#roll = np.roll(gcube, (_zl/2-z, _yl/2-y, _xl/2-x), (0, 1, 2))
									#_cube = roll[_zl/2-zl/2: _zl/2+zl/2+1, _yl/2-yl/2: _yl/2+yl/2+1, _xl/2-xl/2: _xl/2+xl/2+1]
									#_cube[_cube == i] = 0
									#_m = _cube
									_m = fits.getdata(flaes + '/EAGLE/gals/snap%s_x/%d.gmask%d.fits' % (s, j, galrad))
									zz, yy, xx = _m.shape
									_m = _m[zz / 2 - zw0:zz / 2 + zw0 + 1]
									#_mt += _m
									_mt.append(_m)
									gal = (_m > 0)# & (_m != j)
									#_nt += ~gal
									_s[gal] = 0
									del _m
								#_st += _s
								_st.append(_s)
								del _s
							#_st /= _nt
							f = np.mean(_st, 0)
							#_std = np.std(_st)
							m = np.sum(_mt, 0)
							stack.append(f)
							mstack.append(m)
							hdu.data = f
							hdu.writeto(_fout, clobber=True)
							#hdu.data = _std
							#gnamehdu.writeto(_fout.replace('.fits', '.STD.fits'), clobber=True)
						else:
							stack.append(fits.getdata(_fout))
			else:
				stack = fits.getdata(cename)
				mstack = fits.getdata(gname)
						#std.append(fits.getdata(_fout.replace('.fits', '.STD.fits')))
			
			if 0:
				stack = np.array(stack)
				print 'snap', snap, 'stack shape', stack.shape
				_stack = np.nanmean(stack, 0)
				hdu.data = _stack
				hdu.writeto(cename, clobber=True)
				# _std = np.nanmean(std, 0)
				#hdu.writeto(cename.replace('.fits', '.STD.fits'), clobber=True)
				
				if do_gmask:
					mstack = np.array(mstack)
					_mstack = np.nansum(mstack, 0)
					_mstack2 = np.nanmean(mstack > 0, 0)
					hdu.data = _mstack
					hdu.writeto(gname.replace('.galcov.fits', '.gals.fits'), clobber=True)
					hdu.data = _mstack2
					hdu.writeto(gname, clobber=True)

		else:
			if not os.path.isfile(cename) or overwrite:
				# _st, _mt, _nt = np.zeros([zl, yl, xl]), np.zeros([zl, yl, xl]), np.zeros([zl, yl, xl])
				_st, _mt = [], []
				for j, off, x, y, z in zip(cat['ID'], offpix, _x, _y, _z):
					print j, off
					
					# roll = np.roll(data_cube, (_zl/2-z, _yl/2-y, _xl/2-x), (0, 1, 2))
					# _s = roll[_zl/2-zl/2: _zl/2+zl/2+1, _yl/2-yl/2: _yl/2+yl/2+1, _xl/2-xl/2: _xl/2+xl/2+1]
					
					_s = fits.getdata(flaes + '/EAGLE/gals/snap%s_x/%d.%s.fits' % (s, j, type))
					zz, yy, xx = _s.shape
					_s = _s[zz / 2 - zw0:zz / 2 + zw0 + 1, yy/2-yl/2:yy/2+yl/2+1, xx/2-xl/2:xx/2+xl/2+1]
					if do_gmask:
						_m = fits.getdata(flaes + '/EAGLE/gals/snap%s_x/%d.gmask%d.fits' % (s, j, galrad))
						zz, yy, xx = _m.shape
						_m = _m[zz / 2 - zw0:zz / 2 + zw0 + 1, yy/2-yl/2:yy/2+yl/2+1, xx/2-xl/2:xx/2+xl/2+1]
						# _mt += _m
						_mt.append(_m)
						gal = (_m > 0)  # & (_m != j)
						# _nt += ~gal
						_s[gal] = 0
						del _m
					# _st += _s
					_st.append(_s)
					del _s
				# _st /= _nt
				stack = np.mean(_st, 0)
				hdu.data = stack
				hdu.writeto(cename, clobber=True)
				
				if do_gmask:
					_mt = np.array(_mt)
					_mstack = np.nansum(_mt, 0)
					_mstack2 = np.nanmean(_mt > 0, 0)
					hdu.data = _mstack
					hdu.writeto(gname.replace('.galcov.fits', '.gals.fits'), clobber=True)
					hdu.data = _mstack2
					hdu.writeto(gname, clobber=True)
					
			else:
				stack = fits.getdata(cename)
				mstack = fits.getdata(gname)
		
		rads = [[0, 1], [1, 2], [2, 4], [4, 6], [6, 10], [10, 16], [16, 24], [24, 30], [40, 50], [50, 60],
				[60, 70], [70, 80]]
		nr = len(rads)

		rpix_eagle = np.array(rads) * asec2pix_eagle

		stack = fits.getdata(cename)
		gstack = fits.getdata(gname)
		#std = fits.getdata(cename.replace('.fits', '.STD.fits'))
		zle, yle, xle = stack.shape
		yy, xx = np.ogrid[0: yle, 0: xle]
		ce = (xx - xle / 2.) ** 2 + (yy - yle / 2.) ** 2

		y, ym, ystd = {}, {}, {}
		x = [(r[0] + r[1]) / 2. for r in rads]
		_xl = [.9, 300]
		_yl = [5e-3, .3]
		
		if do_p:
			label = [r'$\rm{M_U}$', r'$\rm{SFR\,[M_\odot\,yr^{-1}]}$',
					 r'$\rm{M_{200}\,[10^{10}\,M_\odot]}$', r'$\rm{d_{5th}\,[cMpc]}$']
			colors = ['red', 'orange', 'green', 'pink']#, 'blue', 'purple']
			coords = [(0, 0), (0, 1), (1, 0), (1, 1)]
			
			if do_gmask:
				_yg = []
				_ycum = []
				for r in rpix_eagle:
					cool = (ce >= r[0] ** 2) & (ce < r[1] ** 2)
					_yg.append(np.nanmean(gstack[zle / 2 - zw0:zle / 2 + zw0 + 1, cool]))
				yg[snap] = _yg
	
			plt.figure()
			plt.rc('font', size=fsize)
			plt.rc('xtick', labelsize=fsize-2)  # fontsize of the tick labels
			plt.rc('ytick', labelsize=fsize-2)
			fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(9, 7))
			fig.tight_layout()
			fig.suptitle(r'snapshot %s, $z = %.2f$' % (s, redshifts[snap]), x=.5, y=.98)#, fontsize=fsize)
			plt.gcf().subplots_adjust(left=.1, bottom=.08, top=.95)
			ax[0, 0].set_ylabel(r'$\rm{(f_{LLS}-\langle f_{LLS}\rangle)/\langle f_{LLS}\rangle}$')#, fontsize=fsize)
			ax[1, 0].set_ylabel(r'$\rm{(f_{LLS}-\langle f_{LLS}\rangle)/\langle f_{LLS}\rangle}$')#, fontsize=fsize)
			ax[1, 0].set_xlabel('distance [arcsec]')#, fontsize=fsize)
			ax[1, 1].set_xlabel('distance [arcsec]')#, fontsize=fsize)
			for p, l, c in zip(pnames, label, coords):
				ax[c].text(15, .8, l)
				print p, l, c
				ph = np.array(pmed[p])
				_p = (ph[1:]+ph[:-1])/2.
				mat = []
				mat.append(x)
				header = 'radius'
				prange = np.arange(len(ph) - 1)
				if p == 'Mmean200' or p == 'SFR':
					_colors = colors[::-1]
					prange = prange[::-1]
				else:
					_colors = colors
				dashes = prange + 1
				for i in prange:
					ax[c].plot(x, [0] * len(x), color='lightgrey')
					lab = '%s_%d' % (p, _p[i])
					y[lab] = []
					ym[lab] = []
					#ystd[lab] = []
					st = fits.getdata(cename.replace('.fits', '.%s_%.2f_%.2f.fits' % (p, ph[i], ph[i + 1])))
					if np.nansum(st) > 0:
						for r in rpix_eagle:
							cool = (ce >= r[0] ** 2) & (ce < r[1] ** 2)
							_y = np.nansum(np.nanmean(st[zle/2-zw0:zle/2+zw0+1, cool], 1))
							_ym = np.nansum(np.nanmean(stack[zle/2-zw0:zle/2+zw0+1, cool], 1))
							#_ystd = np.nansum(np.nanmean(std[zle/2-zw0:zle/2+zw0+1, cool], 1))
							if type == 'LLS' or type == 'LLScorr':
								if np.isnan(_y): _y = 0.
								if np.isnan(_ym): _ym = 0.
								if _y > 1: _y = 1
								if _ym > 1: _ym = 1
							y[lab].append(_y)
							ym[lab].append(_ym)
							#ystd[lab].append(_ystd)
						_l = ''
	
						_y = np.array(y[lab])
						_ym = np.array(ym[lab])
						ydiff = (_y-_ym)/_ym
						if p == 'Mmean200': _l = r'$%.1f$' % (10 ** (_p[i]-10))
						elif p == 'n5d': _l = r'$%.1f$' % (_p[i]/1000.)
						else: _l = r'$%.1f$' % _p[i]
						ax[c].plot(x, ydiff, color=_colors[i], label=_l, dashes=(5, dashes[i]))#, fontsize=fsize)
						mat.append(y[lab])
						header += ' '+lab
				flls = np.nanmean(mat[1:], 0)
				sb = sbpeaks[snap]*flls
				mat.append(flls)
				mat.append(sb)
				header += ' fLLS SB_HM12'
				#np.savetxt('%s/fLLS/%s_prof-%s_snap%s%s.dat' %
				#           (folder, type, p, s, smodel), np.array(mat).T, header=header)
	
				ax[c].legend(loc=3, fancybox=False, framealpha=0.2, ncol=2)#, fontsize=fsize)
				ax[c].set_xlim([.5, 25])
				ax[c].set_xticks([5, 10, 15, 20])#, fontsize=fsize)
				ax[c].set_ylim([-1, 1])
				ax[c].set_yticks([-.9, -.6, -.3, 0, .3, .6, .9])#, fontsize=fsize)
				if 0:#p == 'U' or p == 'SFR':
					ax[c].twiny()
					ax[c].set_xlim([.5, 75])
					kpc = np.array([10, 100, 200, 300, 400, 500])
					if snap > 12: kpc = np.array([10, 100, 200, 300, 400, 500, 600])
					_kpc = kpc / asec2kpcs[snap]
					ax[c].set_xticks(_kpc)
					ax[c].set_xticklabels(kpc.astype('int'))
					ax[c].set_xlabel('distance [pkpc]')
	
			plt.savefig('%s/fLLS/%s_prof_snap%s.pdf' % (folder, type, s))
			plt.close()
			_flls[snap] = flls
		
		else:
			flls = []
			fcum = []
			for r in rpix_eagle:
				cool = (ce >= r[0] ** 2) & (ce < r[1] ** 2)
				coolc = (ce <= r[0] ** 2)
				_stack = stack[zle / 2 - zw0:zle / 2 + zw0 + 1, cool]
				_stackc = stack[zle / 2 - zw0:zle / 2 + zw0 + 1, coolc]
				if type == 'NHI':
					_stack = np.nansum(_stack, 0) > 10**17.5
					_stackc = np.nansum(_stackc, 0) > 10**17.5
					flls.append(np.nanmean(_stack))
					fcum.append(np.nanmean(_stackc))
				else:
					flls.append(np.nansum(np.nanmean(_stack, 1)))
					fcum.append(np.nansum(np.nanmean(_stackc, 1)))
			flls = np.array(flls)
			fcum = np.array(fcum)
			flls[flls > 1] = 1
			fcum[fcum > 1] = 1
			_flls[snap] = flls
			_fcum[snap] = fcum

	plt.figure()

	i=0
	for snap in snaps[::-1]:
		s = '%s' % snap
		y = np.copy(_flls[snap])
		if type=='LLS' or type == 'LLScorr': y[y>1] = 1
		if type == 'NHI':
			y[y<10] = 1
			y = np.log10(y)
		plt.plot(x, y, label=r'$z = %.2f$' % redshifts[snap], color=zcolors[i], dashes=(5, i+1))
		i += 1
		#plt.scatter(x, y)
	if type == 'LLS' or type == 'LLScorr':
		plt.semilogy()
		plt.xlim([0.5, 25])
	if type == 'NHI':
		plt.xlim([0.5, 27])

	#plt.xticks([2, 6, 20, 60], ['2', '6', '20', '60'])
	#plt.xticks(xticks, xts)
	plt.ylim([.01, 1.1])
	plt.yticks([.01, .03, .1, .3, 1], ['0.01', '0.03', '0.1', '0.3', '1.0'])
	plt.ylabel(r'f$_{\rm{LLS}}$')

	plt.xlabel('distance [arcsec]')

	#plt.yticks([.1, .3, 1], ['0.1', '0.3', '1.0'])
	plt.legend()
	plt.tight_layout()
	plt.savefig('%s/fLLS/%sprof%s.pdf' % (folder, type, smodel))
	plt.close()
	
	i=0
	
	
	if 0:
		for snap in snaps:
			s = snap
			yc = np.copy(_fcum[snap])
			if type=='LLS' or type == 'LLScorr': y[y>1] = 1
			if type == 'NHI':
				yc[yc<10] = 1
				yc = np.log10(yc)
			plt.plot(x, yc, label=r'$z = %.2f$' % redshifts[snap], color=zcolors[i], dashes=(5, i+1))
			i += 1
			#plt.scatter(x, y)
		plt.semilogy()
		plt.xlim([0.5, 25])
	
	
		#plt.xticks([2, 6, 20, 60], ['2', '6', '20', '60'])
		#plt.xticks(xticks, xts)
		plt.ylim([.01, 1.1])
		plt.yticks([.01, .03, .1, .3, 1], ['0.01', '0.03', '0.1', '0.3', '1.0'])
		plt.ylabel(r'f$_{\rm{LLS}}$')
	
		plt.xlabel('distance [arcsec]')
	
		#plt.yticks([.1, .3, 1], ['0.1', '0.3', '1.0'])
		plt.legend()
		plt.tight_layout()
		plt.savefig('%s/fLLS/%sprofcum%s.pdf' % (folder, type, smodel))
		plt.close()

	
	plt.figure()
	plt.rc('font', size=fsize)
	i = 0
	for snap in snaps[::-1]:
		s = str(snap)
		y = np.copy(_flls[snap])
		y[y>1] = 1
		plt.plot(x, y*sbpeaks[snap], label=r'$z = %.2f$' % redshifts[snap], color=zcolors[i], linestyle='--', dashes=(5, i+1))#, fontsize=fsize)
		i += 1
		#plt.scatter(x, y*sbpeaks[snap])
	
	plt.semilogy()
	xmin, xmax = 0.5, 25
	plt.xlim([xmin, xmax])
	xmax = 7
	#plt.xticks([2, 6, 20, 60], ['2', '6', '20', '60'])
	plt.ylim([.01, 5])
	plt.yticks([.01, .03, .1, .3, 1, 3], ['0.01', '0.03', '0.1', '0.3', '1.0', '3.0'], fontsize=fsize)
	if 0:
		SBs = 1.3, .8, .7*2/np.sqrt(5), .8/np.sqrt(5), 8e3/np.sqrt(800000)
		labels = 'UDF-mosaic', 'UDF-10', 'KCWI 5h', 'MXDF', 'HETDEX'
		for ii in range(len(SBs)):
			plt.plot([xmin, xmax], [SBs[ii], SBs[ii]], color='grey')
			plt.text(xmax, SBs[ii], labels[ii], fontsize=fsize-2)
	plt.xlabel('distance [arcsec]', fontsize=fsize)
	plt.ylabel(r'SB$\rm{_{HM12} [10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize)
	#plt.yticks([.1, .3, 1], ['0.1', '0.3', '1.0'])
	plt.legend()
	plt.tight_layout()
	plt.savefig('%s/fLLS/SBHM12_%sprof%s.pdf' % (folder, type, smodel))
	plt.close()


	if 0:
		plt.figure()
		i = 0
		for snap in snaps[::-1]:
			s = str(snap)
			plt.plot(x, yg[snap], label=r'$z = %.2f$' % redshifts[snap], color=zcolors[i], linestyle='--', dashes=(5, i+1))
			i += 1
			#plt.scatter(x, y*sbpeaks[snap])
		plt.xlabel('distance [arcsec]')
		plt.ylabel(r'$\rm{f_{CGM}}$')
		plt.xlim([.5, 75])
		plt.legend()
		plt.savefig('%s/fLLS/fCGM_%s_prof%s.pdf' % (folder, type, smodel))
		plt.close()

if do_muse:
	snr = True
	galcov = False
	pnames = ['U', 'n5d']
	colors = 'darkorchid', 'blue', 'green', 'gold', 'red'
	colors = 'blue', 'green', 'gold', 'red'
	markers = "H", "d", "s", "^", "p", "o", "h", "D", "d"
	asec2pix_muse = 5  # muse
	#zps = [[2.8, 3.2], [3.2, 3.8], [3.8, 4.4], [4.4, 4.6], [4.6, 5.6]]
	zps = [[2.9, 3.4], [3.4, 4.5], [4.5, 5.5], [5.5, 6.5]]#[[2.8, 3.4], [3.4, 4.4], [4.4, 5.4], [5.4, 6.6]]
	zoffs = [-1, -1, -1, -1]#[-4, -4, -5, -6]
	rads = [[0, .2], [0, .5], [.5, 1], [1, 1.5], [1.5, 2], [2, 3], [3, 4], [4, 6], [6, 9], [9, 12], [12, 15],
			[15, 18], [18, 24], [24, 30]]
	x = [(r[0] + r[1]) / 2. for r in rads]
	plt.figure()
	if not galcov: plt.plot(x, np.zeros(len(x)), color='grey')
	rpix = np.array(rads) * asec2pix_muse
	i = 0
	zrange = [2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
       19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
       36, 37, 38, 39, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
       73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
       90, 91, 92, 93, 94, 95, 96, 97, 98]

	_ys = []

	for zp in zps:
		odat = '../../all/radprof%.1f-%.1f.dat' % (zp[0], zp[1])
		if not os.path.isfile(odat) or overwrite:
			stack = getdata('../../all/simplestacks/stack_UDF_mosaic.csub.g6mask.sclip3.corr.z%.2f-%.2f.fits' % (zp[0], zp[1]))
			gstack = getdata('../../all/simplestacks/gstack_UDF_mosaic.z%.2f-%.2f.fits' % (zp[0], zp[1]))
			zl, yl, xl = stack.shape
			zz, yy, xx = np.ogrid[0: zl, 0: yl, 0: xl]
			c = (xx - xl / 2.) ** 2 + (yy - yl / 2.) ** 2
			y = []
			ys = []
			yr = []
			yg = []
			for r in rpix:
				cool = (c >= r[0] ** 2) & (c < r[1] ** 2)
				_fr = []
				for _zz in zrange:
					zmin = max(0, zl / 2 - zw0 + _zz)
					zmax = min(zl, zl / 2 + zw0 + _zz + 1)
					fm = stack[_zz - zw0: _zz + zw0 + 1, cool[0]]
					if do_sclip: fm[abs(fm) > nsigma * np.nanstd(fm)] = np.nan
					_fr.append(np.nansum(np.nanmean(fm, 1)))
				fr = np.nanmean(_fr)
				fstd = np.nanstd(_fr)
				f = stack[zl/2 - zw0+zoff: zl/2 + zw0 + 1+zoff, cool[0]]
				g = gstack[zl / 2 - zw0+zoff: zl / 2 + zw0 + 1+zoff, cool[0]]
				if do_sclip: f[abs(f) > nsigma * np.nanstd(f)] = np.nan
				y.append((np.nansum(np.nanmean(f, 1)))*flux2sb)
				ys.append(fstd*flux2sb)
				yr.append(fr*flux2sb)
				yg.append(np.nanmean(g))
			m = np.array([x, y, ys, yr, yg])
			np.savetxt(odat, m, header='radius SB SB_std SB_rand gal_cov')

		else:
			m = np.loadtxt(odat)
		y = np.array(m[1])
		ys = np.array(m[2])
		_ys.append(ys)
		yr = np.array(m[3])
		yg = np.array(m[4])
		filt = 5, 1
		if snr:
			plt.plot(x, sp.savgol_filter((y-yr)/ys,filt[0],filt[1]), color=colors[i])
			plt.scatter(x, sp.savgol_filter((y-yr)/ys,filt[0],filt[1]), color=colors[i], label=r'$%.1f<z<%.1f$'% (zp[0], zp[1]), marker=markers[i])
		elif galcov:
			plt.plot(x, sp.savgol_filter(yg,filt[0],filt[1]), color=colors[i])
			plt.scatter(x, sp.savgol_filter(yg,filt[0],filt[1]), color=colors[i], label=r'$%.1f<z<%.1f$'% (zp[0], zp[1]), marker=markers[i])

		else:
			plt.plot(x, sp.savgol_filter(y-yr,filt[0],filt[1]), color=colors[i])
			plt.scatter(x, sp.savgol_filter(y-yr,filt[0],filt[1]), color=colors[i], label=r'$%.1f<z<%.1f$'% (zp[0], zp[1]), marker=markers[i])

		#plt.plot(x, y, label=r'$%.1f<z<%.1f$'% (zp[0], zp[1]), color=colors[i])
		#plt.scatter(x, y , color = colors[i])
		#plt.fill_between(x, y-ys/2., y+ys/2., facecolor=colors[i], alpha=0.25, lw=0, edgecolor='none')
		#plt.errorbar(x, y, yerr=ys, xerr=0, fmt="none", capsize=4, color = colors[i])
		i += 1
	plt.legend()

	if not snr and not galcov:
		ys = np.nanmean(_ys, 0)+np.nanstd(_ys, 0)
		plt.fill_between(x, -2*ys, 2*ys, facecolor='lightgrey', alpha=0.25, lw=0, edgecolor='none', zorder=0)
		plt.fill_between(x, -ys, ys, facecolor='grey', alpha=0.2, lw=0, edgecolor='none', zorder=0)
	#plt.plot(x, ys, color='grey')
	#plt.plot(x, sp.savgol_filter(y - yr, 5, 1), label=r'$%.1f<z<%.1f$' % (zp[0], zp[1]), color=colors[i])
	plt.xlim([4, 20])
	xticks = [6, 9, 12, 15, 18]
	plt.xticks(xticks, ['%d' % xt for xt in xticks], fontsize=fsize)
	plt.xlabel('distance [arcsec]', fontsize=fsize)

	if snr:
		plt.ylabel(r'SNR')
		plt.ylim([-2, 6])
	elif galcov:
		plt.ylabel(r'f$_{\rm{CGM}}$', fontsize=fsize)
		plt.ylim([0, .2])
	else:
		plt.ylabel(r'SB $\rm{[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize)
		plt.ylim([-1, 3])
	#plt.yscale('symlog', linthreshy=0.1)

	#plt.yticks([-.1, 0, .1, 1], ['-0.1', '0', '0.1', '1.0'])
	if snr: plt.savefig('../../Figures/SNRprof_MUSE.pdf')
	elif galcov: plt.savefig('../../Figures/fCGMprof_MUSE.pdf')
	else: plt.savefig('../../Figures/SBprof_MUSE.pdf')
	plt.close()

