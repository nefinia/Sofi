from typing import Any, Union

__author__ = 'gallegos'

from pyfits import getdata, PrimaryHDU
import numpy as np
import os
from tools_sofi import astroim
from sys import argv
import matplotlib.mlab as mlab
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.io import fits
from tools_sofi import rdarg, sclipping
import params

fitcat = rdarg(argv, key='fitcat', type=str, default='all')  # 'all'#'HDFS'
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
variables = ['Luminosity', 'Environment', 'Redshift', "Com_Distance", "Proj_Distance"]
sbplot = rdarg(argv, 'sbprof', bool, False)
pitheta = rdarg(argv, 'pitheta', bool, False)
pithetaEAGLE = rdarg(argv, 'pithetaEAGLE', bool, False)
halfplot = rdarg(argv, 'halfplot', bool, False)
lutz = rdarg(argv, 'lutz', bool, False)
paperims = rdarg(argv, 'paperims', bool, False)
subcubes = rdarg(argv, 'subcubes', bool, False)
dndz = rdarg(argv, 'dndz', bool, False)
hm = rdarg(argv, 'hm', bool, False)
circhist = rdarg(argv, 'circhist', bool, False)

init = ''  # 'not-normalized_masked_cubical_'
n = len(variables)
medians = [91, 9, 3.5, 7, 32]
mins = [1, 1, 2.9, 0, 10]
maxs = [2000, 100, 4.0, 20, 80]
name = 'stack.avs.dat'
extraname = ''  # '.reject.abitmore'#extraname = ''#
offset = 0
zw0 = 2
xmin = 46
xmax = 52
ymin = 35
ymax = 45
zmin = 10
zmax = 12
dpi = 100

envplot = False
ksplot = False
radprof = False
zneig = False
radangplot = False
spectra = False


def histplot(datos, drand, dmix, hname, title=None, gaussian=None):
	n, bins, patches = plt.hist(dmix, 16, normed=1, facecolor='blue', alpha=0.75)
	fsize = 20
	plt.tick_params(labelsize=fsize)
	plt.xlabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize)
	plt.savefig(hname.replace(".png", "_full.png"))
	plt.close()
	n, bins, patches = plt.hist(drand, 16, normed=1, facecolor='grey', alpha=0.75, label='Random Subsets')
	bins = np.append(bins, (np.arange(5) + 1) * abs(bins[1] - bins[0]) + bins[-1])
	print bins
	n, bins, patches = plt.hist(datos, bins, normed=1, facecolor='green', alpha=0.75, label='Oriented Subsets')
	if gaussian:
		(mu, sigma) = norm.fit(datos)
		(mu2, sigma2) = norm.fit(drand)
		y = mlab.normpdf(bins, mu, sigma)
		y2 = mlab.normpdf(bins, mu2, sigma2)
		l = plt.plot(bins, y, 'r--', linewidth=2)
		l2 = plt.plot(bins, y2, 'r--', linewidth=2)
	if title:
		plt.title(title)
	plt.legend(prop={'size': 10})
	print hname
	plt.savefig(hname, dpi=dpi, format='pdf')
	plt.close()


if pitheta:
	cat1 = '%s/%s/cats/lae_pairs.fits' % (folder, 'UDF')
	cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'HDFS')
	data1 = getdata(cat1, 1)
	data2 = getdata(cat2, 1)
	figx, figy = [10, 8]
	ym = 30
	yl = 20
	theta0 = 16
	plt.figure(figsize=(figx, figy))
	plt.ylim(-1, ym)
	colors = [(0, 0.4, 1), (.5, 0, .5), (1, 0.1, 0)]  # B -> P -> G
	cmap_name = 'my_list'
	cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=5)
	xmin = -30
	xmax = 600
	fsize = 23
	vmin = 2.86
	vmax = 4.04
	plt.xlim(xmin, xmax)
	plt.xlabel(r'$r_p\,\rm{[pkpc]}$', fontsize=fsize)
	plt.ylabel(r'$\pi\,\rm{[cMpc]}$', fontsize=fsize)
	pi = np.concatenate([np.array(data1['pi_mpc']), np.array(data2['pi_mpc'])])
	piv = np.abs(np.concatenate([np.array(data1['pi_v']), np.array(data2['pi_v'])]))
	theta = np.concatenate([np.array(data1['theta']), np.array(data2['theta'])])
	pd = np.concatenate([np.array(data1['proj_pdist']), np.array(data2['proj_pdist'])]) * 1000.
	z1 = np.concatenate([np.array(data1['redshift1']), np.array(data2['redshift1'])])
	z2 = np.concatenate([np.array(data1['redshift2']), np.array(data2['redshift2'])])
	zs = (z1 + z2) / 2.
	close = (pi < yl) & (pi > .5) & (z1 < 4.) & (z2 < 4.) & (theta > theta0)
	far = (((pi > yl) | (pi < .5)) & (pi < ym)) | (theta < theta0)
	x = pd[close]
	y = pi[close]
	z = zs[close]
	x2 = pd[far]
	y2 = pi[far]
	plt.scatter(x2, y2, 50, marker='.', lw=0, vmin=vmin, vmax=vmax, c='gray', alpha=.5)
	plt.scatter(x, y, 100, z, '*', lw=0, vmin=vmin, vmax=vmax, cmap=cm, alpha=.8)
	plt.hlines(20, xmin, xmax, linestyles='dashed', lw=2)
	plt.hlines(.5, xmin, xmax, linestyles='dashed', lw=2)
	cbar = plt.colorbar()
	cbar.set_label(r'$\rm{redshift}$', fontsize=fsize, labelpad=0.6)
	plt.savefig("../../Figures/pi_theta.pdf", dpi=dpi, format='pdf')
	plt.show()
	plt.close()
	ymin = -50
	ymax = 1800
	plt.figure(figsize=(figx, figy))
	plt.ylim(ymin, ymax)
	plt.xlim(-5, 80)
	plt.xlabel(r'$\theta\,\rm{[arcsec]}$', fontsize=fsize)
	plt.ylabel(r'$\Delta v\,\rm{[km/s]}$', fontsize=fsize)
	close = (pi < yl) & (pi > .5) & (z1 < 4.) & (z2 < 4.) & (theta > theta0)
	far = (((pi > yl) | (pi < .5)) & (pi < ym)) | (theta < theta0)
	y = piv[close]
	x = theta[close]
	z = zs[close]
	y2 = piv[far]
	x2 = theta[far]
	plt.scatter(x2, y2, 50, marker='.', lw=0, vmin=vmin, vmax=vmax, c='gray', alpha=.5)
	plt.scatter(x, y, 100, z, '*', lw=0, vmin=vmin, vmax=vmax, cmap=cm, alpha=.8)
	plt.vlines(theta0, ymin, ymax, linestyles='dashed', lw=2)
	cbar = plt.colorbar()
	cbar.set_label(r'$\rm{redshift}$', fontsize=fsize, labelpad=0.6)
	plt.savefig("../../Figures/pi_theta_obs.pdf", dpi=dpi, format='pdf')
	plt.show()

if pithetaEAGLE:
	colors = [(0, 0.4, 1), (.5, 0, .5), (1, 0.1, 0)]  # B -> P -> G
	cmap_name = 'my_list'
	cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=3)
	fsize = 15
	figx, figy = [7, 6]
	msize = 60
	h0 = 69.6
	omega_l = 0.714
	coords = ['x', 'y', 'z']
	redshifts = [3.984, 3.528, 3.017]
	colors = ['red', 'green', 'blue']
	theta0 = 16
	lcube = 4096.
	coml = 25.  # cMpc
	pix2com = coml / lcube
	ym = 30
	yl = 20


	def plot(x, y, z, x2, y2, outname, xmin, xmax, ymin, ymax,
			 xlabel, ylabel, hlines=False, vlines=False, show=False):
		plt.figure(figsize=(figx, figy))
		plt.ylim(ymin, ymax)
		plt.xlim(xmin, xmax)
		plt.xlabel(xlabel, fontsize=fsize)
		plt.ylabel(ylabel, fontsize=fsize)
		for red, color in zip(redshifts, colors):
			cool = z == red
			plt.scatter(x[cool], y[cool], msize, marker='*', lw=0, c=color, alpha=.8, label='z=%.3f' % red)
		plt.scatter(x2, y2, msize, marker='.', lw=0, c='gray', alpha=.5)
		if hlines:
			plt.hlines(20, xmin, xmax, linestyles='dashed', lw=2)
			plt.hlines(.5, xmin, xmax, linestyles='dashed', lw=2)
		if vlines:
			plt.vlines(theta0, ymin, ymax, linestyles='dashed', lw=2)
		plt.savefig(outname, dpi=dpi, format='pdf')
		if show: plt.show()
		plt.close()


	cat = '%s/%s/cats/lae_pairs.fits' % (folder, 'EAGLE')
	data = getdata(cat, 1)

	# v1, v2, pi, pi2, piv, piv2, theta, pd, close, sclose, far, sfar = [{}]*12
	zs = data['redshift']
	c = coords[0]
	if 1:
		co = ['x', 'y', 'z']
		co.remove(c)
		p1, p2 = [data['%st' % co[0]], data['%st' % co[1]]]
		pd = np.sqrt(p1 ** 2 + p2 ** 2) * pix2com / (1 + zs)
		v1 = data['v%s1' % c]
		v2 = data['v%s2' % c]
		pi = data['%st' % c] * pix2com
		pi2 = pi + (v2 - v1) / h0 / np.sqrt((1 - omega_l) * (1 + zs) ** 3 + omega_l)
		piv = pi * h0
		piv2 = piv + (v2 - v1)
		theta = data['theta%s' % c]
		close = (pi < yl) & (pi > .5) & (theta > theta0)
		close2 = (pi2 < yl) & (pi2 > .5) & (theta > theta0)
		nclose = np.sum(close)
		sclose = np.random.choice(range(nclose), int(nclose * .02))
		far = (((pi > yl) | (pi < .5)) & (pi < ym)) | (theta < theta0)
		far2 = (((pi2 > yl) | (pi2 < .5)) & (pi2 < ym)) | (theta < theta0)
		nfar = np.sum(far)
		sfar = np.random.choice(range(nfar), int(nfar * .02))

	print nclose, np.sum(close2), nfar, np.sum(far2)

	x = pd[close][sclose]
	y = pi[close][sclose]
	z = zs[close][sclose]
	x2 = pd[far][sfar]
	y2 = pi[far][sfar]
	outname = '../../Figures/pi_theta_EAGLE.pdf'
	xlabel = r'$r_p\,\rm{[pMpc]}$'
	ylabel = r'$\pi\,\rm{[cMpc]}$'
	xmin = 0
	xmax = 8
	ymin = -.5
	ymax = 23
	plot(x, y, z, x2, y2, outname, xmin, xmax, ymin, ymax,
		 xlabel, ylabel, hlines=True, vlines=False, show=False)

	y = pi2[close][sclose]
	y2 = pi2[far][sfar]
	outname = '../../Figures/pi_theta_galv_EAGLE.pdf'
	plot(x, y, z, x2, y2, outname, xmin, xmax, ymin, ymax,
		 xlabel, ylabel, hlines=True, vlines=False, show=False)

	xlabel = r'$\theta\,\rm{[arcsec]}$'
	ylabel = r'$\Delta v\,\rm{[km/s]}$'
	ymin = -200
	ymax = 1600
	xmin = 0
	xmax = 550
	y = piv[close][sclose]
	x = theta[close][sclose]
	y2 = piv[far][sfar]
	x2 = theta[far][sfar]
	outname = '../../Figures/piv_theta_EAGLE.pdf'
	plot(x, y, z, x2, y2, outname, xmin, xmax, ymin, ymax,
		 xlabel, ylabel, hlines=False, vlines=True, show=False)

	y = piv2[close][sclose]
	y2 = piv2[far][sfar]
	outname = '../../Figures/piv_theta_galv_EAGLE.pdf'
	plot(x, y, z, x2, y2, outname, xmin, xmax, ymin, ymax,
		 xlabel, ylabel, hlines=False, vlines=True, show=False)

	if 0:
		xmin = 0
		xmax = 8
		ymin = 0
		ymax = 23

		plot(pd[c][close][sclose])
		y = pi2[close][sclose]
		y2 = pi2[far][sfar]

		plt.ylim(ymin, ymax)
		plt.xlim(xmin, xmax)
		plt.xlabel(r'$r_p\,\rm{[pMpc]}$', fontsize=fsize)
		plt.ylabel(r'$\pi\,\rm{[cMpc]}$', fontsize=fsize)
		plt.scatter(x2, y2, msize, marker='.', lw=0, vmin=vmin, vmax=vmax, c='gray', alpha=.5)
		plt.scatter(x, y, msize, z, '.', lw=0, vmin=vmin, vmax=vmax, cmap=cm, alpha=.8)
		plt.scatter(x, y, msize, z, '.', lw=0, vmin=vmin, vmax=vmax, cmap=cm, alpha=.8)
		plt.hlines(20, xmin, xmax, linestyles='dashed', lw=2)
		plt.hlines(.5, xmin, xmax, linestyles='dashed', lw=2)
		cbar = plt.colorbar(ticks=redshifts)
		cbar.set_label(r'$\rm{redshift}$', fontsize=fsize, labelpad=0.6)
		plt.savefig("../../Figures/pi_theta_EAGLE.pdf", dpi=dpi, format='pdf')
		plt.show()
		plt.close()

		plt.figure(figsize=(figx, figy))
		y = pi2[close][sclose]
		y2 = pi2[far][sfar]
		plt.ylim(ymin, ymax)
		plt.xlim(xmin, xmax)
		plt.xlabel(r'$r_p\,\rm{[pMpc]}$', fontsize=fsize)
		plt.ylabel(r'$\pi\,\rm{[cMpc]}$', fontsize=fsize)
		plt.scatter(x2, y2, msize, marker='.', lw=0, vmin=vmin, vmax=vmax, c='gray', alpha=.5)
		plt.scatter(x, y, msize, z, '.', lw=0, vmin=vmin, vmax=vmax, cmap=cm, alpha=.8)
		plt.hlines(20, xmin, xmax, linestyles='dashed', lw=2)
		plt.hlines(.5, xmin, xmax, linestyles='dashed', lw=2)
		cbar = plt.colorbar(ticks=redshifts)
		cbar.set_label(r'$\rm{redshift}$', fontsize=fsize, labelpad=0.6)
		plt.savefig("../../Figures/pi_theta_galv_EAGLE.pdf", dpi=dpi, format='pdf')
		plt.show()
		plt.close()

	if 0:
		ymin = -50
		ymax = 1500
		xmin = 0
		xmax = 550

		plt.figure(figsize=(figx, figy))
		plt.ylim(ymin, ymax)
		plt.xlim(xmin, xmax)
		plt.xlabel(r'$\theta\,\rm{[arcsec]}$', fontsize=fsize)
		plt.ylabel(r'$\Delta v\,\rm{[km/s]}$', fontsize=fsize)
		y = piv[close][sclose]
		x = theta[close][sclose]
		y2 = piv[far][sfar]
		x2 = theta[far][sfar]
		plt.scatter(x2, y2, msize, marker='.', lw=0, vmin=vmin, vmax=vmax, c='gray', alpha=.5)
		plt.scatter(x, y, msize, z, '.', lw=0, vmin=vmin, vmax=vmax, cmap=cm, alpha=.8)
		plt.vlines(theta0, ymin, ymax, linestyles='dashed', lw=2)
		cbar = plt.colorbar(ticks=redshifts)
		cbar.set_label(r'$\rm{redshift}$', fontsize=fsize, labelpad=0.6)
		plt.savefig("../../Figures/piv_theta_EAGLE.pdf", dpi=dpi, format='pdf')
		plt.show()
		plt.close()

		plt.figure(figsize=(figx, figy))
		plt.ylim(ymin, ymax)
		plt.xlim(xmin, xmax)
		plt.xlabel(r'$\theta\,\rm{[arcsec]}$', fontsize=fsize)
		plt.ylabel(r'$\Delta v\,\rm{[km/s]}$', fontsize=fsize)
		y = piv2[close][sclose]
		x = theta[close][sclose]
		y2 = piv2[far][sfar]
		x2 = theta[far][sfar]
		plt.scatter(x2, y2, msize, marker='.', lw=0, vmin=vmin, vmax=vmax, c='gray', alpha=.5)
		plt.scatter(x, y, msize, z, '.', lw=0, vmin=vmin, vmax=vmax, cmap=cm, alpha=.8)
		plt.vlines(theta0, ymin, ymax, linestyles='dashed', lw=2)
		cbar = plt.colorbar()
		cbar.set_label(r'$\rm{redshift}$', fontsize=fsize, labelpad=0.6)
		plt.savefig("../../Figures/piv_theta_galv_EAGLE.pdf", dpi=dpi, format='pdf')
		plt.show()
		plt.close()

if subcubes:
	print 'subcubesss'

	folders = ['../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100/',
			   '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n8-100/']

	is_sb1 = os.path.isfile('sb1.dat')
	is_sb2 = os.path.isfile('sb2.dat')
	overwrite = True
	luniq = 76634
	npix = 91897

	if is_sb1 and is_sb2 and not overwrite:
		sb1 = np.loadtxt('sb1.dat')
		sb2 = np.loadtxt('sb2.dat')

	else:
		zoff = 1
		zw0 = 2

		sb1 = []
		sb2 = []
		lst = open(folders[0] + "stack.lst", 'r')
		names = [l[:-1] for l in lst]
		ffits = [np.array(getdata(n)) for n in names]
		ffits = np.array(ffits)
		nuns = [np.array(getdata(n.replace('.fits', '.nun.fits'))) for n in names]
		# npixs = [np.array(getdata(n.replace('.fits', '.NPIX.fits'))) for n in names]
		ffits = np.array(ffits)
		nuns = np.array(nuns)
		vars = [getdata(ff + 'stack.prevar.fits')[0] for ff in folders]
		vars = np.array(vars)
		n, zl, yl, xl = ffits.shape
		zf = zl / 2 + zw0 + 1 + zoff
		zi = zl / 2 - zw0 + zoff
		# x1, x2, y1, y2, z1, z2 = [45, 52, 38, 43, zi, zf]
		x1, x2, y1, y2, z1, z2 = [40 + 5, 40 + 14, 38, 43, zi, zf]  # where the cgm emission should be!
		xx1, xx2, yy1, yy2, zz1, zz2 = [40 + 5, 40 + 14, 38, 43, zi, zf]  # where the filament emission should be!
		conv = 7.8125 * (z2 - z1)  # includes arcsec conversion!!!
		conv2 = 7.470  # conversion to physical projected distance
		nsigma = 3
		npix = []
		# for ff, l, nn, npix in zip(ffits, names, nuns, npixs):
		for ff, l, nn in zip(ffits, names, nuns):

			high_sigma = np.abs(ff) > nsigma * vars[0]
			ff[high_sigma] = np.nan
			aaa = conv * np.nanmean(ff[z1: z2, y1: y2, z1: z2])
			sb1.append(aaa)
			if aaa > 24: print l, aaa

			if luniq is None:
				nn[high_sigma] = 0
				ns = nn[z1: z2, y1: y2, x1: x2]
				ns = ns[ns > 0]
				for n in ns:
					sn = str(int(n))
					nnn = [int(sn[i:i + 10]) for i in np.arange(0, len(sn), 10)]
					for n_ in nnn:
						npix.append(n)
					# npix[high_sigma] = np.nan
					# ntot += np.nansum(npix)
					# print 'N unique %d ' % (len(np.unique(nuniq)))

		if luniq is None:
			uniq = np.unique(npix)
			luniq = len(uniq)
			lpix = len(npix)
			print 'Total %d unique %d fraction %f' % (lpix, luniq, luniq / float(lpix))

		sb1 = np.array(sb1)[np.isfinite(sb1)]
		lst = open(folders[1] + "stack.lst", 'r')
		names = [l[:-1] for l in lst]
		ffits = [getdata(n) for n in names]
		ffits = np.array(ffits)
		n, zl, yl, xl = ffits.shape
		for ff, l in zip(ffits, names):
			high_sigma = np.abs(ff) > nsigma * vars[1]
			ff[high_sigma] = np.nan
			aaa = conv * np.nanmean(ff[z1: z2, y1: y2, z1: z2])
			sb2.append(aaa)
			if aaa > 20: print l, aaa
		sb2 = np.array(sb2)[np.isfinite(sb2)]
		np.savetxt('sb1.dat', sb1)
		np.savetxt('sb2.dat', sb2)

	fsize = 32
	plt.figure(figsize=(18, 10))
	bins = np.arange(-46, 48, 4)
	plt.hist(sb1[np.isfinite(sb1)], bins, facecolor='gray', alpha=0.75, label='Full sample')
	plt.hist(sb2[np.isfinite(sb2)], bins, facecolor='purple', alpha=0.75, label='8+ Neighbours')
	plt.legend(prop={'size': fsize - 6})
	plt.xlabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize)
	plt.ylabel(r'$\rm{Counts}$', fontsize=fsize)
	plt.tick_params(labelsize=fsize * .8)
	plt.xlim(-46, 46)
	plt.savefig('../../all/plots/subcubes.pdf', dpi=dpi, format='pdf')
	plt.close()

if halfplot:
	variables = ['Luminosity', 'Environment', 'Redshift', "Com_Distance", "Proj_Distance"]

	mflux = 91.0415050986  # 30 # now is luminosity! 29.53166
	mred = 3.47385  # 3.47268
	mtheta = 32.586  # 32.0235
	mdists = 8.1767  # 7.78095
	mnpairs = 8  # 6
	extraname = ''
	medians = [91, 8, 3.5, 8, 32]
	mins = [1, 1, 2.9, 0, 16]
	maxs = [2000, 100, 4.0, 20, 80]
	# max flux around 1700 (check units...1e-20?)
	ff = [[[]] * n] * 2
	fs = []
	frs = []
	fout = open('%s/%s/cats/halfs.dat' % (folder, fitcat), 'w')
	fout.write('#id name half median SB_stack SB_mean SB_err\n')
	nmin = np.zeros((2, n))
	nmax = np.zeros((2, n))
	mu = np.zeros((2, n))
	mu2 = np.zeros((2, n))
	sigma = np.zeros((2, n))
	# zoff = 1
	# zw0 = 3
	zw = 9
	zw0 = zw / 2
	zoff = 0

	x1, x2, y1, y2, z1, z2 = [45, 52, 38, 43, 10 - zw0 + zoff, 10 + zw0 + 1 + zoff]

	conv = 7.125 * (z2 - z1)  # /((x2-x1)*(y2-y1))
	usedat = False

	for i in np.arange(n):
		for j in range(n):
			if j == i:
				nmin[0][j] = medians[j]
				nmax[0][j] = maxs[j]
				nmin[1][j] = mins[j]
				nmax[1][j] = medians[j]
			# print nmin[0][j], medians[j], nmax[0][j], maxs[j], nmin[1][j], mins[j], nmax[0][j], medians[j], j, i
			else:
				nmin[1][j] = mins[j]
				nmax[1][j] = maxs[j]
				nmin[0][j] = mins[j]
				nmax[0][j] = maxs[j]
			# print nmin[0][j], nmax[0][j], nmin[1][j], nmax[0][j], j, i

		s = ['min', 'max']
		s2 = ['>', '<']
		s3 = ['_ht_', '_lt_']

		for k in range(2):
			if usedat:
				ss = '%s/%s/stacks/d%d-%d_pd%d-%d_z%.1f-%.1f_l%d-%d_n%d-%d%s/%s' % (
					folder, fitcat, nmin[k][3], nmax[k][3], nmin[k][4], nmax[k][4], nmin[k][2], nmax[k][2],
					nmin[k][0], nmax[k][0], nmin[k][1], nmax[k][1], extraname, name)

				print ss
				ff = np.loadtxt(ss).transpose()

			else:
				ss = '%s/%s/stacks/d%d-%d_pd%d-%d_z%.1f-%.1f_l%d-%d_n%d-%d%s/' % (
					folder, fitcat, nmin[k][3], nmax[k][3], nmin[k][4], nmax[k][4], nmin[k][2], nmax[k][2],
					nmin[k][0], nmax[k][0], nmin[k][1], nmax[k][1], extraname)
				print ss
				# hdu = fits.open(ss + 'stack.fits')
				# mu[k][i] = conv*np.nanmean((hdu[0].data)[z1:z2, y1:y2, x1:x2])
				frs = []
				fts = []
				for j in range(200):
					# print i, 'of 200'
					hdu = fits.open(ss + 'randoms/random_stack.%d.fits' % j)
					frs.append(conv * np.nanmean((hdu[0].data)[z1:z2, y1:y2, x1:x2]))
					hdu = fits.open(ss + 'tests/stack.t%d.fits' % j)
					fts.append(conv * np.nanmean((hdu[0].data)[z1:z2, y1:y2, x1:x2]))
				fmix = []
				for f1 in frs:
					for f2 in fts:
						fmix.append(f2 - f1)
				del frs, fts
				fmix = np.array(fmix)
				mu[k][i] = np.nanmean(fmix)
				sigma[k][i] = np.nanstd(fmix)
				del fmix

			if usedat:
				fmix = []
				for ii in ff[2]:
					for jj in ff[1]:
						fmix.append(ii - jj)
				sigma[k][i] = np.nanstd(fmix)
				mu[k][i] = np.nanmean(fmix)
				mu2[k][i] = ff[0][0]
				ss2 = '%s/%s/plots/%s/%s' % (
				folder, fitcat, extraname[1:], '%s%s%d.png' % (variables[i], s3[k], medians[i]))
				histplot(ff[2], ff[1], fmix, ss2, r'$%s%s%.1f$' % (variables[i], s2[k], medians[i]))
				sss = '%d %s %s %.1f %f %f %f\n' % (
					i, variables[i], s2[k], medians[i], mu2[k][i], mu[k][i], sigma[k][i])
				print sss
				fout.write(sss)

	ss = '%s/%s/stacks/d%d-%d_pd%d-%d_z%.1f-%.1f_l%d-%d_n%d-%d%s/%s' % (
		folder, fitcat, mins[3], maxs[3], mins[4], maxs[4], mins[2], maxs[2],
		mins[0], maxs[0], mins[1], maxs[1], extraname, name)

	ss2 = '%s/%s/stacks/d%d-%d_pd%d-%d_z%.1f-%.1f_l%d-%d_n%d-%d%s/' % (
		folder, fitcat, mins[3], maxs[3], mins[4], maxs[4], mins[2], maxs[2],
		mins[0], maxs[0], mins[1], maxs[1], extraname)

	if usedat:
		print ss
		ff = np.loadtxt(ss).transpose()
		for ii in ff[2]:
			for jj in ff[1]:
				fmix.append(ii - jj)
		fsigma = np.nanstd(fmix)
		fmu = np.nanmean(fmix)
		fmu2 = ff[0][0]
		ss2 = '%s/%s/plots/%s/%s' % (folder, fitcat, extraname[1:], name.replace('.dat', '.png'))
		histplot(ff[2], ff[1], fmix, ss2)
		sss = '%d full - -1 %f %f %f\n' % (n, fmu2, fmu, fsigma)
		print sss
		fout.write(sss)
		fout.close()

	else:
		print ss2
		hdu = fits.open(ss2 + 'stack.fits')
		# fmu = conv*np.nanmean((hdu[0].data)[z1:z2, y1:y2, x1:x2])
		frs = []
		fts = []
		for j in range(200):
			# print i, 'of 200'
			hdu = fits.open(ss2 + 'randoms/random_stack.%d.fits' % j)
			frs.append(conv * np.nanmean((hdu[0].data)[z1:z2, y1:y2, x1:x2]))
			hdu = fits.open(ss2 + 'tests/stack.t%d.fits' % j)
			fts.append(conv * np.nanmean((hdu[0].data)[z1:z2, y1:y2, x1:x2]))
		fmix = []
		for f1 in frs:
			for f2 in fts:
				fmix.append(f2 - f1)
		del frs, fts
		fmix = np.array(fmix)
		fmu = np.nanmean(fmix)
		fsigma = np.nanstd(fmix)
		del fmix

	plt.figure(figsize=(7, 4))
	plt.errorbar([fmu], [n + .5], xerr=[fsigma], fmt='o', color='purple', capsize=10)
	plt.scatter([fmu], [n + .5], color='purple')
	labels = [r'$\rm{Luminosity}$', r'$\rm{Number\,\,of\,\,Neighbours}$', r'$\rm{Redshift}$',
			  r"$\rm{Comoving\,\,Distance}$",
			  r"$\rm{Projected\,\,Distance}$", r"$\rm{Full}$"]
	y = np.arange(n + 1) + .5
	plt.ylim(-.2, n + 1)
	plt.xlim(-0.8, 3.5)
	xticks = np.arange(0, 3.5, 1)
	plt.yticks(y, labels, rotation=0, fontsize=16)
	plt.xticks(xticks, xticks)
	y = np.arange(n) + .5
	x1 = mu[0]
	x2 = mu[1]
	xerr1 = sigma[0]
	xerr2 = sigma[1]
	plt.errorbar(x1, y, xerr=xerr1, fmt='o', color='red', capsize=10, label=r'$\rm{Higher}$')
	plt.scatter(x1, y, color='red')
	plt.errorbar(x2, y, xerr=xerr2, fmt='o', color='blue', capsize=10, label=r'$\rm{Lower}$')
	plt.scatter(x2, y, color='blue')
	plt.xlabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=16)
	plt.tick_params(labelright=True, labelleft=False)
	# leg = plt.legend(loc='upper left', numpoints=1)
	# text = leg.get_texts()
	# plt.setp(text[0], color = 'red')
	# plt.setp(text[1], color = 'blue')
	# plt.text(0, 0, r'$\rm{Half\,higher}$', fontsize=20)#, fontcolor='red'
	# plt.text(0, 2, r'$\rm{Half\,lower}$', fontsize=20)
	plt.tight_layout()
	plt.savefig('%s/%s/plots/%s/halfs.pdf' % (folder, fitcat, extraname[1:]), dpi=dpi, format='pdf')
	plt.close()

if envplot:
	n = 8
	x = np.arange(n) + 2
	randmean = np.zeros(n)
	sigma = np.zeros(n)
	mu = np.zeros(n)
	mu2 = np.zeros(n)

	for i in range(n):
		ss = '%s/%s/stacks/not-normalized_masked_cubical_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d-%d_nmin%d_nmax%d/%s' % (
			folder, fitcat, mins[3], maxs[3], mins[4], maxs[4], mins[2], maxs[2],
			mins[0], maxs[0], i + 1, i + 2, name)
		ff = np.loadtxt(ss).transpose()
		fmix = []
		for ii in ff[2]:
			for jj in ff[1]:
				fmix.append(ii - jj)
		sigma[i] = np.nanstd(fmix)
		mu[i] = np.nanmean(fmix)
		mu2[i] = ff[0][0]
		ss2 = '%s/%s/plots/%/%s' % (folder, fitcat, extraname, '%s%d-%d.png' % (variables[1], i + 1, i + 2))
		histplot(ff[2], ff[1], fmix, ss2)  # , r'$%d \leq \rm{N_{neighbours}} \leq%d$' % (i + 1, i + 2))
	y = mu - randmean
	yerr = sigma
	plt.errorbar(x, y, yerr=yerr, xerr=np.zeros(len(x)) + .5, fmt='o', color='red')
	plt.scatter(x, y, color='red')
	plt.xticks(np.arange(n + 2) + 1, rotation=0)
	plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=16)
	plt.xlabel('# of neighbours', fontsize=16)
	plt.legend()
	plt.savefig('%s/%s/plots/%s/environment.png' % (folder, fitcat, extraname))
	plt.close()
	# plt.show()

if sbplot:
	print 'SB prof'

	# run this before starting!!
	# ulimit -S -n 10000

	fitcat = 'all'
	random = True
	nrand = 200
	tests = False

	if 0:
		prename = '.unique'  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
		# folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % prename
		folder = '../../%s/stacks/d0-20_pd16-80_z2.9-4.0_l0-2000_n20-100%s/' % (fitcat, prename)
		extraname = '_unique'  # '_nmin8'

		##define collection of region sizes
		xw = 1
		yw = 0
		zw = 9
		zoff = 1  # 1
		width = 4
		sigmas = 2
		x1 = 2
		x2 = 13.6
		y1 = -1
		y2 = 10
		hmSB = 1.14
		y1b = -6
		y2b = 10

	if 0:
		fitcat = 'mosaic'
		random = True
		nrand = 19
		tests = True
		prename = ''  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
		# folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % prename
		folder = '../../%s/stacks/d0-20_los0-20_pd16-300_z2.9-4.0_l0-2000_n1-100_v0-10000%s/' % (fitcat, prename)
		extraname = '_mosaic'  # '_unique'#'_nmin8'

		##define collection of region sizes
		xw = 1
		yw = 2
		zw = 9
		zoff = 0  # 1
		width = 4
		sigmas = 2
		x1 = 1.7
		x2 = 13.6
		y1 = .1
		y2 = 100
		hmSB = 1.14
		y1b = -1
		y2b = 6

	if 0:
		fitcat = 'mosaic'
		random = True
		nrand = 20
		tests = True
		prename = '.UDF10'  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
		# folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % prename
		folder = '../../%s/stacks/d0-20_los0-20_pd16-300_z2.9-4.0_l0-2000_n1-100_v0-10000%s/' % (fitcat, prename)
		extraname = '_UDF10'  # '_unique'#'_nmin8'

		##define collection of region sizes
		xw = 1
		yw = 2
		zw = 9
		zoff = 0  # 1
		width = 4
		sigmas = 2
		x1 = 1.7
		x2 = 13.6
		y1 = .1
		y2 = 100
		hmSB = 1.14
		y1b = -1
		y2b = 6
		ext = 'png'

	if 0:
		fitcat = 'UDF'
		random = True
		nrand = 20
		tests = True
		prename = '.mosaic'  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
		# folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % prename
		folder = '../../%s/stacks/d0-20_los0-20_pd16-300_z2.9-4.0_l0-2000_n1-100_v0-10000%s/' % (fitcat, prename)
		extraname = '_mosaic2'  # '_unique'#'_nmin8'

		##define collection of region sizes
		xw = 1
		yw = 2
		zw = 9
		zoff = 0  # 1
		width = 4
		sigmas = 2
		x1 = 1.7
		x2 = 13.6
		y1 = .1
		y2 = 100
		hmSB = 1.14
		y1b = -1
		y2b = 6
		ext = 'png'

	if 0:
		fitcat = 'UDF'
		random = True
		nrand = 50
		tests = True
		prename = ''  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
		# folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % prename
		folder = '../../%s/stacks/d0-20_los0-20_pd16-300_z2.9-4.0_l0-2000_n1-100_v0-10000%s/' % (fitcat, prename)
		extraname = '_mosaic'  # '_unique'#'_nmin8'

		##define collection of region sizes
		xw = 1
		yw = 2
		zw = 9
		zoff = 0  # 1
		width = 4
		sigmas = 2
		x1 = 1.7
		x2 = 13.6
		y1 = .1
		y2 = 100
		hmSB = 1.14
		y1b = -1
		y2b = 6
		ext = 'png'

	if 1:
		prename = ''  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
		folder = '../../%s/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % (fitcat, prename)
		extraname = '_nmin8'
		ext = 'pdf'

		##define collection of region sizes
		xw = 1
		yw = 2
		zw = 5
		zoff = 1
		width = 4
		sigmas = 2

		x1 = 1.7
		x2 = 13.6
		y1 = -1.5
		y2 = 10
		hmSB = 1.14
		y1b = -4
		y2b = 6

	if 0:
		prename = ''  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
		folder = '../../%s/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % (fitcat, prename)
		extraname = '_full'
		ext = 'pdf'

		##define collection of region sizes
		xw = 1
		yw = 2
		zw = 5
		zoff = 1
		width = 4
		sigmas = 2

		x1 = 1.7
		x2 = 13.6
		y1 = -1.5
		y2 = 10
		hmSB = 1.14
		y1b = -4
		y2b = 6


	def sbplot(x, y, yr=None, yrstd=None, name='SB.pdf', ystd=None, ycombstd=None, fsize=40,
			   label="Oriented stack towards neighbours", hm=True, randlabel="Random stacks", dir=None):

		plt.figure(figsize=(20, 12))
		plt.semilogx()
		# plt.loglog()
		ax = plt.axes()
		cax = plt.gca()
		plt.xlim([x1, x2])
		xrange = np.arange(2, 11, 2)
		plt.xticks(xrange, xrange, fontsize=fsize)
		plt.xlabel(r"$\theta$ [arcsec]", fontsize=fsize)
		plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize + 2)
		yticks = [0, 2, 4, 6, 8, 10]
		plt.yticks(yticks, yticks, fontsize=fsize)
		plt.minorticks_off()
		plt.twiny(ax=None)
		plt.semilogx()
		# plt.loglog()
		plt.xlim([x1, x2])
		plt.xlabel(r'$r_p\,\rm{[pkpc]}$', fontsize=fsize)
		kpc = np.arange(15, 90, 15) / 7.47
		plt.xticks(kpc, (kpc * 7.47).astype('int'), fontsize=fsize)
		plt.ylim([y1, y2])
		plt.plot((x1, x2), (0, 0), '--', label='', color="gray", linewidth=width + 1)
		if hm: plt.plot((x1, x2), (hmSB, hmSB), '--', label=r"LLS Fluorescence from HM12", color="green",
						linewidth=width + 1)
		print len(dir), 'len dirirrrr'
		if dir is not None: plt.text(x2 * .4, y2 * .5, r'$\rm{%s \%sarrow}$'
									 % (dir.replace('top', 'top\,\,direction').replace('bottom', 'bottom\,\,direction')
										.replace('right', 'towards\,\,neighbours\,')
										.replace('left', 'away\,\,from\,\,neighbours\,'),
										dir.replace('top', 'up').replace('bottom', 'down')), fontsize=fsize)
		plt.plot(x, y, label=label, lw=width, color='black')
		plt.scatter(x, y, s=width * 10, color='black')
		if ystd is not None:
			_yd = y - sigmas * ystd
			# for iii in range(len(_yd)):
			# if _yd[iii] < 0: _yd[iii] = 1.e-10
			plt.fill_between(x, _yd, y + sigmas * ystd, facecolor='gray', alpha=0.3, lw=0, edgecolor='none')

		if yr is not None:
			_yd = yr - sigmas * yrstd
			# for iii in range(len(_yd)):
			# if _yd[iii] < 0: _yd[iii] = 1.e-10
			plt.plot(x, yr, label=randlabel, c='dodgerblue', lw=width)
			plt.scatter(x, yr, c='dodgerblue', s=width * 10)
			plt.fill_between(x, _yd, yr + sigmas * yrstd, facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')
		plt.minorticks_off()
		if hm: plt.legend(fontsize=fsize, loc='best')  # (3.5,2))
		# plt.text(x[5]*.95, yr[5]-sigmas*yrstd[5]*1.5, r'$%d\,\sigma$' % sigmas, size=fsize, color='dodgerblue')
		# mdir = '../../all/plots/%s/' % extraname[1:]
		mdir = '../../%s/plots/%s/' % (fitcat, extraname[1:])
		if not os.path.isdir(mdir):
			os.system('mkdir %s' % mdir)
		plt.savefig(mdir + name.replace('.pdf', '.' + ext), dpi=dpi, format=ext)
		# plt.show()
		plt.close()

		if ycombstd is not None:

			plt.figure(figsize=(30 - xw, 15))
			plt.semilogx()
			ax = plt.axes()
			cax = plt.gca()
			plt.xlim([x1, x2])
			xrange = np.arange(2, 13, 2)
			plt.xticks(xrange, xrange, fontsize=fsize)
			plt.xlabel(r"$\theta$ [arcsec]", fontsize=fsize)
			plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize + 2)
			plt.yticks(fontsize=fsize)
			plt.minorticks_off()
			plt.twiny(ax=None)
			plt.semilogx()
			plt.xlim([x1, x2])
			plt.xlabel(r'$\theta$ [kpc]', fontsize=fsize)
			kpc = np.arange(15, 92, 15) / 7.47
			plt.xticks(kpc, (kpc * 7.47).astype('int'), fontsize=fsize)
			plt.ylim([y1b, y2b])
			plt.plot((x1, x2), (0, 0), '--', label='', color="gray", linewidth=width + 1)
			if hm: plt.plot((x1, x2), (hmSB, hmSB), '--', label=r"LLS Fluorescence from HM12", color="green",
							linewidth=width + 3)
			ydiff = y - yr
			plt.plot(x, ydiff, label="Oriented - Random", lw=width, color='black')
			plt.scatter(x, ydiff, s=width * 10, color='black')
			if random:
				plt.fill_between(x, ydiff - sigmas * ycombstd, ydiff + sigmas * ycombstd, facecolor='gray', alpha=0.3,
								 lw=0, edgecolor='none')
			plt.minorticks_off()
			plt.legend(fontsize=fsize, loc='best')  # (3.5,2))
			mdir = '../../all/plots/%s/' % extraname[1:]
			if not os.path.isdir(mdir):
				os.system('mkdir %s' % mdir)
			plt.savefig(mdir + name.replace('.pdf', '.2.pdf'), dpi=dpi, format='pdf')
			# plt.show()
			plt.close()


	plt.loglog()
	plt.xticks([1, 2, 3, 4], ['1', '2', '3', '4'])
	plt.twiny()
	plt.xlabel(r'$\theta\,\rm{[kpc]}$', fontsize=16)
	plt.loglog()
	plt.xlim(0.5, 5)

	# folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100.nosclip/'
	# "../../all/stacks/not-normalized_masked_cubical_d0-20_pd16-80_z2.9-4.0_f0-99999.nosclip/"
	# folder='not-normalized_masked_cubical_d0-20_pd16-80_z2.9-4.0_f0_nmin6.half/'
	#    folder = '../../all/stacks/d0-20_pd10-80_z3.5-4.0_l1-2000_n9-100%s/' % extraname

	fname = folder + "stack.fits"
	frname = folder + "random_stack.fits"
	israndom = os.path.isfile(frname)
	hdu = fits.open(fname)
	f = hdu[0].data

	binsize = 2
	# conv0 = 1.25 / pow(binsize * .2, 2)

	if random:
		if israndom:
			hdur = fits.open(frname)
			fr = hdur[0].data
		frs = []
		ft = []
		for i in range(nrand):
			hdur = fits.open(folder + 'randoms/random_stack.%d.fits' % i)
			frs.append(hdur[0].data)  # * conv0)
			if tests:
				hdut = fits.open(folder + 'tests/stack.t%d.fits' % i)
				ft.append(hdut[0].data)  # * conv0)
		frs = np.array(frs)
		if tests: ft = np.array(ft)

		if not israndom:
			sclip = 2
			stds = np.nanstd(frs, 0)
			high_sigma = np.abs(frs) > sclip * stds
			frs[high_sigma] = np.nan
			fr = np.nanmean(frs, 0)

	## similar analysis that in tools but more flexible
	zl, yl, xl = f.shape
	xmin, xmax, ymin, ymax, zmin, zmax = [xl / 2, xl - 10, 0, yl - yw, 0, zl - zw]
	s = []
	sr = []
	plt.figure(figsize=(20, 15))
	fsize = 30

	## sb profile
	t = np.concatenate([np.arange(xmin, xmin + 12, 1), [53, 58, 67, 81]]).astype('int')
	# np.concatenate([np.arange(xmin, xmin+10, 1), np.logspace(0, 1.6, num=6, endpoint=True, base=10)[:-2]+xmin+9,
	# [70]]).astype('int')
	dir = ['right', 'top', 'left', 'bottom']
	hm = [True, False, False, False]
	randlabel = ["Random stacks", '', '', '']
	zi = zl / 2 - zw / 2 + zoff
	zf = zl / 2 + zw / 2 + zoff + 1
	sb = []
	sb_std = []
	sbr = []
	sbr_std = []
	sbrcomb_std = []
	yymin = yl / 2 - yw / 2
	yymax = yl / 2 + yw / 2 + 1
	for i in range(len(t) - 1):
		# enterder bien estos limites dependiendo de la direccion!!!!!

		yis = [yymin, t[i] - 1, yymin, xl - t[i + 1]]
		yfs = [yymax, t[i + 1], yymax, xl - t[i] + 1]
		xis = [t[i] - 1, yymin, xl - t[i + 1], yymin]
		xfs = [t[i + 1], yymax, xl - t[i] + 1, yymax]
		for yi, yf, xi, xf, ds in zip(yis, yfs, xis, xfs, dir):
			# 7.8125 comes from 1.25 Angstrom width times 1/(.4)**2 from convertion to the correct aperture
			c = zw * 7.8125  # / (yf-yi) / (xf-xi) / (zf-zi)
			n1, n2, n3 = f[zi: zf, yi: yf, xi:xf].shape
			# print ds, xi, xf, yi, yf, c
			# print (yf-yi) *(xf-xi), n1*n2*n3
			sb.append(np.nanmean(f[zi:zf, yi:yf, xi:xf]) * c)
			if tests: sb_std.append(np.nanstd([np.nanmean(ff[zi:zf, yi:yf, xi:xf]) * c for ff in ft]))
			if random:
				sbr.append(np.nanmean(fr[zi:zf, yi:yf, xi:xf]) * c)
				sbr_std.append(np.nanstd([np.nanmean(ff[zi:zf, yi:yf, xi:xf]) * c for ff in frs]))
			if tests: sbrcomb_std.append(np.nanstd(
				[np.nanmean(f1[zi:zf, yi:yf, xi:xf] - f2[zi:zf, yi:yf, xi:xf]) * c for f1, f2 in zip(ft, frs)]))

	xi, xf = [15, 30]  # from 4 to 12 arcsec
	if random:
		SBcool = c * (np.nanmean(f[zi:zf, yi:yf, xi:xf]) - np.nanmean(frs[:, zi:zf, yi:yf, xi:xf]))
		SBstdcool = np.nanstd(
			[np.nanmean(f1[zi:zf, yi:yf, xi:xf] - f2[zi:zf, yi:yf, xi:xf]) * c for f1, f2 in zip(ft, frs)])
		print 'SB %f 2 sigma SB %f' % (SBcool, 2 * SBstdcool)

	lx = len(t) - 1
	ld = len(dir)
	# labels = ["Oriented stack %s direction" % s for s in dir]
	# labels = ["Oriented towards galaxy neighbours"]*4
	labels = ["Oriented stack"] * 4
	sb = np.array(sb).reshape((lx, ld))
	if tests: sb_std = np.array(sb_std).reshape((lx, ld))
	if random:
		sbr = np.array(sbr).reshape((lx, ld))
		sbr_std = np.array(sbr_std).reshape((lx, ld))
	if tests: sbrcomb_std = np.array(sbrcomb_std).reshape((lx, ld))
	x = (t[:-1] + (t[1:] - t[:-1]) / 2. - xl / 2) * .4
	for i in range(ld):
		if random:
			_sbr = sbr[:, i]
			_sbr_std = sbr_std[:, i]
		else:
			_sbr = None
			_sbr_std = None

		if tests:
			_ystd = None  # sb_std[:, i]
			_ycombstd = sbrcomb_std[:, i]
		else:
			_ycombstd = None
			_ystd = None
		sbplot(x, sb[:, i], _sbr, np.abs(_sbr_std), 'SB_%s.pdf' % (dir[i] + extraname), ystd=_ystd,
			   ycombstd=_ycombstd, label=labels[i], hm=hm[i], randlabel=randlabel[i], dir=dir[i])

	dospectra = False
	if dospectra:
		## Collapsed spectra
		cs = []
		csr = []
		zrange = np.arange(0, zl)
		conv1 = 88.13482042236556  # 300000*1.25/3.5/1215.67
		for i in zrange:
			cs.append(np.nanmean(f[i: i + 1, yymin: yymax, xmin: xmax]))
			csr.append(np.nanmean(fr[i: i + 1, yymin: yymax, xmin: xmax]))

		plt.plot((zrange - zl / 2) * conv1, cs, label="Oriented")
		plt.plot((zrange - zl / 2) * conv1, csr, label="Random")
		plt.ylabel("Mean Flux [1e-20 erg/s/cm^2]")
		plt.xlabel("velocity [km/s]")
		plt.legend()
		# plt.show()
		plt.savefig('../../all/plots/%s/spectra.png' % extraname[1:])
		plt.close()

if lutz:
	print 'Lutz'

	folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100.reject.unique/'
	sbdat = folder + "SBprofile.dat"
	overwrite = True


	def radprof(fits, rads, zw0=2, zmean=3.5, conv=7.8125, sides=False, std=False):
		zl, yl, xl = fits.shape

		y, x = np.ogrid[-yl / 2 + 1: yl / 2 + 1, -xl / 2 + 1: xl / 2 + 1]
		zcen = range(zl / 2 - zw0, zl / 2 + zw0 + 1)
		zblue = range(zl / 2 - 5)
		zred = range(zl / 2 - 5, zl)
		fcen = fits[zcen, :]
		fblue = fits[:zl / 2 - 5, :]
		fred = fits[zl / 2 + 6:, :]
		SB = []
		if std: std0 = []

		if sides:
			SBred = []
			SBblue = []
			stdred = []
			stdblue = []

		# rads = np.concatenate([rads, [xl]])
		nrads = len(rads)
		for i in range(nrads - 1):
			# conv1 = conv# / (np.pi * (rads[i+1] ** 2 - rads[i] ** 2))
			cool = (x ** 2 + y ** 2 >= rads[i] ** 2) & (x ** 2 + y ** 2 < rads[i + 1] ** 2)
			fc = fcen[:, cool]
			if i > nrads / 3:
				sclip = 3
				stds = np.nanstd(fc, 0)
				high_sigma = np.abs(fc) > sclip * stds
				fc[high_sigma] = np.nan

			SB.append(conv * np.nanmean(fc))
			if std: std0.append(conv * np.nanstd(fc))
			# print conv1, sum, rads[i], rads[i+1], sum*conv1

			if sides:
				# print rads[i] + 1, sum, sum * conv
				SBred.append(conv * np.nanmean(fred[:, cool]))
				stdred.append(conv * np.nanstd(fred[:, zl / 2 + 6:, cool]))
				SBblue.append(conv * np.nanmean(fblue[:, cool]))
				stdblue.append(conv * np.nanstd(fblue[:, :zl / 2 - 5, cool]))

		out = []
		out.append(SB)
		if std:
			out.append(std0)
		if sides:
			out.append(SBred)
			out.append(stdred)
			out.append(SBblue)
			out.append(stdblue)

		return np.array(out)


	if 0:  # os.path.isfile(sbdat) and not overwrite:
		print sbdat, "already exists"
		a = np.loadtxt(sbdat).T
		t = a[0]  # in arcsec
		t2 = a[1]  # in kpc
		SB = a[2]
		SBmin = a[3]
		SBmax = a[4]
		SBstd = a[5]

	else:
		# lst = open(folder + "stack.lst", 'r')
		a = open('../../HDFS/LAEs/laes.c.dat', 'r')
		lst = []
		for i in a:
			lst.append(i[:-1])
		a = open('../../UDF/LAEs/laes.c.dat', 'r')
		for i in a:
			lst.append(i[:-1])

		ffits = [getdata(l) for l in lst]
		ffits = np.array(ffits)

		print '%d individual LAEs!!!' % len(ffits)

		simplestack = True
		if simplestack:
			fits, stds = sclipping(ffits, nsigma=3, dim=0)
			f = np.nanmean(fits, 0)
			zw0 = 2
			offset = 1
			hdu = PrimaryHDU()
			name = 'lae_stack.fits'
			zl, xl, yl = f.shape
			print 'length: x %d y %d z %d ' % (xl, yl, zl)
			zw = '%d:%d' % (zl / 2 - zw0 + 1 + offset, zl / 2 + zw0 + 1 + offset)
			hdu.data = f
			print 'cubim', name
			hdu.writeto(name, clobber=True)
			hdu.data = np.nansum(f[zl / 2 - zw0 + offset:zl / 2 + zw0 + 1 + offset, :, :], 0)
			hdu.writeto(name.replace('.fits', '.IM.fits'), clobber=True)
			s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (name, zw, name.replace('.fits', '.IM.fits'))
			print s
			os.system(s)
			smooth = True
			makeim = True
			contours = True
			highq = True
			label = ''
			title = ''
			stdd = 1
			arrow = False
			sb = True
			vmi = -.05
			vma = .5
			xmin, xmax, ymin, ymax = [-xl / 2, xl / 2, -yl / 2, yl / 2]
			astroim(name.replace('.fits', '.IM.fits'), smooth=smooth, saveim=makeim, show=False, cbfrac=.08, pad=.006,
					dfig=(8, 10), contours=contours, scb_label=label, title=title, vmin=vmi, vmax=vma, std=stdd,
					pair=False,
					highq=highq, sb=sb, y0=0, x0=0, gray=True, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, arrow=arrow)

		whitestack = False
		if whitestack:
			fws = [getdata(l.replace('.c.', '_white.')) for l in lst]
			fws = np.array([np.nanmedian(fw, 0) for fw in fws])
			fws, stdsw = sclipping(fws, nsigma=3)
			fw = np.nanmedian(fws, 0)
			name = 'lae_stack_white.fits'
			hdu = PrimaryHDU()
			hdu.data = fw
			hdu.writeto(name, clobber=True)
			smooth = True
			makeim = True
			contours = True
			highq = True
			label = ''
			title = ''
			stdd = 1
			arrow = False
			sb = True
			vmi = -.05
			vma = .5
			xmin, xmax, ymin, ymax = [-xl / 2, xl / 2, -yl / 2, yl / 2]
			astroim(name.replace('_white.fits', '.fits'), smooth=smooth, saveim=makeim, show=False, cbfrac=.08,
					pad=.006,
					dfig=(8, 10), contours=contours, scb_label=label, title=title, vmin=vmi, vmax=vma, std=stdd,
					pair=False,
					highq=highq, sb=sb, y0=0, x0=0, gray=True, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, arrow=arrow)

		zw0 = 2
		binsize = 1
		n, zl, yl, xl = ffits.shape
		zf = zl / 2 + zw0 + 1
		zi = zl / 2 - zw0
		conv = 1.25 * (2 * zw0 + 1) / (0.2 * binsize) ** 2.  # includes arcsec and wavelength conversion!!!!
		conv2 = 7.470  # conversion to physical projected distance (distance must be in arcsec)
		basex = 10
		rads = np.array([0, 2, 4, 8, 16, 32, 64, 100])  # in pixels
		# rads = np.arange(0, 70, 1)
		# rads = np.array([0, .5, 1, 2, 3, 4, 6, 8, 16, 32, 64]) #in pixels

		t = (rads + np.roll(rads, -1))[:-1] * .2 * binsize  # rads in arcsec
		if 0:
			SBs = np.array([radprof(fit, rads, zw0, conv=conv) for fit in ffits])
			SBgood = []
			for ss in SBs:
				if np.sum(ss < 0) == 0:
					SBgood.append(ss)

			pos = 1
			SBgoodmin = min(np.array(SBgood)[:, pos])

			SB2 = []
			for ss in SBgood:
				if ss[pos] > SBgoodmin: SB2.append(ss)

			SB2 = SBs
			SBstd = np.nanstd(SB2, 0)
			good = (SB2) < 3 * SBstd
			SB = np.array([np.nanmean(sb[g]) for sb, g in zip(SB2.T, good.T)])
			SBmax = np.array([np.max(sb[g]) for sb, g in zip(SB2.T, good.T)])
			SBmin = np.array([np.min(sb[g]) for sb, g in zip(SB2.T, good.T)])

			np.savetxt(sbdat, np.array([t, t * conv2, SB, SBmin, SBmax, SBstd]).T,
					   header='theta[arcsec] theta[kpc] SBmean SBmin SBmax SBstd')
		else:
			fit = getdata('lae_stack.fits')
			SB, std0 = radprof(fit, rads, zw0, conv=conv, std=True)
			fitw = getdata('lae_stack_white.fits')
			SBw, stdw0 = radprof(fitw, rads, zw0, conv=conv, std=True)
		# minSB = min(np.where(SB < 0)[0])

	c3 = 1  # for plotting purposes
	fsize = 14
	fig, ax = plt.subplots()
	# [plt.plot(t, Ss*c3, color='gray', lw=2, alpha=.3) for Ss in SBgood]
	# plt.plot([0], [0], color='gray', lw=2, alpha=.3, label='LAE sample (this study)')
	plt.plot(t, SB * c3, label=r'Average profile', color='red', lw=2)
	plt.scatter(t, SB * c3, color='red', lw=2)
	plt.fill_between(t, c3 * (SB - std0), c3 * (SB + std0), facecolor='red', alpha=0.3, lw=0, edgecolor='none')
	plt.plot(t, SBw * c3, label=r'Average profile white', color='purple', lw=2)
	plt.scatter(t, SBw * c3, color='purple', lw=2)
	plt.fill_between(t, c3 * (SBw - stdw0), c3 * (SBw + stdw0), facecolor='purple', alpha=0.3, lw=0, edgecolor='none')
	xticks = [2, 4, 8, 16]
	plt.xlabel(r'$\theta$ [arcsec]', fontsize=fsize)
	exponent = 20 - np.log10(c3)
	plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-%d}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$' % exponent, fontsize=fsize + 2)
	# plt.loglog()

	ax.set_xscale('log', basex=basex)
	# ax.set_yscale('log', basey=10)
	plt.ylim(-5 * c3, 25 * c3)
	plt.xlim(2, 10)
	plt.xticks(xticks, xticks)
	z_wish = 3.5
	plt.minorticks_off()

	zLutz, rhLutz, r19Lutz = np.loadtxt('lutz.dat', usecols=(1, 3, 6), unpack=1)
	sbCoeffLutz = ((1. + zLutz) / (1. + z_wish) / 4.) ** 4  # to scale / 4.
	cMpcCoeffLutz = (1. + zLutz) / (1. + z_wish)

	####plot
	from math import exp


	def lyaW2016(r, rh, r19):
		return 10 * c3 * exp(r19 / rh) * np.exp(-r / rh)  # if c3=1 -> units: 1e-20 cgs


	imin = 20
	imax = 0
	imax2 = 3
	lineW = 3
	sbLutzmean = np.zeros((len(rhLutz), len(np.arange(0., 170., 1))))
	sbLutzmeanSB = np.zeros((len(rhLutz), len(np.arange(0., 170., 1))))

	for i in range(len(rhLutz)):
		b = np.arange(0., 170., 1)
		sb = lyaW2016(b, rhLutz[i], r19Lutz[i])
		sbLutzmean[i] = sb
		sbLutzmeanSB[i] = sb * sbCoeffLutz[i]
	# plt.plot(b/conv2, sbLutzmean.mean(axis=0), c='dodgerblue', lw=lineW, label='Wisotzki+16')
	# plt.fill_between(b/conv2, sbLutzmean[imin], sbLutzmean[[imax, imax2]].max(axis=0), facecolor='dodgerblue', alpha=0.3, lw=0,
	#                 edgecolor='none')

	plt.legend()
	# plt.plot(b, sbLutzmeanSB.mean(axis=0), c='dodgerblue', lw=lineW)
	# plt.fill_between(b,sbLutzmeanSB[imin], sbLutzmeanSB[[imax,imax2]].max(axis=0), facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')
	# plt.plot(b * cMpcCoeffLutz.mean(), sbLutzmeanSB.mean(axis=0), c='dodgerblue', lw=lineW) #, label='Wisotzki+16 (LAEs)')
	# plt.fill_between(b * cMpcCoeffLutz.mean(),sbLutzmeanSB[imin], sbLutzmeanSB[[imax,imax2]].max(axis=0),
	#                 facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')

	twiny = False
	if twiny:
		# not working properly!!!
		plt.twiny(ax=None)
		ax.set_xscale('log', basex=basex)
		plt.xticks(xticks, np.array(xticks) * conv2)
		plt.xlabel(r'$\theta$ [kpc]', fontsize=fsize)

	plt.savefig('../../all/plots/SB_prof.pdf', dpi=dpi, format='pdf')
# plt.show()

if zneig:
	cat = '../../all/cats/laes.fits'
	data = getdata(cat, 1)
	npair = data['npairs']
	good = np.isfinite(npair)
	npair = npair[good]
	z = data['Z'][good]
	am = []
	mm = []
	t = []
	b = []
	bins = []
	bin = .3
	for i in np.arange(2.9, 4.1, bin):
		p = npair[(z > i) & (z <= i + bin)]
		mm.append(np.median(p))
		t.append(np.percentile(p, 75))
		b.append(np.percentile(p, 25))
		bins.append(i + bin / 2.)
		print mm, bin / 2.
	am = np.array([bins, mm, b, t]).T
	yerr = [np.array(mm) - np.array(b), np.array(t) - np.array(mm)]
	np.savetxt('../../all/cats/nvsz.dat', am, header='bins median n_25 n_75')
	plt.scatter(bins, mm)
	plt.errorbar(bins, mm, yerr=yerr, fmt='o', color='red')
	plt.axvline(x=3.5, linestyle='--')
	plt.xlabel("z")
	plt.ylabel(r"$n_{\rm{neighbours}}$", fontsize=18)
	plt.savefig("../../AAAPaper/zn.png")
	plt.show()

if paperims:
	overwrite = True
	title = ''
	makeim = True
	highq = True
	dpi = 120
	hdu = PrimaryHDU()
	arrow = True
	show = True
	zw0 = 9
	offset = 0  # 1


	def cubim(data, name, zw0=3, offset=0):
		zl, xl, yl = data.shape
		zw = '%d:%d' % (zl / 2 - zw0 + 1 + offset, zl / 2 + zw0 + 1 + offset)
		hdu.data = data
		# print 'cubim', name
		# hdu.writeto(name, clobber=True)
		s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (name, zw, name.replace('.fits', '.IM.fits'))
		print s
		os.system(s)


	# astroim(name.replace('.fits', '.IM.fits'), smooth=smooth, saveim=makeim, show=False, cbfrac=.08, pad=.006,
	#        dfig=(8, 10), contours=contours, scb_label=label, title=title, vmin=vmi, vmax=vma, std=stdd,
	#        highq=highq, sb=sb, y0=0, x0=0, gray=True, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, arrow=arrow,
	#        pcolor=pcolor, imout=imout, text=text, textpos=textpos)

	print 'Paper images'

	if 1:
		folders = ['../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100/',
				   '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n8-100/',
				   '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-8/',
				   '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l91-2000_n1-100/',
				   '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-91_n1-100/',
				   '../../all/stacks/d0-20_pd32-80_z2.9-4.0_l1-2000_n1-100/',
				   '../../all/stacks/d0-20_pd16-32_z2.9-4.0_l1-2000_n1-100/',
				   '../../all/stacks/d8-20_pd16-80_z2.9-4.0_l1-2000_n1-100/',
				   '../../all/stacks/d0-8_pd16-80_z2.9-4.0_l1-2000_n1-100/',
				   '../../all/stacks/d0-20_pd16-80_z3.5-4.0_l1-2000_n1-100/',
				   '../../all/stacks/d0-20_pd16-80_z2.9-3.5_l1-2000_n1-100/',
				   '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n8-100/']
		nims = len(folders)
		#        xwidths = [40, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40]
		#        ywidths = [40, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 7]
		xwidths = [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]
		ywidths = [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30]

		stds = np.concatenate(([1.4], [2.] * (nims - 2), [1.4]), 0)
		names = ['../../all/plots/stack-full.pdf',
				 '../../all/plots/stack-nmin8.pdf',
				 '../../all/plots/stack-nmax8.pdf',
				 '../../all/plots/stack-lmin91.pdf',
				 '../../all/plots/stack-lmax91.pdf',
				 '../../all/plots/stack-pdmin32.pdf',
				 '../../all/plots/stack-pdmax32.pdf',
				 '../../all/plots/stack-dmin8.pdf',
				 '../../all/plots/stack-dmax8.pdf',
				 '../../all/plots/stack-rmin35.pdf',
				 '../../all/plots/stack-rmax35.pdf',
				 '../../all/plots/stack-nmin8_2.png',
				 ]
		texts = [[r'$\rm{Full\,\,Sample}$', r'$\rm{390\,subcubes,\,96\,LAEs}$', r'$\rm{Oriented}$'],
				 [r'$\rm{Number\,of\,Neighbours\geqslant8}$', r'$\rm{Oriented}$'],
				 [r'$\rm{Number\,of\,Neighbours<8}$', r'$\rm{Oriented}$'],
				 [r'$\rm{Luminosity}>9.1\times10^{41}\rm{[erg/s]}$', r'$\rm{Oriented}$'],
				 [r'$\rm{Luminosity}<9.1\times10^{41}\rm{[erg/s]}$', r'$\rm{Oriented}$'],
				 [r'$32 < \theta < 80\,\rm{arcsec}$', r'$\rm{Oriented}$'],
				 [r'$16 < \theta < 32\,\rm{arcsec}$', r'$\rm{Oriented}$'],
				 [r'$8<d<20\,\rm{cMpc}$', r'$\rm{Oriented}$'], [r'$0.5<d<8\,\rm{cMpc}$', r'$\rm{Oriented}$'],
				 [r'$3.5<z<4$', r'$\rm{Oriented}$'], [r'$2.9<z<3.5$', r'$\rm{Oriented}$'], None]
		tpos = [[[-xwidths[0] * .95, -ywidths[0] * .9], [-xwidths[0] * .95, ywidths[0] * .95],
				 [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .89], [xwidths[0] * .6, -ywidths[0] * .89]],
				[[-xwidths[0] * .95, -ywidths[0] * .89], [xwidths[0] * .6, -ywidths[0] * .89]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				[[-xwidths[0] * .95, -ywidths[0] * .9], [xwidths[0] * .6, -ywidths[0] * .9]],
				None]

		texts2 = [[r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  [r'$\rm{Oriented - Random}$'],
				  None]
		tpos2 = [[[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 [[xwidths[0] * .2, -ywidths[0] * .9]],
				 None]
		regions = [None, [5, 14, -3, 3], [5, 14, -3, 3], [5, 14, -3, 3], [5, 14, -3, 3], [5, 14, -3, 3],
				   [5, 14, -3, 3], [5, 14, -3, 3], [5, 14, -3, 3], [5, 14, -3, 3], [5, 14, -3, 3],
				   None]  # [5, 14, -3, 3]]
		arrows = np.concatenate(([True], [False] * (nims - 1)), 0)
		scale = np.concatenate(([False], [False] * (nims - 1)), 0)
		scale2 = np.concatenate(([True] * (nims - 1), [False]), 0)

		colors = ['purple', 'red', 'blue', 'red', 'blue', 'red', 'blue', 'red', 'blue', 'red', 'blue', 'purple']
		vmin = -.05
		vmax = .5
		xmin, xmax, ymin, ymax, zmin, zmax, offset = [-20, 20, -20, 20, -3, 3, 1]

		show = False
		zw0 = zmax

		for f, xw, yw, n, s, t, tp, t2, tp2, a, sc, sc2, r, c in zip(folders, xwidths, ywidths, names,
																	 stds, texts, tpos,
																	 texts2, tpos2, arrows, scale, scale2,
																	 regions, colors)[-2:]:
			# a = open(f+'stack.lst', 'r')
			# ss = ''
			# for line in a: ss += line
			# '%d subcubes' % nsc, '%d LAEs' % nlae, #[w*.45, -w*.9], [w*.49, -w*.8],
			# b=re.findall('\d+', ss)
			# nsc = len(b)/2
			# nlae = len(np.unique(b))
			fit = getdata(f + 'stack.fits')
			cubim(fit, f + 'stack.fits', zw0, offset=offset)
			# x1, x2, y1, y2, z1, z2 = [40+11, 40+14, 38, 43]

			# region = None
			astroim(f + 'stack.IM.fits', smooth=True, saveim=True, show=show,
					cbfrac=.08, pad=.006, dfig=(8, 10), contours=True,
					x0=0, y0=0, imout=n, std=s, scale=sc,
					text=t, textpos=tp, regcolor=c,
					scb_label=r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}',
					title='', vmin=vmin, vmax=vmax, gray=True, arrow=a, region=r,
					xmin=-xw, xmax=xw, ymin=-yw, ymax=yw, pcolor='white', highq=False, dpi=dpi, nbins=5)

			if 0:
				fit = getdata(f + 'stack.randtest.fits')
				cubim(fit, f + 'stack.randtest.fits', offset=offset)
				astroim(f + 'stack.randtest.IM.fits', smooth=True, saveim=True, show=show,
						cbfrac=.08, pad=.006, dfig=(8, 10), contours=False,
						x0=0, y0=0, imout=n.replace('.pdf', '_randtest.pdf'),
						scb_label=r'$\rm{SB}\,\rm{[}10^{-%s}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}', regcolor=c,
						title='', vmin=vmin, vmax=vmax, gray=True, arrow=a, region=r,
						xmin=-xw, xmax=xw, ymin=-yw, ymax=yw, pcolor='black', highq=True, dpi=dpi,
						text=t2, textpos=tp2, scale=sc2, nbins=5)

	if 0:
		# os.system('python images.py -single True -id1 139 -id2 585 -prename .pair -angle True')
		# os.system('python images.py -single True -id1 144 -id2 216 -prename .pair -angle True')
		# os.system('python images.py -single True -id1 139 -id2 585 -prename .pair -angle True')

		files = ['../../HDFS/pairs/not-normalized/masked/cubical_pair139-585.pair.fits',
				 '../../HDFS/pairs/not-normalized/masked/cubical_pair422-501.pair.fits',
				 '../../HDFS/pairs/not-normalized/masked/cubical_pair437-501.pair.fits',
				 '../../UDF/pairs/not-normalized/masked/cubical_pair97-6292.pair.fits',
				 '../../UDF/pairs/not-normalized/masked/cubical_pair118-197.pair.fits',
				 '../../UDF/pairs/not-normalized/masked/cubical_pair391-6676.pair.fits',
				 '../../UDF/pairs/not-normalized/masked/cubical_pair6296-6672.pair.fits']

		angles = [-152.17, 170.85, 23.32, 98.08, 73.32, 73.32, 73.32]
		xmins = [20, 2, 10, None, None, None, None]
		xmaxs = [90, 90, 30, None, None, None, None]
		ymins = [10, 10, 10, None, None, None, None]
		ymaxs = [70, 70, 50, None, None, None, None, None]
		ids = [[139, 585, 216], None, None, None, None, None, None]
		idpos = [[[30, 40], [82, 40], [50, 55]], None, None, None, None, None, None]

		names = ['../../all/plots/HDFS-pair_139-585.pdf',
				 '../../all/plots/HDFS-pair_422-501.pdf',
				 '../../all/plots/HDFS-pair_437-501.pdf',
				 '../../all/plots/UDF-pair_97-6292.pdf',
				 '../../all/plots/UDF-pair_118-197.pdf',
				 '../../all/plots/UDF-pair_391-6676.pdf',
				 '../../all/plots/UDF-pair_6296-6672.pdf']

		vmin = -2
		vmax = 12
		zw0 = 2
		offset = 0
		for f, a, xmin, xmax, ymin, ymax, n, ii, ip in zip(files, angles, xmins, xmaxs, ymins, ymaxs, names, ids,
														   idpos):
			fit = getdata(f)
			cubim(fit, f, zw0, offset=offset)
			astroim(f.replace('.fits', '.IM.fits'), smooth=True, saveim=True, show=show,
					cbfrac=.08, pad=.006, dfig=(8, 10), contours=True,
					imout=n, std=None, pair=True,
					text=ii, textpos=ip,
					scb_label=r'Flux [$10^{-20}\,\rm{erg/s/cm^2}$]',
					title='', vmin=vmin, vmax=vmax, gray=True, arrow=False, angle=a,
					xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, pcolor='white', highq=True)

	if 0:

		ids1 = [422, 43, 162, 433, 200, 106, 63, 318, 344, 277, 437, 148, 97, 171]
		ids2 = [501, 162, 449, 520, 503, 6291, 197, 510, 6307, 6297, 501, 400, 6292, 633]
		fitcats = ['HDFS', 'HDFS', 'HDFS', 'HDFS', 'HDFS', 'UDF', 'UDF', 'UDF', 'UDF', 'UDF', 'HDFS', 'UDF', 'UDF',
				   'UDF']
		nims = len(ids1)
		# ids1 = [200, 246, 308, 422, 437, 97, 118, 148, 171, 109]
		# ids2 = [501, 501, 6292, 197, 400, 633]
		xw = 10
		yw = 10
		xmin = -xw
		xmax = xw
		ymin = -yw
		ymax = yw
		data = getdata('../../all/cats/laes_pairs.fits', 1)
		angles = data['angle']
		ii1 = data['id1']
		ii2 = data['id2']
		angs = [angles[(ii1 == i1) & (ii2 == i2)][0] for i1, i2 in
				zip(ids1, ids2)]  # [0.291037, 1.570796, -0.245091, 1.325986, -2.466991, 2.356512]
		angs2 = np.pi + np.array(angs)
		angs = [None] * nims
		angs2 = angs
		stds = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 4, 1.2, 1.2, 1.2, 1.2, 3.5, 1.2, 2.5]
		zw0 = 2
		offset = 0
		vmin = -.1
		vmax = 1
		show = False
		nsigma = 6

		if 0:
			for i in range(nims):
				s = 'python images.py -pair False -single True -xw0 %d -yw0 %d -angle True' % (xw, yw) + \
					' -id1 %d -id2 %d -overwrite True -zw0 %d -offset %d -fitcat %s -vmin %f -vmax %f' \
					% (ids1[i], ids2[i], zw0, offset, fitcats[i], vmin, vmax)
				print s
				os.system(s)

		if 1:
			for i in (range(nims)):
				astroim(
					'../../%s/pairs/not-normalized/masked/cubical_pair%d-%d.IM.fits' % (fitcats[i], ids1[i], ids2[i])
					, smooth=True, saveim=True, show=show,
					cbfrac=.08, pad=.006, dfig=(8, 10), contours=True, highq=True,
					x0=0, y0=0, std=stds[i], imout='../../all/plots/%s_%d-%d.pdf' % (fitcats[i], ids1[i], ids2[i]),
					text=['%s %d' % (fitcats[i], ids1[i])], textpos=[[-xw * .9, -yw * .9]],
					title='', vmin=vmin, vmax=vmax, gray=True, angle=angs[i],
					xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, pcolor='white', nsigma=nsigma)
				astroim(
					'../../%s/pairs/not-normalized/masked/cubical_pair%d-%d.IM.fits' % (fitcats[i], ids2[i], ids1[i])
					, smooth=True, saveim=True, show=show,
					cbfrac=.08, pad=.006, dfig=(8, 10), contours=True, highq=True,
					x0=0, y0=0, std=stds[i], imout='../../all/plots/%s_%d-%d.pdf' % (fitcats[i], ids2[i], ids1[i]),
					text=['%s %d' % (fitcats[i], ids2[i])], textpos=[[-xw * .9, -yw * .9]],
					title='', vmin=vmin, vmax=vmax, gray=True, angle=angs2[i],
					xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, pcolor='white', nsigma=nsigma)

if dndz:
	print 'aa'
	zw0 = 2
	zoff = 1
	zf = zl / 2 + zw0 + 1 + zoff
	zi = zl / 2 - zw0 + zoff
	x1, x2, y1, y2, z1, z2 = [40 + 5, 40 + 14, 38, 43, zi, zf]

	uniq = []
	anun = np.squeeze(nuns[:, z1: z2, y1: y2, x1: x2])
	for an in anun:
		aaaa = [int(an[i:i + 10]) for i in np.arange(0, len(an) / 10 + 1, 10)]
		for aaa in aaaa: uniq.append(aaa)
	if 0:
		def p(r): return 5. / (np.pi * (2 * r + 1))


		N = 5 * 390
		a = [[np.power(p(r), i) * np.power(1 - p(r), 4 - i) / i for i in np.arange(1, 4)] for r in np.arange(5, 11)]
		print np.sum(a) / np.power(N, 4.)

if hm:

	def p_lya(T4=1, case='B'):
		# T4 is Temperature/1e4K, taken from Cantalupo+08
		if case == 'A': return .41 - .165 * np.log10(T4) - .015 * np.power(T4, -.44)
		if case == 'B': return .686 - .106 * np.log10(T4) - .009 * np.power(T4, -.44)

	gamma_bkg = {}

	gamma_bkg['HM12'] = {'7': [3.60E-13, 2.14E-13, 3.23E-18],
						 '8': [4.26E-13, 2.45E-13, 1.33E-17],
						 '9': [4.86E-13, 2.81E-13, 5.60E-17],
						 '10': [5.71E-13, 3.28E-13, 2.31E-16],
						 '11': [6.77E-13, 3.89E-13, 8.70E-16],
						 '12': [7.66E-13, 4.42E-13, 2.10E-15],
						 '13': [9.50E-13, 5.55E-13, 9.06E-15],
						 '14': [9.64E-13, 5.67E-13, 1.13E-14]}

	gamma_bkg['HM01'] = {'7': [3.60E-13, 2.14E-13, 3.23E-18],
						 '8': [.6E-12, 2.45E-13, 1.33E-17],
						 '9': [.7E-12, 2.81E-13, 5.60E-17],
						 '10': [.9E-12, 3.28E-13, 2.31E-16],
						 '11': [1E-12, 3.89E-13, 8.70E-16],
						 '12': [1.5E-12, 4.42E-13, 2.10E-15],
						 '13': [1.7E-12, 5.55E-13, 9.06E-15],
						 '14': [1.9E-12, 5.67E-13, 1.13E-14]}

	def calc(z, snap):
		p0 = max(0, max(np.where(z > reds)[0]))
		p1 = min(len(reds) - 1, min(np.where(z < reds)[0]))
		wmin = 200
		wmax = 912
		cool = np.where((wav >= wmin) & (wav <= wmax))[0]
		ncool = len(cool)
		m = (z - reds[p0]) / (reds[p1] - reds[p0])
		J = data['col%d' % (p0 + 2)] * (1 - m) + data[
			'col%d' % (p1 + 2)] * m  # type: Union[int, Any] # +2 since J columns start from #2
		suma = 0
		suma2 = 0
		def sigma(v):
			x = 3.041e-16*v#h*1Hz/13.6eV * v
			return 1.34/x**2.99-.34/x**3.99
			#return (v*3.041e-16)**-3

		for i in range(ncool - 1):
			k0, k1 = [cool[i], cool[i + 1]]
			_J = (J[k0] + J[k1])/2.
			v0 = nu[k0]
			v1 = nu[k1]
			v = (v0+v1)/2.
			dv = v0-v1
			_suma = _J / v * dv
			suma += _suma
			suma2 += _suma * sigma(v)
		R_HM12 = suma / h
		sigma0 = 6.3e-18 # for HI
		gamma_HM12 = sigma0 * 4 * np.pi * suma2 / h

		#wmin = 100
		#wmax = 1e10


		factor = et * 3.8397e-22 / (1. + z) ** 4
		#div = gamma_bkg['HM12'][snap][0]*1e16/R_HM12
		# factor: result of 4 pi^2*c*h/(1215.67 angstrom)/(360*60*60)^2 in ergs and times et
		#print 'redshift', z, 'div', div, 'Gamma HM12 integral', gamma_HM12, 'Gamma_HM12 paper', gamma_bkg['HM12'][snap][0], 'R_HM12', R_HM12, 'SB', R_HM12 * factor
		print 'redshift', z, 'Gamma HM12 integral', gamma_HM12, 'R_HM12', R_HM12, 'SB', R_HM12 * factor
		#return z, R_HM12, gamma_bkg['HM01'][snap][0]*1e16, gamma_bkg['HM12'][snap][0]*1e16



	# how many photons are coverted into lya
	et = p_lya()
	h = 6.62e-27
	#c = 2.998e10
	# lya_rest = 1215.67
	data = getdata('../../UVB/UVB_spec.fits', 1)
	wav = data['wavelength']
	nu = 2.998e18 / wav
	reds = np.loadtxt('../../UVB/redshift_cuba.dat')
	a = []
	redshifts = params.redshifts
	for s in np.arange(6, 15).astype(str):
		print s
		a.append(calc(redshifts[s], s))
	a = np.array(a).T
	if 0:
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(a[0], a[1], label='rm12')
		plt.plot(a[0], a[2], label='hm01')
		plt.plot(a[0], a[3], label='hm12')
		plt.legend()
		plt.savefig('../../UVB/UVB_models.png')
		plt.close()


	if 0:
		a = np.loadtxt('../../cuba_madau.dat')
		b = a.T
		l = 1 / b[0]
		J3 = b[29]
		J35 = b[31] * .6 + b[32] * .4
		lipos = 229
		lfpos = 913
		sum3 = 0
		sum35 = 0
		for i, j, k in zip(l[lipos:lfpos], J3[lipos:lfpos], J35[lipos:lfpos]):
			sum3 += j * i
			sum35 += k * i

		R3 = et / h * sum3
		R35 = et / h * sum35

		print R3, R35


def toSBthin(N=1e17):
	def sigma(ll):
		s0 = 6.3e-18
		return s0 * np.power(ll / li, -2.75)

	h = 6.62e-27
	c = 2.998e10
	a = np.loadtxt('../../cuba_madau.dat')
	b = a.T
	l = 1 / b[0]
	Lyaen = h * c / 1216.
	arcsec = np.power(360. * 60. * 60., -2.)
	# li = 228
	# lf = 912
	lipos = 226
	lfpos = 377
	li = l[lipos]
	J3 = b[29]
	J35 = b[31] * .6 + b[32] * .4
	sum3 = 0
	sum35 = 0
	et = .42
	for i, j, k in zip(l[lipos:lfpos], J3[lipos:lfpos], J35[lipos:lfpos]):
		s = sigma(i)
		sum3 += j * i * s
		sum35 += k * i * s

	sum3 *= N
	sum35 *= N
	R3 = et / h * sum3 * np.power(4., -4.) * Lyaen * arcsec
	R35 = et / h * sum35 * np.power(4.5, -4.) * Lyaen * arcsec

	print R3, R35


if spectra:
	prename = ''  # '.reject3'#'.test2noreject'#''.test1sclip3'#'.reject-abitless'#
	folder = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100%s/' % prename
	extraname = '_full'  # '_unique'#'_nmin8'

	##define collection of region sizes
	xw = 1
	yw = 2
	zw = 9
	zoff = 0  # 1
	width = 4
	sigmas = 2
	x1 = 1.7
	x2 = 13.6
	y1 = -1.5
	y2 = 10
	hmSB = 1.14
	y1b = -4
	y2b = 6

	fname = folder + "stack.fits"
	frname = folder + "random_stack.fits"
	israndom = os.path.isfile(frname)
	hdu = fits.open(fname)
	f = hdu[0].data
	zl, yl, xl = f.shape
	rads = [5, 7, 10]
	y, x = np.ogrid[-yl / 2 + 1: yl / 2 + 1, -xl / 2 + 1: xl / 2 + 1]
	specs = []
	for r in rads:
		cool = x ** 2 + y ** 2 <= r ** 2

if circhist:
	import glob

	warea = False
	fitcat = 'EAGLE'  # 'UDF'
	pairs = getdata('../../%s/cats/lae_pairs.fits' % fitcat, 1)
	laes = getdata('../../%s/cats/laes.fits' % fitcat, 1)
	outf = '../../' + fitcat
	xl, yl, zl = [4096, 4096, 0]  # [945, 947, 1300]
	nws = laes['nw']
	Nlae = len(nws)
	ids = laes['ID']
	id1 = pairs['id1']
	id2 = pairs['id2']
	x2 = pairs['x2']
	y2 = pairs['y2']
	z2 = pairs['z2']
	x1 = laes['x']
	y1 = laes['y']
	z1 = laes['z']
	if fitcat == 'EAGLE':
		com2pix = 163.84
		x1 *= com2pix
		y1 *= com2pix
		z1 *= com2pix

	angles = pairs['angle'] + np.pi
	cds = pairs['com_dist']
	nwsort = np.argsort(nws)[::-1]
	N = 30
	bottom = yl * .1
	max_height = yl * .2
	max_weight = max(nws) / N
	max_rad = yl
	theta = np.linspace(0., 2 * np.pi, N, endpoint=False)
	width = (2 * np.pi) / N
	delta = np.pi / float(N * 4)
	figfold = outf + '/plots/angle_hist/'
	if os.path.isdir(figfold):
		print 'Folder %s already exists.' % figfold
	else:
		glob.os.makedirs(figfold)
	catfold = outf + '/cats/angle_hist/'
	if os.path.isdir(catfold):
		print 'Folder %s already exists.' % catfold
	else:
		glob.os.makedirs(catfold)

	if fitcat == 'EAGLE':
		com2pix = 163.84
	# bottom /= 10
	# max_height /= 10
	# max_rad /= 10


	def weight(x, y, z, x2, y2, z2):
		r = np.sqrt((x2 - x) ** 2 + (y2 - y) ** 2)
		rmax = max(r)
		xd = np.abs(x2 - x)
		yd = np.abs(y2 - y)
		zd = np.abs(z2 - z)
		wl = (x + xd > xl) * .5 * np.abs(xl + xd - x) / xd
		wr = (x - xd < 0) * .5 * (x + xd) / xd
		wx = ((x + xd < xl) & (x - xd > 0)) + wl + wr
		wl = (y + yd > yl) * .5 * np.abs(yl + yd - y) / yd
		wr = (y - yd < 0) * .5 * (y + yd) / yd
		wy = ((y + yd < yl) & (y - yd > 0)) + wl + wr
		wl = (z + zd > zl) * .5 * np.abs(zl + zd - z) / zd
		wr = (z - zd < 0) * .5 * (z + zd) / zd
		wz = ((z + zd < zl) & (z - zd > 0)) + wl + wr
		w = wx * wy * wz * rmax / r
		return w


	for nj in range(Nlae):
		j = nwsort[nj]
		good = (id1 == ids[j]) & (cds > .5) & (cds < 20.) & (np.abs(x2 - x1[j]) > 3) & (np.abs(y2 - y1[j]) > 3) & (
				np.abs(y2 - y1[j]) > 3)
		ang = angles[good]
		cd = cds[good]
		if warea: wd = weight(x1[j], y1[j], z1[j], x2[good], y2[good], z2[good])
		w = np.zeros(N)
		for i in range(N):
			if i == 0:
				good2 = (ang % (2 * np.pi) > (theta[i] - delta) % (2 * np.pi)) | (
						ang % (2 * np.pi) <= (theta[i] + width + delta) % (2 * np.pi))
			elif i == N - 1:
				good2 = (ang % (2 * np.pi) > (theta[i] - delta) % (2 * np.pi)) | (
						ang % (2 * np.pi) <= (theta[i] + width + delta) % (2 * np.pi))
			else:
				good2 = (ang > theta[i] - delta) & (ang <= theta[i] + width + delta)
			if warea:
				ww = 1 / (cd[good2] * wd[good2])
			else:
				ww = 1 / cd[good2]
			w[i] = np.sum(ww[np.isfinite(ww)])
		radii = w * max_height / max_weight
		maxr = max(radii)
		print 'id', ids[j], 'max weight', max(w), 'max radii', max(radii)
		ax = plt.subplot(111, polar=True)
		ax.set_rmax(max_rad * 2)
		ax.set_yticklabels([])
		bars = ax.bar(theta + np.pi, radii, width=width + 2 * delta, bottom=bottom, edgecolor="none", linewidth=0)
		# Use custom colors and opacity
		for r, bar in zip(radii, bars):
			bar.set_facecolor(plt.cm.viridis(r / maxr))
			bar.set_alpha(0.8)
		dy = yl - y1[j]
		dx = xl - x1[j]
		rs = [dx ** 2 + y1[j] ** 2, dx ** 2 + dy ** 2, x1[j] ** 2 + dy ** 2, x1[j] ** 2 + y1[j] ** 2]
		rs = np.sqrt(np.array(rs))
		xs = np.array([dx, dx, -x1[j], -x1[j]])
		ys = np.array([-y1[j], dy, dy, -y1[j]])
		atans = [np.arctan2(yy, xx) for yy, xx in zip(ys, xs)]
		plt.plot((atans[0], atans[1]), (rs[0], rs[1]), '-', color="green", linewidth=2)
		plt.plot((atans[1], atans[2]), (rs[1], rs[2]), '-', color="green", linewidth=2)
		plt.plot((atans[2], atans[3]), (rs[2], rs[3]), '-', color="green", linewidth=2)
		plt.plot((atans[3], atans[0]), (rs[3], rs[0]), '-', color="green", linewidth=2)
		plt.savefig(figfold + 'anghist%d-id%d-nw%1.2f.png' % (nj, ids[j], nws[j]))
		plt.close()
		np.savetxt(catfold + 'id%d.dat' % ids[j], np.array([theta + np.pi, radii]).T, header='theta radius')
