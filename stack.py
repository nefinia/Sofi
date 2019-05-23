__author__ = 'gallegos'

import os
from os.path import isfile, isdir
from pyfits import getdata, PrimaryHDU
from sys import argv

import numpy as np

from tools_sofi import astroim, rdarg, stack, pdfim, analysis

doanalysis = rdarg(argv, 'analysis', bool, False)
binsize = rdarg(argv, 'binsize', int, 2)
contours = rdarg(argv, 'contours', bool, True)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)
cubexmask = rdarg(argv, 'cubexmask', bool, False)
dosclip = rdarg(argv, 'dosclip', bool, True)
extraname = rdarg(argv, 'extraname', str, '')
fitcat = rdarg(argv, key='fitcat', type=str, default='all')  # 'all'#'HDFS'
flipx = rdarg(argv, 'flipx', bool, False)
flipy = rdarg(argv, 'flipy', bool, False)
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
foldercat = rdarg(argv, 'foldercat', str, '../../')  # '/scratch/gallegos/MUSE/'
folderim = rdarg(argv, 'folderim', str, '/net/galaxy-data/export/galaxydata/gallegos/')  # '/scratch/gallegos/MUSE/'
badlist = rdarg(argv, 'badlist', str, 'badpairs.dat')  # '/scratch/gallegos/MUSE/'
fshear = rdarg(argv, 'flipshear', bool, False)
highq = rdarg(argv, 'highq', bool, False)
imtype = rdarg(argv, 'imtype', str, 'mean')  # 'mean'
jin = rdarg(argv, 'jin', int, 0)
jackniffe = rdarg(argv, 'jackniffe', bool, False)
mask = rdarg(argv, 'mask', bool, True)
makeim = rdarg(argv, 'makeim', bool, True)
makepdf = rdarg(argv, 'makepdf', bool, False)
meancorr = rdarg(argv, 'meancorr', bool, True)
nodist = rdarg(argv, 'nodist', bool, True)
norm = rdarg(argv, 'norm', bool, False)
nrand = rdarg(argv, 'nrand', int, 200)
ntest = rdarg(argv, 'ntest', int, 200)
overwrite = rdarg(argv, 'overwrite', bool, False)
parallel = rdarg(argv, 'parallel', bool, False)
prename = rdarg(argv, 'prename', str, '')
propvar = rdarg(argv, 'pvar', bool, False)
statvar = rdarg(argv, 'svar', bool, False)
random = rdarg(argv, 'random', bool, False)
randflipy = rdarg(argv, 'randflipy', bool, False)
reject = rdarg(argv, 'reject', bool, False)
rotate = rdarg(argv, 'rotate', bool, False)
sb = rdarg(argv, 'sb', bool, False)
scalelims = rdarg(argv, 'scalelims', str, '-0.01 0.01')
sclip = rdarg(argv, 'sclip', int, 3)
# scliptype 0: one single std for the full cube
# scliptype 1: std for each pixel in the sub-cube set
# scliptype 2: std for each z layer with mask in the LAE
scliptype = rdarg(argv, 'scliptype', int, 1)
sconfmin = rdarg(argv, 'sconfmin', int, 1)
simon = rdarg(argv, 'cubical', bool, True)
std = rdarg(argv, 'std', float, 2)  # 0.0707064#
tests = rdarg(argv, 'tests', bool, False)
unique = rdarg(argv, 'unique', bool, False)
unique2 = rdarg(argv, 'unique2', bool, False)
unique3 = rdarg(argv, 'unique3', bool, False)
# using unique pixels
unique4 = rdarg(argv, 'unique4', bool, False)
unique5 = rdarg(argv, 'unique5', bool, False)
verbose = rdarg(argv, 'verbose', int, 1)
weights = rdarg(argv, 'weights', bool, False)

vmin = rdarg(argv, 'vmin', float, -.05)
vmax = rdarg(argv, 'vmax', float, .5)
# if fitcat == 'all': vmin, vmax, std, sconfmin = [-.05, .5, 2, 1]
# if fitcat == 'mosaic': vmin, vmax, std, sconfmin = [-.1, 1, 4, 2]
# if fitcat == 'HDFS': vmin, vmax, std, sconfmin = [-.4, 4, 10, 1]

# if unique4: vmin, vmax, std = [vmin*3, vmax*3, std*3]

xmin = rdarg(argv, 'xmin', float, -40)
xmax = rdarg(argv, 'xmax', float, 40)
ymin = rdarg(argv, 'ymin', float, -40)
ymax = rdarg(argv, 'ymax', float, 40)
zmin = rdarg(argv, 'zmin', float, -20)
zmax = rdarg(argv, 'zmax', float, 20)
yw0 = rdarg(argv, 'yw0', int, 2)
zw0 = rdarg(argv, 'zw0', int, 2)
offset = rdarg(argv, 'offset', int, 1)

if verbose < 2:
	vb = ' > /dev/null 2>&1'
else:
	vb = ''

if cubesmooth:
	smooth = False
	ssmooth = ' -boxsm 3 '
else:
	smooth = True
	ssmooth = ''
if imtype == 'median': extraname = '.' + imtype
if unique: extraname += '.unique'
if not prename: prename = ''
if not extraname: extraname = ''

# d distance
# p proj distance
# r redshift
# los line of sight distance
# n number of neighbours
# vel velocity
dmin = rdarg(argv, 'dmin', float, .5)
dmax = rdarg(argv, 'dmax', float, 20)
pmin = rdarg(argv, 'pmin', float, 16)
pmax = rdarg(argv, 'pmax', float, 2000)
rmin = rdarg(argv, 'rmin', float, 2.9)
rmax = rdarg(argv, 'rmax', float, 4.0)
lmin = rdarg(argv, 'lmin', float, 0)
lmax = rdarg(argv, 'lmax', float, 2000)
losmin = rdarg(argv, 'losmin', float, 0.5)
losmax = rdarg(argv, 'losmax', float, 20)
nmin = rdarg(argv, 'nmin', int, 1)
d5min = rdarg(argv, 'd5min', int, 0)
d5max = rdarg(argv, 'd5max', int, 20)
nwmin = rdarg(argv, 'nwmin', int, 0)
velmin = rdarg(argv, 'velmin', float, 0)
velmax = rdarg(argv, 'velmax', float, 10000)

if unique3:
	unique = True
	xmin, xmax, ymin, ymax, zmin, zmax = [-30, 30, -30, 30, -10, 10]

if cubexmask:
	prename = '.Residuals'
	vmin = -.5
	vmax = .5

if mask:
	sm = 'masked'
else:
	sm = 'unmasked'

if weights:
	sw = 'weighted'
else:
	sw = 'unweighted'

if norm:
	sn = 'normalized'
else:
	sn = 'not-normalized'

if simon:
	st = 'cubical'
else:
	st = 'cylindrical'

if norm: lenx = 40

if fitcat == 'all':
	cat2 = '%s/%s/cats/lae_pairs.fits' % (foldercat, 'HDFS')
	data = getdata(cat2, 1)
	ids1 = data['id1']
	ids2 = data['id2']

	laedata = getdata('%s/%s/cats/laes.fits' % (foldercat, 'HDFS'), 1)
	flae = laedata['LYALPHA_LUM']
	idlae = laedata['ID']
	# fhst = laedata['HST_F606W']
	lum_lae1 = []
	lum_lae2 = []
	for i, j in zip(ids1.astype(int), ids2.astype(int)):
		lum_lae1.append(np.float(flae[idlae == i]))
		lum_lae2.append(np.float(flae[idlae == j]))
	sconf1 = data['sconf1'].astype(int)
	sconf2 = data['sconf2'].astype(int)
	pdists = data['x2'] - data['x1']
	zs1 = data['z1']
	zs2 = data['z2']
	redshift = (data['redshift1'] + data['redshift2']) / 2.
	dists = data['pi_Mpc']
	vel = data['pi_v']
	theta = data['theta']
	fitcats = ['HDFS'] * len(data)
	dclose = np.where((dists > .5) & (dists <= 20))
	npairs1 = np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids1])
	npairs2 = np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids2])
	cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'UDF')
	data = getdata(cat2, 1)
	ids1 = np.concatenate((ids1, data['id1']), 0)
	ids2 = np.concatenate((ids2, data['id2']), 0)
	laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'UDF'), 1)
	flae = laedata['LYALPHA_LUM']
	idlae = laedata['ID']
	for i, j in zip(data['id1'].astype(int), data['id2'].astype(int)):
		lum_lae1.append(np.float(flae[idlae == i]))
		lum_lae2.append(np.float(flae[idlae == j]))
	sconf1 = np.concatenate((sconf1, data['sconf1'].astype(int)), 0)
	sconf2 = np.concatenate((sconf2, data['sconf2'].astype(int)), 0)
	pdists = np.concatenate((pdists, data['x2'] - data['x1']), 0)
	zs1 = np.concatenate((zs1, data['z1']), 0)
	zs2 = np.concatenate((zs2, data['z2']), 0)
	redshift = np.concatenate((redshift, (data['redshift1'] + data['redshift2']) / 2.), 0)
	dists = np.concatenate((dists, data['pi_Mpc']), 0)
	theta = np.concatenate((theta, data['theta']), 0)
	nhdfs = len(fitcats)
	fitcats = np.concatenate((fitcats, ['UDF'] * len(data)), 0)
	lum_lae2 = np.array(lum_lae2)
	lum_lae1 = np.array(lum_lae1)
	dclose = np.where((data['pi_Mpc'] > .5) & (data['pi_Mpc'] <= 20))
	npairs1 = np.concatenate((npairs1, data['npairs1']))
	npairs2 = np.concatenate((npairs2, data['npairs2']))
	vel = np.concatenate((vel, data['pi_v']), 0)

else:
	if fitcat == 'mosaic10':
		cat2 = '%s/%s/cats/lae_pairs.fits' % (foldercat, 'UDF')
	else:
		cat2 = '%s/%s/cats/lae_pairs.fits' % (foldercat, fitcat)
	if fitcat == 'mosaic-b' or fitcat == 'udf10-mosaic': cat2 = '%s/%s/cats/lae_pairs.fits' % (foldercat, 'mosaic')
	data = getdata(cat2, 1)
	ids1 = data['id1']
	ids2 = data['id2']
	sconf1 = data['sconf1'].astype(int)
	sconf2 = data['sconf2'].astype(int)
	if fitcat == 'mosaic10':
		laedata = getdata('%s/%s/cats/laes.fits' % (foldercat, 'UDF'), 1)
		fitcat = 'mosaic'
		extraname += '.UDF10'
	if fitcat == 'mosaic-b' or fitcat == 'udf10-mosaic':
		laedata = getdata('%s/%s/cats/laes.fits' % (foldercat, 'mosaic'), 1)
	else:
		laedata = getdata('%s/%s/cats/laes.fits' % (foldercat, fitcat), 1)
	idlae = laedata['ID']
	if fitcat == 'EAGLE':
		flae = idlae - idlae + 1
	else:
		fhst = laedata['HST_F606W']
		flae = laedata['LYALPHA_LUM']
		fhst1 = []
		fhst2 = []
		for i, j in zip(ids1, ids2):
			fhst1.append(np.float(fhst[idlae == i]))
			fhst2.append(np.float(fhst[idlae == j]))
	lum_lae1 = data['lum1']
	lum_lae2 = data['lum2']
	pdists = data['x2'] - data['x1']
	zs1 = data['z1']
	zs2 = data['z2']
	redshift = (data['redshift1'] + data['redshift2']) / 2.
	dists = data['com_dist']  # data['pi_Mpc']
	los_dist = data['pi_Mpc']
	theta = data['theta']
	vel = data['pi_v']
	if fitcat == 'mosaic-b': fitcats = ['mosaic'] * len(data)
	if fitcat == 'udf10-mosaic':
		fitcats = ['UDF'] * len(data)
	else:
		fitcats = [fitcat] * len(data)
	npairs1 = data['npairs1']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids1])
	npairs2 = data['npairs2']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids2])
	npairs = data['npairs1']
	nw1 = data['nw1']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids1])
	nw2 = data['nw2']  # np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids2])
	# d5th1 = data['d5th1']
	# d5th2 = data['d5th2']
	angles = data['angle'] + np.pi
	xs1 = data['x1']
	xs2 = data['x2']
	ys1 = data['y1']
	ys2 = data['y2']

if fitcat == 'udf10-mosaic':
	extraname += '.' + fitcat
	fitcat = 'UDF'

h = 'theta '
hdu = PrimaryHDU()

if 0:  # nodist:
	temp = np.copy(ids1)
	ids1 = np.concatenate((temp, ids2), 0).astype(int)
	ids2 = np.concatenate((ids2, temp), 0).astype(int)
	temp = np.copy(sconf1)
	sconf1 = np.concatenate((temp, sconf2), 0).astype(int)
	sconf2 = np.concatenate((sconf2, temp), 0).astype(int)
	redshift = np.concatenate((redshift, redshift), 0)
	dists = np.concatenate((dists, dists), 0)
	theta = np.concatenate((theta, theta), 0)
	vel = np.concatenate((vel, vel), 0)
	temp = np.copy(lum_lae1)
	lum_lae1 = np.concatenate((temp, lum_lae2), 0)
	lum_lae2 = np.concatenate((lum_lae2, temp), 0)
	fitcats = np.concatenate((fitcats, fitcats), 0)
	temp = np.copy(zs1)
	zs1 = np.concatenate((temp, zs2), 0)
	zs2 = np.concatenate((zs2, temp), 0)
	npairs = np.concatenate((npairs1, npairs2), 0)

	temp = np.copy(xs1)
	xs1 = np.concatenate((temp, xs2), 0)
	xs2 = np.concatenate((xs2, temp), 0)
	temp = np.copy(ys1)
	ys1 = np.concatenate((temp, ys2), 0)
	ys2 = np.concatenate((ys2, temp), 0)

rejected = False
if reject:
	badpairs = np.loadtxt('%s/%s/cats/%s' % (folder, fitcat, badlist))
	print 'Bad pairs list', badlist

	for bp in badpairs:
		rejected |= (ids1 == bp[0]) & (ids2 == bp[1])

close = np.where((redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax) \
				 & (theta > pmin) & (lum_lae1 >= lmin) & (lum_lae1 < lmax)  # & (d5th1 > d5min) & (d5th1 < d5max) \
				 & (sconf1 >= sconfmin) & (sconf2 >= sconfmin) & (~rejected))[0]  # \
# & (npairs1 >= nmin) & (npairs1 < nmax) & (vel >= velmin) & (vel < velmax) & (nw1 > nwmin)
# & (fhst1>0))[0]

if unique5:
	close = (redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax) \
			& (theta > pmin) & (lum_lae1 >= lmin) & (lum_lae1 < lmax) \
			& (sconf1 >= sconfmin) & (sconf2 >= sconfmin) & (~rejected) \
			& (npairs1 >= nmin) & (npairs1 < nmax) & (vel >= velmin) & (vel < velmax) & (nw1 > nwmin)
	reg = False
	print 'unique 5'
	for i in idlae:
		fff = '../../mosaic/cats/angle_hist/id%d.dat' % i
		if os.path.isfile(fff):
			a = np.loadtxt(fff).T
			a[0] = a[0] - np.pi
			width = a[0][1] - a[0][0]
			# print a[0], a[1]
			wmean = np.nanmean(a[1])
			wstd = np.nanstd(a[1])
			angs = a[0][a[1] > 2 * wstd + wmean]
			# angs = np.array(angs)
			for _a in angs:
				reg |= ((angles > _a) & (angles <= _a + width) & (ids1 == i))
			print i, angs, np.sum(reg)

	close &= reg
	print np.sum(close), 'aaaaaaaaaaaaaaaaaaaaaaaaa'
	close = np.where(close)[0]

if fitcat == 'mosaic-b':
	reg = False
	xl = [1, 316, 631, 946]
	yl = [1, 317, 632, 948]
	close = (redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax) \
			& (theta > pmin) & (lum_lae1 >= lmin) & (lum_lae1 < lmax) \
			& (sconf1 >= sconfmin) & (sconf2 >= sconfmin) & (~rejected) \
			& (npairs1 >= nmin) & (npairs1 < nmax) & (vel >= velmin) & (vel < velmax)
	for i in range(len(xl) - 1):
		for j in range(len(yl) - 1):
			reg |= ((redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax) \
					& (theta > pmin) & (lum_lae1 >= lmin) & (lum_lae1 < lmax) \
					& (sconf1 >= sconfmin) & (sconf2 >= sconfmin) & (~rejected) \
					& (npairs1 >= nmin) & (npairs1 < nmax) & (vel >= velmin) & (vel < velmax) &
					(xs1 > xl[i]) & (xs1 <= xl[i + 1]) & (ys1 > yl[j]) & (ys1 <= yl[j + 1])
					& (xs2 > xl[i]) & (xs2 <= xl[i + 1]) & (ys2 > yl[j]) & (ys2 <= yl[j + 1]))
	close &= reg
	close = np.where(close)[0]
	print '*******************'
	print len(close)
	print '*******************'
	extraname = '.mosaic-b'
	fitcat = 'mosaic'

if len(close) > 1:
	id1 = ids1[close]
	uid = np.unique(id1)
	ngal = len(uid)
	print 'Number of galaxies to be stacked', ngal

	id2 = ids2[close]
	z1 = zs1[close]
	z2 = zs2[close]
	fcat = np.array(fitcats)[close]
	dist = dists[close]
	theta = theta[close]
	red = redshift[close]
	npair = npairs1[close]
	lst = []

	random_lst = [[] for i in range(nrand)]

	if unique:

		ucat = np.unique(fcat)
		dclose = []
		for i in range(len(id1)):
			iis = (id1 == id1[i]) & (fcat == fcat[i])
			if dist[i] == min(dist[iis]):
				dclose.append(i)
		dclose = np.array(dclose)
		id1 = id1[dclose]
		id2 = id2[dclose]
		z1 = z1[dclose]
		z2 = z2[dclose]
		fcat = fcat[dclose]
		dist = dist[dclose]
		thet = theta[dclose]
		red = red[dclose]
		npair = npair[dclose]
		print 'Number of subcubes to be stacked', len(dclose)
		close = dclose

	nsub = len(close)
	print 'Number of subcubes to be stacked', nsub

	#overlap_roi = 2 * (12 - 6) / np.pi / (12 ** 2 - 6. ** 2) # overlap for my region of interest (2" tall and bwt 6" and 12")
	#print 'Effective decrease in noise', nsub - (nsub - ngal)*overlap_roi, 'average overlaping fraction', overlap_roi

if 1:

	for i in range(len(id1)):

		if unique3:
			pair = '%s/%s/LAEs/lae%d.fits' % (folderim, fcat[i], id1[i])
			lst.append(pair)

		else:
			# pair = '%s/%s/pairs/%s/%s/%s_%d-%d%s.fits' % (folderim, fcat[i], sn, sm, st, id1[i], id2[i], prename)
			pair = '%s/%s/pairs/%d-%d%s.fits' % (folderim, fcat[i], id1[i], id2[i], prename)
			# print "pair", pair
			if os.path.isfile(pair): lst.append(pair)
			# else: print '.',#"Pair %s does not exist" % pair

			if random:
				for j in range(nrand):
					pair = '%s/%s/pairs/%s/%s/random%d/random_%s_pair%d-%d%s.%d.fits' % (
						folderim, fcat[i], sn, sm, j, st, id1[i], id2[i], prename, j + jin)
					if os.path.isfile(pair):
						random_lst[j].append(pair)
					else:
						print '',  # "Random pair %s does not exist" % pair
	nstack = len(lst)

	print '\n Number of existing subcubes to be stacked', nstack

	foldername = '%s/%s/stacks/d%d-%d_los%d-%d_pd%d-%d_z%.1f-%.1f_l%d-%d_d5th%d-%d_v%d-%d%s%s' % (
		folder, fitcat, dmin, dmax, losmin, losmax,
		pmin, pmax, rmin, rmax, lmin, lmax, d5min, d5max, velmin, velmax, prename, extraname)

	print "Output files in", foldername

	if not isdir(foldername): os.system('mkdir %s' % foldername)

	lstname = '%s/stack.lst' % (foldername)
	flst = open(lstname, 'w')
	for l in lst: flst.write('%s\n' % l)
	flst.close()

	stackname = lstname.replace('.lst', '.fits')

	title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$,' % (
		nstack, dmin, dmax, pmin, pmax, rmin, rmax) + '\n' + \
			r'$%d<L<%d$, $%d<v<%d$ %s' % (lmin, lmax, velmin, velmax, extraname)
	# title = ''
	f, nf, var = stack(lst, stackname, imtype, sclip, zw0, makeim, title, vmin, vmax, npix=True, var=True,
					   flipshear=fshear, corrmean=meancorr, flipy=flipy, overwrite=overwrite,
					   scliptype=scliptype, dosclip=dosclip, vb=vb, offset=offset, z1=z1, z2=z2,
					   xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, highq=highq, arrow=True, std=std, unique=unique4)

	if scliptype == 0:
		sf = np.sqrt(var)
	else:
		sf = None

	if makepdf:
		imout = lstname.replace('.lst', '.pdf')
		if not os.path.isfile(imout) or overwrite:
			lstIM = [l for l in lst]
			pdfim(lstIM, fcats=fcat, imout=imout, contours=False, zw0=zw0)  # xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax

	zl, yl, xl = f.shape
	zw = '%d:%d' % (zl / 2 - zw0 + 1, zl / 2 + zw0 + 1)
	stdbin = 1  # 6/binsize
	stdbins = stdbin * np.arange(xl / stdbin)

	if propvar:
		proplst = [l.relace('.fits', '.PROPVAR.fits') for l in lst]
		sp, pvar = stack(lst, stackname, imtype='var', sclip=sclip, zw=zw0, makeim=makeim, title=title, vmin=vmin,
						 vmax=vmax, npix=False, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
						 overwrite=overwrite, stds=sf, scliptype=scliptype, dosclip=dosclip, vb=vb, offset=offset,
						 z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, unique=unique4)

	if unique2:
		ulst = []
		rulst = []
		uid = np.unique(id1)
		ucat = np.unique(fcat)
		funique = '%s/unique2/' % foldername
		if not os.path.isdir(funique): os.system('mkdir %s' % funique)
		j = 0
		for uc in ucat:
			for ui in uid:
				set = (id1 != ui) | (fitcats[close] != uc)
				print 'Discarding LAE %d %s' % (ui, uc)
				if len(set) > 0:
					ulst = np.array(lst)[set]
					ulstname = '%s/stack.%d.lst' % (funique, j)
					ustackname = ulstname.replace('.lst', '.fits')
					flst = open(ulstname, 'w')
					for l in ulst: flst.write('%s\n' % l)
					flst.close()
					title = 'Discarded %s-%d N %d' % (ucat[0], ui, len(ulst)) + '\n' + \
							r'$%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $lum>%d$' % \
							(dmin, dmax, pmin, pmax, rmin, rmax, lmin) \
							+ '\n' + extraname

					fu = stack(ulst, ustackname, imtype, sclip, zw0, makeim, title, vmin * 2, vmax * 2, npix=False,
							   var=False, corrmean=meancorr, randflipy=True, overwrite=overwrite, stds=sf,
							   vb=vb, offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
							   scliptype=scliptype, unique=unique4)
					j += 1

	if random:
		full = np.where([len(random_lst[i]) == nstack for i in range(nrand)])[0]
		print "Random lists with full subcubes %d of %d" % (len(full), nrand)
		frs = []

		if doanalysis:
			avr_right = []
			avr_left = []
			avr_top = []
			avr_bottom = []
			avu_right = []
			avu_left = []
			avu_top = []
			avu_bottom = []
			avru_right = []
			avru_left = []
			avru_top = []
			avru_bottom = []
			all = np.array([(stdbins - stdbins[-1] / 2.) * (binsize * .2)])
		h = 'theta'
		if not isdir(foldername + '/randoms'): os.system('mkdir %s/randoms' % foldername)

		for i in range(nrand):
			rlstname = '%s/randoms/random_stack.%d.lst' % (foldername, i)
			rstackname = rlstname.replace('.lst', '.fits')
			title = extraname + '\n' + r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $lum>%d$' % (
				len(random_lst[i]), dmin, dmax, pmin, pmax, rmin, rmax, lmin)
			if len(random_lst[i]) > 0:
				flst = open(rlstname, 'w')
				for l in random_lst[i]: flst.write('%s\n' % l)
				flst.close()

				if doanalysis:
					fr, nfr = stack(random_lst[i], rstackname, imtype, sclip, zw0, makeim, title, vmin, vmax,
									npix=True, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
									overwrite=overwrite, std=std, stds=sf, scliptype=scliptype, dosclip=dosclip, vb=vb,
									offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
									unique=unique4)
					aa, _aa, hh = analysis(fr, nfr, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
					h += ' SB%d' % i
					all = np.concatenate((all, np.array(aa[1:])))
					avr_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]))
					avr_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]))
					avr_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]))
					avr_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]))
				else:
					fr = stack(random_lst[i], rstackname, imtype, sclip, zw0, makeim, title, vmin, vmax,
							   npix=False, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
							   overwrite=overwrite, std=std, stds=sf, scliptype=scliptype, dosclip=dosclip, vb=vb,
							   offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, unique=unique4)

				frs.append(rlstname.replace('.lst',
											'.fits'))  # if len(random_lst[i]) == nstack: frs.append(rlstname.replace('.lst', '.fits'))

				if unique and 0:
					ulst = []
					rulst = []
					funique = '%s/unique/' % foldername
					if not os.path.isdir(funique): os.system('mkdir %s' % funique)
					for uc in ucat:
						for ui in uid:
							set = (id1 == ui) & (fitcats[close] == uc)
							if np.sum(set) > 0:
								lset = np.array(lst)[set]
								choice = np.random.choice(full)
								rlset = np.array(random_lst[choice])[set]
								ulst.append(np.random.choice(lset))
								rulst.append(np.random.choice(rlset))
					ulstname = '%s/stack.u%d.lst' % (funique, i)
					ustackname = ulstname.replace('.lst', '.fits')
					rulstname = '%s/stack.ru%d.lst' % (funique, i)
					rustackname = rulstname.replace('.lst', '.fits')
					flst = open(ulstname, 'w')
					for l in ulst: flst.write('%s\n' % l)
					flst.close()
					title = extraname + '\n' + r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $lum>%d$' % (
						len(ulst), dmin, dmax, pmin, pmax, rmin, rmax, lmin)

					if doanalysis:
						fu, nu = stack(ulst, ustackname, imtype, sclip, zw0, makeim, title, vmin * 2, vmax * 2,
									   npix=True, scliptype=scliptype,
									   var=False, corrmean=meancorr, randflipy=True, overwrite=overwrite, stds=sf,
									   vb=vb, offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
									   unique=unique4)
						rfu, rnu = stack(rulst, rustackname, imtype, sclip, zw0, makeim, title, vmin * 2, vmax * 2,
										 npix=True, scliptype=scliptype,
										 var=False, corrmean=meancorr, randflipy=True, overwrite=overwrite, stds=sf,
										 vb=vb, offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
										 unique=unique4)

						aa, _aa, hh = analysis(fu, nu, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
						h += ' SB%d' % i
						all = np.concatenate((all, np.array(aa[1:])))
						avu_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]))
						avu_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]))
						avu_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]))
						avu_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]))

						aa, _aa, hh = analysis(rfu, rnu, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
						h += ' SB%d' % i
						all = np.concatenate((all, np.array(aa[1:])))
						avru_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]))
						avru_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]))
						avru_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]))
						avru_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]))


					else:
						fu = stack(ulst, ustackname, imtype, sclip, zw0, makeim, title, vmin * 2, vmax * 2,
								   npix=False, scliptype=scliptype,
								   var=False, corrmean=meancorr, randflipy=True, overwrite=overwrite, stds=sf,
								   vb=vb, offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
						rfu = stack(rulst, rustackname, imtype, sclip, zw0, makeim, title, vmin * 2, vmax * 2,
									npix=False, scliptype=scliptype,
									var=False, corrmean=meancorr, randflipy=True, overwrite=overwrite, stds=sf,
									vb=vb, offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

		if ntest > nrand and tests:
			for i in range(nrand, ntest):
				rlstname = '%s/randoms/random_stack.%d.lst' % (foldername, i)
				rstackname = rlstname.replace('.lst', '.fits')
				flst = open(rlstname, 'w')
				lists = [l.replace('.fits', '.%d.fits' %
								   np.random.choice(full)).replace("/cubical", "/random_cubical") for l in lst]
				for l in lists: flst.write('%s\n' % l)
				flst.close()
				title = extraname + '\n' + r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $lum>%d$' % (
					nrand, dmin, dmax, pmin, pmax, rmin, rmax, lmin)

				if doanalysis:
					fr, nfr = stack(lists, rstackname, imtype, sclip, zw0, makeim, title, vmin, vmax,
									npix=True, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
									overwrite=overwrite, std=std, stds=sf, scliptype=scliptype, vb=vb, offset=offset,
									z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, unique=unique4)

					aa, _aa, hh = analysis(fr, nfr, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
					h += ' SB%d' % (i + nrand)
					all = np.concatenate((all, np.array(aa[1:])))

					avr_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]))
					avr_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]))
					avr_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]))
					avr_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]))
				else:
					fr = stack(lists, rstackname, imtype, sclip, zw0, makeim, title, vmin, vmax,
							   npix=False, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
							   overwrite=overwrite, std=std, stds=sf, scliptype=scliptype, vb=vb, offset=offset,
							   z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, unique=unique4)

		superrandom = '%s/random_stack.fits' % (foldername)
		print "Creating super-random file"
		title = ''
		frs = stack(frs, superrandom, imtype, std, zw0, makeim, title, vmin, vmax, highq=highq, contours=False,
					npix=False, var=False, flipshear=fshear, flipy=flipy, overwrite=overwrite,
					std=std, stds=sf, dosclip=False, scliptype=scliptype, vb=vb, offset=offset, z1=z1, z2=z2,
					xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, unique=unique4)

		if doanalysis:
			nfrs = [fff.replace('.fits', '.NPIX.fits') for fff in frs]
			nfrs = stack(nfrs, superrandom.replace('.fits', '.NPIX.fits'), 'flux', std, zw0, makeim, title, vmin=None,
						 vmax=None, npix=False, var=False, corrmean=False, flipy=flipy, overwrite=overwrite,
						 dosclip=False, scliptype=scliptype, vb=vb, offset=offset, z1=z1, z2=z2,
						 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, unique=unique4)

			np.savetxt('%s/random_stack.dat' % foldername, np.matrix(all).T, header=h)

	if tests:
		fts = []
		fjs = []

		if doanalysis:
			avt_right = []
			avt_left = []
			avt_top = []
			avt_bottom = []
			all = np.array([(stdbins - stdbins[-1] / 2.) * (binsize * .2)])
		h = 'theta'
		if not isdir(foldername + '/tests'): os.system('mkdir %s/tests' % foldername)
		if randflipy:
			if not isdir(foldername + '/flipy'): os.system('mkdir %s/flipy' % foldername)

		for i in range(ntest):
			sample = np.random.choice(lst, nstack)
			print "\n\nStacking bootstrapped sample of %d subcubes" % len(sample)
			title = extraname + '\n' + r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $L>%d$' % (
				len(sample), dmin, dmax, pmin, pmax, rmin, rmax, lmin)
			tlstname = '%s/tests/stack.t%d.lst' % (foldername, i)
			tstackname = tlstname.replace('.lst', '.fits')
			flst = open(tlstname, 'w')
			for l in sample: flst.write('%s\n' % l)
			flst.close()

			if doanalysis:
				ftest, nftest = stack(sample, tstackname, imtype, sclip, zw0, makeim, title,
									  vmin, vmax, npix=True, var=False, corrmean=meancorr, randflipy=False,
									  overwrite=overwrite, stds=sf, vb=vb, offset=offset, z1=z1, z2=z2,
									  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, scliptype=scliptype, unique=unique4)

				aa, _aa, hh = analysis(ftest, nftest, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
				h += ' SB%d' % i
				all = np.concatenate((all, np.array(aa[1:])))
				avt_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]))
				avt_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]))
				avt_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]))
				avt_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]))
			else:
				ftest = stack(sample, tstackname, imtype, sclip, zw0, makeim, title,
							  vmin, vmax, npix=False, var=False, corrmean=meancorr, randflipy=False,
							  overwrite=overwrite, stds=sf, vb=vb, offset=offset, z1=z1, z2=z2,
							  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
			fts.append(tlstname.replace('.lst', '.fits'))

			if randflipy:
				ylstname = '%s/flipy/stack.y%d.lst' % (foldername, i)
				ystackname = ylstname.replace('.lst', '.fits')

				fy = stack(lst, ystackname, imtype, sclip, zw0, makeim, title, vmin, vmax, npix=False, var=False,
						   corrmean=meancorr, randflipy=True, overwrite=overwrite, stds=sf, vb=vb,
						   offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
						   scliptype=scliptype, unique=unique4)

		fts = stack(fts, '%s/stack.t.fits' % (foldername), imtype, sclip, zw0, makeim, title, vmin, vmax,
					npix=False, var=False, flipy=flipy, overwrite=overwrite, std=std, dosclip=False, vb=vb,
					offset=offset, z1=z1, z2=z2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, scliptype=scliptype,
					unique=unique4)

		if jackniffe:
			for i in range(nstack):
				sample = np.delete(lst, i)
				print "\n\nStacking Jackniffe sample of %d subcubes" % int(nstack - 1)
				jlstname = '%s/tests/stack.j%d.lst' % (foldername, i)
				jstackname = jlstname.replace('.lst', '.fits')
				fjack, nfjack = stack(sample, jstackname, imtype, sclip, zw0, makeim, title,
									  vmin, vmax, npix=True, var=False, corrmean=meancorr, randflipy=False,
									  overwrite=overwrite, stds=sf, dosclip=dosclip, offset=offset,
									  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, scliptype=scliptype, unique=unique4)
			fjs.append(tlstname.replace('.lst', '.fits'))
			# aa, _aa, hh = analysis(fjack, nfjack, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)

		if 0: fts = stack(fjs, '%s/stack.j.fits' % (foldername), imtype, sclip, zw0, makeim, title, vmin, vmax,
						  npix=False, var=False, flipy=flipy, overwrite=overwrite, std=std, dosclip=False, vb=vb,
						  offset=offset, z1=z1, z2=z2)
		if doanalysis:
			h += '\n'
			np.savetxt('%s/stack.t.dat' % foldername, np.matrix(all).T, header=h)

		if random:
			frt = f - frs
			hdu.data = frt
			hdu.writeto(stackname.replace('.fits', '.randtest.fits'), clobber=True)
			s = 'CubeArit  %s - %s %s %s' % (
				stackname, '%s/random_stack.fits' % (foldername), stackname.replace('.fits', '.randtest.fits'), vb)
			# print s
			# os.system(s)
			s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s %s' % (
				stackname.replace('.fits', '.randtest.fits'), zw, stackname.replace('.fits', '.randtest.IM.fits'),
				vb)
			print s
			os.system(s)
			astroim(stackname.replace('.fits', '.randtest.IM.fits'), smooth=smooth, saveim=makeim,
					show=False,
					cbfrac=.08, pad=.006, dfig=(8, 10), contours=False,
					scb_label=r'Flux [$10^{-20}\,\rm{erg/s/cm^2}$]',
					title='', vmin=vmin, vmax=vmax, gray=True, arrow=True,
					xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, pcolor='black')

	if makeim and os.path.isfile(stackname.replace('.fits', '.randtest.IM.fits')):
		astroim(stackname.replace('.fits', '.randtest.IM.fits'), smooth=smooth, saveim=makeim,
				show=False, highq=highq,
				cbfrac=.08, pad=.006, dfig=(8, 10), contours=False,
				scb_label=r'Flux [$10^{-20}\,\rm{erg/s/cm^2}$]',
				title='', vmin=vmin, vmax=vmax, gray=True,
				xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

	if not random: frs = None
	if tests and doanalysis:
		all, _all, h = analysis(f, nf, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0, frs)
		avf_right = [np.nanmean(all[1][(all[0] <= 6) & (all[0] >= 2)])]
		avf_left = [np.nanmean(all[1][(all[0] <= -2) & (all[0] >= -6)])]
		avf_top = [np.nanmean(_all[1][(_all[0] <= 6) & (_all[0] >= 2)])]
		avf_bottom = [np.nanmean(_all[1][(_all[0] <= -2) & (_all[0] >= -6)])]
		np.savetxt('%s/stack_n%d.dat' % (foldername, nstack), np.matrix(all).T, header=h)
		np.savetxt('%s/stack_n%d_vertical.dat' % (foldername, nstack), np.matrix(_all).T, header=h)

		if random:
			avp_right = []
			avp_left = []
			avp_top = []
			avp_bottom = []
			for l in lst:
				fp = getdata(l)
				if isfile(l.replace('.fits', '.NPIX.fits')):
					nfp = getdata(l.replace('.fits', '.NPIX.fits'))
				else:
					nfp = fp
				aa, _aa, hh = analysis(fp, nfp, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
				avp_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]))
				avp_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]))
				avp_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]))
				avp_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]))
				del fp, nfp
			all = [avp_right, avp_top, avp_left, avp_bottom]
			np.savetxt('%s/subcubes.avs.dat' % foldername, np.matrix(all).T,
					   header='subcubes_r subcubes_t subcubes_l subcubes_b')

			if unique:
				all = [avf_right * ntest, avr_right, avt_right, avru_right, avu_right,
					   avf_top * ntest, avr_top, avt_top, avru_top, avu_top,
					   avf_left * ntest, avr_left, avt_left, avru_left, avu_left,
					   avf_bottom * ntest, avr_bottom, avt_bottom, avru_bottom, avu_bottom]
				np.savetxt('%s/stack.avs.dat' % foldername, np.matrix(all).T,
						   header='full_r rand_r sets_r urand_r usets_r full_t rand_t sets_t urand_t usets_t ' +
								  'full_l rand_l sets_l urand_l usets_l full_b rand_b sets_b urand_b usets_b')

			else:
				all = [avf_right * ntest, avr_right, avt_right, avf_top * ntest, avr_top, avt_top,
					   avf_left * ntest, avr_left, avt_left, avf_bottom * ntest, avr_bottom, avt_bottom]
				np.savetxt('%s/stack.avs.dat' % foldername, np.matrix(all).T,
						   header='full_r random_r sets_r full_t random_t sets_t full_l random_l sets_l full_b random_b sets_b')

	f_snr = f * np.sqrt(1. / var)
	hdu.data = f_snr
	hdu.writeto(stackname.replace('.fits', '.SNR.fits'), clobber=True)
	s = 'Cube2Im -cube %s[*,*,%s] -varcube %s -snrmap %s -imtype flux %s' \
		% (stackname, zw, stackname.replace('.fits', '.VAR.fits'), stackname.replace('.fits', '.SNR.IM.fits'), vb)

	print s
	os.system(s)

	astroim(stackname.replace('.fits', '.SNR.IM.fits'), smooth=smooth, saveim=makeim, show=False,
			cbfrac=.08, gray=True, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
			pad=.006, dfig=(8, 10), contours=contours, scb_label='SNR', title='', nsigma=6, vmin=-1,
			vmax=15, sb=False)

	if tests:
		field = fcat == 'UDF'
		pairs_data = [field, id1, id2, red, dist, theta, npair]
		np.savetxt('%s/stack_pairs.dat' % foldername, np.matrix(pairs_data).T,
				   header="Field id1 id2 redshift comoving_distance theta n_neigbors",
				   fmt='%d %d %d %1.3f %1.3f %1.3f %d')
