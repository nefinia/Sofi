#!/usr/bin/env python
__author__ = 'gallegos'
import glob
import h5py
# import eagleSqlTools
import matplotlib.pyplot as plt
import os
import scipy.interpolate as inter
from matplotlib.colors import LinearSegmentedColormap
from pyfits import getdata, PrimaryHDU
from sys import argv

import numpy as np

from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
snaps = rdarg(argv, 'snap', list, [10, 11, 12, 13, 14], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
halfplot = rdarg(argv, 'halfplot', bool, False)
histover = rdarg(argv, 'histover', bool, False)
LLScubes = rdarg(argv, 'LLScubes', bool, False)
sql = rdarg(argv, 'sql', bool, False)
sphericalLLS = rdarg(argv, 'sphericalLLS', bool, False)
pairLLS = rdarg(argv, 'pairLLS', bool, False)
circlehist = rdarg(argv, 'circlehist', bool, False)
h2d = rdarg(argv, 'h2d', bool, False)
kde = rdarg(argv, 'kde', bool, False)
sbhist = rdarg(argv, 'sbhist', bool, False)
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
radprof = rdarg(argv, 'radprof', bool, False)
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)

cdict1 = {'red': ((0.0, 0.0, 0.0),
				  (0.5, 0.0, 0.1),
				  (1.0, 1.0, 1.0)),

		  'green': ((0.0, 0.0, 0.0),
					(1.0, 0.0, 0.0)),

		  'blue': ((0.0, 0.0, 1.0),
				   (0.5, 0.1, 0.0),
				   (1.0, 0.0, 0.0))
		  }
colors = [(0, 0, 1), (1, 1, 0), (1, 0, 0)]
cm = LinearSegmentedColormap.from_list('sofi', colors, N=50)

snames = {'10': '010_z003p984', '11': '011_z003p528', '12': '012_z003p017', '13': '013_z002p478', '14': '014_z002p237',
		  '15': '015_z002p012'}
reds = {'10': 984, '11': 528, '12': 17}
redshifts = {'10': 3.984, '11': 3.528, '12': 3.017, '13': 2.478, '14': 2.237, '15': 2.012}
asec2kpcs = {'10': 7.842, '11': 7.449, '12': 7.108, '13': 8.241, '14': 8.396, '15': 8.516}
gamma_bkg = {'10': [5.71E-13, 3.28E-13, 2.31E-16], '11': [6.77E-13, 3.89E-13, 8.70E-16],
			 '12': [7.66E-13, 4.42E-13, 2.10E-15], '13': [9.50E-13, 5.55E-13, 9.06E-15],
			 '14': [9.64E-13, 5.67E-13, 1.13E-14]}
heat_bkg = {'10': [2.27E-12, 2.18E-12, 7.24E-15],
			'11': [2.68E-12, 2.62E-12, 2.33E-14],
			'12': [3.02E-12, 3.05E-12, 5.01E-14],
			'13': [3.75E-12, 4.22E-12, 1.78E-13],
			'14': [3.81E-12, 4.42E-12, 2.18E-13]}
dz = {'10': 0.0255, '11': .035, '12': .03}
types = ['NHI', 'NHII']

# sbpeaks = [2.396]
# sbpeaks = {'10': 1.011, '11': 1.484, '12':2.396} # for an aperture of 1 asec^2!!! I am converting to flux later on
sbpeaks = {'10': .739, '11': 1.293, '12': 2.579}  # for an aperture of 1 asec^2!!! I am converting to flux later on
zlens = {'10': 34, '11': 29, '12': 25, '13': 20, '14': 18, '15': 16}
lcube = 4096
coml = 25  # cMpc
com2pix = 163.84  # lcube/coml
kpc2pix = lcube / float(coml * 1e3)
rads = np.array([0, 2, 4, 8, 12, 20, 30, 50, 100, 200])

lognhi = [13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21]
nhi = np.power(10, lognhi)
sb = [0.0001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96, 0.98, 1, 1, 1]

nhi2sb = inter.interp1d(nhi, sb)
sb2nhi = inter.interp1d(sb, nhi)

if temperature:
	res, lmax = 512, 0
	f0 = '/net/galaxy-data/export/galaxydata/saeed/EAGLE/RefL0025N0752/'
	colors = {'10': 'red', '11': 'green', '12': 'blue'}
	for snap in snaps:
		data = getdata('../../EAGLE/cats/gals_snap%d.fits' % snap, 1)
		ids = data['ID']
		if 0:
			xc, yc, zc, u, ids = data['x'], data['y'], data['z'], data['U'], data['ID']
			ngal = len(xc)
			cube = getdata('%s/snapshot_%s/fits_3d/snap_%s_%s_%s_8_var_uniform-temperature_lmax_%d.fits' %
						   (f0, snames[str(snap)], snames[str(snap)], res, res, lmax))
			zl, yl, xl = cube.shape
			conv = zl / float(coml)
			z, y, x = np.ogrid[1: zl + 1, 1: yl + 1, 1: xl + 1]
			T = []  # 20, 40 and 80 kpc
			for i in range(ngal):
				print 'id', ids[i], '%d/%d' % (i, ngal)
				cool = ((.05 * conv) < (z - zc[i] * conv) ** 2 + (y - yc[i] * conv) ** 2 + (x - xc[i] * conv) ** 2) & \
					   ((z - zc[i] * conv) ** 2 + (y - yc[i] * conv) ** 2 + (x - xc[i] * conv) ** 2 < (.08 * conv))
				_t = np.nanmean(cube[cool])
				print 'Temperature between 50 and 80 kpc %f, U %f, number of pixels used %d' % (_t, u[i], np.sum(cool))
				T.append(_t)
			np.savetxt('../../EAGLE/cats/temp_snap%d' % snap, np.array([ids, T]).T, header='ID Temperature')
		if 1:
			i1, i2, fconns, frands = np.loadtxt('../../EAGLE/analysis/avsb_snap%d_x.dat' % snap).T
			fconn = []
			for i in ids:
				cool = i1 == i
				fconn.append(np.nanmean(frands[cool]))

			np.savetxt('../../EAGLE/cats/fconn_snap%d' % snap, np.array([ids, fconn]).T, header='ID fconn')

		else:
			T = data['Temperature']
			U = data['U']
			plt.scatter(U, T, label='z %.2f' % redshifts[str(snap)], color=colors[str(snap)])
	if 0:
		plt.semilogy()
		plt.xlabel(r'U')
		plt.ylabel(r'T[K]')
		plt.legend(scatterpoints=1)
		plt.show()

if LLScubes:
	cubefolder = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/'
	saeed = '/net/galaxy-data/export/galaxydata/saeed/EAGLE/RefL0025N0752/'
	init_evol = False
	col_dens = True
	halos_mask = False
	do_cubes = False
	do_inter = False

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

	if do_inter:
		plt.close()
		nhi2sb2 = inter.InterpolatedUnivariateSpline(nhi, sb)
		plt.plot(lognhi, sb, 'bo', label='Data')
		x = np.arange(13, 21, .01)
		plt.plot(x, nhi2sb2(10 ** x), 'k--', label='Spline fit')
		plt.plot(x, nhi2sb(10 ** x), 'k--', label='1d fit', color='red')
		plt.semilogy()
		plt.legend()
		plt.xlim(15, 21.5)
		plt.hlines(y=1, xmin=15, xmax=21.5)
		plt.ylim(.01, 2.1)
		plt.savefig('interp.png')
		plt.close()

		y = np.arange(.0001, 1, .001)
		plt.plot(sb, nhi, 'bo', label='Data')
		plt.plot(y, sb2nhi(y), 'k--', label='1d fit')
		plt.semilogy()
		plt.legend()
		plt.savefig('interp2.png')
		plt.close()

	for snap in snaps:
		s = str(snap)
		print 'snapshot', snap
		red = redshifts[s]
		# sbpeak = sbpeaks[snap]
		asec2kpc = asec2kpcs[s]
		zlen = zlens[s]
		asec2pix = asec2kpc * (1 + red) * kpc2pix
		sb2flux = asec2pix ** -2.
		rad = 6 * asec2pix
		cat = getdata('../../EAGLE/cats/gals_snap%d.fits' % snap, 1)
		sgamma = '%.2E %.2E %.2E' % tuple(gamma_bkg[s])
		sheat = '%.2E %.2E %.2E' % tuple(heat_bkg[s])
		sname = snames[s]

		for type in ['NHI']:

			fname = '%s/snapshot_%s/' % (saeed, sname)
			if init_evol:
				glob.os.chdir(fname+'init_evol')
				srun = './init_evol.sh 512 4096 8 25. %.3f "%s" "%s"' % (red, sgamma, sheat)
				print srun
				os.system(srun)
				glob.os.chdir(scodes)
			if col_dens:
				glob.os.chdir(fname + 'column_density')
				finit = '../init_evol/SO.snap_%s_512_4096_8_%s_%s.ts0000' % \
						(sname, sgamma.replace(' ', '_'), sheat.replace(' ', '_'))
				if os.path.isfile(finit):
					lmax = 3
					nl = zlens[s]  # number of layers
					nc = 10  # number of cores
					srun = './column_density_layers.sh %s %s %d %d %d' % (finit, type, lmax, nl, nc)
					print srun
					os.system(srun)
				else: print '%s not there?' % finit
				glob.os.chdir(scodes)

			if halos_mask:
				outname = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/snap%s_%s.%s.fits' % (
				snap, coord, type)
				if not os.path.isfile(fname) or overwrite:
					if os.path.isfile(outname) and not overwrite:
						cube = getdata(outname)
					else:
						cube = np.zeros((zlen, lcube, lcube))

						if mask:
							mask = np.zeros((zlen, lcube, lcube))
							c = ['x', 'y', 'z']
							c.remove(coord)
							xc, yc = [cat[c[0]] * com2pix, cat[c[1]] * com2pix]
							zc = cat[coord] * zlen / coml

							for i in range(len(xc)):
								print i
								halos = ((x - xc[i]) ** 2 + (y - yc[i]) ** 2 < rad ** 2) & (abs(z - zc[i]) < zw0)
								_NHI[halos] = 0

						for i in range(zlen):
							print 'Coord', coord, 'snap', snap, 'layer', i + 1
							cubename = 'snap%s_%s_%d.%s.fits' % (sname, coord, i, type)
							data = getdata('%s/%s' % (fname, cubename))
							cube[i, :, :] = data
							hdu.data = cube
							hdu.writeto(outname, clobber=True)
		if do_cubes:
			nl = zlens[s]
			#check lower or upper case numbers!!!
			sgamma = '%.2E_%.2E_%.2E' % tuple(gamma_bkg[s])
			sheat = '%.2E_%.2E_%.2E' % tuple(heat_bkg[s])
			sname = snames[s]
			_s = saeed + '/snapshot_%s/column_density/layers/' % sname

			for coord, proj in zip(coordnames, [1, 2, 3]):
				if 0:
					cubes = []
					for i in range(nl):
							flayer = '%sSO.snap_%s_512_4096_8_%s_%s.ts0000_var_NHI_proj_%d_lmax_3_l_%d_%d.fits' % \
									 (_s, sname, sgamma, sheat, proj, i+1, nl)
							cubes.append(getdata(flayer))
					print 'done loading layers'
					cubes = np.array(cubes)
					hdu.data = cubes
					hdu.writeto('%s/snap%d_%s.NHI.fits' % (cubefolder, snap, coord), clobber=True)
					print 'done NHI cube'
					lls = np.copy(cubes)
					_lls = cubes > 10**17.2 #from prochaska
					lls[_lls] = 1
					lls[~_lls] = 0
					hdu.data = lls
					hdu.writeto('%s/snap%d_%s.LLS.fits' % (cubefolder, snap, coord), clobber=True)
					print 'done LLS cube'
				cubes = getdata('%s/snap%d_%s.NHI.fits' % (cubefolder, snap, coord))
				cubes[cubes < nhi[0]] = nhi[0]
				cubes[cubes > nhi[-1]] = nhi[-1]
				sbcube = nhi2sb(cubes)
				sbcube[sbcube == sb[0]] == 0
				hdu.data = sbcube
				hdu.writeto('%s/snap%d_%s.SB.fits' % (cubefolder, snap, coord), clobber=True)
				print 'done SB cube'

if sql:
	conn = eagleSqlTools.connect('sgallego', 'dpTT852J')

	for snap in snaps:
		sql = "SELECT \
				gal.CentreOfMass_x as x,\
				gal.CentreOfMass_y as y,\
				gal.CentreOfMass_z as z,\
				gal.Velocity_x as vx,\
				gal.Velocity_y as vy,\
				gal.Velocity_z as vz,\
				gal.Redshift as redshift,\
				gal.MassType_Star as stellar_mass,\
				gal.StarFormationRate as SFR,\
				gal.MassType_DM as DM_mass,\
				gal.MassType_Gas gass_mass,\
				sizes.R_halfmass30 as size,\
				gal.SubGroupNumber as SubGroupNumber,\
				mag.u_nodust as U,\
				mag.g_nodust as G,\
				mag.r_nodust as R,\
				mag.i_nodust as I,\
				mag.z_nodust as Z,\
				mag.Y_nodust as Y,\
				mag.J_nodust as J,\
				mag.H_nodust as H,\
				mag.K_nodust as K,\
				flux.Johnson_V as f_V,\
				flux.Johnson_R as f_R,\
				fof.Group_R_Mean200 as Rmean200,\
				fof.Group_M_Mean200 as Mmean200\
			  FROM\
				RefL0025N0752_SubHalo as gal,\
				RefL0025N0752_Magnitudes as mag,\
				RefL0025N0752_Aperture as ape,\
				RefL0025N0752_Sizes as sizes,\
				RefL0025N0752_FOF as fof,\
				RefL0025N0752_DustyFluxes as flux\
			  WHERE\
				gal.SnapNum = %d and\
				ape.Mass_Star > 1.0e8 and\
				ape.ApertureSize = 30 and\
				gal.GalaxyID = mag.GalaxyID and\
				gal.GalaxyID = ape.GalaxyID and\
				gal.GalaxyID = SIZES.GalaxyID and\
				gal.GalaxyID = flux.GalaxyID and\
				gal.GroupID = fof.GroupID" % snap

		data = eagleSqlTools.execute_query(conn, sql)
		np.savetxt('../../EAGLE/cats/gals_%d.dat' % snap, data,
				   header='x y z vx vy vz redshift stellar_mass SFR DM_mass gass_mass size SubGroupNumber '
						  'U G R I Z Y J H K f_V f_R Rmean200 Mmean200')

if sphericalLLS:
	LLSfracs_xyz = []
	LLScums_xyz = []

	cn = coordnames[0]
	c = coords[0]
	fit = fits[0]
	# for cn, c, fit in zip(coordnames, coords, fits):
	for snap in [10, 11, 12]:
		print 'Collapsed axis:', cn
		llsfrac_name = '../../EAGLE/LLS/LLSfracs_snap%d_%s.dat' % (snap, cn)
		llscum_name = '../../EAGLE/LLS/LLScums_snap%d_%s.dat' % (snap, cn)
		isfiles = os.path.isfile(llsfrac_name) & os.path.isfile(llscum_name)

		if not isfiles:
			yl, xl = fit.shapenethz
			_y, _x = np.ogrid[1: yl + 1, 1: xl + 1]
			LLSfracs = []
			LLScums = []

			for i in (ids - 1):
				print 'Galaxy %d of %d' % (i + 1, ngal)
				x0 = c[0][i]
				y0 = c[1][i]
				# cool = ids == i
				nrads = len(rads)
				LLSfrac = []
				LLScum = []

				for j in range(nrads - 1):
					cool = ((_x - x0 - .5) ** 2 + (_y - y0 - .5) ** 2 < rads[j + 1] ** 2)
					LLScum.append(np.nanmean(fit[cool]))
					cool &= ((_x - x0 - .5) ** 2 + (_y - y0 - .5) ** 2 >= rads[j] ** 2)
					LLSfrac.append(np.nanmean(fit[cool]))

				if 0:
					plt.plot(rad_asec, LLSfrac)
					plt.xlabel('distance[arcsec]')
					plt.ylabel('LLS fraction')
					plt.savefig('../../EAGLE/LLS/plots/lls_frac%d.png' % (i + 1))
					# plt.show()
					plt.close()

				LLSfracs.append(LLSfrac)
				LLScums.append(LLScum)

			LLSfracs = np.array(LLSfracs).T
			np.savetxt(llsfrac_name, LLSfracs)
			LLSfracs_xyz = np.concatenate(LLSfracs_xyz, LLSfracs)
			LLScums = np.array(LLScums).T
			np.savetxt(llscum_name, LLScums)
			LLScums_xyz = np.concatenate(LLScums_xyz, LLScums)

		if isfiles: LLSfracs = np.loadtxt(llsfrac_name)
		LLS1 = [np.median(l) for l in LLSfracs]
		LLSstd1 = [np.std(l) for l in LLSfracs]
		plt.plot(rad_asec, LLS1, lw=2)
		plt.ylim(0, 1)
		plt.xlabel('distance[arcsec]')
		plt.ylabel('LLS fraction')
		plt.errorbar(rad_asec, LLS1, yerr=LLSstd1)
		plt.savefig('../../EAGLE/LLS/plots/LLSfracs_snap%d_%s.png' % (snap, cn))
		# plt.show()
		plt.close()
		plt.plot(rad_asec, LLS1, lw=2)
		plt.ylim(0, 1)
		plt.xlim(0, 8)
		plt.xlabel('distance[arcsec]')
		plt.ylabel('LLS fraction')
		plt.errorbar(rad_asec, LLS1, yerr=LLSstd1)
		plt.savefig('../../EAGLE/LLS/plots/LLSfracs-close_snap%d_%s.png' % (snap, cn))
		# plt.show()
		plt.close()

		if isfiles: LLScums = np.loadtxt(llscum_name)
		LLS = [np.median(l) for l in LLScums]
		LLSstd = [np.std(l) for l in LLScums]
		plt.plot(rad_asec, LLS, lw=2)
		plt.semilogy()
		plt.ylim([0.015, 1])
		plt.yticks([.03, 0.1, .3, 1], ('0.03', '0.1', '0.3', '1'))
		plt.xlabel('distance[arcsec]')
		plt.ylabel('f_LLS(r<R)')
		plt.errorbar(rad_asec, LLS, yerr=LLSstd)
		plt.savefig('../../EAGLE/LLS/plots/LLScums_%s_%s.png' % (snap, cn))
		# plt.show()
		plt.close()

		LLS = [np.median(l) for l in LLScums]
		LLSstd = [np.std(l) for l in LLScums]
		x = 0.042 * np.array(rad_asec) ** 2
		dndz = np.array(LLS) * x
		dndzstd = np.array(LLSstd) * x
		plt.plot(rad_asec, dndz, lw=2)
		# plt.ylim([0.1, 1])
		plt.semilogy()
		plt.xlabel('distance[arcsec]')
		plt.ylabel('dn/dz(LLS,r<R)')
		plt.errorbar(rad_asec, dndz, yerr=dndzstd)
		plt.savefig('../../EAGLE/LLS/plots/lls_dndz_%s.png' % cn)
		# plt.show()
		plt.close()

if pairLLS:
	for snap in [10, 11, 12]:
		for coord in ['x', 'y', 'z']:
			pairs = getdata('../../EAGLE/cats/lae_pairs_snap%d.fits' % snap, 1)
			id1 = pairs['id1']
			id2 = pairs['id2']
			x1 = pairs['x1']
			y1 = pairs['y1']
			z1 = pairs['z1']
			x2 = pairs['x2']
			y2 = pairs['y2']
			z2 = pairs['z2']
			nw1 = pairs['nw1']
			nw2 = pairs['nw2']
			theta = pairs['theta%s' % coord]
			dist = pairs['com_dist']
			xs = np.concatenate([np.arange(0, 24, 2), [28, 36, 54, 82, 150]]).astype('int')
			nsteps = len(xs) - 1
			# Approx 1 arcsec height
			ymin = -2.5
			ymax = 2.5
			fracs = []
			fracs_av = []
			f1 = '../../EAGLE/LLS/all_fracs_snap%d_%s.dat' % (snap, coord)
			f2 = '../../EAGLE/LLS/all_fracs_av_snap%d_%s.dat' % (snap, coord)
			pairsname = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/LLScats/cubical_pair'
			good_pairs = [os.path.isfile(pairsname + '%d-%d.h5' % (i1, i2)) for i1, i2 in zip(id1, id2)]
			i1good = id1[good_pairs]
			i2good = id2[good_pairs]
			overwrite = False

			if not os.path.isfile(f1) or not os.path.isfile(f2) or overwrite:
				for i1, i2 in zip(i1good, i2good):
					print 'Pair %d-%d' % (i1, i2)
					cat = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/cats/logNHI_snap%d_%scubical_pair%d-%d.h5' \
						  % (snap, coord, i1, i2)
					data = h5py.File(cat, 'r')
					x = np.array(data['px'])
					y = np.array(data['py'])
					f = np.array(data['flux'])
					frac = []
					frac_av = []
					for i in range(nsteps):
						oriented = (ymin <= y) & (y <= ymax) & (xs[i] <= x) & (x < xs[i + 1])
						average = ((ymin <= x) & (x <= ymax) & (xs[i] <= y) & (y < xs[i + 1])) | \
								  ((ymin <= y) & (y <= ymax) & (-xs[i] > x) & (x >= -xs[i + 1])) | \
								  ((ymin <= x) & (x <= ymax) & (-xs[i] > y) & (y >= -xs[i + 1]))
						frac.append(np.nanmean(f[oriented]))
						frac_av.append(np.nanmean(f[average]))
					fracs.append(frac)
					fracs_av.append(frac_av)
				np.savetxt(f1, fracs)
				np.savetxt(f2, fracs_av)
				fracs = np.array(fracs)
				fracs_av = np.array(fracs_av)
			else:
				fracs = np.loadtxt(f1)
				fracs_av = np.loadtxt(f2)

			good_theta = theta > 20
			# distbins = [1, 3, 5, 10, 15, 20]
			mubins = [[-22, -19], [-19, -16.5]]  # [-17, -19, -20, -21, -23]
			mrbins = [[-23, -20], [-20, -17]]  # [-17, -19, -20, -21, -23]
			nwbins = [[0, 15], [10, 22], [18, 22], [15, 22]]
			rads = xs[1:] * pix2deg * 3600

			means_av = []
			stds_av = []
			for xx in range(nsteps):
				means_av.append(np.nanmean(fracs_av[:, xx]))
				stds_av.append(np.nanstd(fracs_av[:, xx]))

			dmax = 5
			umax = -19
			ugood = ids[mu <= umax]
			others = [(ii1 in ugood) & (ii2 in ugood) for ii1, ii2 in zip(id1, id2)]
			others &= (dist <= dmax)
			_s = 'dmax%d_umax%d' % (dmax, umax)
			# for i in range(len(magbins)-1):
			#    idgood = ids[(magbins[i]>=mu) & (mu>=magbins[i+1])]

			props = [mr, mu, nw]
			bins = [mrbins, mubins, nwbins]
			names = ['mr', 'mu', 'nw']

			for j in range(len(bins)):
				for i in range(len(bins[j])):
					idgood = ids[(bins[j][i][0] <= props[j]) & (props[j] <= bins[j][i][1])]
					cool = []
					for ii1, ii2 in zip(id1, id2):
						cool.append((ii1 in idgood) & (ii2 in idgood))
					cool &= good_theta & others
					cool = cool[good_pairs]
					frac_cool = fracs[cool]
					means = []
					stds = []
					for xx in range(nsteps):
						means.append(np.nanmean(frac_cool[:, xx]))
						stds.append(np.nanstd(frac_cool[:, xx]))

					plt.plot(rads, means, label='Oriented to neighbour', color='blue')
					plt.errorbar(rads, means, yerr=stds, color='blue')
					plt.plot(rads, means_av, label='Random orientations', color='gray')
					# plt.errorbar(rads, means_av, yerr=stds_av, color='gray')
					plt.xlabel('distance [arcsec]')
					plt.ylabel('LLS fraction')
					plt.legend(loc='best')
					# plt.semilogy()
					plt.ylim(0, 0.3)
					plt.xlim(2, 15.5)
					plt.savefig('../../EAGLE/LLS/plots/LLSfrac_pairs+circular_%s_%s%d_%d.png' % (
						_s, names[j], bins[j][i][0], bins[j][i][1]))
					plt.close()

if sbhist:

	def sb_roi(cat, xmin, xmax, ymin, ymax, zmin, zmax):
		data = h5py.File(cat, 'r')
		pz = np.array(data['pz'])
		zgood = (pz >= zmin) & (pz <= zmax)
		f = np.array(data['SB'])[zgood]
		px = np.array(data['px'])[zgood]
		py = np.array(data['py'])[zgood]
		right = (px > xmin) & (px < xmax) & (py > ymin) & (py < ymax)
		left = (px < -xmin) & (px > -xmax) & (py > ymin) & (py < ymax)
		down = (px > ymin) & (px < ymax) & (py < -xmin) & (py > -xmax)
		up = (px > ymin) & (px < ymax) & (py > xmin) & (py < xmax)
		rand = left | down | up
		fright = np.nanmean(f[right])
		frand = np.nanmean(f[rand])
		return fright, frand


	folder = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/cats/'
	props = [['u1', 'd'], ['nw1', 'd']]
	xmin
	for snap in snaps:
		for coord in coordnames:
			catfolder = folder + 'SB_snap%d_%s/' % (snap, coord)
			glob.os.chdir(catfolder)
			cats = glob.glob('*.h5')
			cats.sort()
			id1, id2 = np.array([np.array(c[:-3].split('-')).astype(int) for c in cats]).T

if h2d:

	if 0:
		freexmax = False
		if freexmax:
			extraname = '.freexmax'  # ''_3-7asec'
		else:
			extraname = ''  # ''.d<6'
		from itertools import product

		zr = []  # np.arange(-2,3)

		# trying just u1 and d
		pnameshist = [['u1',
					   'd']]  # [['nw1', 'd'], ['nw1', 'd'], ['nw1', 'u1']]#, ['u2', 'd'], ['nw1', 'd']]  # [pnames[0], pnames[2]]
		# pnameshist = [['u1', 'u2'], ['u1', 'nw1']]#, ['nw1', 'u2']]  # [pnames[0], pnames[2]]

		thetamin = 20
		h_all = []
		n_all = []
		hsb_all = []
		superall = {}
		ssss = 0

		pn0, pn1 = pnameshist[0]
		txt0 = '../../EAGLE/analysis/histograms/hist_XXX_snapYY_%s-%s%s.dat' % (pn0, pn1, extraname)

		for snap in snaps:
			houts = {}
			txt = txt0.replace('YY', str(snap))
			houts['fright'] = txt.replace('XXX', 'fright')
			houts['frand'] = txt.replace('XXX', 'frand')
			houts['nsubs'] = txt.replace('XXX', 'nsubs')

			hbool = True
			for ho in houts: hbool &= os.path.isfile(houts[ho])

			if hbool:  # not hbool:
				for snap in snaps:
					sbpeak = sbpeaks[str(snap)]
					red = redshifts[str(snap)]
					asec2kpc = asec2kpcs[str(snap)]

					sb_dimming = np.power(1 + red, -4)
					all = {}
					thetamin = 16
					# height 2 arcsec
					# width 6 to 12 arcsec -> same as my paper!
					xasecmin = 6
					xasecmax = 12
					yasecmin = -1
					yasecmax = 1

					asec2pix = asec2kpc * (1 + red) * kpc2pix
					sb2flux = asec2pix ** -2.
					# limits in pixels
					xmin = xasecmin * asec2pix
					if not freexmax: xmax = xasecmax * asec2pix
					ymin = yasecmin * asec2pix
					ymax = yasecmax * asec2pix
					zmax = 2
					zmin = -zmax
					rad_asec = np.array([rads[i + 1] for i in range(len(rads) - 1)]) / asec2pix

					fsize = 20
					keys = ['right', 'left', 'top', 'down']

					# Use observable quantities if possible
					# dz = x1-x2
					pairs = getdata('../../EAGLE/cats/lae_pairs_snap%d.fits' % snap, 1)

					pnames = ['u1', 'u2', 'd', 'm1', 'm2', 'nw1']

					params = {'u1': pairs['U1'], 'u2': pairs['U2'], 'd': pairs['com_dist'],
							  'm1': np.log10(pairs['stellar_mass1']), 'm2': np.log10(pairs['stellar_mass2']),
							  'nw1': pairs['nw1']}

					umin = min(params['u1'])
					umax = max(params['u1'])

					urange = [-21.405001, -19.551001, -19.080999, -18.549, -18.066, -17.507999]
					# urange = [-21.405001, -20.125999, -19.551001, -19.306999, -19.080999, -18.822001, -18.549, -18.291,
					#          -18.066, -17.811001, -17.507999]
					pranges = {'u1': urange, 'u2': urange,
							   'd': np.arange(0, 21, 2), 'nw1': np.arange(0.08, .17, .02),
							   'm1': np.arange(10.5, 12.1, .17), 'm2': np.arange(10.5, 12.1, .17)}
					plabels = {'u1': u'$\mathrm{U_1}$', 'u2': u'$\mathrm{U_2}$', 'd': u'$\mathrm{d\,[cMpc]}$',
							   'm1': u'$\mathrm{M_{stars,\,1}}$', 'm2': u'$\mathrm{M_{stars,\,2}}$',
							   'nw1': u'$\mathrm{Number\,of\,neighbors}$'}

					props1 = [[[], []]] * 18
					props2 = [[]] * 4
					pnames = ['id', 'x', 'y', 'z', 'nw', 'U', 'G', 'R', 'I', 'Z', 'Y', 'J', 'H', 'K',
							  'gas_mass', 'DM_mass', 'size', 'stellar_mass']

					for i in range(len(props1)):
						props1[i] = [pairs['%s1' % pnames[i]], pairs['%s2' % pnames[i]]]

					id, x, y, z, nw, u, g, r, i, z, y, j, h, k, mgas, mdm, size, mstars = props1

					for coord in coordnames:

						pairsname = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/cats/SB_snap%d_%s/' \
									% (snap, coord)
						outname2 = '../../EAGLE/analysis/avsb_snap%d_%s%s.dat' % (snap, coord, extraname)

						pnames2 = ['theta%s' % coord, 'com_dist', 'ang%s' % coord, 'shear%s' % coord]

						for i in range(len(props2)): props2[i] = pairs[pnames2[i]]
						theta, dist, angle, shear = props2

						# props = [props1, props2]

						cool = []
						for i1, i2 in zip(id[0], id[1]):
							cool.append(os.path.isfile(pairsname + '%d-%d.h5' % (i1, i2)))
						cool &= (theta > thetamin) & (dist < 20) & (dist > .5)  # (dist<6)
						ncool = np.sum(cool)
						cool = np.where(cool)[0]

						if freexmax:
							xmax = (theta[cool] - 6) * asec2pix
						else:
							xmax = np.zeros(ncool) + xasecmax * asec2pix

						if 1:  # ncool>0:
							condition = (not os.path.isfile(outname2)) or overwrite
							if condition:

								fout = open(outname2, 'w')
								fout.write('#id1 id2 f_right f_rand\n')
								foutz = {}
								for _z in zr:
									foutz[str(_z)] = open(outname2.replace('.dat', '.%d.dat' % _z), 'w')
									foutz[str(_z)].write('#id1 id2 f_right f_rand\n')

								for i in range(ncool):
									print 'Pair %d-%d, snap%d %s. %d of %d' % (
										id[0][cool[i]], id[1][cool[i]], snap, coord, i, ncool)
									cat = '%s%d-%d.h5' % (pairsname, id[0][cool[i]], id[1][cool[i]])
									data = h5py.File(cat, 'r')
									pz = np.array(data['pz'])
									zgood = (pz >= zmin) & (pz <= zmax)
									f = np.array(data['SB'])[zgood]
									px = np.array(data['px'])[zgood]
									py = np.array(data['py'])[zgood]
									right = (px > xmin) & (px < xmax[i]) & (py > ymin) & (py < ymax)
									left = (px < -xmin) & (px > -xmax[i]) & (py > ymin) & (py < ymax)
									down = (px > ymin) & (px < ymax) & (py < -xmin) & (py > -xmax[i])
									up = (px > ymin) & (px < ymax) & (py > xmin) & (py < xmax[i])
									rand = left | down | up
									f_ = np.nanmean(f[right])
									# f_ = np.nanmedian(f[right])
									# f_n = np.nanmean(right)
									f_rand = np.nanmean(f[rand])
									# f_rand = np.nanmedian(f[rand])
									# f_rand_n = np.nanmean(rand)
									fout.write('%d %d %f %f\n' % (id[0][cool[i]], id[1][cool[i]], f_, f_rand))
									for _z in zr:
										# remember f[zgood].....
										rightz = right & (pz == _z)
										randz = rand & (pz == _z)
										f_ = np.nanmean(f[rightz])
										# f_n = np.nanmean(rightz)
										f_rand = np.nanmean(f[randz])
										# f_rand_n = np.nanmean(randz)
										foutz[str(_z)].write(
											'%d %d %f %f\n' % (id[0][cool[i]], id[1][cool[i]], f_, f_rand))

									print 'sb: right %1.3f rand %1.3f' % \
										  (f_, f_rand)
								fout.close()
								for _z in zr: foutz[str(_z)].close()

							else:
								print "Output file(s) already exist(s)", outname2

							_f = np.loadtxt(outname2)
							_fz = {}
							for _z in zr: _fz[str(_z)] = np.loadtxt(outname2.replace('.dat', '.%d.dat' % _z))

							for p in pnameshist:
								pn0, pn1 = p

								if pn0 != pn1:
									r0, r1 = [pranges[pn0], pranges[pn1]]
									n0, n1 = [len(r0) - 1, len(r1) - 1]
									id1 = list(_f[:, 0])
									id2 = list(_f[:, 1])
									cool2 = []
									nid0 = len(id[0])
									nid1 = len(id1)
									for _j in range(len(id[0])):
										end = False
										iii = 0
										while not end:
											if iii >= len(id1):
												end = True
											else:
												if (id1[iii] == id[0][_j]) and (id2[iii] == id[1][_j]):
													id1.remove(id1[iii])
													id2.remove(id2[iii])
													cool2.append(_j)
													end = True
											iii += 1
									_p0, _p1 = [params[pn0][cool2], params[pn1][cool2]]
									f_right = _f[:, 2]
									# f_right_n = _f[:, 3]
									f_rand = _f[:, 3]
									# f_rand_n = _f[:, 5]
									print 'Combining %s with %s' % (pn0, pn1)

									outs = ['../../EAGLE/analysis/histograms/hist_fright_snap%d_%s_%s-%s%s.dat' % (
										snap, coord, pn0, pn1, extraname),
											'../../EAGLE/analysis/histograms/hist_frand_snap%d_%s_%s-%s%s.dat' % (
												snap, coord, pn0, pn1, extraname),
											# '../../EAGLE/analysis/histograms/hist_n_snap%d_%s_%s-%s%s.dat' % (snap, coord, pn0, pn1, extraname),
											'../../EAGLE/analysis/histograms/hist_nsubs_snap%d_%s_%s-%s%s.dat' % (
												snap, coord, pn0, pn1, extraname)]

									isfiles = True
									for out in outs: isfiles &= os.path.isfile(out)

									if isfiles and not overwrite and not histover:
										print 'h2d files exist'
										all[coord, pn0, pn1] = [np.loadtxt(out) for out in outs]
									else:
										if histover and isfiles: 'Overwriting h2d files!'
										if not isfiles: print 'No h2d files!'
										# n = np.zeros((n0, n1))
										nsubs = np.zeros((n0, n1))
										hfright = np.zeros((n0, n1))
										hfrand = np.zeros((n0, n1))

										for i in range(n0):
											for j in range(n1):
												coolh = (_p0 > r0[i]) & (_p0 <= r0[i + 1]) & (_p1 > r1[j]) & (
													_p1 <= r1[j + 1])
												# n[i, j] = np.nansum(f_right_n[coolh])
												nsubs[i, j] = np.nansum(coolh)  # /n[i, j]
												hfright[i, j] = np.nanmean(f_right[coolh])
												hfrand[i, j] = np.nanmean(f_rand[coolh])  # / np.nansum(f_rand_n[coolh])

										all[coord, pn0, pn1] = [hfright, hfrand, nsubs]

										for out, _h in zip(outs, all[coord, pn0, pn1]):
											print 'Saving file', out
											np.savetxt(out, _h)

									if 0:
										def doplot(h, name, labelname, vmin=None, vmax=None):
											# from matplotlib.colors import LogNorm
											if pn0 == 'd': plt.figure(figsize=(7, 5))
											if pn1 == 'd': plt.figure(figsize=(6, 7))
											plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin,
													   vmax=vmax)
											dx = (r0[1] - r0[0])
											dy = (r1[1] - r1[0])
											plt.xticks(2 * np.arange(n0 / 2) + .5,
													   np.arange(r0[0] + dx, r0[-1], 2 * dx))
											plt.yticks(2 * np.arange(n1 / 2) + .5,
													   np.arange(r1[0] + dy, r1[-1], 2 * dy))
											plt.xlabel(plabels[pn0], fontsize=fsize)
											plt.ylabel(plabels[pn1], fontsize=fsize)
											cbar = plt.colorbar()
											cbar.set_label(labelname, size=fsize)
											plt.savefig(
												'../../EAGLE/analysis//plots/coords/hist_%s_snap%d_%s_%s-%s%s.png'
												% (name, snap, coord, pn0, pn1, extraname))
											plt.close()


										# Number of neighbors per bin
										doplot(n, 'N', u'$\mathrm{n}$')  # , vmin=0, vmax=50)

										# sb
										doplot(hfright, 'SB', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=.4, vmax=.7)

										# sb rand
										doplot(hfrand, 'SB_rand', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=.4, vmax=.7)

										# sb - sb rand
										doplot((hfright - hfrand) * 100 / hfrand, 'SBminusSBrand',
											   u'$\mathrm{SB\,\,increase\,[percentage]}$', vmin=0, vmax=20)

										# sb wrt sbpeak
										doplot(hfright / sbpeak, 'SBvsUVB',
											   u'$\mathrm{SB/SB_{UVB}}$', vmin=0, vmax=.5)

						else:
							print 'No subcubes for %s coordinate' % coord

					for p in pnameshist:  # list(product(pnames, pnames)):
						superall['%d' % snap] = all
						ssss += 1
						print 'ssss', ssss


						def doplot(h, name, labelname, vmin=None, vmax=None, title=''):
							# from matplotlib.colors import LogNorm
							if pn0 == 'd': plt.figure(figsize=(7, 5))
							if pn1 == 'd': plt.figure(figsize=(5.5, 6.21))
							plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin, vmax=vmax,
									   cmap=cm)
							dx = (r0[1] - r0[0])
							dy = (r1[1] - r1[0])
							plt.title(title)
							plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
							plt.yticks([-.5, 1.5, 3.5, 5.5, 7.5, 9.5], [0, 4, 8, 12, 16, 20])
							# plt.xticks(2 * np.arange(n0 / 2) + .5, np.arange(r0[0] + dx, r0[-1], 2 * dx))
							# plt.yticks(2 * np.arange(n1 / 2) + .5, np.arange(r1[0] + dy, r1[-1], 2 * dy))
							plt.xlabel(plabels[pn0], fontsize=fsize)
							plt.ylabel(plabels[pn1], fontsize=fsize)
							cbar = plt.colorbar()
							cbar.set_label(labelname, size=fsize)
							plt.savefig(
								'../../EAGLE/analysis/plots/hist_%s_snap%d_%s-%s%s.png'
								% (name, snap, pn0, pn1, extraname))
							plt.close()


						pn0, pn1 = p
						print 'Combining %s with %s for all coords and snap %d' % (pn0, pn1, snap)

						r0, r1 = [pranges[pn0], pranges[pn1]]
						n0, n1 = [len(r0), len(r1)]
						outs = ['../../EAGLE/analysis/histograms/hist_fright_snap%d_%s-%s%s.dat' % (
							snap, pn0, pn1, extraname),
								'../../EAGLE/analysis/histograms/hist_frand_snap%d_%s-%s%s.dat' % (
									snap, pn0, pn1, extraname),
								# '../../EAGLE/analysis/histograms/hist_n_snap%d_%s-%s%s.dat' % (snap, pn0, pn1, extraname),
								'../../EAGLE/analysis/histograms/hist_nsubs_snap%d_%s-%s%s.dat' % (
									snap, pn0, pn1, extraname)]
						isfiles = True
						for out in outs: isfiles &= os.path.isfile(out)

						if isfiles and not overwrite and not histover:
							print 'h2d files exist'
							_all = [np.loadtxt(out) for out in outs]
						else:
							print 'No h2d files!'
							_all = []
							for coord in coordnames: _all.append(all[coord, p[0], p[1]])
							_all = np.nanmean(_all, 0)
							for out, _h in zip(outs, _all):
								print 'Saving file', out
								np.savetxt(out, _h)
							hfright, hfrand, nsubs = _all

							# from matplotlib.colors import LogNorm

							# Number of pixels per bin
							# doplot(np.log10(n), 'N', u'$\mathrm{log(pixels)}$')#, vmin=0, vmax=40)
							vmin, vmax = [.04, .13]  # [None, None]#
							doplot(hfright, 'SBvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=vmin, vmax=vmax,
								   title=u'snapshot %d, z=%.3f' % (snap, red))

							# sb rand
							doplot(hfrand, 'SBrandvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=vmin, vmax=vmax,
								   title=u'snapshot %d, z=%.3f' % (snap, red))

						if 0:
							# Number of neighbors per bin
							doplot(nsubs, 'Nsub', u'$\mathrm{Number\,\,of\,\,orientations}$')  # , vmin=0, vmax=40)
							# sb
							doplot(hfright * sbpeak, 'SB', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=0, vmax=.14)

							doplot(hfright, 'SBvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=0, vmax=.14)

							# sb rand
							doplot(hfrand * sbpeak, 'SB_rand', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=0, vmax=.14)

							# sb - sb rand
							doplot((hfright - hfrand) * 100 / hfrand, 'SBminusSBrand',
								   u'$\mathrm{SB\,\,increase\,[percentage]}$', vmin=-30, vmax=80)

							# SNR
							doplot(hfright * np.sqrt(nsubs), 'SNR', u'$\mathrm{Relative\,\,SNR}$', vmin=0.2, vmax=1.5)

			else:

				houts2 = ['../../EAGLE/analysis/histograms/hist_fright_%s-%s%s.dat' % (pn0, pn1, extraname),
						  '../../EAGLE/analysis/histograms/hist_frand_%s-%s%s.dat' % (pn0, pn1, extraname),
						  '../../EAGLE/analysis/histograms/hist_nsubs_%s-%s%s.dat' % (pn0, pn1, extraname)]
				hbool = True
				for ho in houts2: hbool &= os.path.isfile(ho)

				hright, hrand, nsubs = [[], [], []]

				if hbool:
					print 'aa'
					hists = [hright, hrand, nsubs]
					for ho, h in zip(houts2, hists): h = np.loadtxt(ho)
				else:
					for ho in houts:
						h = np.loadtxt(ho)


	def doplot(h, name, labelname, p0='', p1='', r0=None, r1=None, vmin=None, vmax=None, title='',
			   xlabel='HST_F606W', ylabel='d [cMpc]', cticks=None):

		# from matplotlib.colors import LogNorm
		# if pn0 == 'd': plt.figure(figsize=(7, 5))
		plt.figure(figsize=(5, 6.2))  # if pn1 == 'd':
		plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin, vmax=vmax,
				   cmap=cm)
		if p0 == 'u1':
			plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
		else:
			plt.xticks([0, 1.5, 3, 4.5], [26, 28, 30, 32])
		# plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
		plt.yticks([-.5, 1.5, 3.5, 5.5, 7.5, 9.5], [0, 4, 8, 12, 16, 20])

		if r0 is not None:
			dx = (r0[1] - r0[0])
			plt.xticks(2 * np.arange(n0 / 2) + .5, np.arange(r0[0] + dx, r0[-1], 2 * dx))
		if r1 is not None:
			dy = (r1[1] - r1[0])
			plt.yticks(2 * np.arange(n1 / 2) + .5, np.arange(r1[0] + dy, r1[-1], 2 * dy))
		plt.xlabel(xlabel, fontsize=fsize)
		plt.ylabel(ylabel, fontsize=fsize)
		plt.title(title)
		cbar = plt.colorbar()
		cbar.set_label(labelname, size=fsize)
		if cticks is not None: cbar.set_ticks(cticks)
		plt.savefig(
			'../../EAGLE/analysis/plots/hist_%s_%s-%s%s.png'
			% (name, p0, p1, extraname))
		plt.close()


	# recheck all of this part!!

	if 0:
		print 'Final analysis'
		uvb = None
		huvb = None
		uvbrand = None
		iii = 0
		for k in superall.keys():
			sa = superall[k]
			for kk in sa.keys():
				iii += 1
				_all = sa[kk]
				if huvb is None:
					# n = _all[2]
					nsub = _all[2]
					hfright = _all[0] * sbpeaks[k]  # *nsub
					hfrand = _all[1] * sbpeaks[k]  # *nsub
					huvb = _all[0]  # *nsub
					uvb = np.nanmean(_all[0])
					uvbrand = np.nanmean(_all[1])
					nuvb = _all[2]
					print iii, k, kk, 'oriented %.3f random %.3f' % (uvb, uvbrand)

				else:
					# n += _all[2]
					_nsub = _all[2]
					nsub += _nsub
					_hright = _all[0]  # *_all[2]
					hfright += _hright * sbpeaks[k]
					_hfrand = _all[1]  # *_all[2]
					hfrand += _hfrand * sbpeaks[k]
					huvb += _all[0]
					_uvb = np.nanmean(_hright)
					_uvbrand = np.nanmean(_hfrand)
					uvb += _uvb
					uvbrand += _uvbrand
					nuvb += _nsub
					print iii, k, kk, 'oriented %.4f random %.4f' % (_uvb, _uvbrand)
		ntot = 9.  # np.sum(n)
		# uvb = uvb/ntot/sbpeaks['11']
		# uvbrand = uvbrand/ntot/sbpeaks['11']

		hfright /= ntot
		hfrand /= ntot
		huvb /= ntot

	else:
		urange = [-21.405001, -19.551001, -19.080999, -18.549, -18.066, -17.507999]
		drange = np.arange(0, 21, 2)
		n0 = len(urange)
		n1 = len(drange)
		fsize = 14
		extraname = ''
		hright = {}
		hrand = {}
		huvb = {}
		for snap in snaps:
			ss = str(snap)
			vmin, vmax = [.04, .16]
			_f = getdata('../../EAGLE/cats/lae_pairs_snap%d.fits' % snap, 1)
			fright = np.concatenate([_f['f_right_x'], _f['f_right_y'], _f['f_right_z']])
			fright[np.isnan(fright)] = 0
			frand = np.concatenate([_f['f_rand_x'], _f['f_rand_y'], _f['f_rand_z']])
			frand[np.isnan(frand)] = 0
			u1 = np.concatenate([_f['U1'], _f['U1'], _f['U1']])
			dist = np.concatenate([_f['com_dist'], _f['com_dist'], _f['com_dist']])
			n = np.histogram2d(u1, dist, [urange, drange], weights=fright > 0)[0]
			h = np.histogram2d(u1, dist, [urange, drange], weights=fright)[0]
			huvb[ss] = h / n
			hright[ss] = huvb[ss] * sbpeaks[ss]
			n = np.histogram2d(u1, dist, [urange, drange], weights=frand > 0)[0]
			h = np.histogram2d(u1, dist, [urange, drange], weights=frand)[0]
			hrand[ss] = h * sbpeaks[ss] / n
			doplot(huvb[ss], 'SBvsUVB_snap%d' % snap, u'$\mathrm{SB/SB_{UVB}}$', 'u1', 'd',  # urange, drange,
				   xlabel=u'$\mathrm{U}}$', ylabel=u'$\mathrm{d\,[cMpc]}}$',
				   title='snapshot %d, z=%.3f' % (snap, redshifts[ss]), vmin=vmin, vmax=vmax)
			doplot(hright[ss], 'SB_snap%d' % snap, u'$\mathrm{SB\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', 'u1', 'd',
				   # urange, drange,
				   xlabel=u'$\mathrm{U}}$', ylabel=u'$\mathrm{d\,[cMpc]}}$',
				   title='snapshot %d, z=%.3f' % (snap, redshifts[ss]), vmin=.04, vmax=.26,
				   cticks=[.05, .1, .15, .2, .25])
		hfright = np.nanmean([hright['10'], hright['11'], hright['12']], 0)
		hfrand = np.nanmean([hrand['10'], hrand['11'], hrand['12']], 0)
		huvb = np.nanmean([huvb['10'], huvb['11'], huvb['12']], 0)

	# print 'total SB/SB_UVB: oriented %.4f random %.4f' % (uvb, uvbrand)

	noise1 = 1.  # assuming 1e-20 e.g. UDF -> but 0.44 with the combination of 390 subcubes
	noise2 = noise1 * np.sqrt(3)  # assuming 1e-20 e.g. UDF

	udf = getdata('../../UDF/cats/lae_pairs.fits', 1)
	mosaic = getdata('../../mosaic/cats/lae_pairs-2018.fits', 1)

	d1 = udf['com_dist']
	t1 = udf['theta']
	f1 = udf['HST_F606W']
	d2 = mosaic['com_dist']
	f2 = mosaic['HST_F606W']
	t2 = mosaic['theta']

	g1 = (d1 > .5) & (d1 < 20) & (t1 > 16) & (f1 > 0)
	g2 = (d2 > .5) & (d2 < 20) & (t2 > 16) & (f2 > 0)

	dbins = np.arange(0, 21, 2)
	# fbins = [25.5921, 27.3024, 27.9615, 28.5004, 28.9428, 29.343,
	#         29.724,  30.0098,  30.5226,  31.3483,  32.046]
	fbins = [25.5921, 27.9615, 28.9428, 29.724, 30.5226, 32.046]

	p0 = 'u1'
	p1 = 'd'
	# n0 = len(pranges[p0])
	# n1 = len(pranges[p1])
	r0 = fbins
	r1 = dbins

	h1 = np.histogram2d(f1, d1, [fbins, dbins])[0]
	h2 = np.histogram2d(f2, d2, [fbins, dbins])[0]
	doplot(h1, 'Nsub_UDF', u'$\mathrm{Number\,\,of\,\,orientations}$', 'HST_F606W', p1, r0, r1)  # , vmin=.2, vmax=1.5)
	doplot(h2, 'Nsub_mosaic', u'$\mathrm{Number\,\,of\,\,orientations}$', 'HST_F606W', p1, r0,
		   r1)  # , vmin=.2, vmax=1.5)#, vmin=.2, vmax=.5)

	SNR1 = hfright * np.sqrt(h1) / noise1
	SNRrand1 = hfrand * np.sqrt(h1) / noise1
	SNR2 = hfright * np.sqrt(h2) / noise2
	SNRrand2 = hfrand * np.sqrt(h2) / noise2
	vmin, vmax = [None, None]  # [.3, 1.6]
	doplot(huvb, 'SBvsUVBtot', u'$\mathrm{SB/SB_{UVB}}$', 'u1', p1, vmin=.06, vmax=.14,
		   xlabel=r'$\mathrm{U}$')  # r0, r1,
	doplot(SNR1, 'SNR_UDF', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)
	doplot(SNRrand1, 'SNRrand_UDF', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin,
		   vmax=vmax)  # , vmin=.2, vmax=.5)
	doplot((SNR1 - SNRrand1) * 100 / SNRrand1, 'OrientedvsRandom_UDF', u'$\mathrm{SB\,\,increase\,[percentage]}$',
		   'HST_F606W', p1, r0, r1, vmin=0, vmax=40)  # , vmin=.2, vmax=.5)
	doplot(SNR2, 'SNR_mosaic', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)  # , vmin=.2, vmax=.5)
	doplot(SNRrand2, 'SNRrand_mosaic', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin,
		   vmax=vmax)  # , vmin=.2, vmax=.5)
	doplot((SNR2 - SNRrand2) * 100 / SNRrand2, 'OrientedvsRandom_mosaic', u'$\mathrm{SB\,\,increase\,[percentage]}$',
		   'HST_F606W', p1, r0, r1, vmin=0, vmax=40)  # , vmin=.2, vmax=.5)

	nh1 = np.nansum(h1)
	fconn1 = np.nansum((huvb * h1).reshape(-1)) / nh1
	SNR_tot1 = np.nansum((SNR1 * h1).reshape(-1)) / nh1 / noise1
	SNRrand_tot1 = np.nansum((SNRrand1 * h1).reshape(-1)) / nh1 / noise1

	nh2 = np.nansum(h2)
	fconn2 = np.nansum((huvb * h2).reshape(-1)) / nh2
	SNR_tot2 = np.nansum((SNR2 * h2).reshape(-1)) / nh2 / noise2
	SNRrand_tot2 = np.nansum((SNRrand2 * h2).reshape(-1)) / nh2 / noise2

	print 'Total SNR UDF %.3f rand %.3f' % (SNR_tot1, SNRrand_tot1)
	print 'f_conn UDF %.3f' % (fconn1)
	print 'Total SNR MOSAIC %.3f rand %.3f' % (SNR_tot2, SNRrand_tot2)
	print 'f_conn MOSAIC %.3f' % (fconn2)

if halfplot:
	print 'halfplot'

if circlehist:
	ax = plt.subplot(111, polar=True)
	max_rad = 2
	bottom = 1
	delta = .5
	# ax.set_rmax(max_rad * 2)
	ax.set_yticklabels([])
	pairs = getdata('../../EAGLE/cats/lae_pairs.fits', 1)
	id1 = pairs['id1']
	id2 = pairs['id2']
	x1 = pairs['x1']
	y1 = pairs['y1']
	z1 = pairs['z1']
	x2 = pairs['x2']
	y2 = pairs['y2']
	z2 = pairs['z2']
	nw1 = pairs['nw1']
	nw2 = pairs['nw2']
	theta = pairs['theta']
	dist = pairs['com_dist']
	umax = -19
	idgood = ids[mu <= umax]
	nwmin = 15
	cool = []
	pairsname = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/LLScats/cubical_pair'
	for i1, i2 in zip(id1, id2):
		cool.append((i1 in idgood) & (i2 in idgood) & (os.path.isfile(pairsname + '%d-%d.h5' % (i1, i2))))
	cool &= (id1 != id2) & (theta > 16) & (dist < 5) & (nw1 > nwmin) & (nw2 > nwmin)
	print 'Cools ids', np.sum(cool), 'of', len(id1)
	LLSrads = []
	angbins = []
	nrad = 0
	nrads = len(rads)
	nbins = []
	angbins = []
	delta_angs = []
	nangs = []
	llsangs = []
	rads = [0, 2, 6, 12, 20, 30, 60]
	rad2 = [.5, 1, 2, 2, 1, 1, 1]
	rad3 = [0, .5, 1.5, 3.5, 5.5, 6.5, 7.5]
	rad4 = [1, 6, 12, 20, 30, 60]
	nrads = len(rads)
	for j in range(nrads - 1):
		delta_ang = np.pi / float(rad4[j])
		angbin = np.arange(0., 2 * np.pi + 1e-10, 2 * delta_ang)
		nbin = len(angbin) - 1
		nbins.append(nbin)
		angbins.append(angbin)
		LLSrads.append(np.zeros(nbin))
		delta_angs.append(delta_ang)(px > xmin) & (px < xmax) & (py > ymin) & (py < ymax)

	filedat = '../../EAGLE/LLS/llscumang.dat'
	overwrite = True

	if not os.path.isfile(filedat) or overwrite:
		x = np.array([])
		y = np.array([])
		f = np.array([])
		for i in np.where(cool)[0]:
			print 'Pair %d-%d' % (id1[i], id2[i])
			cat = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/LLScats/cubical_pair%d-%d.h5' % (id1[i], id2[i])
			data = h5py.File(cat, 'r')
			x = np.concatenate((x, np.array(data['px'])))
			y = np.concatenate((y, np.array(data['py'])))
			f = np.concatenate((f, np.array(data['flux'])))

		for j in range(nrads - 1):
			d2 = x ** 2 + y ** 2
			ang = (np.arctan2(y, x) + np.pi - delta_angs[j]) % (2 * np.pi)

			for a in range(nbins[j]):
				inside = (ang > angbins[j][a]) & (ang <= angbins[j][a + 1]) & (d2 < rads[j + 1] ** 2) & (
					d2 < rads[j] ** 2)
				LLSrads[j][a] = np.nanmean(f[inside])

		fout = open(filedat, 'w')
		for line in LLSrads:
			for row in line:
				if np.isnan(row):
					fout.write('0. ')
				else:
					fout.write('%f ' % float(row))
			fout.write('\n')
		fout.close()

	else:
		fin = open(filedat, 'r')
		j = 0
		for line in fin:
			a = 0
			for row in line.split():
				LLSrads[j][a] = float(row)
				a += 1
			j += 1
		fin.close()

	for j in range(nrads - 1):
		if j == 0:
			ec = "none"
		else:
			ec = "black"
		bars = ax.bar(angbins[j][:-1] + delta_angs[j], np.zeros(len(angbins[j]) - 1) + rad2[j], width=2 * delta_angs[j],
					  bottom=rad3[j], edgecolor=ec, linewidth=.5)
		mlr = max(LLSrads[j])
		for bar, lr in zip(bars, LLSrads[j]): bar.set_facecolor(plt.cm.jet(1.6 * lr))

	plt.savefig('../../EAGLE/LLS/plots/circhist.png')
	plt.close()

if lutzmodel:

	for snap in [10, 11, 12]:
		data = getdata('../../EAGLE/cats/gals_snap%s.fits' % coord, 1)
		comx = data['x']
		comy = data['y']
		comz = data['z']
		x = comx * com2pix
		y = comy * com2pix
		z = comz * com2pix
		pos = [[y, z], [x, z], [x, y]]

		for coord, xp, yp in zip(['x', 'y', 'z'], ):

			outfile = '../../EAGLE/LLS/LLSfracs_%s_%s.dat' % (snap, coord)
			overwrite = False
			if os.path.isfile(outfile) and not overwrite:
				mat = np.loadtxt(outfile).T
				rad_asec = mat[0] / 5.
				llsin = mat[1]
				stdin = mat[2]
				llsout = mat[3]
				stdout = mat[4]
			else:
				fits = getdata(
					'/net/galaxy-data/export/galaxydata/gallegos/EAGLE/logNHI_snap%d_%s.fits' % (snap, coord))
				inside = False
				xpos = y
				ypos = z
				l = 4096
				# 1 spaxel is (0.2 arcsec)^2
				rads = [5, 10, 20, 40, 60, 80, 100, 120, 160, 200]
				rad_asec = np.array(rads) / 5
				_y, _x = np.ogrid[0: l, 0: l]
				llsin = []
				llsout = []
				stdin = []
				stdout = []
				for r in rads:
					print 'Inside %d arcsec' % (r / 5)
					for i, j, gal in zip(xpos, ypos, ids):
						# print 'Gal %d x %d y %d' % (gal, i, j)
						inside |= (_x - i) ** 2 + (_y - j) ** 2 < r ** 2

					_llsin = np.nanmean(fits[inside])
					_llsout = np.nanmean(fits[~inside])
					_stdin = np.nanstd(fits[inside])
					_stdout = np.nanstd(fits[~inside])
					print 'f_LLS inside %f, outside %f' % (_llsin, _llsout)
					llsin.append(_llsin)
					llsout.append(_llsout)
					stdin.append(_stdin)
					stdout.append(_stdout)

				mat = np.array([rads, llsin, stdin, llsout, stdout]).T
				np.savetxt(outfile, mat, header='distance_pix LLS_in std_LLS_in LLS_out std_LLS_out')

			plt.plot(rad_asec, llsin, lw=2, color='blue', label='Inside')
			plt.plot(rad_asec, llsout, lw=2, color='red', label='Outside')
			plt.ylim(0, .4)
			# plt.semilogy()
			plt.xlabel(u'$\mathrm{distance\,[arcsec]}$', fontsize=18)
			plt.ylabel(u'$f_{\mathrm{LLS}}(\mathrm{r<R})}$', fontsize=18)
			plt.errorbar(rad_asec, llsin, yerr=stdin, color='blue')
			plt.errorbar(rad_asec, llsout, yerr=stdout, color='red')
			plt.savefig('../../EAGLE/LLS/plots/LLSfracs_snap%d_%s.png' % (snap, coord))
			# plt.show()
			plt.close()

if 1:
	zw = 1
	feagle = '/net/galaxy-data/export/galaxydata/gallegos/EAGLE/'
	rads = np.array([[0, 2], [2, 4], [4, 6], [6, 8], [8, 12], [12, 16],
			[16, 20], [20, 24], [24, 30], [30, 40], [40, 60], [60, 80], [80, 100], [100, 150], [150, 200]])
	feagles = []
	l = 601
	hdu = PrimaryHDU()
	yo, xo = np.ogrid[0:lcube, 0:lcube]
	for s in snaps:
		asec2kpc, red = [asec2kpcs[str(s)], redshifts[str(s)]]
		data = getdata('../../EAGLE/cats/gals_snap%d.fits' % s, 1)
		asec2pix = asec2kpc*(1+red)*kpc2pix
		radpix = rads*asec2pix
		for c in coordnames:
			cube = getdata('%s/snap%d_%s.LLS.fits' % (feagle, s, c))
			zl, yl, xl = cube.shape
			ids = data['ID']
			cs = ['x', 'y', 'z']
			cs.remove(c)
			xs = np.round(data[cs[0]] * xl / coml).astype(int)
			ys = np.round(data[cs[1]] * yl / coml).astype(int)
			zs = np.round(data[c] * zl / coml).astype(int)
			for i, x, y, z in zip(ids, xs, ys, zs):
				print i, x, y, z
				pos = (xo - x) ** 2 + (yo - y) ** 2
				frad = []
				for r in radpix:
					cool = (pos > (r[0])**2) & \
						   (pos < (r[1])**2)
					fcool = np.nanmean(cube[:, cool], 1)
					frad.append(np.roll(fcool, zl/2-z))
				hdu.data = np.array(frad)
				hdu.writeto('%s/gals/snap%d_%s/%d.RADPROF.fits' % (feagle, s, c, i), clobber=True)

if snr:
	def sclipping(fits, nsigma, dim=None, mask=None):
		# print 'Asign nan to values > nsigma in a fits array'
		if dim is None:
			stds = np.nanstd(fits[:, mask])
		else:
			stds = np.nanstd(fits[:, mask], dim)
		high_sigma = np.abs(fits) > nsigma * stds
		fits[high_sigma] = np.nan
		return fits, stds


	folder = '/net/abnoba/scratch2/gallegos/Research/MUSE/'
	flaes = '/net/galaxy-data/export/galaxydata/gallegos/'
	fstack = folder + '/all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100/'

	ncores = 10
	zw = 2
	SB = 1.27  # SB at z=3.5
	type = 'LLS'  # 'NHI'

	# oriented
	if 0:

		feagles = []
		for s in [10, 11, 12]:
			for c in ['x', 'y', 'z']:
				feagles.append(getdata(
					folder + '/EAGLE/stacks/snap%d_%s_d0-20_pd16-2000_d5th0.5-20.0_u0.0-99.0/stack.fits' % (s, c)))
		feagle = np.nanmean(feagles, 0)
		fits = getdata(fstack + 'stack.fits')
		zl, xl, yl = fits.shape
		frands = []
		for i in range(200):
			frands.append(
				np.nansum(getdata(fstack + 'randoms/random_stack.%d.fits' % i)[zl / 2 - zw:zl / 2 + zw + 1, :, :], 0))
		# 1 spaxel is (0.4 arcsec)^2
		# 7.8125 comes from 1.25 Angstrom wavelength width times 1/(0.4)**2 from conversion to the correct aperture
		conv = 7.8125
		fits = np.nansum(fits[zl / 2 - zw:zl / 2 + zw + 1, :, :], 0)  # stack is in flux!
		zle, xle, yle = feagle.shape
		feagle = np.nanmean(feagle[zle / 2 - zw:zle / 2 + zw + 1, :, :], 0)
	# end oriented

	# 1 spaxel is (0.2 arcsec)^2
	# conv comes from 1.25 Angstrom wavelength width times 1/(0.2)**2 from conversion to the correct aperture
	conv = 31.25

	# laes
	do_laes=False
	if do_laes:
		fall = []
		fstds = []
		hdu = PrimaryHDU()
		for s in snaps:
			cat = getdata('../../EAGLE/cats/gals_snap%d.fits' % s, 1)
			pairs = getdata('../../EAGLE/cats/lae_pairs_snap%d.fits' % s, 1)
			ids = cat['ID']
			us = cat['U']
			mstar = cat['stellar_mass']
			x, y, z = [cat['x'], cat['y'], cat['z']]
			n5d = []
			#fit = np.zeros([lcube, lcube, lcube])
			for i, _x, _y, _z in zip(ids, x, y, z):
				d = np.sqrt(((x - _x + 12.5) % 25 - 12.5) ** 2 + ((y - _y + 12.5) % 25 - 12.5) ** 2 + (
						(z - _z + 12.5) % 25 - 12.5) ** 2)
				dsort = np.argsort(d)
				d.sort()
				n5d.append(d[5])
				cool = pairs['id1'] == i
				id2, x2, y2, z2, cd2 = np.array([pairs['id2'][cool], pairs['x2'][cool], pairs['y2'][cool],
												 pairs['z2'][cool], pairs['com_dist'][cool]])
				#cls = cd2 < 2
				#for cl in cls:
			#		fit[x2[cl], y2[cl], z2[cl]] += id2[cl]

			# n5d = cat['n5d']
			med_u = np.nanmedian(us)
			med_n5d = np.nanmedian(n5d)
			med_mstar = np.nanmedian(mstar)

			feagles = []
			feagle_med = {'n5d_high': [], 'n5d_low': [],  'u_high': [], 'u_low': []}
			try:
				glob.os.makedirs('../../EAGLE/simplestacks/')
			except:
				pass
			if 1:

				for c in coordnames:
					for i, u, m in zip(ids, us, n5d):
						ffits = flaes + '/EAGLE/gals/snap%d_%s/%d.fits' % (s, c, i)
						print ffits
						if os.path.isfile(ffits):
							fit = getdata(ffits)
							feagles.append(fit)
							if u < med_u: feagle_med['u_low'].append(fit)
							if u > med_u: feagle_med['u_high'].append(fit)
							if m < med_n5d: feagle_med['n5d_low'].append(fit)
							if m > med_n5d: feagle_med['n5d_high'].append(fit)

			if 1:
				feagle = np.nanmean(feagles, 0)
				fstd = np.nanstd(feagles, 0)
				fall.append(feagle)
				fstds.append(fstd)

				hdu.data = feagle
				hdu.writeto('../../EAGLE/simplestacks/snap%d.%s.fits' % (s, type), clobber=True)
				hdu.data = fstd
				hdu.writeto('../../EAGLE/simplestacks/snap%d.STD.%s.fits' % (s, type), clobber=True)

				for k in feagle_med.keys():
					feagle_med[k] = np.nanmean(feagle_med[k], 0)
					hdu.data = feagle_med[k]
					hdu.writeto('../../EAGLE/simplestacks/snap%d_%s.%s.fits' % (s, k, type), clobber=True)
		fcomb = np.nanmean(fall, 0)
		hdu.data = fcomb
		hdu.writeto('../../EAGLE/simplestacks/snaps10-11-12.%s.fits' % type, clobber=True)

	if 0:

		fdat = 'muse-vs-eagle.dat'
		if os.path.isfile(fdat):
			_r0, _r1, _fin, _fe, _fstd, _fe_std, _fmean_out, _fstd_out = np.loadtxt(fdat)
		else:
			_fin = None
		if overwrite: _fin = None
		# _fin = None
		feagle = [getdata('../../EAGLE/simplestacks/snap%d.%s.fits' % (s, type)) for s in [10, 11, 12]]
		feagle = np.nanmean(feagle, 0)
		zl, yl, xl = feagle.shape
		feagle_out = (feagle[0, :, :] + feagle[-1, :, :]) / 2.
		feagle = feagle[zl / 2, :, :]

		# feagle = np.nansum(feagle[zl / 2 - zw:zl / 2 + zw + 1, :, :], 0)
		# fnhi = sb2nhi(feagle[zl / 2 - zwe:zl / 2 + zwe + 1, :, :])
		# fnhi[fnhi<15] = 0
		# fnhi = np.power(10, fnhi)
		# fnhisum = np.log10(np.nansum(fnhi, 0))
		# feagle = nhi2sb(np.log10(np.nansum(fnhi, 0)))


		asec2pix_eagle = .1836  # average btw snapshots
		asec2pix_muse = .2  # muse
		fall = []
		fouts = []
		n5d = []
		try:
			glob.os.makedirs('../../all/simplestacks/')
		except:
			pass

		if _fin is None:
			for fitcat in ['UDF', 'HDFS']:
				cat = getdata('../../%s/cats/laes.fits' % fitcat, 1)
				pairs = getdata('../../%s/cats/lae_pairs.fits' % fitcat, 1)
				ids = cat['ID']
				zs = cat['redshift']
				sconf = cat['sconf']
				cool = (sconf >= 2) & (zs < 4)
				d = pairs['com_dist']
				id1 = pairs['id1']
				for i in ids[cool]:
					ffits = flaes + '/%s/LAEs/%d.fits' % (fitcat, i)
					if os.path.isfile(ffits):
						fit = getdata(ffits)
						zlm, ylm, xlm = fit.shape
						fitin = fit[zlm / 2 - zw:zlm / 2 + zw + 1, :, :]
						fitout = fit[:2 * zw + 1, :, :]
						fall.append(np.nanmean(fitin, 0) * (2 * zw + 1))
						fouts.append(np.nanmean(fitout, 0))
						n5d.append(np.sort(d[id1 == i])[5])
			fall = np.array(fall)
			fouts = np.array(fouts)
			fall = sclipping(fall, 3, 0)[0]
			fouts = sclipping(fouts, 3, 0)[0]
			fmuse = np.nanmean(fall, 0)
			fout = np.nanmean(fouts, 0)
		# hdu = PrimaryHDU()
		# hdu.data = fmuse
		# hdu.writeto('../../all/simplestacks/stack.fits', clobber=True)
		# fmuse = np.nansum(fmuse[zl/2-zw:zl/2+zw+1, :, :], 0)
		rads = [[0, 2], [2, 4], [4, 6], [6, 8], [8, 12], [12, 16], [16, 20], [20, 24], [24, 30]]
		# rads = [[4, 6], [4, 8], [4, 12], [6, 14], [8, 16], [10, 20], [14, 20], [18, 30]]
		rpix_muse = np.array(rads) / asec2pix_muse
		rpix_eagle = np.array(rads) / asec2pix_eagle
		y, x = np.ogrid[0: yl, 0: xl]
		c = (x - xl / 2.) ** 2 + (y - yl / 2.) ** 2
		fconn = []
		fstds = []
		fig, ax = plt.subplots(figsize=(7, 4))
		cmap = get_cmap('tab10')
		color = 0
		all = []
		if _fin is None:
			for k in range(len(rads)):
				rm = rpix_muse[k]
				re = rpix_eagle[k]
				r = rads[k]
				print 'Between %.1f and %.1f arcsec' % (r[0], r[1])
				if _fin is None:
					inside_muse = (c > rm[0] ** 2) & (c < rm[1] ** 2)
					inside_eagle = (c > re[0] ** 2) & (c < re[1] ** 2)
					nin = np.sum(inside_muse)
					# print 'Number of pixels', nin
					fin = np.nanmean(fmuse[inside_muse])

					if 0:
						# jackniffe
						fjack = []
						fout = []
						fstds = []
						nall = len(fall)
						# for i in range(nin):
						#	fjack.append(np.nanmean(np.delete(fmuse[inside_muse], i)))
						for i in range(nall):
							fstds.append(np.nanstd(fall[i, inside_muse]))
							_f = np.delete(fall, i, 0)
							_fo = np.delete(fouts, i, 0)
							fjack.append(np.nanmean(_f[:, inside_muse]))
							fout.append(np.nanmean(_fo[:, inside_muse]))
						fstd = np.nanstd(fjack)
						fmean_out = np.nanmean(fout)
						fstd_out = np.nanstd(fout)
						fstd_stack = np.nanstd(fmuse[inside_muse])
						fstd_ind = np.nanmean(fstds)

						print 'individual std %f, overall %f, individual/sqrt(nall) %f' % (
						fstd_ind, fstd_stack, fstd_ind / np.sqrt(nall))
					else:
						nall = len(fall)
						fstds = [np.nanstd(fall[i, inside_muse]) for i in range(nall)]
						fstds_out = [np.nanstd(fouts[i, inside_muse]) for i in range(nall)]
						fix = np.nanstd(fmuse[inside_muse]) * np.sqrt(nall) / np.nanmean(fstds)
						fix_out = np.nanstd(fout[inside_muse]) * np.sqrt(nall) / np.nanmean(fstds_out)
						# print np.nanstd(fmuse[inside_muse]), np.nanmean(fstds)/np.sqrt(nall), fix, fix_out
						fstd = np.nanstd(fmuse[inside_muse]) * fix / np.sqrt(nin)
						fmean_out = np.nanmean(fout[inside_muse])
						fstd_out = np.nanstd(fout[inside_muse]) * fix_out / np.sqrt(nin)
					fe = np.nanmean(feagle[inside_eagle])
					fe_out = np.nanmean(feagle_out[inside_eagle])
					# fjack = []
					# for i in range(len(feagles)): fjack.append(np.nanmean(np.delete(feagles, i, 0)[:, inside_eagle]))
					# del feagles
					fe_std = np.nanstd(feagle[inside_eagle]) / np.sqrt(nin)
				else:
					fin = _fin[k]
					fe = _fe[k]
					fstd = _fstd[k]
					fe_std = _fe_std[k]
					fmean_out = _fmean_out[k]
					fstd_out = _fstd_out[k]
				tot_std = np.sqrt((fe_std / fe) ** 2 + (fstd / fin) ** 2)
				fy = fin * conv / fe * (1 + 2 * tot_std)
				print 'SB %.3f std %.3f SBout %.3f stdout %.3f. fconn %.3f out %.3f. SB EAGLE %.3f. 2sigma SB uplim %.3f' % \
					  (fin * conv, fstd * conv, fmean_out * conv, fstd_out * conv, fe, fe_out, fe * SB, fy)
				rm = (r[0] + r[1]) / 2
				if color == 0:
					lbs = ['Measured', 'Expected', 'Upper limit']
				else:
					lbs = [None, None, None]
				# ax.scatter(rm, fin*conv, c=cmap(color), marker='o', label=lbs[0])
				# ax.errorbar(rm, fin*conv, xerr=(r[1]-r[0])/2., yerr=2*fstd*conv, c=cmap(color), capsize=4, alpha=.4)
				# ax.scatter(rm, fe*SB, c=cmap(color), marker='*', label=lbs[1])
				# ax.errorbar(rm, fe*SB, yerr=2*fe_std*SB, c=cmap(color), capsize=4, alpha=.4)
				# ax.errorbar(rm, fe * SB, yerr=2 * fe_std * SB, c='red', capsize=4, alpha=.4)
				# ax.errorbar(rm, fy, xerr=.3, yerr=.2*fy, uplims=fy, c=cmap(color), capsize=None)
				# ax.errorbar(rm, fy, xerr=.3, yerr=.04, uplims=fy, c=cmap(color), capsize=None)
				# ax.set_yticklabels([.03, .05, .1, .3, .5, 1, 3])
				color += 1
				all.append([r[0], r[1], fin, fe, fstd, fe_std, fmean_out, fstd_out, fe_out])
				fconn.append(fin)
				fstds.append(fstd)
		if _fin is None:
			all = np.array(all).T
			_r0, _r1, _fin, _fe, _fstd, _fe_std, _fmean_out, _fstd_out, _fe_out = all
			np.savetxt(fdat, all)
		rm = (_r0 + _r1) / 2
		ax.scatter(rm, _fin * conv, c='blue', marker='o', label=r'MUSE')
		ax.plot(rm, _fin * conv, c='blue', marker='o')
		plt.fill_between(rm, (_fin - 2 * _fstd) * conv, (_fin + 2 * _fstd) * conv, facecolor='blue', alpha=0.3, lw=0,
						 edgecolor='none')
		# ax.errorbar(rm, _fin*conv, xerr=(_r1-_r0)/2., c='grey', capsize=4, alpha=.4)
		ax.scatter(rm, _fe * SB, c='red', marker='*', label=r'HM12 EAGLE')
		ax.plot(rm, _fe * SB, c='red')
		ax.plot(rm, _fmean_out*SB, c='pink')
		ax.scatter(rm, _fmean_out*SB, c='pink', marker='.', label=r'MUSE outside')
		ax.plot(rm, _fe_out*SB, c='gray')
		ax.scatter(rm, _fe_out*SB, c='gray', marker='.', label=r'HM12 EAGLE outside')
		plt.fill_between(rm, (_fe - 2 * _fe_std) * SB, (_fe + 2 * _fe_std) * SB, facecolor='red', alpha=0.3, lw=0,
						 edgecolor='none')
		plt.semilogy()
		# plt.ylim(0, 1)
		plt.xlabel('distance [arcsec]')
		plt.ylabel(r'$\mathrm{SB_{Ly\alpha}\,[10^{-20}erg/s/cm^2/arcsec^2]}}$')
		plt.legend()
		# leg = ax.get_legend()
		# for ls in leg.legendHandles: ls.set_color('black')
		plt.savefig('../../Figures/SB_upper-limits.pdf', ext='pdf', pdi=200)
		# plt.show()
		plt.close()

if superstack:

	def sclipping(fits, nsigma, dim=None, mask=None):
		# print 'Asign nan to values > nsigma in a fits array'
		if dim is None:
			stds = np.nanstd(fits[:, mask])
		else:
			stds = np.nanstd(fits[:, mask], dim)
		high_sigma = np.abs(fits) > nsigma * stds
		fits[high_sigma] = np.nan
		return fits, stds


	ds = [[0, 2], [1, 3], [2, 4], [3, 6], [5, 8], [7, 10], [9, 15], [14, 20]]
	spmin = {}
	spmax = {}
	sbmin = {}
	sbmax = {}
	srand = {}
	stdmin = {}
	stdmax = {}
	zw0 = 5
	folder = '../../EAGLE/stacks/'
	f2 = '../..all/stacks/'
	f3 = '../../EAGLE/cats/'
	hdu = PrimaryHDU()
	imtype = 'mean'  # 'median'#
	mask = False
	smask = '.mask' * mask
	stype = ('.' + imtype) * (imtype == 'median')
	for dd in ds:
		dmin = dd[0]
		dmax = dd[1]
		all = []
		ncores = dmax
		for snap in snaps:
			red = redshifts[str(snap)]
			asec2kpc = asec2kpcs[str(snap)]
			thetamin = 16
			# height 2 arcsec, width 6 to 12 arcsec -> same as my paper!
			xasecmin = 6
			xasecmax = 12
			yasecmin = -1
			yasecmax = 1
			bin = 2.  # binning after creating subcubes
			asec2pix = asec2kpc * (1 + red) * kpc2pix / bin
			sb2flux = asec2pix ** -2.
			# limits in pixels
			xmin = int(xasecmin * asec2pix)
			xmax = int(xasecmax * asec2pix)
			ymin = int(yasecmin * asec2pix)
			ymax = int(yasecmax * asec2pix)
			# limits in pixels
			# xmin =

			rad_asec = np.array([rads[i + 1] for i in range(len(rads) - 1)]) / asec2pix
			data = getdata(f3 + 'lae_pairs_snap%d.fits' % snap, 1)
			dist = data['com_dist']

			props = ['U1', 'n5d1']
			pr = props[1]
			doall = True
			prop = data[pr]
			cool = (dist > dmin) & (dist < dmax)
			pmedian = abs(np.median(prop[cool]))
			print '%s median' % pr, pmedian
			if pr == props[0]:
				smax = 'd%s-%d_pd16-2000_d5th0.5-20.0_u0.0-%.1f' % (dmin, dmax, pmedian)
				smin = 'd%s-%d_pd16-2000_d5th0.5-20.0_u%.1f-99.0' % (dmin, dmax, pmedian)
				pmin = 'umin'
				pmax = 'umax'
				sprop = 'U'
			if pr == props[1]:
				smax = 'd%s-%d_pd16-2000_d5th0.5-%.1f_u0.0-99.0' % (dmin, dmax, pmedian)
				smin = 'd%s-%d_pd16-2000_d5th%.1f-20.0_u0.0-99.0' % (dmin, dmax, pmedian)
				pmin = 'd5min'
				pmax = 'd5max'
				sprop = 'd5th'

			stacks_pmin = []
			stacks_pmax = []
			std_pmin = []
			std_pmax = []
			for coord in coordnames:  # ['x']:#
				# os.system(
				#		'mpirun -n %d python stackEAGLE.py -snap %d -coord %s -overwrite True -imtype %s -mask %s'
				#		% (ncores, snap, coord, imtype, mask))
				print 'd %d-%d snap %d coord %s' % (dmin, dmax, snap, coord)
				sname = folder + 'snap%d_%s_%s%s%s/stack.fits' % (
					snap, coord, smax, stype, smask)
				if not os.path.isfile(sname) or overwrite:
					os.system(
						'mpirun -n %d python stackEAGLE.py -snap %d -coord %s -dmin %d -dmax %d -%s %f -overwrite True -imtype %s -mask %s'
						% (ncores, snap, coord, dmin, dmax, pmax, pmedian, imtype, mask))
				_stack = getdata(sname)
				zl, yl, xl = _stack.shape
				stacks_pmax.append(_stack[:, :, :yl])
				if doall:
					_std = getdata(sname.replace('.fits', '.STD.fits'))
					std_pmax.append(_std[:, :, :yl])
				sname = folder + 'snap%d_%s_%s%s%s/stack.fits' % (
					snap, coord, smin, stype, smask)
				if not os.path.isfile(sname) or overwrite:
					os.system(
						'mpirun -n %d python stackEAGLE.py -snap %d -coord %s -dmin %d -dmax %d -%s %f -overwrite True -imtype %s -mask %s' %
						(ncores, snap, coord, dmin, dmax, pmin, pmedian, imtype, mask))

				_stack = getdata(sname)
				zl, yl, xl = _stack.shape
				stacks_pmin.append(_stack[:, :, :yl])
				if doall:
					_std = getdata(sname.replace('.fits', '.STD.fits'))
					std_pmin.append(_std[:, :, :yl])
			stack_pmin = np.nanmean(stacks_pmin, 0)
			stack_pmax = np.nanmean(stacks_pmax, 0)
			if doall:
				std_pmin = np.nanmean(std_pmin, 0)
				std_pmax = np.nanmean(std_pmax, 0)
			zl, yl, xl = stack_pmin.shape
			hdu.data = stack_pmax
			fs = folder + 'snap%d_d%d-%d_pd16-2000_d5th0.5-4.0_u0.0-%.1f%s%s/' % (
				snap, dmin, dmax, pmedian, stype, smask)
			if not os.path.isdir(fs): glob.os.makedirs(fs)
			if not os.path.isfile(fs + 'stack.IM.fits'):
				hdu.writeto(fs + 'stack.fits', clobber=True)
				stackim = np.nanmean(stack_pmax[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
				hdu.data = stackim
				hdu.writeto(fs + 'stack.IM.fits', clobber=True)
			hdu.data = stack_pmin
			fs = folder + 'snap%d_d%d-%d_pd16-2000_d5th0.5-4.0_u%.1f-99.0%s%s/' % (
				snap, dmin, dmax, pmedian, stype, smask)
			if not os.path.isdir(fs): glob.os.makedirs(fs)
			if not os.path.isfile(fs + 'stack.IM.fits'):
				hdu.writeto(fs + 'stack.fits', clobber=True)
				stackim = np.nanmean(stack_pmin[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
				hdu.data = stackim
				hdu.writeto(fs + 'stack.IM.fits', clobber=True)

			x0, x1 = xl / 2 + xmin, xl / 2 + xmax + 1
			_x0, _x1 = xl / 2 - xmax, xl / 2 - xmin + 1
			y0, y1 = yl / 2 + ymin, yl / 2 + ymax + 1
			_y0, _y1 = yl / 2 - ymax, yl / 2 - ymin + 1
			z0, z1 = zl / 2 - zw0, zl / 2 + zw0 + 1
			spnhi = np.power(10, sb2nhi(stack_pmin[z0: z1, y0: y1, x0: x1]))
			spmin[('%d' % snap, '%d' % dmax)] = np.nanmean(nhi2sb(np.log10(np.nansum(spnhi, 0))))
			spnhi = np.power(10, sb2nhi(stack_pmax[z0: z1, y0: y1, x0: x1]))
			spmax[('%d' % snap, '%d' % dmax)] = np.nanmean(nhi2sb(np.log10(np.nansum(spnhi, 0))))
			if doall:
				# npix to properly propagate the std
				npix = 1. / np.sqrt(3 * (xmax - xmin) * (ymax - ymin) * zw0)

				stdmin[('%d' % snap, '%d' % dmax)] = np.nanmean(std_pmin[z0: z1, y0: y1, x0: x1]) * npix
				stdmax[('%d' % snap, '%d' % dmax)] = np.nanmean(std_pmax[z0: z1, y0: y1, x0: x1]) * npix

				# srand, random means average over other directions (up, left, bottom)
				sr = np.nanmean(stack_pmin[z0: z1, y0: y1, _x0: _x1])
				sr += np.nanmean(stack_pmin[z0: z1, _y0: _y1, x0: x1])
				sr += np.nanmean(stack_pmin[z0: z1, _y0: _y1, _x0: _x1])
				sr += np.nanmean(stack_pmax[z0: z1, y0: y1, _x0: _x1])
				sr += np.nanmean(stack_pmax[z0: z1, _y0: _y1, x0: x1])
				sr += np.nanmean(stack_pmax[z0: z1, _y0: _y1, _x0: _x1])
				srand[('%d' % snap, '%d' % dmax)] = sr * (2 * zw0 + 1) / 6.

			if 0:
				all.append(stack * sbpeaks[str(snap)])
				stack = np.nanmean(all, 0)
				zl, yl, xl = stack.shape
				hdu.data = stack
				fs = folder + 'd%d-%d_pd16-2000_nw0-1000%s%s/' % (d, d2, stype, smask)
				if not os.path.isdir(fs): glob.os.makedirs(fs)
				hdu.writeto(fs + '/stack.fits', clobber=True)
				stackim = np.nansum(stack[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
				hdu.data = stackim
				hdu.writeto(fs + 'stack.IM.fits', clobber=True)
		if 0:

			fn = f2 + 'd0-%d_los0-20_pd16-2000_z2.9-4.0_l0-2000_n1-1000_v0-10000/' % d
			nname = fn + 'stack.NPIX.IM.fits'
			if not os.path.isfile(nname): os.system('python stack.py -dmax %d -overwrite True' % d)
			npix = getdata(nname)
			yl, xl = npix.shape
			noise = 10.  # assuming 1e-19 cgs noise
			snr = stackim[:, :xl] * np.sqrt(npix) / noise / (
				2 * zw0 + 1.)  # for the image stack I am doing a sum therefore the noise does not decrease, so I added the factor 1/(2*zw0+1.)
			hdu.data = snr
			hdu.writeto(fs + 'stack.SNR.fits', clobber=True)

	types = ['o', '^', 's', 'D', 'p']
	colors = ['red', 'green', 'blue', 'orange', 'magenta']
	SB = [.73, 1.27, 2.49, 5.14, 6.97]
	fLLS = [0.01329666667, 0.008964333333, 0.005651, 0.004793716667, 0.003494351]
	fconn = [0.01999333333, 0.01287333333, 0.007092333333, 0.004155633333, 0.0030469]
	dpi = 200
	ext = 'pdf'
	zs = [redshifts['%d' % s] for s in snaps]
	plt.figure()
	plt.ylim([0, .021])
	plt.yticks([0.000, .005, .01, .015, .02])
	plt.xlabel(r'$\mathrm{redshift}$')
	# plt.xticks(zs)
	plt.scatter(zs, fLLS, color='red', label=r'$\mathrm{f_{LLS}}$')
	plt.plot(zs, fLLS, color='red')
	plt.scatter(zs, fconn, color='blue', label=r'$\mathrm{SB/SB_{UVB}}$')
	plt.plot(zs, fconn, color='blue')
	plt.legend()
	plt.savefig('../../EAGLE/plots/fLLS_fconn-vs-redshift.%s' % ext, dpi=dpi, format=ext)
	plt.close()
	sb2cMpc = []

	for i in range(len(snaps)):
		snap = snaps[i]
		t = types[i]
		c = colors[i]
		sb = SB[i]
		z = redshifts['%d' % snap]
		print snap
		fmin = np.array([spmin[('%d' % snap, '%d' % dd[1])] for dd in ds])
		fmax = np.array([spmax[('%d' % snap, '%d' % dd[1])] for dd in ds])
		f = (fmin + fmax) / 2.
		x = np.array([(dd[0] + dd[1]) / 2. for dd in ds])
		f0 = np.zeros(len(fmax))
		xlabel = r'neighbour distance [cMpc]'
		xlim = [1, 17]

		plt.figure(1)
		plt.ylim([-.25, .25])
		plt.xlim(xlim)
		plt.xlabel(xlabel)
		plt.ylabel(r'$\Delta\mathrm{f_{conn,\,%s}}$' % sprop, fontsize=13)
		plt.plot(x, (fmin - fmax) / f, color=c)
		plt.scatter(x, (fmin - fmax) / f, color=c, marker=t, label=r'$z= %.1f$' % z)  # , label='Snap %d' % snap)

		plt.figure(5)
		plt.xlim()
		plt.xlabel(r'$\mathrm{redshift}$')
		plt.ylabel(r'$\mathrm{SB\,[10^{-20}erg/s/cm^2/arcsec^2]}}$', fontsize=13)
		plt.scatter(z, f[0] * sb, color=c, marker=t)
		sb2cMpc.append(f[0] * sb)

		if doall:
			frand = np.array([srand[('%d' % snap, '%d' % dd[1])] for dd in ds])
			smin = np.array([stdmin[('%d' % snap, '%d' % dd[1])] for dd in ds])
			smax = np.array([stdmax[('%d' % snap, '%d' % dd[1])] for dd in ds])
			s = (smin + smax) / 2.
			plt.figure(2)
			# plt.ylim([0, .5])
			# plt.yticks([0, .03, .06, .09, .12, .15, .18, .21])
			plt.xlim(xlim)
			plt.xlabel(xlabel)
			plt.ylabel(r'$\mathrm{f_{conn}}$', fontsize=13)
			plt.plot(x, f, color=c)
			plt.scatter(x, f, color=c, marker=t, label=r'$z=%.1f$' % z)
			plt.fill_between(x, np.max([f0, f - s], 0), f + s, facecolor=c, alpha=0.3, lw=0, edgecolor='none')

			plt.figure(3, figsize=(7, 5))
			plt.ylim([-.07, .17])
			plt.xlim(xlim)
			plt.xlabel(xlabel)
			plt.ylabel(r'$\mathrm{(f_{conn}-f_{rand})/f_{rand}}$', fontsize=13)
			plt.plot(x, (f - frand) / frand, color=c)
			plt.scatter(x, (f - frand) / frand, color=c, marker=t, label=r'$z=%.1f$' % z)
			# plt.fill_between(x, np.max([f0, f-frand-s], 0), f-frand+s, facecolor=c, alpha=0.3, lw=0, edgecolor='none')

			plt.figure(4)
			# plt.ylim([0, .2])
			# plt.yticks([0, .03, .06, .09, .12, .15, .18, .21])
			plt.xlim(xlim)
			plt.xlabel(xlabel)
			plt.ylabel(r'$\mathrm{SB\,[10^{-20}erg/s/cm^2/arcsec^2]}}$', fontsize=13)
			plt.plot(x, f * sb, color=c)
			plt.scatter(x, f * sb, color=c, marker=t, label=r'$z=%.1f$' % z)
		# plt.fill_between(x, np.max([f0, (f-s)*sb], 0), (f+s)*sb, facecolor=c, alpha=0.3, lw=0, edgecolor='none')

	plt.figure(1)
	plt.legend()
	plt.savefig('../../EAGLE/plots/%s_diff-vs-d.%s' % (sprop, ext), dpi=dpi, format=ext)
	plt.close()
	plt.figure(5)
	plt.plot(zs, sb2cMpc, color='black', zorder=-1)
	plt.savefig('../../EAGLE/plots/SB-vs-z.%s' % ext, dpi=dpi, format=ext)
	plt.close()

	if doall:
		plt.figure(2)
		plt.legend()
		plt.savefig('../../EAGLE/plots/fconn-vs-d.%s' % ext, dpi=dpi, format=ext)
		plt.close()
		plt.figure(3)
		plt.legend()
		plt.savefig('../../EAGLE/plots/fconn-frand.%s' % ext, dpi=dpi, format=ext)
		plt.close()
		plt.figure(4)
		plt.legend()
		plt.savefig('../../EAGLE/plots/SB-vs-d.%s' % ext, dpi=dpi, format=ext)
		plt.close()

if 0:
	fout = '../../Figures/'
	highq = True
	if highq:
		ext = 'pdf'
	else:
		ext = 'png'
	if 1:
		f = '../../EAGLE/stacks/d0-2_pd16-2000_nw0-1000/'
		fits = getdata(f + 'stack.fits')
		zl, yl, xl = fits.shape
		ff = np.mean(fits[zl / 2 - 1: zl / 2 + 1, :, :], 0)
		hdu = PrimaryHDU()
		hdu.data = ff

		hdu.writeto(f + 'stack.IM.2.fits', clobber=True)
		astroim(f + 'stack.IM.2.fits', smooth=True, saveim=True, show=False,
				cbfrac=.08, pad=.006, dfig=(8, 8), contours=True,
				x0=0, y0=0, imout=fout + 'stack-d2-snap11.%s' % ext, std=None, scale=False, regcolor=None,
				scb_label=r'$\rm{f_{conn}}$', xticks=[-16, -12, -8, -4, 0, 4, 8, 12, 16],
				yticks=[-12, -8, -4, 0, 4, 8, 12],
				title='', vmin=0, vmax=.25, gray=False, sb=False, cmap='rainbow',
				xmin=-25, xmax=25, ymin=-20, ymax=20, pcolor='white', highq=highq, dpi=200, nbins=5)

	f = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n8-100/'
	astroim(f + 'stack.IM.fits', smooth=True, saveim=True, show=False,
			cbfrac=.08, pad=.006, dfig=(8, 8), contours=True,
			x0=0, y0=0, imout=fout + 'stack-muse-nmin8.%s' % ext, std=2, scale=False, regcolor=None,
			scb_label=r'$\rm{SB\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}',
			xticks=[-16, -12, -8, -4, 0, 4, 8, 12, 16], yticks=[-12, -8, -4, 0, 4, 8, 12],
			title='', vmin=-0.3, vmax=2, gray=False, cmap='rainbow', sb=True,
			xmin=-25, xmax=25, ymin=-20, ymax=20, pcolor='white', highq=highq, dpi=200, nbins=5)

if radprof:
	def sclipping(fits, nsigma, dim=None, mask=None):
		# print 'Asign nan to values > nsigma in a fits array'
		if dim is None:
			stds = np.nanstd(fits[:, mask])
		else:
			stds = np.nanstd(fits[:, mask], dim)
		high_sigma = np.abs(fits) > nsigma * stds
		fits[high_sigma] = np.nan
		return fits, stds


	zw = 1
	yw = 2
	zoff = 0
	dw = 1
	c = (2 * zw + 1) * 7.8125
	fin = '/net/galaxy-data/export/galaxydata/gallegos/'
	folder = '../../EAGLE/stacks/snap11_x_galstack/'
	fstacks = [getdata(f) for f in glob.glob('%s/EAGLE/gals/[!IM]*.fits' % fin)]

	folder2 = '../../all/stacks/laestack/'
	fits = [getdata(f) for f in glob.glob('%s/UDF/LAEs/lae[!IM]*.fits' % fin)]
	fits, stds = sclipping(fits, 3, 0)

	stack = getdata(folder + 'stack.fits')
	stack2 = getdata(folder2 + 'stack.fits')
	frs = []
	for i in range(200):
		hdur = getdata(folder2 + 'randoms/random_stack.%d.fits' % i)
		frs.append(hdur)
	zl, yl, xl = stack.shape
	zl2, yl2, xl2 = stack2.shape
	print xl, xl2, yl, yl2, zl, zl2
	xmin = yl / 2
	t = np.concatenate([np.arange(xmin, xmin + 12, 1), [xmin + 13, xmin + 18, xmin + 27, xmin + 41]]).astype('int')
	zi = zl / 2 - zw / 2
	zf = zl / 2 + zw / 2 + 1
	zi2 = zl2 / 2 - zw / 2 + zoff
	zf2 = zl2 / 2 + zw / 2 + 1 + zoff
	y, x = np.ogrid[-yl / 2 + 1: yl / 2 + 1, -xl / 2 + 1: xl / 2 + 1]
	y2, x2 = np.ogrid[-yl2 / 2 + 1: yl2 / 2 + 1, -xl2 / 2 + 1: xl2 / 2 + 1]
	for i in range(len(t) - 1):
		sb.append((t[i] + t[i + 1]) * .5)
		cool = (x ** 2 + y ** 2 >= t[i] ** 2) & (x ** 2 + y ** 2 < t[i + 1] ** 2)
		cool2 = (x2 ** 2 + y2 ** 2 >= t[i] ** 2) & (x2 ** 2 + y2 ** 2 < t[i + 1] ** 2)
		sb.append(np.nanmean(stack[zi:zf, cool]))
		sb.append(np.nanmean(stack2[zi2:zf2, cool2]) * c)
		sb.append(np.nanstd([np.nanmean(ff[zi:zf, cool2]) * c for ff in frs]))
	lx = len(t) - 1
	# ld = len(dir)*3 + 1
	ld = 5
	sb = np.array(sb).reshape((lx, ld))
	np.savetxt('../../EAGLE/plots/radprof.txt', sb,
			   header='dpix sim muse std')

	x1 = 2
	x2 = 9
	y1 = 0.05
	y2 = 12
	hm = False
	hmSB = 1.14
	y1b = -4
	y2b = 6
	fsize = 30
	width = 4
	sigmas = 2

	plt.figure(figsize=(20, 12))
	plt.loglog()
	ax = plt.axes()
	cax = plt.gca()
	plt.xlim([x1, x2])
	xrange = np.arange(2, 9, 2)
	plt.xticks(xrange, xrange, fontsize=fsize)
	plt.xlabel(r"$\theta$ [arcsec]", fontsize=fsize)
	plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize + 2)
	yticks = [0, 2, 4, 6, 8, 10]
	plt.yticks(yticks, yticks, fontsize=fsize)
	plt.minorticks_off()
	plt.twiny(ax=None)
	plt.loglog()
	plt.xlim([x1, x2])
	plt.xlabel(r'$r_p\,\rm{[pkpc]}$', fontsize=fsize)
	kpc = np.arange(15, 90, 15) / 7.47
	plt.xticks(kpc, (kpc * 7.47).astype('int'), fontsize=fsize)
	plt.ylim([y1, y2])
	plt.plot((2, x2), (0, 0), '--', label='', color="gray", linewidth=width + 1)
	if hm: plt.plot((x1, x2), (hmSB, hmSB), '--', label=r"LLS Fluorescence from HM12", color="green",
					linewidth=width + 1)

	x = (t[:-1] + (t[1:] - t[:-1]) / 2. - yl / 2) * .4

	# plt.plot(x, sb[:, 2], label='MUSE oriented', lw=width, color='black')
	# plt.scatter(x, sb[:, 2], s=width * 10, color='black')
	sb = np.array(sb)
	sb[sb < 0] = 0.001
	plt.plot(x, sb[:, 3], label='MUSE', lw=width, color='dodgerblue')
	plt.scatter(x, sb[:, 3], s=width * 10, color='dodgerblue')
	a = sb[:, 3] - sigmas * sb[:, 4]
	a[a < 0] = 0.001
	plt.fill_between(x, a, sb[:, 3] + sigmas * sb[:, 4], facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')

	sbuvmean = 0.8
	sbuvup = 1.2
	sbuvdown = 0.6
	plt.plot(x, sb[:, 1] * sbuvmean, label='EAGLE+HM12', lw=width, color='red')
	plt.scatter(x, sb[:, 1] * sbuvmean, s=width * 10, color='red')
	plt.fill_between(x, sb[:, 1] * sbuvdown, sb[:, 1] * sbuvup, facecolor='red', alpha=0.3, lw=0, edgecolor='none')
	plt.minorticks_off()
	plt.legend(fontsize=fsize, loc='best')  # (3.5,2))
	plt.savefig('../../analysis/muse-vs-eagle_radprof.png')

if sbprof:
	zw = 5
	yw = 2
	zoff = -1
	dw = 1
	c = (2 * zw + 1) * 7.8125
	asec2kpc = 7.47
	folder = '../../EAGLE/stacks/snap11_x_d0-20_pd16-2000_nw0.000-1.000/'
	folder2 = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100/'
	stack = getdata(folder + 'stack.fits')
	stack2 = getdata(folder2 + 'stack.fits')
	srand = getdata(folder2 + 'random_stack.fits')
	frs = []
	for i in range(200):
		hdur = getdata(folder2 + 'randoms/random_stack.%d.fits' % i)
		frs.append(hdur)
	zl, yl, xl = stack.shape
	zl2, yl2, xl2 = stack2.shape
	print xl, xl2, yl, yl2, zl, zl2
	xmin = yl / 2
	t = np.concatenate([np.arange(xmin, xmin + 12, 1), [xmin + 13, xmin + 18, xmin + 27, xmin + 41]]).astype('int')
	zi = zl / 2 - zw / 2
	zf = zl / 2 + zw / 2 + 1
	zi2 = zl2 / 2 - zw / 2 + zoff
	zf2 = zl2 / 2 + zw / 2 + 1 + zoff
	if 1:
		dir = ['right', 'top', 'left', 'bottom']
		yymin = yl / 2 - yw / 2
		yymax = yl / 2 + yw / 2 + 1
		sb = []

		for i in range(len(t) - 1):
			yis = [yymin, t[i] - 1 - dw, yymin, yl - t[i + 1] - dw]
			yfs = [yymax, t[i + 1] + dw, yymax, yl - t[i] + 1 + dw]
			xis = [t[i] - 1 - dw, yymin, yl - t[i + 1] - dw, yymin]
			xfs = [t[i + 1] + dw, yymax, yl - t[i] + 1 + dw, yymax]
			sb.append((t[i] + t[i + 1]) * .5)
			# for yi, yf, xi, xf, ds in zip(yis, yfs, xis, xfs, dir):
			yi, yf, xi, xf, ds = zip(yis, yfs, xis, xfs, dir)[0]
			if 1:
				sb.append(np.nanmean(stack[zi:zf, yi:yf, xi:xf]))
				sb.append(np.nanmean(stack2[zi2:zf2, yi:yf, xi:xf]) * c)
				sb.append(np.nanmean(srand[zi2:zf2, yi:yf, xi:xf]) * c)
				sb.append(np.nanstd([np.nanmean(ff[zi:zf, yi:yf, xi:xf]) * c for ff in frs]))
	lx = len(t) - 1
	# ld = len(dir)*3 + 1
	ld = 5
	sb = np.array(sb).reshape((lx, ld))
	np.savetxt('../../EAGLE/plots/fconn.txt', sb,
			   header='dpix sim muse rand rand_std')

	x1 = 2
	x2 = 9
	y1 = 0.0001
	y2 = 12
	hm = False
	hmSB = 1.14
	y1b = -4
	y2b = 6
	fsize = 30
	width = 4
	sigmas = 2

	plt.figure(figsize=(20, 12))
	plt.loglog()
	ax = plt.axes()
	cax = plt.gca()
	plt.xlim([x1, x2])
	xrange = np.arange(2, 9, 2)
	plt.xticks(xrange, xrange, fontsize=fsize)
	plt.xlabel(r"$\theta$ [arcsec]", fontsize=fsize)
	# plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize + 2)
	yticks = [0, 2, 4, 6, 8, 10]
	plt.yticks(yticks, yticks, fontsize=fsize)
	plt.minorticks_off()
	plt.twiny(ax=None)
	plt.loglog()
	plt.xlim([x1, x2])
	plt.xlabel(r'$r_p\,\rm{[pkpc]}$', fontsize=fsize)
	kpc = np.arange(15, 90, 15)
	plt.xticks(kpc / asec2kpc, kpc.astype('int'), fontsize=fsize)
	plt.ylim([y1, y2])
	plt.plot((2, x2), (0, 0), '--', label='', color="gray", linewidth=width + 1)
	if hm: plt.plot((x1, x2), (hmSB, hmSB), '--', label=r"LLS Fluorescence from HM12", color="green",
					linewidth=width + 1)

	x = (t[:-1] + (t[1:] - t[:-1]) / 2. - yl / 2) * .4  # 0.4 arcsec per pixel

	# plt.plot(x, sb[:, 2], label='MUSE oriented', lw=width, color='black')
	# plt.scatter(x, sb[:, 2], s=width * 10, color='black')
	sb = np.array(sb)
	sb[sb < 0] = 1e-10
	plt.plot(x, sb[:, 3], label='SB MUSE', lw=width, color='dodgerblue')
	plt.scatter(x, sb[:, 3], s=width * 10, color='dodgerblue')
	ymin = sb[:, 3] - sigmas * sb[:, 4]
	ymin[ymin < 0] = 1e-10
	plt.fill_between(x, ymin, sb[:, 3] + sigmas * sb[:, 4], facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')

	sbuvmean = 0.8
	sbuvup = 1.2
	sbuvdown = 0.6
	plt.plot(x, sb[:, 1] * sbuvmean, label='SB EAGLE+HM12', lw=width, color='red')
	plt.scatter(x, sb[:, 1] * sbuvmean, s=width * 10, color='red')
	plt.fill_between(x, sb[:, 1] * sbuvdown, sb[:, 1] * sbuvup, facecolor='red', alpha=0.3, lw=0, edgecolor='none')
	plt.minorticks_off()
	gamma = 3.16e13  # ionizing photon rate from the galaxy /1e40
	kpc2cm = 30.86  # kpc2cm /1e20
	zmean = 3.5
	sb_diff = sb[:, 3] - sb[:, 1] * sbuvmean
	E_lya = 1.64e-11
	asec2rad_sqd = 2.35e-11
	fesc = 4 * np.pi * (1 + zmean) ** 4 * (
		x * asec2kpc * kpc2cm) ** 2 * sb_diff * 1e-20 / E_lya / .6 / gamma / asec2rad_sqd
	gamma_std = 3e13
	obs_std = sb[:, 4]
	sim_std = 0.3
	sb_std = np.sqrt(np.power(obs_std, 2) + sim_std ** 2)
	tot_std = np.sqrt((gamma_std / gamma) ** 2 + (sb_std / sb_diff) ** 2) * fesc
	print 'fesccccc', fesc
	print 'sb std', sb_std
	print 'tot std', tot_std
	fesc[fesc < 0] = 1e-10
	plt.plot(x, fesc, label='Escape fraction', lw=width, color='black')
	ymin = fesc - tot_std
	ymin[ymin < 0] = 1e-10
	plt.fill_between(x, ymin, fesc + tot_std, facecolor='gray', alpha=0.3, lw=0, edgecolor='none')

	plt.legend(fontsize=fsize, loc='best')  # (3.5,2))

	plt.savefig('../../analysis/muse-vs-eagle.png')

if unique:
	# I'll use pixelssss
	yw = 5
	yl = 141
	n = 10000
	ngal = 303
	nsub = 2016
	nfrac = nsub / ngal
	tot = 0

	# in asec
	rads = [2, 4, 6, 10, 20]

	# for z 3
	fconn = [.17, .11, .105, .10, .11]
	# for z 3.5
	fconn = [.14, .10, .9, .85, .9]
	frand = [.89, .96, .99, 1.02, 1]
	noise = np.sqrt(3)
	SB = [.73, 1.27, 2.49, 5.14, 6.97]
	asec2pix = 5
	Nsub_udf = [46, 80, 140, 208, 352]
	Ngal_udf = [30, 48, 61, 65, 75]

	Nsub_mosaic = [375, 1413, 2571, 4120, 7091]
	Ngal_mosaic = [207, 292, 311, 324, 328]

	for d, f, fr, nsubm, ngalm, nsubu, ngalu in zip(rads, fconn, frand, Nsub_mosaic, Ngal_mosaic, Nsub_udf, Ngal_udf):
		tot = 0
		npix = int(2 * np.pi * d)
		nfrac = nsubm / ngalm
		print 'Mosaic d', d, 'nsub', nsubm, 'ngal', ngalm
		print npix, 'pixels at rad', d
		for j in range(n):
			xrand = np.random.randint(0, npix, nfrac)
			tot += len(np.unique(xrand))
		lun = tot / float(n * nfrac)
		print 'less unique', lun
		print 'Efective signal', SB[1] * f * np.sqrt(lun * nsubm) / np.sqrt(3) / 10.
		print 'Random prof signal', SB[1] * f * fr * np.sqrt(lun * nsubm * npix) / np.sqrt(3) / 10.

		tot = 0
		nfrac = nsubu / ngalu
		print 'UDF d', d, 'nsub', nsubu, 'ngal', ngalu
		print npix, 'pixels at rad', d
		for j in range(n):
			xrand = np.random.randint(0, npix, nfrac)
			tot += len(np.unique(xrand))
		lun = tot / float(n * nfrac)
		print 'less unique', lun
		print 'Efective signal', SB[1] * f * np.sqrt(lun * nsubu) / 10.
		print 'Random prof signal', SB[1] * f * fr * np.sqrt(lun * nsubu * npix) / 10.
