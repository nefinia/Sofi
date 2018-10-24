#!/usr/bin/env python
__author__ = 'gallegos'
from tools_sofi import rdarg, astroim  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,
from pyfits import getdata, PrimaryHDU
import h5py
import numpy as np
import os, sys
from math import sqrt
import eagleSqlTools
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sys import argv
import glob
from copy import copy

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
snaps = rdarg(argv, 'snap', list, [10, 11, 12], int)
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
superstack = rdarg(argv, 'superstack', bool, False)
radprof = rdarg(argv, 'radprof', bool, False)
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
mask = rdarg(argv, 'mask', bool, False)
do_delaunay = rdarg(argv, 'overwrite', bool, False)

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


snames = {'10': '010_z003p984', '11': '011_z003p528', '12': '012_z003p017', '13': '013_z002p478', '14': '014_z002p237', '15': '015_z002p012'}
reds = {'10': 984, '11': 528, '12': 17}
redshifts = {'10': 3.984, '11': 3.528, '12': 3.017, '13': 2.478, '14': 2.237, '15': 2.012}
asec2kpcs = {'10': 7.842, '11': 7.449, '12': 7.108, '13': 8.241, '14': 8.396, '15': 8.516}
dz = {'10': 0.0255, '11': .035, '12': .03}
types = ['NHI', 'NHII']

#sbpeaks = [2.396]
#sbpeaks = {'10': 1.011, '11': 1.484, '12':2.396} # for an aperture of 1 asec^2!!! I am converting to flux later on
sbpeaks = {'10': .739, '11': 1.293, '12': 2.579} # for an aperture of 1 asec^2!!! I am converting to flux later on
zlens = {'10': 34, '11': 29, '12': 25, '13': 20, '14': 18, '15': 16}
lcube = 4096
coml = 25  # cMpc
com2pix = 163.84 # lcube/coml
kpc2pix = lcube / float(coml * 1e3)
rads = np.array([0, 2, 4, 8, 12, 20, 30, 50, 100, 200])


if LLScubes:
	cubefolder = '/net/astrogate/export/astrodata/gallegos/EAGLE/'
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
	import scipy.interpolate as inter
	nhi = [16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.5,    19,     19.4, 19.6, 19.8, 20.1, 20.5, 21]
	sb = [0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, 0.82,  0.885,  0.938, 0.96, 0.98, 1, 1, 1]
	s = inter.InterpolatedUnivariateSpline(nhi, sb)
	if 0:
		plt.close()
		import pylab as plt
		plt.plot(nhi, sb, 'bo', label='Data')
		x=np.arange(16, 20, .01)
		plt.plot(x,s(x), 'k--', label='Spline fit')
		plt.semilogy()
		#plt.legend()
		plt.show()

	for j in snaps:
		snap = str(j)
		red = redshifts[snap]
		sbpeak = sbpeaks[snap]
		asec2kpc = asec2kpcs[snap]
		zlen = zlens[snap]
		asec2pix = asec2kpc * (1 + red) * kpc2pix
		sb2flux = asec2pix ** -2.
		rad = 6 * asec2pix
		cat = getdata('../../EAGLE/cats/gals_snap%d.fits' % j, 1)

		for coord in coordnames:
			print 'Coord', coord, 'snap', snap
			sname = snames[snap]

			for type in types:
				fname = '/net/astrogate/export/astrodata/saeed/EAGLE/RefL0025N0752/snapshot_%s/layers/' % sname
				outname = '/net/astrogate/export/astrodata/gallegos/EAGLE/snap%s_%s.%s.fits' % (snap, coord, type)
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
							

	if 0:
		nprefs = [8, 12]

		cubefolder = '../../EAGLE/'
		cubenames = []
		for npr in nprefs:
			for c in coordnames:
				cubenames.append('%ssnap_012_z003p017_column_density_%s_npref%d.fits' % (cubefolder, c, npr))
		cubes = [getdata(cn) for cn in cubenames]
		cools = [np.log10(cb)>17.5 for cb in cubes]
		hdu = PrimaryHDU()
		for cool, cb, cn in zip(cools, cubes, cubenames):
			cb[cool] = 1
			cb[~cool] = 0
			cubeout = cn.replace('snap_012_z003p017_column_density', 'LLS')
			print cubeout
			hdu.data = cb
			hdu.writeto(cubeout, clobber=True)


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
			_y, _x = np.ogrid[1: yl+1, 1: xl+1]
			LLSfracs = []
			LLScums = []

			for i in (ids-1):
				print 'Galaxy %d of %d' % (i+1, ngal)
				x0 = c[0][i]
				y0 = c[1][i]
				#cool = ids == i
				nrads = len(rads)
				LLSfrac = []
				LLScum = []

				for j in range(nrads-1):
					cool = ((_x-x0-.5) ** 2 + (_y-y0-.5) ** 2 < rads[j+1] ** 2)
					LLScum.append(np.nanmean(fit[cool]))
					cool &= ((_x-x0-.5) ** 2 + (_y-y0-.5) ** 2 >= rads[j] ** 2)
					LLSfrac.append(np.nanmean(fit[cool]))

				if 0:
					plt.plot(rad_asec, LLSfrac)
					plt.xlabel('distance[arcsec]')
					plt.ylabel('LLS fraction')
					plt.savefig('../../EAGLE/LLS/plots/lls_frac%d.png' % (i+1))
					#plt.show()
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
		#plt.show()
		plt.close()
		plt.plot(rad_asec, LLS1, lw=2)
		plt.ylim(0, 1)
		plt.xlim(0, 8)
		plt.xlabel('distance[arcsec]')
		plt.ylabel('LLS fraction')
		plt.errorbar(rad_asec, LLS1, yerr=LLSstd1)
		plt.savefig('../../EAGLE/LLS/plots/LLSfracs-close_snap%d_%s.png' % (snap, cn))
		#plt.show()
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
		#plt.show()
		plt.close()

		LLS = [np.median(l) for l in LLScums]
		LLSstd = [np.std(l) for l in LLScums]
		x = 0.042 * np.array(rad_asec)**2
		dndz = np.array(LLS)*x
		dndzstd = np.array(LLSstd)*x
		plt.plot(rad_asec, dndz, lw=2)
		#plt.ylim([0.1, 1])
		plt.semilogy()
		plt.xlabel('distance[arcsec]')
		plt.ylabel('dn/dz(LLS,r<R)')
		plt.errorbar(rad_asec, dndz, yerr=dndzstd)
		plt.savefig('../../EAGLE/LLS/plots/lls_dndz_%s.png' % cn)
		#plt.show()
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
			nsteps = len(xs)-1
			# Approx 1 arcsec height
			ymin = -2.5
			ymax = 2.5
			fracs = []
			fracs_av = []
			f1 = '../../EAGLE/LLS/all_fracs_snap%d_%s.dat' % (snap, coord)
			f2 = '../../EAGLE/LLS/all_fracs_av_snap%d_%s.dat' % (snap, coord)
			pairsname = '/net/astrogate/export/astrodata/gallegos/EAGLE/LLScats/cubical_pair'
			good_pairs = [os.path.isfile(pairsname+'%d-%d.h5' % (i1, i2)) for i1, i2 in zip(id1, id2)]
			i1good = id1[good_pairs]
			i2good = id2[good_pairs]
			overwrite = False

			if not os.path.isfile(f1) or not os.path.isfile(f2) or overwrite:
				for i1, i2 in zip(i1good, i2good):
					print 'Pair %d-%d' % (i1, i2)
					cat = '/net/astrogate/export/astrodata/gallegos/EAGLE/cats/logNHI_snap%d_%scubical_pair%d-%d.h5' \
						  % (snap, coord, i1, i2)
					data = h5py.File(cat, 'r')
					x = np.array(data['px'])
					y = np.array(data['py'])
					f = np.array(data['flux'])
					frac = []
					frac_av = []
					for i in range(nsteps):
						oriented = (ymin<=y)&(y<=ymax)&(xs[i]<=x)&(x<xs[i+1])
						average = ((ymin<=x)&(x<=ymax)&(xs[i]<=y)&(y<xs[i+1]))|\
								  ((ymin<=y)&(y<=ymax)&(-xs[i]>x)&(x>=-xs[i+1]))|\
								  ((ymin<=x)&(x<=ymax)&(-xs[i]>y)&(y>=-xs[i+1]))
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

			good_theta = theta>20
			#distbins = [1, 3, 5, 10, 15, 20]
			mubins = [[-22, -19], [-19, -16.5]]#[-17, -19, -20, -21, -23]
			mrbins = [[-23, -20], [-20, -17]]#[-17, -19, -20, -21, -23]
			nwbins = [[0, 15], [10, 22], [18, 22], [15, 22]]
			rads = xs[1:]*pix2deg*3600

			means_av = []
			stds_av = []
			for xx in range(nsteps):
				means_av.append(np.nanmean(fracs_av[:, xx]))
				stds_av.append(np.nanstd(fracs_av[:, xx]))

			dmax = 5
			umax = -19
			ugood = ids[mu<=umax]
			others = [(ii1 in ugood) & (ii2 in ugood) for ii1, ii2 in zip(id1, id2)]
			others &= (dist<=dmax)
			_s = 'dmax%d_umax%d' % (dmax, umax)
			#for i in range(len(magbins)-1):
			#    idgood = ids[(magbins[i]>=mu) & (mu>=magbins[i+1])]

			props = [mr, mu, nw]
			bins = [mrbins, mubins, nwbins]
			names = ['mr', 'mu', 'nw']

			for j in range(len(bins)):
				for i in range(len(bins[j])):
					idgood = ids[(bins[j][i][0]<=props[j]) & (props[j]<=bins[j][i][1])]
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
					#plt.errorbar(rads, means_av, yerr=stds_av, color='gray')
					plt.xlabel('distance [arcsec]')
					plt.ylabel('LLS fraction')
					plt.legend(loc='best')
					#plt.semilogy()
					plt.ylim(0, 0.3)
					plt.xlim(2, 15.5)
					plt.savefig('../../EAGLE/LLS/plots/LLSfrac_pairs+circular_%s_%s%d_%d.png' % (_s, names[j], bins[j][i][0], bins[j][i][1]))
					plt.close()


if sbhist:

	def sb_roi(cat, xmin, xmax, ymin, ymax, zmin, zmax):
		data = h5py.File(cat, 'r')
		pz = np.array(data['pz'])
		zgood = (pz>=zmin) & (pz<=zmax)
		f = np.array(data['SB'])[zgood]
		px = np.array(data['px'])[zgood]
		py = np.array(data['py'])[zgood]
		right = (px>xmin) & (px<xmax) & (py>ymin) & (py<ymax)
		left = (px<-xmin) & (px>-xmax) & (py>ymin) & (py<ymax)
		down = (px>ymin) & (px<ymax) & (py<-xmin) & (py>-xmax)
		up = (px>ymin) & (px<ymax) & (py>xmin) & (py<xmax)
		rand = left | down | up
		fright = np.nanmean(f[right])
		frand = np.nanmean(f[rand])
		return fright, frand

	folder = '/net/astrogate/export/astrodata/gallegos/EAGLE/cats/'
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
		if freexmax: extraname = '.freexmax'#''_3-7asec'
		else: extraname = ''#''.d<6'
		from itertools import product
		zr = []#np.arange(-2,3)

		# trying just u1 and d
		pnameshist = [['u1', 'd']]#[['nw1', 'd'], ['nw1', 'd'], ['nw1', 'u1']]#, ['u2', 'd'], ['nw1', 'd']]  # [pnames[0], pnames[2]]
		#pnameshist = [['u1', 'u2'], ['u1', 'nw1']]#, ['nw1', 'u2']]  # [pnames[0], pnames[2]]

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

			if hbool:#not hbool:
				for snap in snaps:
					sbpeak = sbpeaks[str(snap)]
					red = redshifts[str(snap)]
					asec2kpc = asec2kpcs[str(snap)]

					sb_dimming = np.power(1+red, -4)
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
					#urange = [-21.405001, -20.125999, -19.551001, -19.306999, -19.080999, -18.822001, -18.549, -18.291,
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

						pairsname = '/net/astrogate/export/astrodata/gallegos/EAGLE/cats/SB_snap%d_%s/' \
									% (snap, coord)
						outname2 = '../../EAGLE/analysis/avsb_snap%d_%s%s.dat' % (snap, coord, extraname)

						pnames2 = ['theta%s' % coord, 'com_dist', 'ang%s' % coord, 'shear%s' % coord]

						for i in range(len(props2)): props2[i] = pairs[pnames2[i]]
						theta, dist, angle, shear = props2

						#props = [props1, props2]

						cool = []
						for i1, i2 in zip(id[0], id[1]):
							cool.append(os.path.isfile(pairsname + '%d-%d.h5' % (i1, i2)))
						cool &= (theta > thetamin) & (dist<20) & (dist>.5)#(dist<6)
						ncool = np.sum(cool)
						cool = np.where(cool)[0]

						if freexmax: xmax = (theta[cool] - 6) * asec2pix
						else: xmax = np.zeros(ncool) + xasecmax * asec2pix

						if 1:#ncool>0:
							condition = (not os.path.isfile(outname2)) or overwrite
							if condition:

								fout = open(outname2, 'w')
								fout.write('#id1 id2 f_right f_rand\n')
								foutz = {}
								for _z in zr:
									foutz[str(_z)] = open(outname2.replace('.dat', '.%d.dat' % _z), 'w')
									foutz[str(_z)].write('#id1 id2 f_right f_rand\n')

								for i in range(ncool):
									print 'Pair %d-%d, snap%d %s. %d of %d' % (id[0][cool[i]], id[1][cool[i]], snap, coord, i, ncool)
									cat = '%s%d-%d.h5' % (pairsname, id[0][cool[i]], id[1][cool[i]])
									data = h5py.File(cat, 'r')
									pz = np.array(data['pz'])
									zgood = (pz>=zmin) & (pz<=zmax)
									f = np.array(data['SB'])[zgood]
									px = np.array(data['px'])[zgood]
									py = np.array(data['py'])[zgood]
									right = (px>xmin) & (px<xmax[i]) & (py>ymin) & (py<ymax)
									left = (px<-xmin) & (px>-xmax[i]) & (py>ymin) & (py<ymax)
									down = (px>ymin) & (px<ymax) & (py<-xmin) & (py>-xmax[i])
									up = (px>ymin) & (px<ymax) & (py>xmin) & (py<xmax[i])
									rand = left | down | up
									f_ = np.nanmean(f[right])
									#f_ = np.nanmedian(f[right])
									#f_n = np.nanmean(right)
									f_rand = np.nanmean(f[rand])
									#f_rand = np.nanmedian(f[rand])
									#f_rand_n = np.nanmean(rand)
									fout.write('%d %d %f %f\n' % (id[0][cool[i]], id[1][cool[i]], f_, f_rand))
									for _z in zr:
										#remember f[zgood].....
										rightz = right & (pz==_z)
										randz = rand & (pz==_z)
										f_ = np.nanmean(f[rightz])
										#f_n = np.nanmean(rightz)
										f_rand = np.nanmean(f[randz])
										#f_rand_n = np.nanmean(randz)
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
										iii=0
										while not end:
											if iii >= len(id1): end = True
											else:
												if (id1[iii]==id[0][_j]) and (id2[iii]==id[1][_j]):
													id1.remove(id1[iii])
													id2.remove(id2[iii])
													cool2.append(_j)
													end = True
											iii += 1
									_p0, _p1 = [params[pn0][cool2], params[pn1][cool2]]
									f_right = _f[:, 2]
									#f_right_n = _f[:, 3]
									f_rand = _f[:, 3]
									#f_rand_n = _f[:, 5]
									print 'Combining %s with %s' % (pn0, pn1)

									outs = ['../../EAGLE/analysis/histograms/hist_fright_snap%d_%s_%s-%s%s.dat' % (snap, coord, pn0, pn1, extraname),
											'../../EAGLE/analysis/histograms/hist_frand_snap%d_%s_%s-%s%s.dat' % (snap, coord, pn0, pn1, extraname),
											#'../../EAGLE/analysis/histograms/hist_n_snap%d_%s_%s-%s%s.dat' % (snap, coord, pn0, pn1, extraname),
											'../../EAGLE/analysis/histograms/hist_nsubs_snap%d_%s_%s-%s%s.dat' % (snap, coord, pn0, pn1, extraname)]

									isfiles = True
									for out in outs: isfiles &= os.path.isfile(out)

									if isfiles and not overwrite and not histover:
										print 'h2d files exist'
										all[coord, pn0, pn1] = [np.loadtxt(out) for out in outs]
									else:
										if histover and isfiles: 'Overwriting h2d files!'
										if not isfiles: print 'No h2d files!'
										#n = np.zeros((n0, n1))
										nsubs = np.zeros((n0, n1))
										hfright = np.zeros((n0, n1))
										hfrand = np.zeros((n0, n1))

										for i in range(n0):
											for j in range(n1):
												coolh = (_p0 > r0[i]) & (_p0 <= r0[i + 1]) & (_p1 > r1[j]) & (_p1 <= r1[j + 1])
												#n[i, j] = np.nansum(f_right_n[coolh])
												nsubs[i, j] = np.nansum(coolh)#/n[i, j]
												hfright[i, j] = np.nanmean(f_right[coolh])
												hfrand[i, j] = np.nanmean(f_rand[coolh])# / np.nansum(f_rand_n[coolh])

										all[coord, pn0, pn1] = [hfright, hfrand, nsubs]

										for out, _h in zip(outs, all[coord, pn0, pn1]):
											print 'Saving file', out
											np.savetxt(out, _h)

									if 0:
										def doplot(h, name, labelname, vmin=None, vmax=None):
											# from matplotlib.colors import LogNorm
											if pn0 == 'd': plt.figure(figsize=(7, 5))
											if pn1 == 'd': plt.figure(figsize=(6, 7))
											plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin, vmax=vmax)
											dx = (r0[1] - r0[0])
											dy = (r1[1] - r1[0])
											plt.xticks(2 * np.arange(n0 / 2) + .5, np.arange(r0[0] + dx, r0[-1], 2 * dx))
											plt.yticks(2 * np.arange(n1 / 2) + .5, np.arange(r1[0] + dy, r1[-1], 2 * dy))
											plt.xlabel(plabels[pn0], fontsize=fsize)
											plt.ylabel(plabels[pn1], fontsize=fsize)
											cbar = plt.colorbar()
											cbar.set_label(labelname, size=fsize)
											plt.savefig(
												'../../EAGLE/analysis//plots/coords/hist_%s_snap%d_%s_%s-%s%s.png'
												% (name, snap, coord, pn0, pn1, extraname))
											plt.close()

										# Number of neighbors per bin
										doplot(n, 'N', u'$\mathrm{n}$')#, vmin=0, vmax=50)

										# sb
										doplot(hfright, 'SB', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=.4, vmax=.7)

										# sb rand
										doplot(hfrand, 'SB_rand', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=.4, vmax=.7)

										# sb - sb rand
										doplot((hfright - hfrand) * 100 / hfrand, 'SBminusSBrand',
											   u'$\mathrm{SB\,\,increase\,[percentage]}$', vmin=0, vmax=20)

										# sb wrt sbpeak
										doplot(hfright/sbpeak, 'SBvsUVB',
											   u'$\mathrm{SB/SB_{UVB}}$', vmin=0, vmax=.5)

						else: print 'No subcubes for %s coordinate' % coord

					for p in pnameshist:#list(product(pnames, pnames)):
						superall['%d'%snap] = all
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
							#plt.xticks(2 * np.arange(n0 / 2) + .5, np.arange(r0[0] + dx, r0[-1], 2 * dx))
							#plt.yticks(2 * np.arange(n1 / 2) + .5, np.arange(r1[0] + dy, r1[-1], 2 * dy))
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
						outs = ['../../EAGLE/analysis/histograms/hist_fright_snap%d_%s-%s%s.dat' % (snap, pn0, pn1, extraname),
								'../../EAGLE/analysis/histograms/hist_frand_snap%d_%s-%s%s.dat' % (snap, pn0, pn1, extraname),
								#'../../EAGLE/analysis/histograms/hist_n_snap%d_%s-%s%s.dat' % (snap, pn0, pn1, extraname),
								'../../EAGLE/analysis/histograms/hist_nsubs_snap%d_%s-%s%s.dat' % (snap, pn0, pn1, extraname)]
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
							#doplot(np.log10(n), 'N', u'$\mathrm{log(pixels)}$')#, vmin=0, vmax=40)
							vmin, vmax = [.04, .13]#[None, None]#
							doplot(hfright, 'SBvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=vmin, vmax=vmax, title=u'snapshot %d, z=%.3f' %(snap, red))

							# sb rand
							doplot(hfrand, 'SBrandvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=vmin, vmax=vmax, title=u'snapshot %d, z=%.3f' %(snap, red))

						if 0:
							# Number of neighbors per bin
							doplot(nsubs, 'Nsub', u'$\mathrm{Number\,\,of\,\,orientations}$')#, vmin=0, vmax=40)
							# sb
							doplot(hfright*sbpeak, 'SB', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=0, vmax=.14)

							doplot(hfright, 'SBvsUVB', u'$\mathrm{SB/SB_{UVB}}$', vmin=0, vmax=.14)

							# sb rand
							doplot(hfrand*sbpeak, 'SB_rand', u'$\mathrm{SB\,[10^{-20}cgs]}$', vmin=0, vmax=.14)

							# sb - sb rand
							doplot((hfright-hfrand)*100/hfrand, 'SBminusSBrand',
								   u'$\mathrm{SB\,\,increase\,[percentage]}$', vmin=-30, vmax=80)


							# SNR
							doplot(hfright*np.sqrt(nsubs), 'SNR', u'$\mathrm{Relative\,\,SNR}$', vmin=0.2, vmax=1.5)

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
		#if pn0 == 'd': plt.figure(figsize=(7, 5))
		plt.figure(figsize=(5, 6.2))#if pn1 == 'd':
		plt.imshow(np.array(h).T, interpolation='nearest', origin='low', vmin=vmin, vmax=vmax,
				   cmap=cm)
		if p0=='u1': plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
		else: plt.xticks([0, 1.5, 3, 4.5], [26, 28, 30, 32])
		#plt.xticks([0, 1.5, 3, 4.5], [-21, -19, -18, -17])
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

	#recheck all of this part!!

	if 0:
		print 'Final analysis'
		uvb = None
		huvb = None
		uvbrand = None
		iii=0
		for k in superall.keys():
			sa = superall[k]
			for kk in sa.keys():
				iii += 1
				_all = sa[kk]
				if huvb is None:
					#n = _all[2]
					nsub = _all[2]
					hfright = _all[0]*sbpeaks[k]#*nsub
					hfrand = _all[1]*sbpeaks[k]#*nsub
					huvb = _all[0]#*nsub
					uvb = np.nanmean(_all[0])
					uvbrand = np.nanmean(_all[1])
					nuvb = _all[2]
					print iii, k, kk, 'oriented %.3f random %.3f' % (uvb, uvbrand)

				else:
					#n += _all[2]
					_nsub = _all[2]
					nsub += _nsub
					_hright = _all[0]#*_all[2]
					hfright += _hright*sbpeaks[k]
					_hfrand = _all[1]#*_all[2]
					hfrand += _hfrand*sbpeaks[k]
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
			n = np.histogram2d(u1, dist, [urange, drange], weights=fright>0)[0]
			h = np.histogram2d(u1, dist, [urange, drange], weights=fright)[0]
			huvb[ss] = h/n
			hright[ss] = huvb[ss]*sbpeaks[ss]
			n = np.histogram2d(u1, dist, [urange, drange], weights=frand>0)[0]
			h = np.histogram2d(u1, dist, [urange, drange], weights=frand)[0]
			hrand[ss] = h*sbpeaks[ss]/n
			doplot(huvb[ss], 'SBvsUVB_snap%d' % snap, u'$\mathrm{SB/SB_{UVB}}$', 'u1', 'd', #urange, drange,
					xlabel=u'$\mathrm{U}}$', ylabel=u'$\mathrm{d\,[cMpc]}}$',
				   title='snapshot %d, z=%.3f' % (snap, redshifts[ss]), vmin=vmin, vmax=vmax)
			doplot(hright[ss], 'SB_snap%d' % snap, u'$\mathrm{SB\,[10^{-20}erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', 'u1', 'd', #urange, drange,
					xlabel=u'$\mathrm{U}}$', ylabel=u'$\mathrm{d\,[cMpc]}}$',
				   title='snapshot %d, z=%.3f' % (snap, redshifts[ss]), vmin=.04, vmax=.26, cticks=[.05, .1, .15, .2, .25])
		hfright = np.nanmean([hright['10'], hright['11'], hright['12']], 0)
		hfrand = np.nanmean([hrand['10'], hrand['11'], hrand['12']], 0)
		huvb = np.nanmean([huvb['10'], huvb['11'], huvb['12']], 0)

	#print 'total SB/SB_UVB: oriented %.4f random %.4f' % (uvb, uvbrand)

	noise1 = 1.# assuming 1e-20 e.g. UDF -> but 0.44 with the combination of 390 subcubes
	noise2 = noise1*np.sqrt(3)# assuming 1e-20 e.g. UDF

	udf = getdata('../../UDF/cats/lae_pairs.fits', 1)
	mosaic = getdata('../../mosaic/cats/lae_pairs-2018.fits', 1)

	d1 = udf['com_dist']
	t1 = udf['theta']
	f1 = udf['HST_F606W']
	d2 = mosaic['com_dist']
	f2 = mosaic['HST_F606W']
	t2 = mosaic['theta']

	g1 = (d1 > .5) & (d1 < 20) & (t1 > 16) & (f1>0)
	g2 = (d2 > .5) & (d2 < 20) & (t2 > 16) & (f2>0)

	dbins = np.arange(0, 21, 2)
	#fbins = [25.5921, 27.3024, 27.9615, 28.5004, 28.9428, 29.343,
	#         29.724,  30.0098,  30.5226,  31.3483,  32.046]
	fbins = [25.5921, 27.9615, 28.9428, 29.724,  30.5226,  32.046]

	p0 = 'u1'
	p1 = 'd'
	#n0 = len(pranges[p0])
	#n1 = len(pranges[p1])
	r0 = fbins
	r1 = dbins

	h1 = np.histogram2d(f1, d1, [fbins, dbins])[0]
	h2 = np.histogram2d(f2, d2, [fbins, dbins])[0]
	doplot(h1, 'Nsub_UDF', u'$\mathrm{Number\,\,of\,\,orientations}$', 'HST_F606W', p1, r0, r1)#, vmin=.2, vmax=1.5)
	doplot(h2, 'Nsub_mosaic', u'$\mathrm{Number\,\,of\,\,orientations}$', 'HST_F606W', p1, r0, r1)#, vmin=.2, vmax=1.5)#, vmin=.2, vmax=.5)

	SNR1=hfright*np.sqrt(h1)/noise1
	SNRrand1=hfrand*np.sqrt(h1)/noise1
	SNR2=hfright*np.sqrt(h2)/noise2
	SNRrand2=hfrand*np.sqrt(h2)/noise2
	vmin, vmax = [None, None]#[.3, 1.6]
	doplot(huvb, 'SBvsUVBtot', u'$\mathrm{SB/SB_{UVB}}$', 'u1', p1, vmin=.06, vmax=.14, xlabel=r'$\mathrm{U}$')#r0, r1,
	doplot(SNR1, 'SNR_UDF', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)
	doplot(SNRrand1, 'SNRrand_UDF', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)#, vmin=.2, vmax=.5)
	doplot((SNR1-SNRrand1)*100/SNRrand1, 'OrientedvsRandom_UDF', u'$\mathrm{SB\,\,increase\,[percentage]}$', 'HST_F606W', p1, r0, r1, vmin=0, vmax=40)#, vmin=.2, vmax=.5)
	doplot(SNR2, 'SNR_mosaic', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)#, vmin=.2, vmax=.5)
	doplot(SNRrand2, 'SNRrand_mosaic', u'$\mathrm{SNR}$', 'HST_F606W', p1, r0, r1, vmin=vmin, vmax=vmax)#, vmin=.2, vmax=.5)
	doplot((SNR2-SNRrand2)*100/SNRrand2, 'OrientedvsRandom_mosaic', u'$\mathrm{SB\,\,increase\,[percentage]}$', 'HST_F606W', p1, r0, r1, vmin=0, vmax=40)#, vmin=.2, vmax=.5)

	nh1 = np.nansum(h1)
	fconn1 = np.nansum((huvb*h1).reshape(-1))/nh1
	SNR_tot1 = np.nansum((SNR1*h1).reshape(-1))/nh1/noise1
	SNRrand_tot1 = np.nansum((SNRrand1*h1).reshape(-1))/nh1/noise1

	nh2 = np.nansum(h2)
	fconn2 = np.nansum((huvb*h2).reshape(-1))/nh2
	SNR_tot2 = np.nansum((SNR2*h2).reshape(-1))/nh2/noise2
	SNRrand_tot2 = np.nansum((SNRrand2*h2).reshape(-1))/nh2/noise2

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
	#ax.set_rmax(max_rad * 2)
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
	pairsname = '/net/astrogate/export/astrodata/gallegos/EAGLE/LLScats/cubical_pair'
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
		angbin = np.arange(0., 2*np.pi+1e-10, 2*delta_ang)
		nbin = len(angbin) - 1
		nbins.append(nbin)
		angbins.append(angbin)
		LLSrads.append(np.zeros(nbin))
		delta_angs.append(delta_ang)(px>xmin) & (px<xmax) & (py>ymin) & (py<ymax)

	filedat = '../../EAGLE/LLS/llscumang.dat'
	overwrite = True

	if not os.path.isfile(filedat) or overwrite:
		x = np.array([])
		y = np.array([])
		f = np.array([])
		for i in np.where(cool)[0]:
			print 'Pair %d-%d' % (id1[i], id2[i])
			cat = '/net/astrogate/export/astrodata/gallegos/EAGLE/LLScats/cubical_pair%d-%d.h5' % (id1[i], id2[i])
			data = h5py.File(cat, 'r')
			x = np.concatenate((x, np.array(data['px'])))
			y = np.concatenate((y, np.array(data['py'])))
			f = np.concatenate((f, np.array(data['flux'])))

		for j in range(nrads - 1):
			d2 = x**2+y**2
			ang = (np.arctan2(y, x)+np.pi-delta_angs[j]) % (2*np.pi)

			for a in range(nbins[j]):
				inside = (ang > angbins[j][a]) & (ang <= angbins[j][a + 1]) & (d2 < rads[j+1] ** 2) & (d2 < rads[j] ** 2)
				LLSrads[j][a] = np.nanmean(f[inside])

		fout = open(filedat, 'w')
		for line in LLSrads:
			for row in line:
				if np.isnan(row): fout.write('0. ')
				else: fout.write('%f ' % float(row))
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
		if j == 0: ec = "none"
		else: ec = "black"
		bars = ax.bar(angbins[j][:-1]+delta_angs[j], np.zeros(len(angbins[j])-1)+rad2[j], width=2*delta_angs[j],
					  bottom=rad3[j], edgecolor=ec, linewidth=.5)
		mlr = max(LLSrads[j])
		for bar, lr in zip(bars, LLSrads[j]): bar.set_facecolor(plt.cm.jet(1.6*lr))

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
				rad_asec = mat[0]/5.
				llsin = mat[1]
				stdin = mat[2]
				llsout = mat[3]
				stdout = mat[4]
			else:
				fits = getdata('/net/astrogate/export/astrodata/gallegos/EAGLE/logNHI_snap%d_%s.fits' % (snap, coord))
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
					print 'Inside %d arcsec' % (r/5)
					for i, j, gal in zip(xpos, ypos, ids):
						#print 'Gal %d x %d y %d' % (gal, i, j)
						inside |= (_x-i)**2 + (_y-j)**2 < r**2

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
			#plt.semilogy()
			plt.xlabel(u'$\mathrm{distance\,[arcsec]}$', fontsize=18)
			plt.ylabel(u'$f_{\mathrm{LLS}}(\mathrm{r<R})}$', fontsize=18)
			plt.errorbar(rad_asec, llsin, yerr=stdin, color='blue')
			plt.errorbar(rad_asec, llsout, yerr=stdout, color='red')
			plt.savefig('../../EAGLE/LLS/plots/LLSfracs_snap%d_%s.png' % (snap, coord))
			#plt.show()
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

	ds = [2, 5, 10, 20]
	zw0 = 2
	folder = '/net/eriu/export/data1/Dropbox/MUSE/EAGLE/stacks/'
	f2 = '/net/eriu/export/data1/Dropbox/MUSE/all/stacks/'
	hdu = PrimaryHDU()
	imtype = 'mean'#'median'#
	mask = False
	smask = '.mask' * mask
	stype = ('.'+imtype)*(imtype=='median')
	for d in ds:
		all = []
		ncores = d
		for snap in snaps:
			stacks = []
			for coord in coordnames:#['x']:#
				print 'd %d snap %d coord %s' % (d, snap, coord)
				sname = folder+'snap%d_%s_d0-%d_pd16-2000_d5th0.5-20.0%s%s/stack.fits' % (snap, coord, d, stype, smask)
				if not os.path.isfile(sname) or overwrite: os.system('mpirun -n %d python stackEAGLE.py -snap %d '
						'-coord %s -dmax %d -overwrite True -imtype %s -mask %s' %
				                                                     (ncores, snap, coord, d, imtype, mask))
				_stack = getdata(sname)
				zl, yl, xl = _stack.shape
				stacks.append(_stack[:, :, :yl])
			stack = np.nanmean(stacks, 0)
			zl, yl, xl = stack.shape
			hdu.data = stack
			fs = folder+'snap%d_d0-%d_pd16-2000_d5th0.5-4.0%s%s/' % (snap, d, stype, smask)
			if not os.path.isdir(fs): glob.os.makedirs(fs)
			hdu.writeto(fs+'stack.fits', clobber=True)
			stackim = np.nanmean(stack[zl/2-zw0:zl/2+zw0+1, :, :], 0)
			hdu.data = stackim
			hdu.writeto(fs+'stack.IM.fits', clobber=True)

			all.append(stack*sbpeaks[str(snap)])

		stack = np.nanmean(all, 0)
		zl, yl, xl = stack.shape
		hdu.data = stack
		fs = folder+'d0-%d_pd16-2000_nw0-1000%s%s/' % (d, stype, smask)
		if not os.path.isdir(fs): glob.os.makedirs(fs)
		hdu.writeto(fs+'/stack.fits', clobber=True)
		stackim = np.nansum(stack[zl/2-zw0:zl/2+zw0+1, :, :], 0)
		hdu.data = stackim
		hdu.writeto(fs+'stack.IM.fits', clobber=True)
		fn = f2+'d0-%d_los0-20_pd16-2000_z2.9-4.0_l0-2000_n1-1000_v0-10000/' % d
		nname = fn+'stack.NPIX.IM.fits'
		if not os.path.isfile(nname): os.system('python stack.py -dmax %d -overwrite True' % d)
		npix = getdata(nname)
		yl, xl = npix.shape
		noise = 10. # assuming 1e-19 cgs noise
		snr = stackim[:, :xl]*np.sqrt(npix)/noise/(2*zw0+1.) # for the image stack I am doing a sum therefore the noise does not decrease, so I added the factor 1/(2*zw0+1.)
		hdu.data = snr
		hdu.writeto(fs+'stack.SNR.fits', clobber=True)

if 0:
	f = '../../EAGLE/stacks/d0-2_pd16-2000_nw0-1000/'
	fits = getdata(f+'stack.fits')
	zl, yl, xl = fits.shape
	ff = np.mean(fits[zl/2-1: zl/2, :, :], 0)
	hdu = PrimaryHDU()
	hdu.data = ff
	hdu.writeto(f+'stack.IM.2.fits', clobber=True)
	astroim(f + 'stack.IM.2.fits', smooth=True, saveim=True, show=False,
			cbfrac=.08, pad=.006, dfig=(8, 8), contours=True,
			x0=-25, y0=0, imout='../../EAGLE/plots/stack-d2-snap11.png', std=None, scale=False, regcolor=None,
			scb_label=r'$\rm{f_{conn}}$', xticks=[-16, -12, -8, -4, 0, 4, 8, 12, 16], yticks=[-12, -8, -4, 0, 4, 8, 12],
			title='', vmin=0, vmax=.25, gray=True, sb=False,
			xmin=-50, xmax=0, ymin=-20, ymax=20, pcolor='white', highq=False, dpi=200, nbins=5)
	
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
	fin = '/net/astrogate/export/astrodata/gallegos/'
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
	c = (2*zw+1) * 7.8125
	asec2kpc = 7.47
	folder = '../../EAGLE/stacks/snap11_x_d0-20_pd16-2000_nw0.000-1.000/'
	folder2 = '../../all/stacks/d0-20_pd16-80_z2.9-4.0_l1-2000_n1-100/'
	stack = getdata(folder+'stack.fits')
	stack2 = getdata(folder2+'stack.fits')
	srand = getdata(folder2+'random_stack.fits')
	frs = []
	for i in range(200):
		hdur = getdata(folder2 + 'randoms/random_stack.%d.fits' % i)
		frs.append(hdur)
	zl, yl, xl = stack.shape
	zl2, yl2, xl2 = stack2.shape
	print xl, xl2, yl, yl2, zl, zl2
	xmin = yl/2
	t = np.concatenate([np.arange(xmin, xmin + 12, 1), [xmin+13, xmin+18, xmin+27, xmin+41]]).astype('int')
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
			yis = [yymin, t[i] - 1 - dw, yymin, yl - t[i + 1]-dw]
			yfs = [yymax, t[i + 1] + dw, yymax, yl - t[i] + 1+dw]
			xis = [t[i] - 1 - dw, yymin, yl - t[i + 1]-dw, yymin]
			xfs = [t[i + 1] + dw, yymax, yl - t[i] + 1+dw, yymax]
			sb.append((t[i]+t[i + 1])*.5)
			#for yi, yf, xi, xf, ds in zip(yis, yfs, xis, xfs, dir):
			yi, yf, xi, xf, ds = zip(yis, yfs, xis, xfs, dir)[0]
			if 1:
				sb.append(np.nanmean(stack[zi:zf, yi:yf, xi:xf]))
				sb.append(np.nanmean(stack2[zi2:zf2, yi:yf, xi:xf])*c)
				sb.append(np.nanmean(srand[zi2:zf2, yi:yf, xi:xf])*c)
				sb.append(np.nanstd([np.nanmean(ff[zi:zf, yi:yf, xi:xf])*c for ff in frs]))
	lx = len(t) - 1
	#ld = len(dir)*3 + 1
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
	#plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=fsize + 2)
	yticks = [0, 2, 4, 6, 8, 10]
	plt.yticks(yticks, yticks, fontsize=fsize)
	plt.minorticks_off()
	plt.twiny(ax=None)
	plt.loglog()
	plt.xlim([x1, x2])
	plt.xlabel(r'$r_p\,\rm{[pkpc]}$', fontsize=fsize)
	kpc = np.arange(15, 90, 15)
	plt.xticks(kpc/asec2kpc, kpc.astype('int'), fontsize=fsize)
	plt.ylim([y1, y2])
	plt.plot((2, x2), (0, 0), '--', label='', color="gray", linewidth=width + 1)
	if hm: plt.plot((x1, x2), (hmSB, hmSB), '--', label=r"LLS Fluorescence from HM12", color="green", linewidth=width + 1)
	
	x = (t[:-1] + (t[1:] - t[:-1]) / 2. - yl / 2) * .4 # 0.4 arcsec per pixel
	
	#plt.plot(x, sb[:, 2], label='MUSE oriented', lw=width, color='black')
	#plt.scatter(x, sb[:, 2], s=width * 10, color='black')
	sb = np.array(sb)
	sb[sb<0] = 1e-10
	plt.plot(x, sb[:, 3], label='SB MUSE', lw=width, color='dodgerblue')
	plt.scatter(x, sb[:, 3], s=width * 10, color='dodgerblue')
	ymin = sb[:, 3]-sigmas*sb[:, 4]
	ymin[ymin < 0] = 1e-10
	plt.fill_between(x, ymin, sb[:, 3]+sigmas*sb[:, 4], facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')

	sbuvmean = 0.8
	sbuvup = 1.2
	sbuvdown = 0.6
	plt.plot(x, sb[:, 1]*sbuvmean, label='SB EAGLE+HM12', lw=width, color='red')
	plt.scatter(x, sb[:, 1]*sbuvmean, s=width * 10, color='red')
	plt.fill_between(x, sb[:, 1]*sbuvdown, sb[:, 1]*sbuvup, facecolor='red', alpha=0.3, lw=0, edgecolor='none')
	plt.minorticks_off()
	gamma = 3.16e13#ionizing photon rate from the galaxy /1e40
	kpc2cm = 30.86 #kpc2cm /1e20
	zmean = 3.5
	sb_diff = sb[:, 3]-sb[:, 1]*sbuvmean
	E_lya = 1.64e-11
	asec2rad_sqd = 2.35e-11
	fesc = 4*np.pi*(1+zmean)**4*(x*asec2kpc*kpc2cm)**2*sb_diff*1e-20/E_lya/.6/gamma/asec2rad_sqd
	gamma_std = 3e13
	obs_std = sb[:, 4]
	sim_std = 0.3
	sb_std = np.sqrt(np.power(obs_std, 2)+sim_std**2)
	tot_std = np.sqrt((gamma_std/gamma)**2+(sb_std/sb_diff)**2)*fesc
	print 'fesccccc', fesc
	print 'sb std', sb_std
	print 'tot std', tot_std
	fesc[fesc<0]= 1e-10
	plt.plot(x, fesc, label='Escape fraction', lw=width, color='black')
	ymin = fesc-tot_std
	ymin[ymin < 0] = 1e-10
	plt.fill_between(x, ymin, fesc+tot_std, facecolor='gray', alpha=0.3, lw=0, edgecolor='none')

	plt.legend(fontsize=fsize, loc='best')  # (3.5,2))

	plt.savefig('../../analysis/muse-vs-eagle.png')

	