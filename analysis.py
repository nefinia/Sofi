#!/usr/bin/env python
__author__ = 'gallegos'
import glob
# import h5py
# import eagleSqlTools
# import matplotlib.pyplot as plt
from sys import argv

import numpy as np
import scipy.interpolate as inter
# from matplotlib.cm import get_cmap
# from matplotlib.colors import LinearSegmentedColormap
from pyfits import getdata, PrimaryHDU

import params
from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

coordnames = rdarg(argv, 'coord', list, ['x', 'y', 'z'], str)
fitcat = rdarg(argv, 'fitcat', list, ['HDFS', 'UDF', 'mosaic'], str)
snaps = rdarg(argv, 'snap', list, [8, 9, 10, 11, 12, 13, 14], int)
scodes = rdarg(argv, 'scodes', str, '/net/abnoba/scratch2/gallegos/Research/MUSE/codes/Sofi/')
overwrite = rdarg(argv, 'overwrite', bool, False)
cubecorr = rdarg(argv, 'cubecorr', bool, False)
caseA = rdarg(argv, 'caseA', bool, True)
scase = '_CaseA' * caseA
halfplot = rdarg(argv, 'halfplot', bool, False)
extraname = rdarg(argv, 'extraname', str, '')#'HM12')#
galcov = rdarg(argv, 'galcov', bool, False)
laemask = rdarg(argv, 'laemask', bool, False)
histover = rdarg(argv, 'histover', bool, False)
LLScubes = rdarg(argv, 'LLS', bool, False)
sql = rdarg(argv, 'sql', bool, False)
sphericalLLS = rdarg(argv, 'sphericalLLS', bool, False)
pairLLS = rdarg(argv, 'pairLLS', bool, False)
circlehist = rdarg(argv, 'circlehist', bool, False)
h2d = rdarg(argv, 'h2d', bool, False)
kde = rdarg(argv, 'kde', bool, False)
nhiprof = rdarg(argv, 'nhi', bool, False)
minres = rdarg(argv, 'minres', int, 512)
maxres = rdarg(argv, 'maxres', int, 4096)
model = rdarg(argv, 'model', str, 'HM01')#'HM12')#
npref = rdarg(argv, 'npref', int, 12)
rad = rdarg(argv, 'rad', int, 3)  # arcsec
radprof = rdarg(argv, 'radprof', bool, False)
sbhist = rdarg(argv, 'sbhist', bool, False)
_ssthr = rdarg(argv, 'ssthr', float, None)#6.73e-3 or 1e10
sbprof = rdarg(argv, 'sbprof', bool, False)
snr = rdarg(argv, 'snr', bool, False)
superstack = rdarg(argv, 'superstack', bool, False)
type = rdarg(argv, 'type', str, 'NHI') #NHtot, f_NHI, NHII
lutzmodel = rdarg(argv, 'lutzmodel', bool, False)
unique = rdarg(argv, 'unique', bool, False)
mask = rdarg(argv, 'mask', bool, False)
do_delaunay = rdarg(argv, 'delaunay', bool, False)
temperature = rdarg(argv, 'temperature', bool, False)
zw = rdarg(argv, 'zw', int, 5)  # arcsec
do_not = 0
if do_not: pass


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
#cm = LinearSegmentedColormap.from_list('sofi', colors, N=50)


gamma_bkg = params.gamma_bkg[model]
heat_bkg = params.heat_bkg[model]
snames = params.snames
redshifts = params.redshifts
asec2kpcs = params.asec2kpcs
ssthr = params.nhSS[model]
if _ssthr is not None:
	ssthr = {}
	for snap in snaps:
		ssthr[str(snap)] = _ssthr

dz = params.dz
sbpeaks = params.sbpeaks
zlens = params.zlens
dndz = params.dndz
fLLScorr = params.fLLScorr

nhi_fit = params.nhi_fit[model]

lcube = maxres # 4096
coml = 25  # cMpc
com2pix = 163.84  # lcube/coml
kpc2pix = lcube / float(coml * 1e3)
rads = np.array([0, 2, 4, 8, 12, 20, 30, 50, 100, 200])

lognhi = [10, 11, 12, 13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21, 22, 23, 30]
nhi = np.power(10, lognhi)
sb = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96, 0.98, 1, 1, 1, 1, 1, 1]

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
			flux2sb = zl / float(coml)
			z, y, x = np.ogrid[1: zl + 1, 1: yl + 1, 1: xl + 1]
			T = []  # 20, 40 and 80 kpc
			for i in range(ngal):
				print 'id', ids[i], '%d/%d' % (i, ngal)
				cool = ((.05 * flux2sb) < (z - zc[i] * flux2sb) ** 2 + (y - yc[i] * flux2sb) ** 2 + (x - xc[i] * flux2sb) ** 2) & \
					   ((z - zc[i] * flux2sb) ** 2 + (y - yc[i] * flux2sb) ** 2 + (x - xc[i] * flux2sb) ** 2 < (.08 * flux2sb))
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

if galcov:
	zw0 = 2
	folder = '/net/abnoba/scratch2/gallegos/Research/MUSE/'

	cov = {}
	rads = [[0.1, 2], [2, 4], [4, 8], [8, 12], [12, 16], [16, 20], [20, 30]]#, [30, 60], [60, 100]]
	nr = len(rads)
	for ff in ['HDFS']:#, 'mosaic']:
		pairs = getdata('../../%s/cats/lae_pairs.fits' % ff, 1)
		theta = pairs['theta']
		cool = (abs(pairs['shear']) <= zw0)# & (pairs['redshift1']<=3) & (pairs['redshift2']<=3)
		theta = theta[cool]
		cov[ff] = np.zeros(nr)
		for i in range(nr):
			area = np.pi*(rads[i][1]**2-rads[i][0]**2)
			inside = (theta > rads[i][0]) & (theta < rads[i][1])
			_cov = np.sum(inside)/area
			print 'catalog', ff, 'rads', rads[i], 'fgal %.3f area %.3e arcsec^2' % (_cov, area)
			cov[ff][i] = _cov

	ce = {}
	for coord in coordnames:
		for snap in snaps:
			k = '%d_%s' % (snap, coord)
			ce[k] = np.zeros(nr)
			cs = ['x', 'y', 'z']
			cs.remove(coord)
			pairs = getdata('../../EAGLE/cats/lae_pairs_snap11.fits', 1)
			theta = pairs['theta%s' % coord]
			cool = abs(pairs['shear%s' % coord]) <= zw0
			for i in range(nr):
				area = np.pi*(rads[i][1]**2-rads[i][0]**2)
				inside = (theta > rads[i][0]) & (theta < rads[i][1])
				_cov = np.sum(inside)/area
				print k, 'r', rads[i], 'fgal %.3f area %.3e' % (_cov, area)
				ce[k][i] = _cov
	cov['EAGLE'] = np.nanmean([ce[_ce] for _ce in ce.iterkeys()], 0)
	plt.figure()
	rm = np.mean(rads, 1)
	import itertools
	marker = itertools.cycle(('*', 'X', '^', 'h', 'd'))

	for k in cov.iterkeys():
		plt.plot(rm, cov[k])
		plt.scatter(rm, cov[k], label=k, marker=marker.next(), s=30)
	plt.legend()
	plt.semilogy()
	plt.xlabel(r'distance [arcsec]')
	plt.ylabel(r'$\mathrm{f_{gal}\,[arcsec^{-2}\,\Delta v^{-1}}]$')
	plt.savefig('../../Figures/Gal_covfrac.jpg')

if cubecorr:
	r = 20
	zw0 = 3
	folder = '/net/galaxy-data/export/galaxydata/gallegos/'
	fitcat = ['mosaic', 'UDF', 'HDFS']
	names = ['DATACUBE_UDF-MOSAIC', 'DATACUBE_UDF-10', 'DATACUBE-HDFS-1.35-PROPVAR']
	mnames = ['DATACUBE_UDF-MOSAIC.IM.Objects_Id', 'DATACUBE_UDF-10.IM.Objects_Id',
			  'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id']
	gnames = ['mosaic.galmask_%darcsec' % r, 'UDF.galmask_%darcsec' % r,
			  'HDFS.galmask_%darcsec' % r]
	ext = '.fits'
	do_csub = True
	cut = 20
	dim = (1, 2)
	for ft, name, mname, gname in [zip(fitcat, names, mnames, gnames)[0]]:
		print ft
		spec = open('%sspec.dat' % ft, 'w')
		spec.write('#layer mean_1l mean_%dl\n' % (zw0*2+1))
		fin = '%s%s/%s' % (folder, ft, name)
		if do_csub: fin += '.csub'
		fout = '%s.corr' % fin
		fmask = '%s%s/%s' % (folder, ft, mname)
		fgmask = '%s%s/%s' % (folder, ft, gname)

		print 'get original %scube' % ('csub ' * do_csub)
		cube = getdata(fin+ext)
		cube[cube == -999] = np.nan
		z, y, x = cube.shape
		_f = np.copy(cube)

		print 'mask continuum objects'
		mask = getdata(fmask+ext)
		bad = mask > 0
		_f[:, bad[0]] = np.nan

		print 'mask 3d objects'
		mask = getdata(fmask.replace('.IM', '')+ext)
		bad = mask > 0
		_f[bad] = np.nan

		print 'mask LAE halos'
		gmask = getdata(fgmask + ext)
		#gmask2 = np.nansum(gmask, 0)
		halos = gmask > 0
		_f[halos] = np.nan

		print 'find high bkg layers'

		smooth = False
		print 'Calculate mean bkg for non high layers, not smoothed!!'
		m = np.nanmean(_f[:, cut: y - cut, cut: x - cut], dim)
		for zz in range(z):
			if smooth:
				zmin = max(0, zz-zw0)
				zmax = min(z, zz+zw0)
				m5 = np.nanmean(_f[zmin: zmax+1, cut: y - cut, cut: x - cut])
			else:
				m5 = np.nanmean(_f[zz, cut: y - cut, cut: x - cut])
			print '%d %f %f' % (zz, m[zz], m5)
			spec.write('%d %f %f\n' % (zz, m[zz], m5))
			cube[zz, :, :] -= m5
		mm = np.nanmean(m)
		print 'mean sb value', mm, 'sum of all layers with their non masked pixels', mm*z
		std = np.nanstd(m)
		high = abs(m-mm) > 4*std
		print np.sum(high), 'high bkg layers of', z, 'std', std
		cube[high, :, :] = np.nan

		hdu = PrimaryHDU()
		hdu.data = cube
		hdu.writeto(fout+ext, clobber=True)

if laemask:
	import astropy.io.fits as fits
	from astropy import wcs

	def red2pix(z):
		"""Convert redshift to pixel"""
		lya_rest = 1215.67
		l0 = 4750.
		pixsize = 1.25  # 1 pixel = 1.25 Angstrom
		return (lya_rest * (1 + z) - l0) / pixsize + 1
	fin = '../../'
	fout = '/net/eos/scratch/gallegos/'#'/net/galaxy-data/export/galaxydata/gallegos/'
	ext = '.fits'
	fitcat = ['mosaic', 'UDF', 'HDFS', 'MXDF']
	names = ['DATACUBE_UDF-MOSAIC', 'DATACUBE_UDF-10', 'DATACUBE-HDFS-1.35-PROPVAR',
			 'DATACUBE_MXDF_ZAP_COR']
	hnames = ['DATACUBE_UDF-MOSAIC', 'UDF10.z1300', 'DATACUBE-HDFS-1.35-PROPVAR',
			  'DATACUBE_MXDF_ZAP_COR']
	mnames = ['DATACUBE_UDF-MOSAIC.IM.Objects_Id', 'UDF10.z1300.IM.Objects_Id',
			  'DATACUBE-HDFS-1.35-PROPVAR.IM.Objects_Id',
			  'DATACUBE_MXDF_ZAP_COR.IM.Objects_Id']
	gnames = ['%s.galmask_%darcsec' % (fc, rad) for fc in fitcat]
	for i in range(len(fitcat)):
		ft = fitcat[i]
		name = names[i]
		cat = getdata(fin + '%s/cats/laes.fits' % ft, 1)
		cubename = '/net/galaxy-data/export/galaxydata/gallegos/%s/%s%s' % (ft, name, ext)
		_cube = getdata(cubename)
		zl, yl, xl = _cube.shape
		try:
			print 'Reading x,y,x coordinates!'
			xs = np.round(cat['x'][good]).astype(int)
			ys = np.round(cat['y'][good]).astype(int)
			zs = np.round(cat['z'][good]).astype(int)
		# & (xs<xlim) & (ys>0) & (ys<ylim) (zs>0) & (zs<zlim)
		except:
			cubename = '/net/galaxy-data/export/galaxydata/gallegos/%s/%s%s' % (ft, name, ext)
			_cube = getdata(cubename)
			zl, yl, xl = _cube.shape
			hname = '/net/galaxy-data/export/galaxydata/gallegos/%s/%s%s' % (ft, hnames[i], ext)
			ttt, header_data_cube = getdata(hname, 0, header=True)
			ra = cat['RA']
			dec = cat['DEC']
			try: redshift = cat['z_muse']
			except: redshift = cat['redshift']
			ids = cat['ID']
			# Removing COMMENT key to avoid problems reading non-ascii characters
			cards = header_data_cube.cards
			bad = ['COMMENT' == b[0] for b in cards]
			for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
			hdulist = fits.open(cubename)
			w = wcs.WCS(header_data_cube, hdulist)
			x, y, z = np.round(w.all_world2pix(ra, dec, [1] * len(ra), 1)).astype(int)
			z = np.round(red2pix(redshift)).astype(int)
		_y, _x = np.ogrid[0:yl, 0:xl]
		cube = np.zeros((zl, yl, xl))

		rpix = rad*5
		for i in range(len(z)):
			gal = (x[i]-_x)**2+(y[i]-_y)**2 < rpix**2
			print i, np.sum(gal)
			zmin = max(0, z[i] - zw)
			zmax = min(zl, z[i] + zw + 1)
			cube[zmin: zmax, gal] += ids[i]
		hdu = PrimaryHDU()
		hdu.data = cube
		hdu.writeto(fout+'%s/%s.galmask_%darcsec_zw%d.fits' % (ft, ft, rad, zw), clobber=True)
