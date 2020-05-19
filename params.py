#!/usr/bin/env python
__author__ = 'gallegos'

lya_rest = 1215.67

redshifts = {6: 5.971, 7: 5.487, 8: 5.037, 9: 4.485, 10: 3.984, 11: 3.528, 12: 3.017, 13: 2.478,
             14: 2.237, 15: 2.012}

gamma_bkg = {}
heat_bkg = {}
sgamma = {}
sheat = {}
nhSS = {}

gamma_bkg['HM12'] = {6: [2.50E-13, 1.69E-13, 8.15E-19],
                     7: [3.50E-13, 2.14E-13, 3.23E-18],
                     8: [4.26E-13, 2.45E-13, 1.33E-17],
                     9: [4.86E-13, 2.81E-13, 5.60E-17],
                     10: [5.71E-13, 3.28E-13, 2.31E-16],
                     11: [6.77E-13, 3.89E-13, 8.70E-16],
                     12: [7.66E-13, 4.42E-13, 2.10E-15],
                     13: [9.50E-13, 5.55E-13, 9.06E-15],
                     14: [9.64E-13, 5.67E-13, 1.13E-14],
                     15: [0.944E-12, 0.557E-12, 0.128E-13]}


heat_bkg['HM12'] = {6: [1.11E-12, 1.21E-12, 4.51E-17],
                    7: [3.60E-13, 2.14E-13, 3.23E-18],
                    8: [4.26E-13, 2.45E-13, 1.33E-17],
                    9: [4.86E-13, 2.81E-13, 5.60E-17],
                    10: [2.27E-12, 2.18E-12, 7.24E-15],
                    11: [2.68E-12, 2.62E-12, 2.33E-14],
                    12: [3.02E-12, 3.05E-12, 5.01E-14],
                    13: [3.75E-12, 4.22E-12, 1.78E-13],
                    14: [3.81E-12, 4.42E-12, 2.18E-13],
                    15: [0.374E-11, 0.444E-11, 0.244E-12]}

# from rahmati+13a
gamma_bkg['HM01'] = {6: [2.62E-13, 1.69E-13, 8.15E-19],
                     7: [3.60E-13, 2.14E-13, 3.23E-18],
                     8: [5.43E-13, 2.45E-13, 1.33E-17],
                     9: [6.7E-13, 2.81E-13, 5.60E-17],
                     10: [7.92E-13, 3.28E-13, 2.31E-16],
                     11: [.96E-12, 3.89E-13, 8.70E-16],
                     12: [1.16E-12, 4.42E-13, 2.10E-15],
                     13: [1.32E-12, 5.55E-13, 9.06E-15],
                     14: [1.41E-12, 5.67E-13, 1.13E-14]}


if 0:
	"old"
	gamma_bkg['HM01'] = {6: [2.62E-13, 1.69E-13, 8.15E-19],
	                     7: [3.60E-13, 2.14E-13, 3.23E-18],
	                     8: [.6E-12, 2.45E-13, 1.33E-17],
	                     9: [.7E-12, 2.81E-13, 5.60E-17],
	                     10: [.9E-12, 3.28E-13, 2.31E-16],
	                     11: [1E-12, 3.89E-13, 8.70E-16],
	                     12: [1.5E-12, 4.42E-13, 2.10E-15],
	                     13: [1.7E-12, 5.55E-13, 9.06E-15],
	                     14: [1.9E-12, 5.67E-13, 1.13E-14]}

heat_bkg['HM01'] = {6: [1.11E-12, 1.21E-12, 4.51E-17],
                    7: [1.48E-12, 1.46E-12, 1.56E-16],
                    8: [2.8E-12, 2.45E-13, 1.33E-17],
                    9: [3E-12, 2.81E-13, 5.60E-17],
                    10: [3.5E-12, 2.18E-12, 7.24E-15],
                    11: [4E-12, 2.62E-12, 2.33E-14],
                    12: [5E-12, 3.05E-12, 5.01E-14],
                    13: [5.5E-12, 4.22E-12, 1.78E-13],
                    14: [6E-12, 4.42E-12, 2.18E-13]}


nhSS['HM12'] = {6: 3.06e-3, 7: 3.75e-3, 8: 4.5e-3, 9: 5.2e-3, 10: 5.8e-3, 11: 6.59e-3, 12: 7.4e-3,
                13: 8.08e-3,
                14: 8.42e-3, 15: 8.7e-3}  # snap 6 and 7 using fit -1.43*z+11.6

nhSS['HM01'] = {6: 3.06e-3, 7: 3.75e-3, 8: 4.5e-3, 9: 5.2e-3, 10: 5.8e-3, 11: 6.59e-3, 12: 7.4e-3,
                13: 8.08e-3,
                14: 8.42e-3, 15: 8.7e-3}  # snap 6 and 7 using fit -1.43*z+11.6

nhSS_fixed = 6.73e-3

snames = {6: '006_z005p971',
          7: '007_z005p487', 8: '008_z005p037', 9: '009_z004p485', 10: '010_z003p984', 11: '011_z003p528',
          12: '012_z003p017', 13: '013_z002p478', 14: '014_z002p237', 15: '015_z002p012'}

asec2kpcs = {6: 5.851, 7: 6.092, 8: 6.393, 9: 6.754, 10: 7.842, 11: 7.449, 12: 7.108, 13: 8.241,
             14: 8.396, 15: 8.516}
dz = {6: .0577, 7: 0.0518, 8: 0.0462, 9: 0.0406, 10: 0.035, 11: .03, 12: .0255, 13: .0208, 14: .01695}
dz0 = 1.028e-3
# sbpeaks = [2.396]
# sbpeaks = {10: 1.011, 11: 1.484, 12:2.396} # for an aperture of 1 asec^2!!! I am converting to flux later on
# sbpeaks calculated with python plots.py -hm True
sbpeaks = {6: .0933, 7: .165, 8: .26, 9: .43, 10: .73, 11: 1.27, 12: 2.49, 13: 5.14, 14: 6.97,
           15: 9.12}  # for an aperture of 1 asec^2!!! I am converting to flux later on
zlens = {6: 56, 7: 50, 8: 45, 9: 39, 10: 34, 11: 29, 12: 25, 13: 20, 14: 18, 15: 16}
# dndz = {7: 6., 8: 5., 9: 4.,10: 3., 11: 2.3, 12: 2., 13: 1.3, 14: 1., 15: .5} #old version, high extrapolation w/r to Prochaska+10
# dndz = {7: 10., 8: 7., 9: 4.4,10: 3.2, 11: 2, 12: 1.5, 13: 1, 14: .9, 15: .8} #Prochaska+10 figure 17 + Prochaska+10 fit
dndz_prochaska = []
dndz_prochaska = {7: 9.83, 8: 6.81, 9: 4.18, 10: 2.56, 11: 1.57, 12: .85, 13: .41, 14: .28,
                  15: .2}  # from Prochaska+10 fit 1.9*((1+z)/(1+3.7))**5.1 for tau>2

dndz = {6: 3.75, 7: 3.65, 8: 3.23, 9: 2.75, 10: 2.33, 11: 1.98, 12: 1.62, 13: 1.27, 14: 1.12,
        15: .99}  # from Crighton+19 Figure 8 fit for tau>2 and adjusted to tau>1 (10% more)

dndz_eagle = {6: 5.52, 7: 4.18, 8: 3.34, 9: 2.52, 10: 2.1, 11: 1.64, 12: 1.42, 13: 0.94, 14: 0.86}

def l(z, ls=1.46, zs=3, a=1.7):
	# Crighton+19 dndz fit
	return ls * ((1 + z) / (1. + zs)) ** a


nhi_fit = {}
# threshold correction to obtain observed dndz --> by eye
# nhi_fit['HM01'] = {8: 18.17, 9: 18.13,10: 18.41, 11: 17.89, 12: 17.39, 13: 17.65, 14: 17.69}#for HM01 model
# nhi_fit['HM12'] = {11: 17.985}#for HM12 model

# threshold correction to obtain observed dndz --> from Prochaska+10 Figure 17 fit
nhi_fit_prochaska = {}
nhi_fit_prochaska['HM01'] = {8: 17.92, 9: 18.04, 10: 17.89, 11: 17.9, 12: 18.03, 13: 17.86,
                             14: 17.77}  # for HM01 model

# threshold correction to obtain observed dndz --> Crighton+19 Figure 8 fit
nhi_fit['HM01'] = {8: 18.5, 9: 18.46, 10: 18.21, 11: 18.05, 12: 18.04, 13: 17.73,
                   14: 17.64}  # for HM01 model

# threshold correction to obtain observed dndz --> Crighton+19 Figure 8 fit,  new values with self shielding
nhi_fit['HM01'] = {8: 19.72, 9: 19.56, 10: 19.53, 11: 18.84, 12: 18.66, 13: 18.43,
                   14: 18.33}  # for HM01 model

nhi_fit['HM12'] = {11: 18.18}  # for HM12 model

fLLScorr = {8: .24, 9: .27, 10: .42, 11: .49, 12: .53, 13: .75,
            14: .83}  # correction to the LLS fraction given the Crighton+19 fit
fLLScorr_hi = {8: .24, 9: .27, 10: .42, 11: .49, 12: .53, 13: .75,
               14: .83}  # correction to the LLS fraction given the Crighton+19 fit
fLLScorr_lo = {8: .24, 9: .27, 10: .42, 11: .49, 12: .53, 13: .75,
               14: .83}  # correction to the LLS fraction given the Crighton+19 fit

lcube = 4096
coml = 25  # cMpc
com2pix = 163.84  # lcube/coml
kpc2pix = lcube / float(coml * 1e3)
rads = [0, 2, 4, 8, 12, 20, 30, 50, 100, 200]

lognhi = [13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21]
sb = [0.0001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96, 0.98, 1, 1, 1]


def fcorr(z):
	import numpy as np
	rval = np.array(redshifts.values())
	rkey = np.array(redshifts.keys())
	zlo = np.amax(rval[rval < z])
	zhi = np.amin(rval[rval > z])
	m = (z - zlo) / (zhi - zlo)
	s0, s1 = [rkey[rval == _ze][0] for _ze in [zlo, zhi]]
	return fLLScorr[s0] * (1 - m) + fLLScorr[s1] * m


omega_l = .714
omega_m = 1 - omega_l
H0 = 2.2e-18  # in 1/s
planck_H0 = 2.1963e-18
planck_omega_l = .693
G = 6.6725e-8  # cgs
proton_mass = 1.67e-24
meanbardens = 5.305653871142694e-07  # omega_m*(H0)**2/ proton_mass / (8. * 3.141592 * G)
Mpc2cm = 3.086e24
h = .6777
