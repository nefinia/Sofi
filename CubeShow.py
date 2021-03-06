__author__ = 'gallegos'

from sys import argv

import numpy as np
from pyfits import getdata


def rdarg(argv, key, type=None, default=None, listtype=int):
	if len(argv) > 1:
		opt = np.where([a == '-%s' % key for a in argv])[0] + 1
		if len(opt) > 0:
			name = argv[int(opt)]
			if type is list:
				name = name.split(',')
				if listtype == int:
					name = [int(i) for i in name]
				elif listtype == float:
					name = [float(i) for i in name]
			elif type is bool:
				name = eval(str(name))
			elif type is int:
				name = int(name)
			elif type is float:
				name = float(name)
			elif type is str:
				name = str(name)
			return name
	if default is not None: return default


cube = rdarg(argv, 'name', str, '')
std = rdarg(argv, 'std', float, None)
bin = rdarg(argv, 'bin', int, 1)
dim = rdarg(argv, 'dim', int, 0)
colors = rdarg(argv, 'color', bool, False)
dmin = rdarg(argv, 'min', int, 0)
dmax = rdarg(argv, 'max', int, None)

if len(argv) < 3:
	print "CubeShow\n\
	Print fits files in the terminal.\n\
    each pixel will be displayed as a number corresponding to its\n\
    standard deviation with respect to the mean.\n\
    Negative values will be displayed with '-'.\n\
    Values above 10 sigma will be displayed with '*'.\n\
    3d fits will be collapsed into two axis.\n\n\
    Options:\n\
            -name:  fits file name.\n\
            -bin:   display numbers in multiples of the std (default=1).\n\
            -dim:   dimension to be collapsed (0, 1 or 2, default=0).\n\
            -min:   minimum value of the collapsed dimension.\n\
            -max:   maximum value of the collapsed dimension.\n\
            -std:   define your own standard deviation.\n\
            -color: display numbers with colors (default=False)."

if cube != '':
	fit = getdata(cube)
	shape = fit.shape
	print 'Image shape', shape
	if dmax is None: dmax = shape[dim]
	if len(shape) > 2:
		if dim == 0: fit = np.nanmean(fit[dmin:dmax+1, :, :], dim)
		if dim == 1: fit = np.nanmean(fit[:, dmin:dmax+1, :], dim)
		if dim == 2: fit = np.nanmean(fit[:, :, dmin:dmax+1], dim)
	yl, xl = fit.shape
	if std is None: std = np.nanstd(fit)
	print 'Generating image with std', std
	std = bin * std
	s = ''
	for y in range(yl):
		for x in range(xl):
			d = fit[y, x]
			if d < 0: s += '\033[0m-'
			for i in range(10):
				if (d >= i * std) & (d < (i + 1) * std):
					if colors: s += '\033[1;3%dm%d' % (i, i)
					else: s += '%d' % i
			if d >= 10 * std: s += '\033[0m*'
		s += '\n'
	print s, '\033[0m '
