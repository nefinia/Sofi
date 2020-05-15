import os
from sys import argv

import numpy as np
import scipy.interpolate as inter
from pyfits import getdata, PrimaryHDU

from params import sbpeaks, zlens
from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,

snaps = rdarg(argv, 'snap', list, [6, 7, 8, 9, 10, 11, 12, 13, 14], int)
overwrite = rdarg(argv, 'overwrite', bool, False)

hdu = PrimaryHDU()
f = '../../zap/'
fin = '/net/eos/scratch/gallegos/EAGLE/'
s = '_x_HM12_512_4096_12_56_6.73E-03_CaseA.LLS.fits'

lognhi = [0, 12, 13, 14, 15, 16, 16.5, 17, 17.5, 17.7, 17.8, 18, 18.4, 18.7, 19, 19.4, 19.6, 19.8, 20.1, 20.5, 21,
		  22, 23, 30]
nhi = np.power(10, lognhi)
sb = [0, 0, 0.00001, 0.001, 0.003, 0.03, 0.08, 0.2, 0.45, 0.55, 0.6, 0.7, .8, 0.85, 0.9, 0.938, 0.96,
	  0.98, 1, 1, 1, 1, 1, 1]


nhi2sb = inter.interp1d(nhi, sb)
lognhi2sb = inter.interp1d(lognhi, sb)
sb2nhi = inter.interp1d(sb, nhi)
nname = f+'UDF10.norm.fits'
fnoise = getdata(nname)
zl, yl, xl = fnoise.shape
#h['NAXIS1'], h['NAXIS2'] = 4096, 4096
for snap in snaps:
	print snap
	erand = f+'snap%d.noise.fits' % snap
	sb = f+'snap%d.sb.fits' % snap
	if not os.path.isfile(sb) or overwrite:
		data = getdata(fin + 'snap%d_x_HM12_512_4096_12_%d_6.73E-03_CaseA.NHI.fits' % (snap, zlens[snap]))
		zs, ys, xs = data.shape
		flog = np.log10(data)
		flog[flog > 23] = 23
		hdu.data = lognhi2sb(flog)*sbpeaks[snap]
		hdu.writeto(sb, clobber=True)
	else: print '%s exists' % sb
	#h['NAXIS3'] = zlens[snap]
	#h['CRVAL3'] = lya_rest*(1+redshifts[snap])
	if 0:
		if not os.path.isfile(erand) or overwrite:
			data = getdata(fin + 'snap%d_x_HM12_512_4096_12_%d_6.73E-03_CaseA.NHI.fits' % (snap, zlens[snap]))
			zs, ys, xs = data.shape
			frand = np.copy(fnoise)
			flog = np.log10(data[:, :yl, :xl])
			flog[flog>23] = 23
			frand[:zs, :, :] = lognhi2sb(flog)*sbpeaks[snap]#+np.random.normal(size=data.shape)*10.
			hdu.data = frand
			hdu.writeto(erand, clobber=True)
		else:
			print '%s exists' % erand
			#frand = getdata(erand)
		urand = f+'snap%d.udf.fits' % snap
		if not os.path.isfile(urand) or overwrite: os.system('CubeArit %s - %s %s' % (f+'UDF10.zeros.fits', erand, urand))
		else: print '%s exists' % urand
		zrand = f+'snap%d.zap.fits' % snap
		if not os.path.isfile(zrand): os.system('zap %s -o %s' % (urand, zrand))
