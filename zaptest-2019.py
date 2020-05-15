import os

import numpy as np
from pyfits import getdata, PrimaryHDU

hdu = PrimaryHDU()
f = '../../zap/'
nrand = 10
rs = 10.
sigma_wav = 2.
halov = True
nhalo = 100
fmuse = 'UDF10.zeros.fits'
fm = getdata(f + fmuse)
zl, yl, xl = fm.shape
z, y, x = np.ogrid[0: zl, 0: yl, 0: xl]

for nr in range(nrand):
	fout0 = 'norm.%d.fits' % nr
	fout1 = 'UDF10.norm.%d.fits' % nr
	if os.path.isfile(f+fout0):
		data = getdata(f+fout0)
	else:
		data = np.random.normal(size=fm.shape)
		hdu.data = data
		hdu.writeto(f+fout0, clobber=True)

	if not os.path.isfile(f+fout1):
		os.system('CubeArit %s%s + %s%s %s%s' % (f, fmuse, f, fout0, f, fout1))

	fnorm = getdata(f + fout1)

	fhalos = 'halos.%d.fits' % nr

	if os.path.isfile(f+fhalos) and not halov:
		xs, ys, zs = np.loadtxt('%shalos_rs%d_swav%d.%d.dat' % (f, rs, sigma_wav, nr), dtype=int).T
		halos = getdata(f + fhalos)
	else:
		def prof(p, c, s):
			return np.exp(-((p[2] - c[2]) / s[2]) ** 2) * \
				   np.exp(-((p[1] - c[1]) / s[1]) ** 2) * \
				   np.exp(-(p[0] - c[0]) ** 2 / 2. / s[0] ** 2) / s[0]


		xs = np.random.randint(0, xl, size=nhalo)
		ys = np.random.randint(0, yl, size=nhalo)
		zs = np.random.randint(0, zl, size=nhalo)
		np.savetxt('%shalos_rs%d_swav%d.%d.dat' % (f, rs, sigma_wav, nr), zip(xs, ys, zs), header='x y z', fmt='%d')
		n = 0
		_h = np.zeros((zl, yl, xl))

		for i, j, k in zip(xs, ys, zs):
			print n, i, j, k
			n += 1
			_h += 10*prof([z, y, x], [k, j, i], [sigma_wav, rs, rs])

		hdu.data = _h
		hdu.writeto(f+fhalos, clobber=True)

	fout2 = 'UDF10.norm_halos.%d.fits' % nr
	if not os.path.isfile(f+fout2):
		s = 'CubeArit %s%s + %s%s %s%s' % (f, fout1, f, fhalos, f, fout2)
		print s
		os.system(s)
	halos = getdata(f+fout2)

	fout3 = 'ZAP.norm_halos.%d.fits' % nr
	if not os.path.isfile(f+fout3):
		os.system('zap %s%s -o %s%s' % (f, fout2, f, fout3))
	fzap = getdata(f+fout3)
	frandzap = 'ZAP.norm.%d.fits' % nr
	fnormzap = getdata(f+frandzap)

	###radprof!!
	rw = 2
	zw = 3
	rads = np.array([0, .5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40, 45, 50, 55, 60])
	asec2pix = 5.
	rpix = rads*asec2pix
	npix = len(rpix)
	n = 0
	yhalos = []
	yzap = []
	ynorm = []
	ynormzap = []
	ygal = []
	for r in range(npix-1):
		yh = []
		yz = []
		yn = []
		ynz = []
		for i, j, k in zip(xs, ys, zs):
			print n, i, j, k
			n += 1
			cool = (((x-i)**2+(y-j)**2 >= rpix[r]**2) & ((x-i)**2+(y-j)**2 < rpix[r+1]**2))[0]
			yh.append(np.nanmean(np.nansum(halos[k-zw:k+zw+1, cool], 1)))
			yz.append(np.nanmean(np.nansum(fzap[k-zw:k+zw+1, cool], 1)))
			yn.append(np.nanmean(np.nansum(fnorm[k-zw:k+zw+1, cool], 1)))
			ynz.append(np.nanmean(np.nansum(fnormzap[k-zw:k+zw+1, cool], 1)))
		yhalos.append(np.nanmean(yh))
		yzap.append(np.nanmean(yz))
		ynorm.append(np.nanmean(yn))
		ynormzap.append(np.nanmean(ynz))

	r = (rads[1:]+rads[:-1])/2.
	np.savetxt(f+'sbprof.%d.dat' % nr, zip(r, yhalos, yzap, ynorm, ynormzap), header='radius SB SB_ZAP SB_rand SB_rand_ZAP')

	doplot= True

	if doplot:
		import matplotlib.pyplot as plt
		plt.plot(r, yhalos, label='Initial halos')
		plt.plot(r, yzap, label='ZAP halos')
		plt.plot(r, ynorm, label='Random')
		plt.plot(r, ynormzap, label='Random ZAP')
		plt.plot(r, np.zeros(nr-1), linestyle=':')
		plt.xlabel('distance [arcsec]')
		plt.ylabel('flux [arbitrary units]')
		plt.xlim(3, 60)
		plt.ylim(-.015, .015)
		plt.legend()
		#plt.show()
		plt.savefig(f+'sbprof.%d.png' % nr)
		plt.close()

