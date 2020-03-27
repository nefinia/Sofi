import matplotlib.pyplot as plt
import numpy as np
from pyfits import getdata
from scipy import signal


def pix2l(p,l0=4750.,pixsize=.125):
	return l0+p*pixsize

def l2z(l):
	lya_rest = 1215.67
	return l/lya_rest-1

bins, zmin, zmax, fmin, fmax = 150, 2.9, 6.6, 25, 30
fitcats = 'mosaic', 'UDF'
labels = 'UDF-mosaic', 'UDF-10'
colors = 'gold', 'red'
data0 = getdata('../../all/cats/rafelski2015.fits', 1)
n, l, spec = np.loadtxt('../../UDF/skyspec.txt').T
if 0:
	cube = getdata('/net/galaxy-data/export/galaxydata/gallegos/UDF/DATACUBE_UDF-10.fits')
	obj = getdata('/net/galaxy-data/export/galaxydata/gallegos/UDF/DATACUBE_UDF-10.Objects_Id.fits')
	mask = obj>0
	cube[mask] = 0
	zl, yl, xl = cube.shape

red = l2z(l)
#spec = np.nansum(cube, (1, 2))
#f, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(4,2))

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.12, 0.3, 0.8, 0.6], yticklabels=[20,40,60,80,120,160], yticks=[20,40,60,80,120,160],
                   xticklabels=[], ylim=(0, 200))
ax2 = fig.add_axes([0.12, 0.1, 0.8, 0.2], ylim=(0, 30000),
				   yticklabels=[])

plt.tight_layout()
plt.xlabel('redshift')

#plt.yscale('symlog',linthreshy=1000.)
#plt.plot(red, signal.savgol_filter(np.abs(spec), 101, 8))
spec[662]=30000
ax2.plot(red, np.abs(spec),color='gray')
ax2.set_yticks([])
ax2.set_xlim(zmin, zmax)
ax2.set_ylim(0, 30000)
ax2.set_ylabel('Normalized \nSky Flux')
delta = .02
ax1.set_xlim(zmin-delta, zmax+delta)
#plt.savefig('../../Figures/spec.pdf')
#plt.close()

if 1:
	F775W = data0['F775W']
	zbpz = data0['zph1']
	#cool0 = (F775W >= fmin) & (F775W <= fmax) & (np.isfinite(F775W)) & (np.isfinite(zbpz)) & (zbpz>=zmin) & (zbpz<=zmax)
	cool0 = (zbpz>=zmin) & (zbpz<=zmax)
	#plt.figure(figsize=(5,4))
	ax1.hist(zbpz[cool0 & data0['in_UDF-mosaic']], bins=bins, label='Rafelski+15', histtype='barstacked')
	for fitcat, label, color in zip(fitcats, labels, colors):
		cool = (data0['in_%s' % label]) & cool0
		zh = np.histogram(zbpz[cool], bins=bins)
		zh1 = np.histogram(zbpz[cool], bins=50)
		data = getdata('../../%s/cats/laes.fits' % fitcat, 1)
		z = data['Z_MUSE']
		conf = data['CONFID']
		HST_F775W = data['HST_F775W']
		#cool = (np.isfinite(HST_F775W)) & (conf>=2) & (HST_F775W >= 25) & (HST_F775W <= 30) & (z>=zmin) & (z<=zmax)
		cool = (conf>=2) & (z>=zmin) & (z<=zmax)
		zhmuse = np.histogram(z[cool], range=[zmin, zmax], bins=bins)
		zhmuse1 = np.histogram(z[cool], range=[zmin, zmax], bins=50)
		comp = zhmuse1[0]/zh1[0].astype(float)

		#comp[comp>1] = 1.
		comp[comp<0] = 0.
		x = (zh1[1][:-1]+zh1[1][1:])/2.
		y = comp*100
		ax1.plot(x, signal.savgol_filter(y, 15, 3), color=color, linestyle='--')
		#plt.ylabel('F775W')
		ax1.hist(z[cool], bins=bins, label=label, color=color, histtype='barstacked')
		ax1.set_ylim((0,80))
		#plt.show()
		#plt.close()
		#plt.hist2d(z, HST_F775W, label=label)
	plt.xlim(red[0]-delta,red[-1]+delta)
	ax1.legend()
	plt.xlabel('z')
	plt.savefig('../../Figures/histogram_redshift_UDF_mosaic.pdf')
	plt.close()
