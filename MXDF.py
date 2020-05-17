import numpy as np
from pyfits import getdata, PrimaryHDU
overwrite = False
folder = '/net/galaxy-data/export/galaxydata/gallegos/MXDF/'
cube = getdata(folder+'MXDF_DATACUBE_FINAL_PROPVAR.csub.fits')
bcorr = np.copy(cube)
bcube = np.empty(cube.shape)
if overwrite:
	mask2d = getdata(folder+'MXDF_DATACUBE_FINAL_PROPVAR.IM.Objects_Id.fits')
	mask3d = getdata(folder+'MXDF_DATACUBE_FINAL_PROPVAR.Objects_Id.fits')
	cube[:, mask2d[0]>0] = np.nan
	cube[mask3d>0] = np.nan
zl, yl, xl = cube.shape
y, x = np.ogrid[0:yl, 0:xl]
r2 = (y-yl/2)**2+(x-xl/2)**2
sigma = 3
rmax = int(np.sqrt(2)*xl/2)
rads = np.arange(0, rmax+1, 2)
nr = len(rads)-1
if overwrite: rspec = np.empty((nr, zl))
else: rspec = getdata(folder+'rspec.fits')


for r0, r1, i in zip(rads[:-1], rads[1:], range(nr)):
	print r0, r1
	inside = (r2 <= r1 ** 2) & (r2 > r0 ** 2)
	for z in range(zl):
		if overwrite:
			cin = cube[z, inside]
			std = np.nanstd(cin)
			cin[abs(cin)>sigma*std] = np.nan
			rspec[i,z] = np.nanmean(cin)
		bcube[z, inside] = rspec[i, z]

hdu = PrimaryHDU()
hdu.data = bcube
hdu.writeto(folder+'bkgcube.fits', clobber=True)
bcube[np.isnan(bcube)] = 0
hdu.data = bcorr-bcube
hdu.writeto(folder+'MXDF.bkgcorr.csub.fits', clobber=True)
if overwrite:
	hdu.data = rspec
	hdu.writeto(folder+'rspec.fits', clobber=True)
