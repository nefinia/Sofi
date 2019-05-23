#!/usr/bin/env python
import astropy.io.fits as fits
import os
from astropy import wcs
from pyfits import getdata, PrimaryHDU
from sys import argv

import numpy as np

from tools_sofi import rdarg


def l2pix(l):
    """Convert wavelength to pixels"""
    l0 = 4750.
    pixsize = 1.25  # 1 pixel = 1.25 Angstrom
    return (l - l0) / pixsize + 1


def lya_ob(z):
    """Convert lya redshift to observed wavelength"""
    lya_rest = 1215.67
    return lya_rest * (1 + z)


simplestack = rdarg(argv, 'stack', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
hdu = PrimaryHDU()
pairs = []
makeim = rdarg(argv, 'makeim', bool, True)
mask = rdarg(argv, 'mask', bool, False)
centerfix = rdarg(argv, 'centerfix', bool, False)
cubex = rdarg(argv, 'cubex', bool, False)
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
folder = rdarg(argv, key='folder', type=str, default='/net/astrogate/export/astrodata/gallegos/')
foldercat = rdarg(argv, key='folder', type=str, default='../../')
extraname = rdarg(argv, key='extraname', type=str, default='')
csub = rdarg(argv, 'csub', bool, True)
single = rdarg(argv, 'single', bool, False)
id0 = rdarg(argv, 'id', int, 1)
smooth = rdarg(argv, 'smooth', bool, True)
line = rdarg(argv, 'line', str, 'lya')
white = rdarg(argv, 'white', bool, False)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)

coord = rdarg(argv, 'coord', str, 'x')
snap = rdarg(argv, 'snap', int, 12)

verbose = rdarg(argv, 'verbose', int, 1)

if verbose < 2:
    vb = ' > /dev/null 2>&1'
else:
    vb = ''

dovar = True
xw = rdarg(argv, 'xw', int, 200)#for EAGLE it was last time 500 (15.10.18)
yw = rdarg(argv, 'yw', int, 200)#for EAGLE it was last time 500 (15.10.18)
zw = rdarg(argv, 'zw', int, 20)
vmin = rdarg(argv, 'vmin', float, -.5)
vmax = rdarg(argv, 'vmax', float, 1)
std = rdarg(argv, 'std', float, None)
smask = '.mask' * mask
binsize = 1
xmin = -xw / 2
xmax = xw / 2
ymin = -yw / 2
ymax = yw / 2
zmin = -zw / 2  # min(z)
zmax = zw / 2  # max(z)
lenx = int(xw / binsize) + 1
leny = int(yw / binsize) + 1

cubesmooth = False
if cubesmooth:
    ssmooth = ' -boxsm 3 '
    smooth = False
else:
    ssmooth = ''

if fitcat == 'EAGLE':
    csub = False
    smooth = False
    cubename = 'SB_snap%s_%s%s.fits' % (snap, coord, smask)
    nconf = 1

if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s.fits' % '.csub'*csub
    nconf = 3
    if dovar:
        varname = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
if fitcat == 'UDF':
    if line == 'ovi':
        cubename = 'UDF10.z1300%s.fits' % '.csub' * csub
    cubename = 'UDF10.z1300%s.fits' % '.csub'*csub#'DATACUBEFINALuser_20141021T055145_72f74684%s.fits' % '.csub'*csub  # last version --> not nowwwww Aug2018
    # cubename = 'DATACUBEFINALuser_20141021T055145_212058ad%s.fits' % scsub
    nconf = 1

if fitcat == 'mosaic' or fitcat == 'mosaic10':
    # cubename = '/net/astrogate/export/astrodata/gallegos/DATACUBEFINALuser_20140921T071255_0f2b159e%s.fits' % scsub
    # cubename = 'DATACUBEFINALuser_20141021T055145_212058ad%s.fits' % scsub
    #cubename = '/net/astrogate/export/astrodata/gallegos/mosaic.z1300%s.fits' % scsub  ########
    cubename = '/net/astrogate/export/astrodata/gallegos/mosaic/DATACUBE_UDF-MOSAIC.z1300%s.fits' % '.csub'*csub
    nconf = 1
    filename = cubename
    vmin = -3
    vmax = 10
    if std is None: std = 100

else:
    filename = '%s/%s/%s' % (folder, fitcat, cubename)

data_cube = getdata(filename, 0)
zlim, ylim, xlim = data_cube.shape

if fitcat == 'HDFS' or fitcat == 'mosaic' or fitcat == 'mosaic10':
    cut = 12
else:
    cut = 0

xpix = []
ypix = []
zpix = []
xr = np.arange(xw + 1)
yr = np.arange(yw + 1)
zr = np.arange(zw + 1)

if fitcat == 'mosaic10':
    cat = '%s/%s/cats/laes.fits' % (foldercat, 'UDF')
elif fitcat == 'EAGLE':
    cat = '%sEAGLE/cats/gals_snap%d.fits' % (foldercat, snap)
else:
    cat = '%s/%s/cats/laes.fits' % (foldercat, fitcat)
if centerfix:
    catout = '%s/%s/cats/laesfixed.txt' % (foldercat, fitcat)
    fout = open(catout, 'w')
    fout.write(
        '#id ra dec x y z flux_lya ra_CubEx dec_CubEx lambda_CubEx x_CubEx y_CubEx z_CubEx flux_lya_CubEx redshift_CubEx lya_lum_CubEx diff com_dist\n')
data = getdata(cat, 1)

if fitcat == 'HDFS':
    conf = data['confidence']
if fitcat == 'UDF' or fitcat == 'mosaic':
    conf = data['CONFID']

hdulist = fits.open(filename)

if fitcat == 'EAGLE':
    lcom = 25.  # comoving length of the simulation
    lcube = 4096
    conv = lcube / lcom
    ids = data['id']
    wavelength = np.zeros(len(ids)) + 1
    flya = np.zeros(len(ids)) + 1
    ra = np.zeros(len(ids)) + 1
    dec = np.zeros(len(ids)) + 1
    cs = ['x', 'y', 'z']
    cs.remove(coord)
    xs = np.round(data[cs[0]] * xlim / lcom).astype(int)
    ys = np.round(data[cs[1]] * ylim / lcom).astype(int)
    zs = np.round(data[coord] * zlim / lcom).astype(int)

    #sbpeaks = {'10': .739, '11': 1.293, '12': 2.579}
    #reds = {'10': 984, '11': 528, '12': 17}
    #redshifts = {'10': 3.984, '11': 3.528, '12': 3.017}
    #asec2kpcs = {'10': 7.842, '11': 7.449, '12': 7.108}
    com2pix = 163.84  # lcube/coml
    kpc2pix = lcube / float(lcom * 1e3)
    #vmin = sbpeaks['%d' % snap] * .03
    #vmax = sbpeaks['%d' % snap] * .5
    #asec2pix = asec2kpcs['%d' % snap] * (1 + redshifts['%d' % snap]) * kpc2pix

else:
    good = np.where(conf >= nconf)[0]

if fitcat == 'HDFS':
    m300 = data['m300'][good]
    redshift = data['z'][good]
    wavelength = lya_ob(redshift)
    flya = data['LYALPHA_FLUX'][good]
    ids = data['id'][good]
    ra = data['ra'][good]
    dec = data['dec'][good]
# xs = np.round(data['X'][good]).astype('int')
# ys = np.round(data['Y'][good]).astype('int')
# zs = l2pix(wavelength).astype('int')

if fitcat == 'UDF' or fitcat == 'mosaic':
    ids = data['ID'][good]
    ra = data['RA'][good]
    dec = data['DEC'][good]
    redshift = data['Z_MUSE'][good]
    if fitcat == 'mosaic':
        flya = data['LYALPHA_FLUX_SUM'][good]
    else:
        flya = data['LYALPHA_FLUX'][good]
    wavelength = data['LYALPHA_LBDA_OBS'][good]

ttt, header_data_cube = getdata(filename, 0, header=True)  # .replace('.csub', '')
# Removing COMMENT key to avoid problems reading non-ascii characters
cards = header_data_cube.cards
bad = ['COMMENT' == b[0] for b in cards]
for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')

if fitcat != 'EAGLE':
    w = wcs.WCS(header_data_cube, hdulist)
    #xs, ys, zs = np.round(w.all_world2pix(ra, dec, wavelength, 1))
    #zs = l2pix(wavelength)
    xs, ys, zs = [data['x'][good], data['y'][good], data['z'][good]]
    if line == 'ovi':
        wav = (1+redshift)*1035
        zs = l2pix(wav)


# xs = np.round(xs).astype('int')
# ys = np.round(ys).astype('int')
# zs = np.round(zs).astype('int')

if fitcat == 'HDFS':
    # xs = data['x_sofi'][good]
    # ys = data['y_sofi'][good]
    # zs = data['z_sofi'][good]
    xs = data['x'][good]
    ys = data['y'][good]
    zs = data['z'][good]

fits = []
if fitcat == 'EAGLE': stack_names = open(folder + "%s/gals/gal_stack.lst" % fitcat, 'w')
else: stack_names = open(folder + "%s/gals/laes_stack.lst" % fitcat, 'w')
lst = []

hdu = PrimaryHDU()
tn = 10
zw0 = 3

if fitcat == 'EAGLE':
    cool = ids > 0
    tn = 5
else:
    cool = (zs < (zlim - zw / 2)) & (zs > zw / 2)  # range(len(ids))#cool = ids == 138#
if single: cool = ids == id0
notcool = zs > (zlim - zw / 2)

for r, d, l, x, y, z, f, i in \
    zip(ra[cool], dec[cool], wavelength[cool], xs[cool], ys[cool], zs[cool], flya[cool], ids[cool]):

    print "--------------------"
    print "ID", i
    # print 'Catalogue values', x, y, z, data_cube[
    #    int(z - 1), int(y - 1), int(x - 1)]  # , data_cube[z1 - d:z1 + d, y1 - d:y1 + d, x1 - d:x1 + d]

    if fitcat == 'EAGLE':
        name = "%s/%s/gals/snap%d_%s/%d%s%s.fits" % (folder, fitcat, snap, coord, i, extraname, smask)
    else:
        name = "%s/%s/gals/%d%s%s.fits" % (folder, fitcat, i, extraname, ('.%s' % line) * (line != 'lya'))
    lst.append(name)

    if overwrite or not os.path.isfile(name):
        flux = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
        if line == 'ovi':
                   flux2 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
        if fitcat == 'EAGLE':
            xt = (xr + xmin + x) % xlim
            yt = (yr + ymin + y) % ylim
            zt = (zr + zmin + z) % zlim
            for i in xr:
                for j in yr:
                    for k in zr:
                        flux[k, j, i] = data_cube[zt[k], yt[j], xt[i]]

        else:
            _xmin = max(cut, x - xw / 2)
            _xmax = min(xlim - cut, x + xw / 2)
            _ymin = max(cut, y - yw / 2)
            _ymax = min(ylim - cut, y + yw / 2)
            _zmin = max(1, z - zw / 2)
            _zmax = min(zlim, z + zw / 2)
            xt = xr + xmin + x
            yt = yr + ymin + y
            zt = zr + zmin + z

            for i in xr:
                print i
                for j in yr:
                    for k in zr:
                        zc, yc, xc = np.round([zt[k], yt[j], xt[i]]).astype(int)
                        cond = (xc>=_xmin) and (xc<=_xmax) and (yc>=_ymin) and (yc<=_ymax) and (zc>=_zmin) and (zc<=_zmax)
                        if cond: flux[k, j, i] = data_cube[zc-1, yc-1, xc-1]

        hdu.data = flux
        hdu.writeto(name, clobber=True)
        if fitcat == 'EAGLE': hdu.data = flux[zw / 2, :, :]
        else:
            hdu.data = np.nansum(flux, 0)
        hdu.writeto(name.replace(".fits", ".IM.fits"), clobber=True)
    else:
        print 'File already exists'

if simplestack:
    print 'Simple stack', len(fits)
    from tools_sofi import stack

    zw0 = 1
    if fitcat == 'EAGLE': stack(lst, folder + "%s/gals/stack.fits" % fitcat, zw0=zw0)
    else: stack(lst, folder + "%s/gals/stack.fits" % fitcat, zw0=zw0)
