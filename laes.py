#!/usr/bin/env python
import numpy as np
from pyfits import getdata, PrimaryHDU, getheader
from math import pi
import os, subprocess, sys
from tools_sofi import distance  # , hubble
import astropy.io.fits as fits
from astropy import wcs
from tools_sofi import astroim, rdarg
from sys import argv

def l2pix(l):
    """Convert wavelength to pixels"""
    l0 = 4750.
    pixsize = 1.25  # 1 pixel = 1.25 Angstrom
    return (l - l0) / pixsize + 1


def lya_ob(z):
    """Convert lya to observed wavelength"""
    lya_rest = 1215.67
    return lya_rest * (1 + z)

simplestack = rdarg(argv, 'stack', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
hdu = PrimaryHDU()
pairs = []
makeim = rdarg(argv, 'makeim', bool, False)
centerfix = rdarg(argv, 'centerfix', bool, True)
theia = False
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
csub = rdarg(argv, 'csub', bool, False)
if csub:
    scsub = '.csub'
else:
    scsub = ''
dovar = True
xw = rdarg(argv, 'xw', int, 20)
yw = rdarg(argv, 'yw', int, 20)
zw = rdarg(argv, 'zw', int, 100)

if theia:
    folder = '/scratch/gallegos/MUSE/'
else:
    folder = '../../'

cubesmooth = False
if cubesmooth:
    smooth = False
    ssmooth = ' -boxsm 3 '
else:
    smooth = True
    ssmooth = ''

if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s.fits' % scsub
    nconf = 3
    if dovar:
        varname = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
if fitcat == 'UDF':
    cubename = 'DATACUBEFINALuser_20141021T055145_212058ad%s.fits' % scsub
    nconf = 1
filename = '%s/%s/%s' % (folder, fitcat, cubename)

data_cube = getdata(filename, 0)

zlim, ylim, xlim = data_cube.shape

xpix = []
ypix = []
zpix = []


cat = '%s/%s/cats/laes.fits' % (folder, fitcat)
if centerfix:
    catout = '%s/%s/cats/laesfixed.txt' % (folder, fitcat)
    fout = open(catout, 'w')
    fout.write('#id ra dec x y z flux_lya ra_fixed dec_fixed lambda_fixed x_fixed y_fixed z_fixed flux_lya_CubEx diff\n')
data = getdata(cat, 1)

if fitcat == 'HDFS':
    conf = data['confidence']
if fitcat == 'UDF':
    conf = data['Confidence']

good = np.where(conf >= nconf)[0]

if fitcat == 'HDFS':
    m300 = data['m300'][good]
    redshift = data['z'][good]
    wavelength = lya_ob(redshift)
    flya = data['LYALPHA_FLUX'][good]
    ids = data['id'][good]
    ra = data['ra'][good]
    dec = data['dec'][good]
    xs = data['X'][good]
    ys = data['Y'][good]
    zs = l2pix(wavelength)

if fitcat == 'UDF':
    hdulist = fits.open(filename)
    ids = data['id'][good]
    ra = data['ra'][good]
    dec = data['dec'][good]
    redshift = data['Z_MUSE'][good]
    flya = data['LYALPHA_FLUX'][good]
    wavelength = data['LYALPHA_LBDA_OBS'][good]
    ttt, header_data_cube = getdata(filename.replace('.csub', ''), 0, header=True)
    w = wcs.WCS(header_data_cube, hdulist)
    xs, ys, zs = w.all_world2pix(ra, dec, wavelength, 1)
    zs = l2pix(wavelength)

fits = []
stack_names = open(folder + "%s/LAEs/laes_stack.lst" % fitcat, 'w')

if 1:
    raf = []
    decf = []
    lamf = []
    xn = []
    yn = []
    zn = []
    fn = []

    for r, d, x, y, z, f, i in zip(ra, dec, xs, ys, zs, flya, ids):
        print "--------------------"
        print "ID", i
        print 'Catalogue values', x, y, z, data_cube[int(z-1), int(y-1), int(x-1)]  # , data_cube[z1 - d:z1 + d, y1 - d:y1 + d, x1 - d:x1 + d]

        name = folder + "%s/LAEs/lae%d.fits" % (fitcat, i)
        xmin = max(1, x-xw)
        xmax = min(xlim, x+xw)
        ymin = max(1, y-yw)
        ymax = min(ylim, y+yw)
        zmin = max(1, z-zw)
        zmax = min(zlim, z+zw)

        if os.path.isfile(name) and not overwrite:
            print "File %s already exists" % name
        else:
            os.system("CubeSel -cube %s[%d:%d,%d:%d,%d:%d] -out %s" % (filename, xmin, xmax, ymin, ymax, zmin, zmax, name))
            os.system("Cube2Im -cube %s[*,*,%d:%d]" % (name, zw-1, zw+1))

        if simplestack:
            fl = getdata(name)
            if fl.shape == (zw*2+1, yw*2+1, xw*2+1):
                fits.append(fl)
                stack_names.write('%s\n' % name)
        if makeim:
            astroim(name.replace('.fits', ".IM.fits"), smooth=smooth, saveim=True, show=False, dfig=(7, 7),
                    contours=True, scb_label='Flux [1e-20 cgs]', vmin=-1, vmax=1)
        if centerfix:
            sn = 2
            os.system('rm %s/%s/LAEs/lae*.cat' % (folder, fitcat))
            while not os.path.isfile(name.replace('.fits', '.cat')):
                print 'Using SN threshold', sn
                s = 'CubEx -cube %s -sn %d -m .false.' % (name, sn)
                os.system(s)
                sn /= 2
            fin = open(name.replace('.fits', '.cat'), 'r')
            rm = []
            dm = []
            lm = []
            xm = []
            ym = []
            zm = []
            im = []
            fm = []

            for l in fin:
                a=l.split()
                if a[0].isdigit():
                    a[0] = int(a[0])
                    if a[0] < 10:
                        im.append(a[0])
                        rm.append(a[-4])
                        dm.append(a[-3])
                        lm.append(a[-2])
                        zm.append(float(a[7])+zmin-1)
                        ym.append(float(a[6])+ymin-1)
                        xm.append(float(a[5])+xmin-1)
                        fm.append(float(a[18]))
            dd = np.sqrt((z-zm)**2+(y-ym)**2+(x-xm)**2)
            a = np.where(dd == min(dd))[0]
            idf = im[a]
            rf = rm[a]
            df = dm[a]
            lf = lm[a]
            zf = zm[a]
            yf = ym[a]
            xf = xm[a]
            ff = fm[a]
            xn.append(xf)
            yn.append(yf)
            zn.append(zf)
            fn.append(ff)
            fout.write("%d %s %s %1.4f %1.4f %1.4f %1.4f %s %s %s %1.4f %1.4f %1.4f %1.4f %1.4f\n"%(i, r, d, x, y, z,
                        f, rf, df, lf, xf, yf, zf, ff, dd[a]))

            print 'Fixed values', xf, yf, zf, data_cube[int(zf-1), int(yf-1), int(xf-1)], 'CubEx id selected', idf
            print 'Difference', x-xf, y-yf, z-zf, data_cube[int(z-1),
                    int(y-1), int(x-1)] - data_cube[int(zf-1), int(yf-1), int(xf-1)]

            xmin = max(1, xf-xw)
            xmax = min(xlim, xf+xw)
            ymin = max(1, yf-yw)
            ymax = min(ylim, yf+yw)
            zmin = max(1, zf-zw)
            zmax = min(zlim, zf+zw)

            namefix = folder + "%s/LAEs/laefix%d.fits" % (fitcat, i)
            os.system("CubeSel -cube %s[%d:%d,%d:%d,%d:%d] -out %s" % (
            filename, xmin, xmax, ymin, ymax, zmin, zmax, namefix))
            os.system("Cube2Im -cube %s[*,*,%d:%d] %s" % (namefix, zw-1, zw+1, ssmooth))

            if makeim:
                astroim(namefix.replace('.fits', ".IM.fits"), smooth=smooth, saveim=True, show=False, dfig=(7, 7),
                        contours=True, scb_label='Flux [1e-20 cgs]', vmin=-1, vmax=1)
if centerfix:
    fout.close()

if simplestack:
    print 'Simple stack', len(fits)

    fits2 = np.array(fits)
    #fits2[abs(fits2) > 5*np.nanstd(fits2)] = np.nan
    hdu = PrimaryHDU()
    hdu.data = np.nanmean(fits2, axis=0)
    hdu.writeto(folder + "%s/LAEs/laes_stack.fits" % fitcat, clobber=True)