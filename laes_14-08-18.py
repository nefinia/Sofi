#!/usr/bin/env python
import numpy as np
from pyfits import getdata, PrimaryHDU, getheader, delval
from math import pi
import os, subprocess, sys
from tools_sofi import distance  # , hubble
import astropy.io.fits as fits
from astropy import wcs
from tools_sofi import astroim, rdarg, hms2deg, comdist
from sys import argv
import itertools as it


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
centerfix = rdarg(argv, 'centerfix', bool, False)
cubex = rdarg(argv, 'cubex', bool, False)
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
folder = rdarg(argv, key='folder', type=str, default='/net/astrogate/export/astrodata/gallegos/')
extraname = rdarg(argv, key='extraname', type=str, default='')
csub = rdarg(argv, 'csub', bool, True)
single = rdarg(argv, 'single', bool, False)
id0 = rdarg(argv, 'id', int, 1)
smooth = rdarg(argv, 'smooth', bool, True)
white = rdarg(argv, 'white', bool, False)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)

coord = rdarg(argv, 'coord', str, 'x')
snap = rdarg(argv, 'snap', int, 12)


verbose = rdarg(argv, 'verbose', int, 1)

if verbose < 2: vb = ' > /dev/null 2>&1'
else: vb = ''

if csub: scsub = '.csub'
else: scsub = ''
dovar = True
xw = rdarg(argv, 'xw', int, 500)
yw = rdarg(argv, 'yw', int, 500)
zw = rdarg(argv, 'zw', int, 12)
vmin = rdarg(argv, 'vmin', float, -.5)
vmax = rdarg(argv, 'vmax', float, 1)
std = rdarg(argv, 'std', float, None)

binsize = 1
xmin = -xw/2
xmax = xw/2
ymin = -yw/2
ymax = yw/2
zmin = -zw/2#min(z)
zmax = zw/2#max(z)
lenx = int(xw / binsize) + 1
leny = int(yw / binsize) + 1


cubesmooth = False
if cubesmooth:
    ssmooth = ' -boxsm 3 '
    smooth = False
else: ssmooth = ''


if fitcat == 'EAGLE':
    csub = False
    smooth = False
    cubename = 'SB_snap%s_%s.fits' % (snap, coord)
    nconf = 1

if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR%s.fits' % scsub
    nconf = 3
    if dovar:
        varname = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
if fitcat == 'UDF':
    cubename = 'DATACUBEFINALuser_20141021T055145_72f74684%s.fits' % scsub # last version
    #cubename = 'DATACUBEFINALuser_20141021T055145_212058ad%s.fits' % scsub
    nconf = 1

if fitcat == 'mosaic' or fitcat == 'mosaic10':
    #cubename = '/net/astrogate/export/astrodata/gallegos/DATACUBEFINALuser_20140921T071255_0f2b159e%s.fits' % scsub
    #cubename = 'DATACUBEFINALuser_20141021T055145_212058ad%s.fits' % scsub
    cubename = '/net/astrogate/export/astrodata/gallegos/mosaic.z1300%s.fits' % scsub########
    nconf = 1
    filename = cubename
    vmin = -3
    vmax = 10
    if std is None: std = 100

else: filename = '%s/%s/%s' % (folder, fitcat, cubename)

data_cube = getdata(filename, 0)
zlim, ylim, xlim = data_cube.shape

if fitcat == 'HDFS' or fitcat == 'mosaic' or fitcat == 'mosaic10':
    cut = 12
else:
    cut = 0

xpix = []
ypix = []
zpix = []

if fitcat == 'mosaic10': cat = '%s/%s/cats/laes.fits' % (folder, 'UDF')
elif fitcat == 'EAGLE': cat = '%sEAGLE/cats/gals_snap%d.fits' % (folder, snap)
else: cat = '%s/%s/cats/laes.fits' % (folder, fitcat)
if centerfix:
    catout = '%s/%s/cats/laesfixed.txt' % (folder, fitcat)
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
    lcom = 25. # comoving length of the simulation
    lcube = 4096
    conv = lcube/lcom
    ids = data['id']
    wavelength = np.zeros(len(ids))+1
    flya = np.zeros(len(ids))+1
    ra = np.zeros(len(ids))+1
    dec = np.zeros(len(ids))+1
    cs = ['x', 'y', 'z']
    cs.remove(coord)
    xs = np.round(data[cs[0]]*xlim/lcom).astype(int)
    ys = np.round(data[cs[1]]*ylim/lcom).astype(int)
    zs = np.round(data[coord]*zlim/lcom).astype(int)
    xr = np.arange(xw + 1)
    yr = np.arange(yw + 1)
    zr = np.arange(zw + 1)
    sbpeaks = {'10': .739, '11': 1.293, '12': 2.579}
    reds = {'10': 984, '11': 528, '12': 17}
    redshifts = {'10': 3.984, '11': 3.528, '12': 3.017}
    asec2kpcs = {'10': 7.842, '11': 7.449, '12': 7.108}
    com2pix = 163.84  # lcube/coml
    kpc2pix = lcube / float(lcom * 1e3)
    vmin = sbpeaks['%d' % snap]*.03
    vmax = sbpeaks['%d' % snap]*.5
    asec2pix = asec2kpcs['%d' % snap] * (1 + redshifts['%d' % snap]) * kpc2pix

else: good = np.where(conf >= nconf)[0]

if fitcat == 'HDFS':
    m300 = data['m300'][good]
    redshift = data['z'][good]
    wavelength = lya_ob(redshift)
    flya = data['LYALPHA_FLUX'][good]
    ids = data['id'][good]
    ra = data['ra'][good]
    dec = data['dec'][good]
    #xs = np.round(data['X'][good]).astype('int')
    #ys = np.round(data['Y'][good]).astype('int')
    #zs = l2pix(wavelength).astype('int')

if fitcat == 'UDF' or fitcat == 'mosaic':
    ids = data['ID'][good]
    ra = data['RA'][good]
    dec = data['DEC'][good]
    redshift = data['Z_MUSE'][good]
    if fitcat == 'mosaic': flya = data['LYALPHA_FLUX_SUM'][good]
    else: flya = data['LYALPHA_FLUX'][good]
    wavelength = data['LYALPHA_LBDA_OBS'][good]

ttt, header_data_cube = getdata(filename, 0, header=True)#.replace('.csub', '')
# Removing COMMENT key to avoid problems reading non-ascii characters
cards = header_data_cube.cards
bad = ['COMMENT' == b[0] for b in cards]
for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')

if fitcat != 'EAGLE':
    w = wcs.WCS(header_data_cube, hdulist)
    xs, ys, zs = np.round(w.all_world2pix(ra, dec, wavelength, 1))
    zs = l2pix(wavelength)

#xs = np.round(xs).astype('int')
#ys = np.round(ys).astype('int')
#zs = np.round(zs).astype('int')

if fitcat == 'HDFS':
    #xs = data['x_sofi'][good]
    #ys = data['y_sofi'][good]
    #zs = data['z_sofi'][good]
    xs = data['x'][good]
    ys = data['y'][good]
    zs = data['z'][good]


fits = []
stack_names = open(folder + "%s/LAEs/laes_stack.lst" % fitcat, 'w')
lst = []
outf = open(folder + "%s/LAEs/laes.c.dat" % fitcat, 'w')

raf = []
decf = []
lamf = []
xn = []
yn = []
zn = []
fn = []
hdu = PrimaryHDU()
tn = 10
zw0 = 3

if fitcat == 'EAGLE':
    cool = ids>0
    tn = 5
else: cool = (zs < (zlim-zw/2)) & (zs > zw/2)#range(len(ids))#cool = ids == 138#
if single: cool = ids == id0
notcool = zs > (zlim-zw/2)

if white:
    outw = open(folder + "%s/LAEs/laes_white.dat" % fitcat, 'w')


for r, d, l, x, y, z, f, i in \
        zip(ra[cool], dec[cool], wavelength[cool], xs[cool], ys[cool], zs[cool], flya[cool], ids[cool]):

    print "--------------------"
    print "ID", i
    #print 'Catalogue values', x, y, z, data_cube[
    #    int(z - 1), int(y - 1), int(x - 1)]  # , data_cube[z1 - d:z1 + d, y1 - d:y1 + d, x1 - d:x1 + d]

    if fitcat == 'EAGLE': name = "%s/%s/LAEs/snap%d_%s/lae%d%s.fits" % (folder, fitcat, snap, coord, i, extraname)
    else: name = "%s/%s/LAEs/lae%d%s.fits" % (folder, fitcat, i, extraname)
    lst.append(name)
    outf.write(name+'\n')

    if white:
        outf.write(name.replace(".fits", "_white.fits")+'\n')

    if fitcat == 'EAGLE':

        if overwrite or not os.path.isfile(name):
            xt = (xr+xmin+x) % xlim
            yt = (yr+ymin+y) % ylim
            zt = (zr+zmin+z) % zlim
            flux = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
            for i in xr:
                for j in yr:
                    for k in zr:
                        flux[k, j, i] = data_cube[zt[k], yt[j], xt[i]]
            hdu.data = flux
            hdu.writeto(name, clobber=True)
            hdu.data = flux[zw/2, :, :]
            hdu.writeto(name.replace(".fits", ".IM.fits"), clobber=True)
        else: print 'File already exists'
    else:
        xmin = max(1 + cut, x - xw/2)
        xmax = min(xlim - cut, x + xw/2)
        ymin = max(1 + cut, y - yw/2)
        ymax = min(ylim - cut, y + yw/2)
        zmin = max(1, z - zw/2)
        zmax = min(zlim, z + zw/2)

        print "x %d y %d x %d x[%d,%d] y[%d,%d] z[%d,%d]" % (
        x, y, z, xmin, xmax, ymin, ymax, zmin, zmax)

        if 1:
            xc = np.arange(np.round(xmin-1+0.0001).astype('int'), np.round(xmax+0.0001).astype('int'))
            yc = np.arange(np.round(ymin-1+0.0001).astype('int'), np.round(ymax+0.0001).astype('int'))
            zc = np.arange(np.round(zmin-1+0.0001).astype('int'), np.round(zmax+0.0001).astype('int'))
            xx = np.round(xc-xc[0]+xmin-x+xw/2).astype('int')
            yy = np.round(yc-yc[0]+ymin-y+yw/2).astype('int')
            zz = np.round(zc-zc[0]+zmin-z+zw/2).astype('int')
            #a = it.product(xx, yy, zz)

        if os.path.isfile(name) and not overwrite:
            print "File %s already exists" % name
        else:
            flux = np.zeros([zw+1, yw+1, xw+1])+float('NaN')
            if white: fluxw = np.zeros([yw+1, xw+1])+float('NaN')
            for ii, iii in zip(xx, xc):
                for jj, jjj in zip(yy, yc):
                    if white: fluxw[jj, ii] = np.nansum(data_cube[:, jjj, iii])
                    for kk, kkk in zip(zz, zc):
                        flux[kk, jj, ii] = data_cube[kkk, jjj, iii]
            hdu.data = flux
            #hdu.writeto(name.replace(".fits", ".c.fits"), clobber=True)
            hdu.writeto(name, clobber=True)
            if cubex:
                os.system("Cube2Im -cube %s[*,*,%d:%d] -out %s %s %s" % (name.replace(".fits", ".c.fits"), zw - zw0, zw + zw0,
                                                                  name.replace(".fits", ".IM.fits"), ssmooth, vb))
            else:
                hdu.data = np.nansum(flux[zw/2-zw0:zw/2+zw0+1, :, :], 0)
                hdu.writeto(name.replace(".fits", ".IM.fits"), clobber=True)

            if white:
                hdu.data = fluxw
                hdu.writeto(name.replace(".fits", "_white.fits"), clobber=True)

            #os.system("CubeSel -cube %s[%d:%d,%d:%d,*] -zmin %d -zmax %d -out %s -multiext .false. -OutVarCube .false."
            #              % (filename, xmin, xmax, ymin, ymax, zmin, zmax, name))

    if makeim:
        # needed to create a new fits file with the right size :/
        scb_label = r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$'
        if fitcat == 'EAGLE':
            if snap == 11:
                #[-45, -27, -9, 0, 9, 27, 45]
                vmin = 0.007
                vmax = 0.7
                #np.round(np.arange(xmin/asec2pix, (xmax+1)/asec2pix, 2)).astype('int')#[-12, -9, -6, -3, 0, 3, 6, 9, 12]#
                #yticks = np.round(np.arange(ymin/com2pix, (ymax+1)/com2pix, 2)).astype('int')#[-1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5]
            if snap == 10:
                #xticks = [-39, -26, -13, 0, 13, 26, 39]
                vmin = 0.01
                vmax = 1
                # np.round(np.arange(xmin/asec2pix, (xmax+1)/asec2pix+2, 2)).astype('int')
            if snap == 12:
                #xticks = [-53, -39, -26, -13, 0, 13, 26, 39, 53]
                vmin = 0.002
                vmax = 2
                # np.round(np.arange(xmin/asec2pix, (xmax+1)/asec2pix+2, 2)).astype('int')
            ntx = 7
            nty = 9
            cbticks = vmax/np.power(10., np.arange(3, -1, -1))
            xticks = np.round(np.linspace(xmin / asec2pix, xmax / asec2pix, ntx), 1)
            yticks = np.round(np.linspace(xmin / com2pix, xmax / com2pix, nty), 1)
            ylabel = r'$\rm{\theta\,[cMpc]}$'
            sb = False
            logscale = True

            #vmin=None
            #vmax=None
        else:
            xticks = np.round(np.arange(x-xw, x+xw+1, tn)).astype('int')
            yticks = np.round(np.arange(y-yw, y+yw+1, tn)).astype('int')
            ylabel = None
            logscale = False
            cbticks = None
            sb = True
        astroim(name.replace('.fits', ".IM.fits"), smooth=smooth, saveim=True, show=False, cbfrac=.08, pad=.006,
                dfig=(7, 7), contours=False, binsize=1, std=std,
                xticks=xticks, yticks=yticks, x0=0, y0=0, vmin=vmin, vmax=vmax, gray=True,
                ylabel=ylabel, sb=sb, logscale=logscale, ntx=ntx, nty=nty, scb_label=scb_label,
                cbticks=cbticks)#xw0=xw/2+.5, y0=yw/2+.5
        if white:
            astroim(name.replace('.fits', "_white.IM.fits"), smooth=smooth, saveim=True, show=False, cbfrac=.08, pad=.006,
                    dfig=(7, 7), contours=True, scb_label='Flux [1e-20 cgs]', binsize=1,  std=std,
                    xticks=xticks, yticks=yticks, x0=xw/2+.5, y0=yw/2+.5, vmin=vmin, vmax=vmax, gray=True)
    if centerfix:
        sn = 4
        catname = name.replace('.fits', '.cat')
        os.system('rm %s' % catname)
        end = False

        while not os.path.isfile(catname) and not end:
            print 'Using SN threshold', sn
            s = 'CubEx -cube %s -sn %d -m .false. -EstVarOutFile %s -n 40 -f .true. -fv .true. -fsr 2 %s' % \
                (name, sn, name.replace('.fits', '.dat'), vb)
            os.system(s)
            sn /= 2
            if sn < .1:
                print "SN threshold too low"
                end = True

        if os.path.isfile(catname):
            fin = open(catname, 'r')
            rrm = []
            ddm = []
            lm = []
            #xm = []
            #ym = []
            #zm = []
            im = []
            fm = []

            for line in fin:
                a = line.split()
                if a[0].isdigit():
                    a[0] = int(a[0])
                    im.append(a[0])
                    rrm.append(a[-4])
                    ddm.append(a[-3])
                    lm.append(a[-2])
                    #zm.append(int(np.round(float(a[7]) - zmin)))
                    #ym.append(int(np.round(float(a[6]) - ymin)))
                    #xm.append(int(np.round(float(a[5]) - xmin)))
                    fm.append(float(a[18]))
            rm, dm = hms2deg(rrm, ddm)
            lm = np.array(lm).astype('float')
            xm, ym, zm = w.all_world2pix(rm, dm, lm, 1)
            zm = l2pix(lm)
            #dd = np.sqrt(np.array(rm-r) ** 2 + np.array(dm-d) ** 2 + np.array(lm-l) ** 2)
            dd = np.sqrt(np.array(zm-z) ** 2 + np.array(ym-y) ** 2 + np.array(xm-x) ** 2)
            a = np.where(dd == min(dd))[0][0]
            idf = im[a]
            xf = xm[a]
            yf = ym[a]
            zf = zm[a]
            rf = rm[a]
            df = dm[a]
            lf = lm[a]
            ff = fm[a]
            xn.append(xf)
            yn.append(yf)
            zn.append(zf)
            fn.append(ff)
            zz = lf/1215.67-1
            cd = comdist(zz)
            ll = 4*np.pi*ff*pow(cd*(1+zz),2)*9.521e28
            fout.write("%d %s %s %1.4f %1.4f %1.4f %1.4f %s %s %s %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n" %
                       (i, r, d, x, y, z, f, rf, df, lf, xf, yf, zf, ff, zz, ll, dd[a], cd))
            print 'Fixed values', xf, yf, zf, data_cube[
                int(zf-1), int(yf-1), int(xf-1)], 'CubEx id selected', idf
            print 'Difference x %1.2f y %1.1f z %1.1f f %1.1f' % (x-xf, y-yf, z-zf,
                data_cube[int(z-1), int(y-1), int(x-1)] - data_cube[int(zf-1), int(yf-1), int(xf-1)])

            xmin = int(max(1 + cut, xf - xw/2))
            xmax = int(min(xlim - cut, xf + xw/2))
            ymin = int(max(1 + cut, yf - yw/2))
            ymax = int(min(ylim - cut, yf + yw/2))
            zmin = int(max(1, zf - zw/2))
            zmax = int(min(zlim, zf + zw/2))
            namefix = folder + "%s/LAEs/lae%d.fixed.fits" % (fitcat, i)
            if 1:
                xc = np.arange(np.round(xmin-1).astype('int'), np.round(xmax).astype('int'))
                yc = np.arange(np.round(ymin-1).astype('int'), np.round(ymax).astype('int'))
                zc = np.arange(np.round(zmin-1).astype('int'), np.round(zmax).astype('int'))
                xx = np.round(xc-xc[0]+xmin-xf+xw/2).astype('int')
                yy = np.round(yc-yc[0]+ymin-yf+yw/2).astype('int')
                zz = np.round(zc-zc[0]+zmin-zf+zw/2).astype('int')
                #a = it.product(xx, yy, zz)
                flux = np.zeros([zw+1, yw+1, xw+1])
                for ii, iii in zip(xx, xc):
                    for jj, jjj in zip(yy, yc):
                        for kk, kkk in zip(zz, zc):
                            flux[kk, jj, ii] = data_cube[kkk, jjj, iii]
                hdu.data = flux
                hdu.writeto(namefix.replace(".fits", ".c.fits"), clobber=True)
            #os.system("CubeSel -cube %s[%d:%d,%d:%d,*] -zmin %d -zmax %d -out %s -multiext .false. -OutVarCube .false." %
            #          (filename, xmin, xmax, ymin, ymax, zmin, zmax, namefix))
            os.system("Cube2Im -cube %s[*,*,%d:%d] -out %s %s %s" % (namefix.replace(".fits", ".c.fits"),
                        zw/2 - zw0, zw/2 + zw0, namefix.replace(".fits", ".IM.fits"), ssmooth, vb))
            if makeim:
                xticks = np.round(np.arange(xf-xw, xf+xw+1, tn)).astype('int')
                yticks = np.round(np.arange(yf-yw, yf+yw+1, tn)).astype('int')
                astroim(namefix.replace('.fits', ".IM.fits"), smooth=smooth, saveim=True, show=False, cbfrac=.08, pad=.006,
                        dfig=(7, 7), contours=True, scb_label='Flux [1e-20 cgs]', binsize=1, std=std,
                        xticks=xticks, yticks=yticks, x0=[xw/2+.5, x+.5], y0=[yw/2+.5, y+.5], vmin=vmin, vmax=vmax)

if centerfix: fout.close()
outf.close()
if white: outw.close()

if simplestack:
    print 'Simple stack', len(fits)
    from tools_sofi import stack
    zw0 = 2
    stack(lst, folder + "%s/LAEs/laes_stack.fits" % fitcat, zw0=zw0)
