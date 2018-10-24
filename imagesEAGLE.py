    #!/usr/bin/env python
__author__ = 'nefinia'
from pyfits import getdata, PrimaryHDU, getheader
import h5py
import numpy as np
import os, sys, time
from math import sqrt, atan2
from tools_sofi import cic, cubex, makejpeg, astroim, rdarg
from sys import argv

binsize = rdarg(argv, 'binsize', int, 2)
coord = rdarg(argv, 'coord', str, 'x')
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)
angle = rdarg(argv, 'angle', bool, False)
csub = rdarg(argv, 'csub', bool, True)
days = rdarg(argv, 'days', float, 7)#1/100000 #1/6 == 4hours
dovar = rdarg(argv, 'dovar', bool, False)
extraname = rdarg(argv, 'extraname', str, '')
fitcat = rdarg(argv, key='fitcat', type=str, default='EAGLE')  # 'UDF
cubetype = rdarg(argv, key='cubetype', type=str, default='SB')  # SB cube
flipx = rdarg(argv, 'flipx', bool, False)
flipy = rdarg(argv, 'flipy', bool, False)
folder = rdarg(argv, 'folder', str, '/net/astrogate/export/astrodata/gallegos/')# '../../')  # '/scratch/gallegos/MUSE/'# '/net/theia/scratch/cantal/public_bin/'
foldercat = rdarg(argv, 'foldercat', str, '/net/astrogate/export/astrodata/gallegos/')  # '/scratch/gallegos/MUSE/'# '/net/theia/scratch/cantal/public_bin/'
half = rdarg(argv, 'half', bool, False)
highq = rdarg(argv, 'highq', bool, False)
i1 = rdarg(argv, 'id1', list, [144], intlist=True)
i2 = rdarg(argv, 'id2', list, [216], intlist=True)
idmin = rdarg(argv, 'idmin', int, 0)
idmax = rdarg(argv, 'idmax', int, 999999)
imtype = rdarg(argv, 'imtype', str, 'flux')  # 'mean'
cubetype = rdarg(argv, 'imtype', str, 'flux')  # 'mean'
jin = rdarg(argv, 'jin', int, 0)
jpeg = rdarg(argv, 'image', bool, True)
laecenter = rdarg(argv, 'laecenter', bool, True)
npixim = rdarg(argv, 'npix', bool, False)
nrand = rdarg(argv, 'nrand', int, 1)
norm = rdarg(argv, 'norm', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
pair = rdarg(argv, 'pair', bool, False)
prename = rdarg(argv, 'prename', str, '')
random = rdarg(argv, 'random', bool, False)
resample = rdarg(argv, 'resample', bool, False)
sb = rdarg(argv, 'sb', bool, False)
scalelims = rdarg(argv, 'scalelims', str, '-0.01 0.01')
sclip = rdarg(argv, 'sclip', bool, False)
show = rdarg(argv, 'show', bool, False)
simon = rdarg(argv, 'cubical', bool, True)
single = rdarg(argv, 'single', bool, False)
snap = rdarg(argv, 'snap', int, 12)
std = rdarg(argv, 'std', float, .00012)#1.25)
unipix = rdarg(argv, 'unipix', bool, False)

xw = rdarg(argv, 'xw', int, 80)
yw = rdarg(argv, 'yw', int, 80)
zw = rdarg(argv, 'zw', int, 10)
xw0 = rdarg(argv, 'xw0', int, 20)
yw0 = rdarg(argv, 'yw0', int, 20)
zw0 = rdarg(argv, 'zw0', int, 3)
vmin = rdarg(argv, 'vmin', float, -.00004)# -.1)
vmax = rdarg(argv, 'vmax', float, .0004)# 1)

#do the x-y collapse to have a spectra in the middle for each pair!

if not extraname: extraname = ''
if not prename: prename = ''
if resample: extraname += '.resamp'
if half: prename += '.half'
if pair: prename += '.pair'
if cubesmooth:
    smooth = False
    ssmooth = ' -boxsm 3 '
else:
    smooth = True
    ssmooth = ''
if sclip: nsclip = 5

if norm: ssnorm = 'normalized'
else: ssnorm = 'not-normalized'

if random:
    sstype = 'random_'
    sstype2 = 'random_'
else:
    sstype = ''
    sstype2 = ''

if simon:
    sstype += 'cubical'
    sstype2 += 'cubical'
else:
    sstype += 'cylindrical'
    sstype2 += 'cylindrical'

data = getdata('../../%s/cats/lae_pairs.fits' % (fitcat), 1)


id1 = data['id1'].astype(int)
id2 = data['id2'].astype(int)
ra1 = data['ra1'].astype(float)
ra2 = data['ra2'].astype(float)
dec1 = data['dec1'].astype(float)
dec2 = data['dec2'].astype(float)
xp1 = data['x1'].astype(int)
yp1 = data['y1'].astype(int)
zp1 = data['z1'].astype(int)
xp2 = data['x2'].astype(int)
yp2 = data['y2'].astype(int)
zp2 = data['z2'].astype(int)
red1 = data['redshift1']
red2 = data['redshift2']
dist = data['com_dist']
theta = data['theta']


ap_conv = 3.125 # *3*1.25/.4 conversion to SB for 1" aperture for a 3 z pix wide (3*1.25 Angstrom)

cubename = '%s_snap%d_%s.fits' % (cubetype, snap, coord)

if single:
    close = False
    for ii1, ii2 in zip(i1, i2):
        close |= ((id1 == int(ii1)) & (id2 == int(ii2)))
else:
    close = np.where((dist <= 20) & (dist > .5) & (theta > 6) & (idmax > id1) & (id1 >= idmin))[0]

angs = np.arctan2(dec2-dec1, ra2-ra1)
lcube = 4096
coml = 25  # cMpc
com2pix = 163.84 # lcube/coml
kpc2pix = lcube / float(coml * 1e3)
red = {'10':3.984, '11':3.528, '12':3.017}
asec2kpcs = {'10':7.842, '11':7.449, '12':7.108}
sb_dimming = np.power(1 + red['%s' % snap], -4)
asec2pix = asec2kpcs['%s' % snap] * (1 + red['%s' % snap]) * kpc2pix
sb2flux = asec2pix ** -2.
sbpeaks = {'10':1.011, '11':1.484, '12':2.396} # for an aperture of 1 asec^2!!! I am converting to flux later on

xmin = -xw
xmax = xw
ymin = -yw
ymax = yw
lenx = int((xw * 2) / binsize) + 1
leny = int((yw * 2) / binsize) + 1
now = time.time()

if norm:
    xmin = -.2
    xmax = 1.2
    lenx = 100

ii1, ii2, dd, zz1, zz2, rr1, rr2, angs = [id1[close], id2[close], dist[close], zp1[close], zp2[close], red1[close], red2[close], angs[close]]

for i1, i2, d, r1, r2, a in zip(ii1, ii2, dd, rr1, rr2, angs):
    print '----------------------------------------------------------------------------------------------'
    print 'id1 %d id2 %d ang %f distance %f' % (i1, i2, a, d)
    if angle: ang = a
    else: ang = None

    def doastroim():
        print 'Generating picture'
        if pair:
            figsize = (int(d)+1, 9)
            arrow = False
            text=[r'$\rm{%s}\,\,%d-%d$' % (fitcat, i1, i2)]
            textpos=[[-xw0*.9, -yw0*.9]]
            xmin, xmax, ymin, ymax = [None]*4
            x0, y0 = [None]*2
        else:
            figsize = (8, 9)
            arrow = False
            text=[r'$\rm{%s}\,\,%d\rightarrow%d$' % (fitcat, i1, i2)]
            textpos=[[-xw0*2, -yw0*2]]
            xmin, xmax, ymin, ymax = [-xw0, xw0, -yw0, yw0]
            x0, y0 = [0, 0]

        if highq: ext = '.pdf'
        else: ext = '.png'
        astroim(outname, smooth=smooth, saveim=True, show=show, cbfrac=.08,
                pad=.006, contours=True, scb_label='Flux', pcolor='white', psize=300, x0=x0, y0=y0,
                pair=pair, gray=True, vmin=vmin, vmax=vmax,
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, arrow=arrow, std=std,
                text=text, textpos=textpos, dfig=figsize, highq=highq, imout=outname.replace('.fits', ext),
                fsize=29, contoursw=1.5, dpi=80, rasterize=False)


    cat = '%s/%s/cats/%s_snap%d_%s/%s_pair%d-%d%s.h5' % (foldercat, fitcat, cubetype, snap, coord, sstype, i1, i2, prename)
    outfolder = '%s/%s/pairs/%s_snap%d_%s/' % (folder, fitcat, cubetype, snap, coord)
    #cat = '%s/%s/cats/%/old_LLScats/%s_pair%d-%d%s.h5' % (foldercat, fitcat, cubetype, snap, coord, sstype, i1, i2, prename)
    #outfolder = '%s/%s/pairs/%s_snap%d_%s/' % (folder, fitcat, cubetype, snap, coord)

    if not os.path.isdir(outfolder): os.system('mkdir %s' % outfolder)

    outname = '%s/%s_%d-%d%s%s.fits' % (outfolder, sstype, i1, i2, prename, extraname)

    print 'Catalog name', cat
    print 'Output name', outname

    if os.path.isfile(cat):
        if os.path.isfile(outname):
            ftime = os.stat(outname).st_mtime
            if now - ftime > days * 86400: print "\n\nOuput file older than %1.2f days!!!!!!!!!!!!!!!!!\n\n" % days
            else: "\nOuput file younger than %1.2f days!!!\n" % days
        else:
            print "Output file does not exist"
            ftime = now + (days-1) * 86400

        if (not os.path.isfile(outname) or overwrite):# and (now - ftime > days * 86400):
            data = h5py.File(cat, 'r')

            xpix = np.array(data['x'])
            ypix = np.array(data['y'])
            y = np.array(data['py'])
            x = np.array(data['px'])
            f = np.array(data['flux'])

            dpix = np.array(data['dpix'])

            inside = (y >= ymin) & (y <= ymax)
            x = x[inside]
            y = y[inside]
            xpix = xpix[inside]
            ypix = ypix[inside]
            f = f[inside]
            if dovar: v = v[inside]
            fmean = np.nanmean(f)
            fstd = np.nanstd(f)

            if sclip:
                good = abs(f - fmean) < nsclip * fstd
                x = x[good]
                y = y[good]
                xpix = xpix[good]
                ypix = ypix[good]
                f = f[good]
                if dovar: v = v[good]

            # check position of LAEs
            if norm:
                if half:
                    _d = .5
                else:
                    _d = 1
            else:
                if half:
                    _d = dpix * .5
                else:
                    _d = dpix

            if dovar: v[v < 0] = float('NaN')
            flux = np.zeros([lenx, leny])
            npix = np.zeros([lenx, leny])
            if unipix: nun = np.zeros([lenx, leny])

            if resample:
                print '%d-%d:' % (i1, i2), 'resample %s' % resample
                xin = x
                yin = y
                fin = f
                a = np.isfinite(fin)

                xn = (xin[a] - xmin) * xmax / float(xmax - xmin)
                yn = (yin[a] - ymin) * ymax / float(ymax - ymin)

                fk, wk, ck = cic(fin[a], xn, lenx, yn, leny, weighting=resample, wraparound=False)

                if imtype == 'mean': flux[:, :] = fk / wk.astype(float)
                if imtype == 'flux': flux[:, :] = fk * np.sqrt(ck) / wk.astype(float)
                npix[:, :] = ck

            else:
                #a = np.isfinite(f)
                #b = np.isfinite(v)
                good = np.where((np.abs(x) < xw) & (np.abs(y) < yw))[0]
                xn = np.round((x[good] - xmin) * (lenx-1) / float(xmax - xmin)).astype(int)#lenx=bug!!! lenx-1 correct
                yn = np.round((y[good] - ymin) * (leny-1) / float(ymax - ymin)).astype(int)
                f = f[good]
                if dovar: v = v[good]
                for m in range(len(xn)):
                    flux[xn[m], yn[m]] += f[m]
                    npix[xn[m], yn[m]] += 1
                    if unipix:
                        aa = int('%04d%04d' % (xpix[m], ypix[m]))
                        nun[xn[m], yn[m]] += aa

            flux[np.isinf(flux)] = float('NaN')
            flux[npix == 0] = float('NaN')

            print '\nPair %d %d done' % (i1, i2)

            flux[flux == -999] = float('NaN')

            hdu = PrimaryHDU()
            hdu.header['ID1'] = (int(i1), 'ID of LAE 1 in the catalogue')
            hdu.header['ID2'] = (int(i2), 'ID of LAE 2 in the catalogue')
            hdu.header['COMDIST'] = (round(d, 4), 'Comoving distance in Mpc')
            hdu.header['COMB'] = (imtype, 'Combination method')
            hdu.header['TYPE'] = ('Projected+Shear', 'Type of coordinate transformation')

            if overwrite:
                print "Warning: Overwriting output fits file"

            # flux fits file
            hdu.data = flux.T
            hdu.writeto(outname, clobber=overwrite)

            # Npixel cube
            if npixim:
                hdu.data = npix.T
                hdu.header['COMMENT'] = 'Pixel map'
                if overwrite: hdu.writeto(outname.replace('.fits', '.NPIX.fits'), clobber=overwrite)

            if unipix:
                hdu.data = nun.T
                hdu.header['COMMENT'] = 'Unique pixels map'
                if overwrite: hdu.writeto(outname.replace('.fits', '.nun.fits'), clobber=overwrite)

            cont = False
            if jpeg:
                if os.path.isfile(outname):
                    doastroim()
                else:
                    print outname, "does not exists!!!!!!!!!!!!!!!"

            if sb:
                hdu = PrimaryHDU()
                fit = getdata(outname)
                hdu.data = fit * ap_conv
                if overwrite: hdu.writeto(outname.replace('.fits', '.SB.fits'), clobber=True)


            del data
            del hdu

        else:
            print 'Computation already done.'

            if jpeg:
                if os.path.isfile(outname):
                    doastroim()

                else:
                    print outname, "does not exists!!!!!!!!!!!!!!!"


            if 0:#jpeg:
                title = r'LAE%d oriented to LAE%d, z=%1.3f, shear $%1.3f\AA$.' % (id1, id2, red1, (red2-red1)*1.25)
                astroim(outname.replace('.fits', '.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                        pad=.006, contours=True, scb_label='Flux', title=title, vmin=-3, vmax=3, x0=xw, y0=yw, gray=True, angle=ang)  # , vmin=-.5, vmax=.5)#,dfig=(8,7))
                # makejpeg(outname.replace('.fits', '.IM.fits'), smooth, outname.replace('.fits', '.jpeg'))

    else:
        print 'No initial catalogue', cat

