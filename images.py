#!/usr/bin/env python
__author__ = 'nefinia'
import os
import time
from sys import argv

import h5py
import numpy as np
from pyfits import getdata, PrimaryHDU

from tools_sofi import cic, astroim, rdarg

binsize = rdarg(argv, 'binsize', int, 2)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)
angle = rdarg(argv, 'angle', bool, False)
coord = rdarg(argv, 'coord', str, 'x')
csub = rdarg(argv, 'csub', bool, True)
cubexmask = rdarg(argv, 'cubexmask', bool, False)
days = rdarg(argv, 'days', float, 1000000000)#1/100000 #1/6 == 4hours
dovar = rdarg(argv, 'dovar', bool, False)
extraname = rdarg(argv, 'extraname', str, '')
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
flipx = rdarg(argv, 'flipx', bool, False)
flipy = rdarg(argv, 'flipy', bool, False)
folder = rdarg(argv, 'folder', str, '/net/astrogate/export/astrodata/gallegos/')#'../../')
if fitcat == 'mosaic' or fitcat == 'mosaic10': folder = '/net/astrogate/export/astrodata/gallegos/'
foldercat = rdarg(argv, 'foldercat', str, '/net/astrogate/export/astrodata/gallegos/')  # '/scratch/gallegos/MUSE/'# '/net/theia/scratch/cantal/public_bin/'
half = rdarg(argv, 'half', bool, False)
hdf5 = rdarg(argv, 'hdf5', bool, True)
highq = rdarg(argv, 'highq', bool, False)
i1 = rdarg(argv, 'id1', list, [144], int)
i2 = rdarg(argv, 'id2', list, [216], int)
idmin = rdarg(argv, 'idmin', int, 0)
idmax = rdarg(argv, 'idmax', int, 999999)
invert = rdarg(argv, 'invert', bool, False)
imtype = rdarg(argv, 'imtype', str, 'flux')  # 'mean'
jin = rdarg(argv, 'jin', int, 0)
jpeg = rdarg(argv, 'image', bool, True)
laecenter = rdarg(argv, 'laecenter', bool, True)
laeunmask = rdarg(argv, 'laeunmask', bool, True) #if False LAEs are not masked
mask = rdarg(argv, 'mask', bool, True)
npixim = rdarg(argv, 'npix', bool, False)
nrand = rdarg(argv, 'nrand', int, 1)
norm = rdarg(argv, 'norm', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
pair = rdarg(argv, 'pair', bool, False)
prename = rdarg(argv, 'prename', str, '')
random = rdarg(argv, 'random', bool, False)
redmin = rdarg(argv, 'rmin', list, [2.8])
redmax = rdarg(argv, 'rmax', list, [4.0])
resample = rdarg(argv, 'resample', bool, False)
sb = rdarg(argv, 'sb', bool, False)
scalelims = rdarg(argv, 'scalelims', str, '-0.01 0.01')
sclip = rdarg(argv, 'sclip', bool, False)
show = rdarg(argv, 'show', bool, False)
simon = rdarg(argv, 'cubical', bool, True)
single = rdarg(argv, 'single', bool, False)
snap = rdarg(argv, 'snap', int, 12)
sncubex = rdarg(argv, 'sncubex', int, 3)
snrim = rdarg(argv, 'snr', bool, False)
std = rdarg(argv, 'std', float, 1.25)
unique = rdarg(argv, 'unique', bool, False)
vartype = rdarg(argv, 'vartype', str, 'PROPVAR')
xw = rdarg(argv, 'xw', int, 80)#x size of the original cube
yw = rdarg(argv, 'yw', int, 80)#y size of the original cube
zw = rdarg(argv, 'zw', int, 10)#z size of the original cube
xw0 = rdarg(argv, 'xw0', int, 40)#x size of the output cube in new coordinates (given the bin size)
yw0 = rdarg(argv, 'yw0', int, 40)
zw0 = rdarg(argv, 'zw0', int, 3)
offset = rdarg(argv, 'offset', int, 0)
vmin = rdarg(argv, 'vmin', float, -.1)
vmax = rdarg(argv, 'vmax', float, 1)

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

if mask: ssmask = 'masked'
else: ssmask = 'unmasked'

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

if fitcat == 'mosaic10': table = open('%s/mosaic/cats/LAE_pairs_info_%s_%s_%s_mosaic10.dat' % (folder, ssnorm, ssmask, sstype), 'w')
else: table = open('%s/%s/cats/LAE_pairs_info_%s_%s_%s.dat' % (folder, fitcat, ssnorm, ssmask, sstype), 'w')
table.write('#ID1 ID2 x1 x2 z1 z2 comdist mean std SB[erg/s/cm**2/arscec**2] 1sigma_SB combtype masked coord_type\n')

if fitcat == 'mosaic10': data = getdata('../../mosaic/cats/lae_pairs_UDF10.fits', 1)
elif fitcat=='EAGLE': data = getdata('../../%s/cats/lae_pairs_snap%d.fits' % (fitcat, snap), 1)
else: data = getdata('../../%s/cats/lae_pairs.fits' % (fitcat), 1)

if 0:#$half or laecenter:
    id1 = np.concatenate((data['id1'], data['id2']), 0).astype(int)
    id2 = np.concatenate((data['id2'], data['id1']), 0).astype(int)
    ra1 = np.concatenate((data['ra1'], data['ra2']), 0).astype(float)
    ra2 = np.concatenate((data['ra2'], data['ra1']), 0).astype(float)
    dec1 = np.concatenate((data['dec1'], data['dec2']), 0).astype(float)
    dec2 = np.concatenate((data['dec2'], data['dec1']), 0).astype(float)
    sconf1 = np.concatenate((data['sconf1'], data['sconf2']), 0).astype(int)
    sconf2 = np.concatenate((data['sconf2'], data['sconf1']), 0).astype(int)
    xp1 = np.concatenate((data['x1'], data['x2']), 0).astype(int)
    yp1 = np.concatenate((data['y1'], data['y2']), 0).astype(int)
    zp1 = np.concatenate((data['z1'], data['z2']), 0).astype(int)
    xp2 = np.concatenate((data['x2'], data['x1']), 0).astype(int)
    yp2 = np.concatenate((data['y2'], data['y1']), 0).astype(int)
    zp2 = np.concatenate((data['z2'], data['z1']), 0).astype(int)
    red1 = np.concatenate((data['redshift1'], data['redshift2']), 0)
    red2 = np.concatenate((data['redshift2'], data['redshift1']), 0)
    dist = np.concatenate((data['com_dist'], data['com_dist']), 0)
    theta = np.concatenate((data['theta'], data['theta']), 0)

id1 = data['id1'].astype(int)
id2 = data['id2'].astype(int)
if fitcat == 'EAGLE':
    ra1 = data['x1']
    ra2 = data['x2']
    dec1 = data['y1']
    dec2 = data['y2']
    sconf1 = np.zeros(len(id1))+1
    sconf2 = np.zeros(len(id1))+1
    red1 = np.zeros(len(id1))+3.5
    red2 = np.zeros(len(id1))+3.5
    theta = data['theta%s' % coord]

else:
    ra1 = data['ra1'].astype(float)
    ra2 = data['ra2'].astype(float)
    dec1 = data['dec1'].astype(float)
    dec2 = data['dec2'].astype(float)
    sconf1 = data['sconf1'].astype(int)
    sconf2 = data['sconf2'].astype(int)
    red1 = data['redshift1']
    red2 = data['redshift2']
    theta = data['theta']
xp1 = data['x1'].astype(int)
yp1 = data['y1'].astype(int)
zp1 = data['z1'].astype(int)
xp2 = data['x2'].astype(int)
yp2 = data['y2'].astype(int)
zp2 = data['z2'].astype(int)
dist = data['com_dist']


if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
    #filevar = '%s/%s/%s' % (folder, fitcat, cubename)
    if dovar: filevar = '%s/%s/DATACUBE-HDFS-1.35-%s.fits' % (folder, fitcat, vartype)
if fitcat == 'UDF':
    cubename = 'DATACUBEFINALuser_20141021T055145_72f74684.fits'
    #cubename = 'DATACUBEFINALuser_20141021T055145_212058ad.fits'
    #filevar = '%s/%s/%s' % (folder, fitcat, cubename)
    if dovar: filevar = '%s/%s/%s' % (folder, fitcat, 'DATACUBEFINALuser_20141021T055145_212058ad.fits')

if fitcat=='EAGLE':
    print 'aaaa'
    cubename = 'SBflux.fits'
    mask = False

if fitcat == 'mosaic' or fitcat == 'mosaic10':
    cubename = '/net/astrogate/export/astrodata/gallegos/mosaic.z1300.fits'
    if dovar: filevar = '/net/astrogate/export/astrodata/gallegos/mosaic.z1300.VAR.fits'
    fitcat = 'mosaic'
if csub:
    cubename = cubename.replace('.fits', '.csub.fits')

if dovar and fitcat != 'mosaic': datavar = getdata(filevar, 2)
if single:
    close = np.where((id1 == i1) & (id2 == i2))[0]
else:
    close = np.where((dist <= 20) & (dist > .5) & (theta > 6) & (sconf1 >= 1) & (sconf2 >= 1)
                     & (red1 < redmax) & (red2 < redmax) & (red1 > redmin) & (red2 > redmin)
                     & (idmax > id1) & (id1 >= idmin))[0]

angs = np.arctan2(dec2-dec1, ra2-ra1)

if mask:
    # Open mask file
    if fitcat == 'HDFS':
        fits = '%s/%s/DATACUBE-HDFS-1.35-STATVAR.IM.Objects_Id.fits' % ('../../', fitcat)
    if fitcat == 'UDF':
        fits = '%s/%s/DATACUBEFINALuser_20141021T055145_212058ad.IM.Objects_Id.fits' % ('../../', fitcat)
    if fitcat == 'mosaic':
        fits = '/net/astrogate/export/astrodata/gallegos/DATACUBE_UDF-MOSAIC.z1300.IM.Objects_Id.fits'
    data_mask, header_data_mask = getdata(fits, 0, header=True)
    data_mask = data_mask[0]
    # unmasked = np.where((data_mask == 0))

    if laeunmask:
        # remove LAE from the mask catalogue
        laeid = {}
        xps = np.concatenate((xp1, [xp2[-1]]), 0)
        yps = np.concatenate((yp1, [yp2[-1]]), 0)
        ids = np.concatenate((id1, [id2[-1]]), 0)

        for x, y, ii in zip(xps, yps, ids):
            yyl, xxl = data_mask.shape
            if (x < xxl) & (y < yyl)& (0 < x) & (0 < y):
                if data_mask[y - 1, x - 1] > 0:
                    laeid['%d' % ii] = data_mask[y - 1, x - 1]
                else:
                    laeid['%d' % ii] = 999

    # Open sky file
    sky = np.loadtxt('%s/%s/skyspec.txt' % ('../../', fitcat)).T
    # z = sky[0].astype('int')
    # wav = sky[1].astype('float')
    skyflux = sky[2].astype('float')
    if len(skyflux) > 1300:
        skyflux = skyflux[:1300]
    skymean = np.mean(skyflux)
    skystd = np.std(skyflux)
    zmask = np.where(abs(skyflux - skymean) < skystd)[0] + 1

xmin = -xw
xmax = xw
ymin = -yw
ymax = yw
zmin = -zw#min(z)
zmax = zw#max(z)
lenx = int(xw0 * 2) + 1
leny = int(yw0* 2) + 1
now = time.time()

if norm:
    xmin = -.2
    xmax = 1.2
    lenx = 100

if simon:
    lenz = int((zw0 * 2)) + 1
    #zw = '%d:%d' % (lenz / 2 - zw0 + 1, lenz / 2 + zw0 + 1)

if invert:
    ii1, ii2, dd, zz1, zz2, angs = [id1[close][::-1], id2[close][::-1], dist[close][::-1], zp1[close][::-1], zp2[close][::-1], angs[close][::-1]]
else:
    ii1, ii2, dd, zz1, zz2, angs = [id1[close], id2[close], dist[close], zp1[close], zp2[close], angs[close]]

for i1, i2, d, z1, z2, a in zip(ii1, ii2, dd, zz1, zz2, angs):
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
            textpos=[[-xw0*.9, -yw0*.87]]
            xmin, xmax, ymin, ymax = [-xw0, xw0, -yw0, yw0]
            x0, y0 = [0, 0]

        if highq: ext = '.pdf'
        else: ext = '.png'
        astroim(outname.replace('.fits', '.IM.fits'), smooth=smooth, saveim=True, show=show, cbfrac=.08,
                pad=.006, contours=True, scb_label='Flux', x0=x0, y0=y0, pcolor='white', psize=300,
                pair=pair, gray=True, angle=ang, vmin=vmin, vmax=vmax,
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, arrow=arrow, std=std,
                text=text, textpos=textpos, dfig=figsize, highq=highq, imout=outname.replace('.fits', ext),
                fsize=29, contoursw=1.5, dpi=80, rasterize=False)

    for j in range(nrand):
        print '----------------------------------------------------------------------------------------------'

        if random:
            j += jin
            cat = '%s/%s/cats/%s_pair%d-%d%s.%s.h5' % (foldercat, fitcat, sstype, i1, i2, prename, j)
            fold = '%s/%s/pairs/%s/%s/random%d' % (folder, fitcat, ssnorm, ssmask, j)
            if not os.path.isdir(fold):
                os.system('mkdir %s' % fold)
            outname = '%s/%s_pair%d-%d%s%s.%s.fits' % (fold, sstype, i1, i2, prename, extraname, j)
        else:
            if fitcat == 'EAGLE':
                cat = '%s/%s/cats/SB_snap%d_%s/%s_pair%d-%d%s.h5' % (foldercat, fitcat, snap, coord, sstype, i1, i2, prename)
                outname = '%s/%s/pairs/SB_snap%d_%s/%s_%d-%d%s%s.fits' % (folder, fitcat, snap, coord, sstype, i1, i2, prename, extraname)
            else:
                cat = '%s/%s/cats/%s_pair%d-%d%s.h5' % (foldercat, fitcat, sstype, i1, i2, prename)
                outname = '%s/%s/pairs/%s/%s/%s_%d-%d%s%s.fits' % (folder, fitcat, ssnorm, ssmask, sstype, i1, i2, prename, extraname)

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

            if (not os.path.isfile(outname) or overwrite) or (now - ftime > days * 86400):
                if hdf5: data = h5py.File(cat, 'r')
                else: data = np.loadtxt(cat.replace('.h5', '.dat')).T

                if mask:
                    print 'Defining good pixels'
                    xpix = np.array(data['x'])
                    ypix = np.array(data['y'])

                    a = (data_mask[ypix - 1, xpix - 1] == 0)
                    if laeunmask: a |= (data_mask[ypix - 1, xpix - 1] == laeid['%d' % i1])
                                 #| (data_mask[ypix - 1, xpix - 1] == laeid['%d' % i2]))[0]
                    if len(a) == 0:
                        continue
                    y = np.array(data['py'])[a]
                    x = np.array(data['px'])[a]
                    xpix = xpix[a]
                    ypix = ypix[a]
                    zpix = np.array(data['z'])[a]
                    f = np.array(data['flux'])[a]
                    if dovar: v = np.array(data['variance'])[a]

                    if simon:
                        z = np.array(data['pz'])[a]
                    good = np.array([])
                    for i in np.unique(zpix):
                        if i in zmask:
                            zgood = np.where(zpix == i)[0]
                            good = np.concatenate((good, zgood), axis=0)
                    good = good.astype(int)
                    x = x[good]
                    y = y[good]
                    if simon:
                        z = z[good]
                    xpix = xpix[good]
                    ypix = ypix[good]
                    zpix = zpix[good]
                    f = f[good]
                    if dovar: v = v[good]
                    good = np.isfinite(f)
                    print 'End defining good pixels'

                else:
                    if hdf5:
                        if unique:
                            xpix = np.array(data['x'])
                            ypix = np.array(data['y'])
                            zpix = np.array(data['z'])
                        x = np.array(data['px'])
                        y = np.array(data['py'])
                        if simon:
                            z = np.array(data['pz'])
                        f = np.array(data['flux'])
                        if dovar: v = np.array(data['variance'])
                    else:
                        x, y, z, f = data

               # dpix = np.array(data['dpix'])

                #Check this part again!!!!!!!!!!
                if 0:#simon:
                    inside = (z >= zmin) & (z <= zmax) & (y >= ymin) & (y <= ymax) & (f > -20)
                    x = x[inside]
                    y = y[inside]
                    z = z[inside]
                    if unique:
                        xpix = xpix[inside]
                        ypix = ypix[inside]
                        zpix = zpix[inside]
                    f = f[inside]
                    if dovar: v = v[inside]
                    fmean = np.nanmean(f)
                    fstd = np.nanstd(f)

                    if sclip:
                        good = abs(f - fmean) < nsclip * fstd
                        x = x[good]
                        y = y[good]
                        z = z[good]
                        if unique:
                            xpix = xpix[good]
                            ypix = ypix[good]
                            zpix = zpix[good]
                        f = f[good]
                        if dovar: v = v[good]


                if dovar: v[v < 0] = float('NaN')
                if simon:  # cubical coords
                    pz = range(lenz)
                    flux = np.zeros([lenz, leny, lenx])
                    if dovar: variance = np.zeros([lenz, leny, lenx])
                    npix = np.zeros([lenz, leny, lenx])
                    if unique: nun = np.zeros([lenz, leny, lenx])


                    if resample:
                        for k in pz:
                            print '%d-%d:' % (
                            i1, i2), 'z', k + 1, 'of', lenz, 'mask', mask, 'norm', norm, 'resample %s' % resample
                            zgood = z == (k + zmin)
                            lzgood = zgood.sum()
                            xin = x[zgood]
                            yin = y[zgood]
                            zin = z[zgood]
                            fin = f[zgood]
                            a = np.isfinite(fin)
                            if dovar:
                                vin = v[zgood]
                                b = np.isfinite(vin)

                            if lzgood > 0:
                                # remember 4 = 2 ** rank -> in this case rank = 2 (2 dims)
                                xn = (xin[a] - xmin) * xmax / float(xmax - xmin)
                                yn = (yin[a] - ymin) * ymax / float(ymax - ymin)

                                fk, wk, ck = cic(fin[a], xn, lenx, yn, leny, weighting=resample, wraparound=False)
                                if dovar: vk, wk, ck = cic(vin[a], xn, lenx, yn, leny, weighting=resample, wraparound=False)

                                if imtype == 'mean':
                                    flux[:, :, k] = fk / wk.astype(float)
                                    if dovar: variance[:, :, k] = vk / wk.astype(float) / ck.astype(float)
                                    #if resample:
                                    #    flux *= 2
                                    #    variance *= 2

                                if imtype == 'flux':
                                    flux[:, :, k] = fk * np.sqrt(ck) / wk.astype(float)
                                    if dovar: variance[:, :, k] = vk * ck.astype(float) / wk.astype(float)

                                npix[:, :, k] = ck

                            x = x[~zgood]
                            y = y[~zgood]
                            z = z[~zgood]
                            f = f[~zgood]
                            if dovar: v = v[~zgood]
                            #xpix = xpix[~zgood]
                            #ypix = ypix[~zgood]
                            #zpix = zpix[~zgood]
                    else:

                        #a = np.isfinite(f)
                        #b = np.isfinite(v)
                        good =(np.abs(x)/binsize <= xw0) & (np.abs(y)/binsize <= yw0) & (np.abs(z) <= zw0)
                        xn = np.round(x[good]/binsize+xw0).astype(int)
                        yn = np.round(y[good]/binsize+yw0).astype(int)
                        zn = np.round(z[good]+zw0).astype(int)
                        #xn = np.round((x[good] - xmin) * (lenx-1) / float(xmax - xmin)).astype(int)#lenx=bug!!! lenx-1 correct
                        #yn = np.round((y[good] - ymin) * (leny-1) / float(ymax - ymin)).astype(int)
                        #zn = (z[good] + zw0).astype(int)
                        f = f[good]
                        if dovar: v = v[good]
                        #xn[xn > xmax] = xmax
                        #yn[yn > ymax] = ymax
                        #zn[zn > 2*zmax] = 2*zmax

                        for m in range(len(xn)):
                            flux[zn[m], yn[m], xn[m]] += f[m]
                            if dovar: variance[zn[m], yn[m], xn[m]] += v[m]
                            npix[zn[m], yn[m], xn[m]] += 1
                            if unique:
                                aa = int('%03d%03d%04d' % (xpix[m], ypix[m], zpix[m]))
                                nun[zn[m], yn[m], xn[m]] += aa

                else:  # cylindrical coords
                    fk, ck = cic(f, (x - xmin) / (xmax - xmin), lenx, (y - ymin) / (ymax - ymin), leny)
                    flux[:, k] = fk[:, ::-1] / ck[:, ::-1].astype(float)
                    if dovar:
                        vk, ck = cic(v, (x - xmin) / (xmax - xmin), lenx, (y - ymin) / (ymax - ymin), leny)
                        variance[:, k] = vk[:, ::-1] / ck[:, ::-1].astype(float) ** 2
                    npix[:, :] = ck[:, ::-1]

                flux[np.isinf(flux)] = float('NaN')
                flux[npix == 0] = float('NaN')
                if dovar:
                    variance[np.isinf(variance)] = float('NaN')
                    variance[npix == 0] = float('NaN')

                print '\nPair %d %d done' % (i1, i2)

                flux[flux == -999] = float('NaN')
                if dovar: variance[variance < 0] = float('NaN')

                hdu = PrimaryHDU()
                hdu.header['ID1'] = (int(i1), 'ID of LAE 1 in the catalogue')
                hdu.header['ID2'] = (int(i2), 'ID of LAE 2 in the catalogue')
                #hdu.header['z1'] = (int(z1), 'LAE 1 z position in the original cube')
                #hdu.header['z2'] = (int(z2), 'LAE 2 z position in the original cube')
                hdu.header['MASKED'] = (mask, 'If true, a mask was applied to the continuum sources')
                hdu.header['COMDIST'] = (round(d, 4), 'Comoving distance in Mpc')
                hdu.header['COMB'] = (imtype, 'Combination method')

                if simon:
                    hdu.header['TYPE'] = ('Projected+Shear', 'Type of coordinate transformation')
                else:
                    hdu.header['TYPE'] = ('Cylindrical', 'Type of coordinate transformation')

                if overwrite:
                    print "Warning: Overwriting output fits file"

                # flux fits file
                hdu.data = flux
                hdu.writeto(outname, clobber=True)
                hdu.data = flux[zw0/2, :, :]
                hdu.writeto(outname.replace('.fits', '.IM.fits'), clobber=True)
                # variance cube
                if dovar:
                    hdu.data = variance.T
                    hdu.header['COMMENT'] = 'Propagated variance in new coordinates'
                    if overwrite: hdu.writeto(outname.replace('.fits', '.%s.fits' % vartype), clobber=overwrite)

                # Npixel cube
                if npixim:
                    hdu.data = npix.T
                    hdu.header['COMMENT'] = 'Pixel map'
                    if overwrite: hdu.writeto(outname.replace('.fits', '.NPIX.fits'), clobber=overwrite)

                # unique pixels cube
                if unique:
                    hdu.data = nun.T
                    hdu.header['COMMENT'] = 'Unique pixels map'
                    if overwrite: hdu.writeto(outname.replace('.fits', '.nun.fits'), clobber=overwrite)

                # snr cube
                if snrim and dovar:
                    hdu.data = (flux / np.sqrt(variance)).T
                    hdu.header['COMMENT'] = 'SNR with %s in new coordinates' % vartype
                    if overwrite: hdu.writeto(outname.replace('.fits', '.SNR.fits'), clobber=overwrite)

                    if simon:
                        # flux image and snr image
                        s = 'Cube2Im -cube %s[*,*,%s] -snrmap %s -varcube %s %s -imtype %s -writeNaN .true.' % (
                            outname, zw, outname.replace('.fits', '.SNR.IM.fits'),
                            outname.replace('.fits', '.%s.fits' % vartype), ssmooth, imtype)
                #else:
                #    s = 'Cube2Im -cube %s %s' % (outname, ssmooth)
                #    os.system(s)

                cont = False
                if jpeg:
                    if os.path.isfile(outname.replace('.fits', '.IM.fits')):
                        doastroim()

                        if snrim:
                            astroim(outname.replace('.fits', '.SNR.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                                    pad=.006, contours=True, scb_label='SNR', vmin=-1, vmax=1,
                                    x0=lenx/2+.5, y0=leny/2+.5, pcolor='white', pair=pair)
                    else:
                        print outname.replace('.fits', '.IM.fits'), "does not exists!!!!!!!!!!!!!!!"

                if sb:
                    hdu = PrimaryHDU()
                    fit = getdata(outname.replace('.fits', '.IM.fits'))
                    hdu.data = fit * (
                    zw0 * 2 + 1) * 3.125  # *3*1.25/.4 conversion to SB for 1" aperture for a 3 z pix wide (3*1.25 Angstrom)
                    if overwrite: hdu.writeto(outname.replace('.fits', '.SB.fits'), clobber=True)

                if cubexmask:  # sb:
                    fit = flux.T
                    fit[np.isnan(fit)] = -999
                    hdu.data = fit
                    if overwrite: hdu.writeto(outname.replace('.fits', '.nonan.fits'), clobber=True)
                    if dovar: fit = variance.T
                    fit[np.isnan(fit)] = -999
                    hdu.data = fit
                    if overwrite: hdu.writeto(outname.replace('.fits', '.PROPVAR.nonan.fits'), clobber=True)
                    #s = 'CubEx -cube %s -var %s -sn %d -m .false. -cct Residuals' % (
                    #outname, outname.replace('.fits', '.%s.fits' % vartype), sncubex)
                    s = 'CubEx -cube %s -n 20 -sn %d -m .false. -cct Residuals' % (
                    outname.replace('.fits', '.nonan.fits'), sncubex)
                    print s
                    os.system(s)
                    s = 'Cube2Im -cube %s %s' % (outname.replace('.fits', '.Residuals.fits'), ssmooth)
                    print s
                    os.system(s)
                    title = r'LAE%d oriented to LAE%d, z=%1.3f, shear $%1.3f\AA$.' % (id1, id2, red1, (red2-red1)*1.25)
                    astroim(outname.replace('.fits', '.Residuals.IM.fits'), smooth=smooth, saveim=True, show=show,
                            cbfrac=.08, pad=.006, contours=True, scb_label='Flux', title=title, vmin=vmin, vmax=vmax,
                            angle=ang, pair=pair, gray=True)#, x0=lenx/2, y0=leny/2)  # , vmin=-.5, vmax=.5)#,dfig=(8,7))
                    # hdu.data = fit*6.25/binsize#zw*1.25/(binsize*0.2) conversion to SB for 1" aperture (originally .2") for a zw z pix wide (zw*1.25 Angstrom)
                    # hdu.writeto(outname.replace('.fits', '.SB.fits'), clobber=True)


                del data
                del hdu

            else:
                print 'Computation already done.'

                if 0:#not os.path.isfile(outname.replace('.fits', '.SNR.IM.fits')):
                    if simon:
                        # flux image and snr image
                        s = 'Cube2Im -cube %s[*,*,%s] -snrmap %s -varcube %s %s -imtype %s -writeNaN .true.' % (
                            outname, zw, outname.replace('.fits', '.SNR.IM.fits'),
                            outname.replace('.fits', '.%s.fits' % vartype), ssmooth, imtype)
                    else:
                        s = 'Cube2Im -cube %s %s' % (outname.replace('.fits', '.SNR.fits'), ssmooth)

                    print s
                    os.system(s)

                if 0:#jpeg:

                    if os.path.isfile(outname.replace('.fits', '.IM.fits')):
                        doastroim()

                        if snrim:
                            astroim(outname.replace('.fits', '.SNR.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                                    pad=.006, contours=True, scb_label='SNR', vmin=-1, vmax=1, angle=ang, pair=pair,
                                    gray=True)#, x0=lenx/2, y0=leny/2)
                    else:
                        print outname.replace('.fits', '.IM.fits'), "does not exists!!!!!!!!!!!!!!!"
                        s = 'Cube2Im -cube %s %s' % (outname, ssmooth)
                        os.system(s)
                        doastroim()
                if cubexmask:  # sb:
                    if 1:  # not os.path.isfile(outname.replace('.fits', '.Residuals.fits')):
                        hdu = PrimaryHDU()
                        fit = getdata(outname.replace('.fits', '.IM.fits'))
                        fit[np.isnan(fit)] = -999
                        hdu.data = fit[::-1, :]
                        hdu.writeto(outname.replace('.fits', '.nonan.fits'), clobber=True)
                        fit = getdata(outname.replace('.fits', '.PROPVAR.fits'))
                        fit[np.isnan(fit)] = -999
                        if simon:
                            hdu.data = fit[:, ::-1, :]
                        else:
                            hdu.data = fit[::-1, :]
                        hdu.writeto(outname.replace('.fits', '.PROPVAR.nonan.fits'), clobber=True)
                        s = 'CubEx -cube %s -var %s -n 20 -sn %d -m .false. -cct Residuals -out %s' % (
                        outname.replace('.fits', '.nonan.fits'),
                        outname.replace('.fits', '.%s.nonan.fits' % vartype),
                        sncubex, outname.replace('.fits', '.Residuals.fits'))
                        print s
                        os.system(s)
                    s = 'Cube2Im -cube %s %s' % (outname.replace('.fits', '.Residuals.fits'), ssmooth)
                    print s
                    os.system(s)
                    astroim(outname.replace('.fits', '.Residuals.IM.fits'), smooth=smooth, saveim=True, show=False,
                            cbfrac=.08, pad=.06, contours=True, scb_label='Flux', vmin=-.05, vmax=.5, gray=True, angle=ang)
                    # , vmin=-.5, vmax=.5)#,dfig=(8,7))
                    # hdu.data = fit*6.25/binsize#zw*1.25/(binsize*0.2) conversion to SB for 1" aperture (originally .2") for a zw z pix wide (zw*1.25 Angstrom)
                    # hdu.writeto(outname.replace('.fits', '.SB.fits'), clobber=True)

                if 0:#jpeg:
                    title = r'LAE%d oriented to LAE%d, z=%1.3f, shear $%1.3f\AA$.' % (id1, id2, red1, (red2-red1)*1.25)
                    astroim(outname.replace('.fits', '.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                            pad=.006, contours=True, scb_label='Flux', title=title, vmin=-3, vmax=3, x0=xw, y0=yw, gray=True, angle=ang)  # , vmin=-.5, vmax=.5)#,dfig=(8,7))
                    # makejpeg(outname.replace('.fits', '.IM.fits'), smooth, outname.replace('.fits', '.jpeg'))

        else:
            print 'No initial catalogue', cat

table.close()
