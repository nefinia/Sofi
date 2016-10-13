#!/usr/bin/env python
__author__ = 'nefinia'
from pyfits import getdata, PrimaryHDU, getheader
import h5py
import numpy as np
import os, sys
from math import sqrt
from tools_sofi import cic, cubex, makejpeg, astroim, rdarg
from sys import argv

half = rdarg(argv, 'half', bool, True)
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
nodist = rdarg(argv, 'nodist', bool, True)
xw = rdarg(argv, 'xw', int, 80)
yw = rdarg(argv, 'yw', int, 80)
zw = rdarg(argv, 'zw', int, 10)
extraname = rdarg(argv, 'extraname', str, '')
if not extraname:
    extraname = ''
scalelims = rdarg(argv, 'scalelims', str, '-0.01 0.01')
cubexmask = rdarg(argv, 'cubexmask', bool, False)
sncubex = rdarg(argv, 'sncubex', int, 3)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
simon = rdarg(argv, 'cubical', bool, True)
random = rdarg(argv, 'random', bool, False)
nrand = rdarg(argv, 'nrand', int, 1)
jin = rdarg(argv, 'jin', int, 0)
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'# '/net/theia/scratch/cantal/public_bin/'
foldercat = rdarg(argv, 'foldercat', str, '/net/astrogate/export/astrodata/gallegos/cats/')  # '/scratch/gallegos/MUSE/'# '/net/theia/scratch/cantal/public_bin/'
flipx = rdarg(argv, 'flipx', bool, False)
flipy = rdarg(argv, 'flipy', bool, False)
mask = rdarg(argv, 'mask', bool, True)
single = rdarg(argv, 'single', bool, False)
i1 = rdarg(argv, 'id1', int, 144)
i2 = rdarg(argv, 'id2', int, 216)
jpeg = rdarg(argv, 'image', bool, True)
norm = rdarg(argv, 'norm', bool, False)
imtype = rdarg(argv, 'imtype', str, 'flux')  # 'mean'
sb = rdarg(argv, 'sb', bool, False)
binsize = rdarg(argv, 'binsize', int, 2)
zw0 = rdarg(argv, 'zw0', int, 1)
vartype = rdarg(argv, 'vartype', str, 'PROPVAR')
resample = rdarg(argv, 'resample', bool, False)
npixim = rdarg(argv, 'npix', bool, False)
snrim = rdarg(argv, 'snr', bool, False)
invert = rdarg(argv, 'invert', bool, False)
prename = rdarg(argv, 'prename', str, '')
if not prename:
    prename = ''


csub = True
sclip = False  # sigma clipping

if cubesmooth:
    smooth = False
    ssmooth = ' -boxsm 3 '
else:
    smooth = True
    ssmooth = ''
if sclip:
    nsclip = 5

if norm:
    ssnorm = 'normalized'
else:
    ssnorm = 'not-normalized'

if mask:
    ssmask = 'masked'
else:
    ssmask = 'unmasked'

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

table = open('%s/%s/cats/LAE_pairs_info_%s_%s_%s.dat' % (folder, fitcat, ssnorm, ssmask, sstype), 'w')
table.write('#ID1 ID2 x1 x2 z1 z2 comdist mean std SB[erg/s/cm**2/arscec**2] 1sigma_SB combtype masked coord_type\n')
data = getdata('%s/%s/cats/lae_pairs.fits' % (folder, fitcat), 1)


if half:
    id1 = np.concatenate((data['id1'], data['id2']), 0).astype(int)
    id2 = np.concatenate((data['id2'], data['id1']), 0).astype(int)
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
    dist = np.concatenate((data['pi_Mpc'], data['pi_Mpc']), 0)
    theta = np.concatenate((data['theta'], data['theta']), 0)
else:
    id1 = data['id1'].astype(int)
    id2 = data['id2'].astype(int)
    sconf1 = data['sconf1'].astype(int)
    sconf2 = data['sconf2'].astype(int)
    xp1 = data['x1'].astype(int)
    yp1 = data['y1'].astype(int)
    zp1 = data['z1'].astype(int)
    xp2 = data['x2'].astype(int)
    yp2 = data['y2'].astype(int)
    zp2 = data['z2'].astype(int)
    red1 = data['redshift1']
    red2 = data['redshift2']
    dist = data['pi_Mpc']
    theta = data['theta']

if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR.fit'
    filevar = '%s/%s/%s' % (folder, fitcat, cubename)
    filevar = '%s/%s/DATACUBE-HDFS-1.35-%s.fits' % (folder, fitcat, vartype)
if fitcat == 'UDF':
    cubename = 'DATACUBEFINALuser_20141021T055145_212058ad.fits'
    filevar = '%s/%s/%s' % (folder, fitcat, cubename)
if csub:
    cubename = cubename.replace('.fits', '.csub.fits')

datavar = getdata(filevar, 2)
# close = np.where(([id1[i] in good for i in range(len(id1))]) & (dist <= 20) & (dist > .5))[0]
if single:
    close = np.where((id1 == i1) & (id2 == i2))[0]
else:
    close = np.where((dist <= 20) & (dist > .5) & (theta > 6))[0]# & (red1 < 4) & (red2 < 4))[0]  # & (sconf1 >= 1) & (sconf2 >= 1))[0]
# close = np.where((id1 >= 499))[0]

if mask:
    # Open mask file
    if fitcat == 'HDFS':
        fits = '%s/%s/DATACUBE-HDFS-1.35-STATVAR.IM.Objects_Id.fits' % (folder, fitcat)
    if fitcat == 'UDF':
        fits = '%s/%s/DATACUBEFINALuser_20141021T055145_212058ad.IM.Objects_Id.fits' % (folder, fitcat)
    data_cube, header_data_cube = getdata(fits, 0, header=True)
    data_cube = data_cube[0]
    # unmasked = np.where((data_cube == 0))

    if 1:
        # remove LAE from the mask catalogue
        laeid = {}
        xps = np.concatenate((xp1, [xp2[-1]]), 0)
        yps = np.concatenate((yp1, [yp2[-1]]), 0)
        ids = np.concatenate((id1, [id2[-1]]), 0)
        for x, y, ii in zip(xps, yps, ids):
            if data_cube[y - 1, x - 1] > 0:
                laeid['%d' % ii] = data_cube[y - 1, x - 1]
            else:
                laeid['%d' % ii] = -1
                # print 'LAE found in mask! id', ii
                # laeids = np.where((data_cube == 0) | ([data_cube == i for i in laeid]))
                # unmasked = np.where([data_cube == i for i in laeid])
                # data_cube[unmasked[1], unmasked[2]] = 0

    # Open sky file
    sky = np.loadtxt('%s/%s/skyspec.txt' % (folder, fitcat)).T
    # z = sky[0].astype('int')
    # wav = sky[1].astype('float')
    skyflux = sky[2].astype('float')
    skymean = np.mean(skyflux)
    skystd = np.std(skyflux)
    zmask = np.where(abs(skyflux - skymean) < skystd)[0] + 1

if norm:
    xmin = -.2
    xmax = 1.2
    lenx = 100
else:
    xmin = -20
    ymax = 80
if simon:
    ymin = -ymax
else:
    ymin = 0

xmin = -xw
xmax = xw
ymin = -yw
ymax = yw
zmin = -zw#min(z)
zmax = zw#max(z)
lenx = int((xw * 2 ) / binsize) + 1
leny = int((yw * 2 ) / binsize) + 1
if simon:
    lenz = int((zw * 2)) + 1
    zw = '%d:%d' % (lenz / 2 - zw0 + 1, lenz / 2 + zw0 + 1)

if invert:
    ii1, ii2, dd, zz1, zz2, rr1, rr2 = [id1[close][::-1], id2[close][::-1], dist[close][::-1], zp1[close][::-1],
                                        zp2[close][::-1], red1[close][::-1], red2[close][::-1]]
else:
    ii1, ii2, dd, zz1, zz2, rr1, rr2 = [id1[close], id2[close], dist[close], zp1[close], zp2[close], red1[close], red2[close]]

for i1, i2, d, z1, z2, r1, r2 in zip(ii1, ii2, dd, zz1, zz2, rr1, rr2):
    print '----------------------------------------------------------------------------------------------'

    for j in range(nrand):
        print '----------------------------------------------------------------------------------------------'

        if random:
            j += jin
            cat = '%s/%s/%s_pair%d-%d%s.%s.h5' % (foldercat, fitcat, sstype, i1, i2, prename, j)
            outname = '%s/%s/pairs/%s/%s/%s_pair%d-%d%s%s.%s.fits' % (folder, fitcat, ssnorm, ssmask, sstype, i1, i2, prename, extraname, j)
        else:
            cat = '%s/%s/%s_pair%d-%d%s.h5' % (foldercat, fitcat, sstype, i1, i2, prename)
            outname = '%s/%s/pairs/%s/%s/%s_pair%d-%d%s%s.fits' % (folder, fitcat, ssnorm, ssmask, sstype, i1, i2, prename, extraname)

        print 'Catalog name', cat
        print 'Output name', outname

        if os.path.isfile(cat):
            if not os.path.isfile(outname) or overwrite:
                data = h5py.File(cat, 'r')

                if mask:
                    print 'Defining good pixels'
                    xpix = np.array(data['x'])
                    ypix = np.array(data['y'])
                    a = np.where((data_cube[ypix - 1, xpix - 1] == 0) | (data_cube[ypix - 1, xpix - 1] == laeid['%d' % i1]))[0]
                                 #| (data_cube[ypix - 1, xpix - 1] == laeid['%d' % i2]))[0]
                    if len(a) == 0:
                        continue
                    y = np.array(data['py'])[a]
                    x = np.array(data['px'])[a]
                    xpix = xpix[a]
                    ypix = ypix[a]
                    zpix = np.array(data['z'])[a]
                    f = np.array(data['flux'])[a]
                    v = np.array(data['variance'])[a]

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
                    v = v[good]
                    good = np.isfinite(f)
                    print 'End defining good pixels'

                else:
                    xpix = np.array(data['x'])
                    ypix = np.array(data['y'])
                    zpix = np.array(data['z'])
                    y = np.array(data['py'])
                    x = np.array(data['px'])
                    if simon:
                        z = np.array(data['pz'])
                    f = np.array(data['flux'])
                    v = np.array(data['variance'])

                dpix = np.array(data['dpix'])

                if not nodist:
                    xmin = min(x) * .99  # int(-.2*dpix)
                    xmax = max(x) * 1.01  # int(1.2*dpix)
                    ymin = min(y)
                    ymax = max(y)
                    lenx = int((xmax - xmin) / binsize)
                    leny = int((ymax - ymin + 1) / binsize)
                    #px = (xmax - xmin) * np.arange(lenx + 1) / lenx + xmin
                    #py = (ymax - ymin) * np.arange(leny + 1) / leny + ymin
                    #leny = len(py)

                if simon:
                    inside = (z >= zmin) & (z <= zmax) & (y >= ymin) & (y <= ymax) & (f > -20)
                    x = x[inside]
                    y = y[inside]
                    z = z[inside]
                    f = f[inside]
                    v = v[inside]
                    fmean = np.nanmean(f)
                    fstd = np.nanstd(f)

                    if sclip:
                        good = abs(f - fmean) < nsclip * fstd
                        x = x[good]
                        y = y[good]
                        z = z[good]
                        f = f[good]
                        v = v[good]

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


                v[v < 0] = float('NaN')
                if simon:  # cubical coords
                    pz = range(lenz)
                    flux = np.zeros([lenx, leny, lenz])
                    variance = np.zeros([lenx, leny, lenz])
                    npix = np.zeros([lenx, leny, lenz])

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
                            vin = v[zgood]
                            a = np.isfinite(fin)
                            b = np.isfinite(vin)

                            if lzgood > 0:
                                # remember 4 = 2 ** rank -> in this case rank = 2 (2 dims)
                                xn = (xin[a] - xmin) * lenx / float(xmax - xmin)
                                yn = (yin[a] - ymin) * leny / float(ymax - ymin)

                                fk, wk, ck = cic(fin[a], xn, lenx, yn, leny, weighting=resample, wraparound=False)
                                vk, wk, ck = cic(vin[a], xn, lenx, yn, leny, weighting=resample, wraparound=False)

                                if imtype == 'mean':
                                    flux[:, :, k] = fk[:, ::-1] / wk[:, ::-1].astype(float)
                                    variance[:, :, k] = vk[:, ::-1] / wk[:, ::-1].astype(float) / ck[:, ::-1].astype(float)
                                    #if resample:
                                    #    flux *= 2
                                    #    variance *= 2

                                if imtype == 'flux':
                                    flux[:, :, k] = fk[:, ::-1] * np.sqrt(ck[:, ::-1]) / wk[:, ::-1].astype(float)
                                    variance[:, :, k] = vk[:, ::-1] * ck[:, ::-1].astype(float) / wk[:, ::-1].astype(float)
                                    if 0:
                                        flux *= 2
                                        variance *= 4
                                npix[:, :, k] = ck[:, ::-1]

                            x = x[~zgood]
                            y = y[~zgood]
                            z = z[~zgood]
                            f = f[~zgood]
                            v = v[~zgood]
                    else:
                        #a = np.isfinite(f)
                        #b = np.isfinite(v)
                        xn = ((x - xmin) * lenx / float(xmax - xmin)).astype(int)
                        yn = ((y - ymin) * leny / float(ymax - ymin)).astype(int)
                        zn = (z - zmin).astype(int)
                        xn[xn > xmax] = xmax
                        yn[yn > ymax] = ymax
                        zn[zn > 2*zmax] = 2*zmax
                        #aaa = np.where((xn<lenx) & (yn<leny) & (zn<lenz))[0]

                        #print 'aaa len', len(aaa)
                        #xn = xn[aaa]
                        #yn = yn[aaa]
                        #zn = zn[aaa]

                        for m in range(len(xn)):
                            flux[xn[m], yn[m], zn[m]] += f[m]
                            variance[xn[m], yn[m], zn[m]] += v[m]
                            npix[xn[m], yn[m], zn[m]] += 1


                else:  # cylindrical coords
                    fk, ck = cic(f, (x - xmin) / (xmax - xmin), lenx, (y - ymin) / (ymax - ymin), leny)
                    flux[:, k] = fk[:, ::-1]  / ck[:, ::-1].astype(float)
                    vk, ck = cic(v, (x - xmin) / (xmax - xmin), lenx, (y - ymin) / (ymax - ymin), leny)
                    variance[:, k] = vk[:, ::-1] / ck[:, ::-1].astype(float) ** 2
                    npix[:, :] = ck[:, ::-1]

                flux[np.isinf(flux)] = float('NaN')
                flux[npix == 0] = float('NaN')

                variance[np.isinf(variance)] = float('NaN')
                variance[npix == 0] = float('NaN')

                print '\nPair %d %d done' % (i1, i2)

                flux[flux == -999] = float('NaN')
                variance[variance < 0] = float('NaN')

                hdu = PrimaryHDU()
                hdu.header['ID1'] = (int(i1), 'ID of LAE 1 in the catalogue')
                hdu.header['ID2'] = (int(i2), 'ID of LAE 2 in the catalogue')
                hdu.header['z1'] = (int(z1), 'LAE 1 z position in the original cube')
                hdu.header['z2'] = (int(z2), 'LAE 2 z position in the original cube')
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
                hdu.data = flux.T
                hdu.writeto(outname, clobber=overwrite)

                # variance cube
                hdu.data = variance.T
                hdu.header['COMMENT'] = 'Propagated variance in new coordinates'
                if overwrite: hdu.writeto(outname.replace('.fits', '.%s.fits' % vartype), clobber=overwrite)

                # Npixel cube
                if npixim:
                    hdu.data = npix.T
                    hdu.header['COMMENT'] = 'Pixel map'
                    if overwrite: hdu.writeto(outname.replace('.fits', '.NPIX.fits'), clobber=overwrite)

                # snr cube
                if snrim:
                    hdu.data = (flux / np.sqrt(variance)).T
                    hdu.header['COMMENT'] = 'SNR with %s in new coordinates' % vartype
                    if overwrite: hdu.writeto(outname.replace('.fits', '.SNR.fits'), clobber=overwrite)

                    if simon:
                        # flux image and snr image
                        s = 'Cube2Im -cube %s[*,*,%s] -snrmap %s -varcube %s %s -imtype %s -writeNaN .true.' % (
                            outname, zw, outname.replace('.fits', '.SNR.IM.fits'),
                            outname.replace('.fits', '.%s.fits' % vartype), ssmooth, imtype)
                else:
                    s = 'Cube2Im -cube %s %s' % (outname, ssmooth)
                    os.system(s)

                cont = False
                if jpeg:
                    if os.path.isfile(outname.replace('.fits', '.IM.fits')):
                        astroim(outname.replace('.fits', '.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                                pad=.006, contours=True, scb_label='Flux', vmin=-3, vmax=3, x0=lenx/2, y0=leny/2, gray=True)

                        if snrim:
                            astroim(outname.replace('.fits', '.SNR.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                                    pad=.006, contours=True, scb_label='SNR', vmin=-1, vmax=1, x0=lenx/2, y0=leny/2, gray=True)
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
                    fit = variance.T
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
                    astroim(outname.replace('.fits', '.Residuals.IM.fits'), smooth=smooth, saveim=True, show=False,
                            cbfrac=.08, pad=.006, contours=True, scb_label='Flux', title=title, vmin=-.25, vmax=.25, x0=lenx/2, y0=leny/2, gray=True)  # , vmin=-.5, vmax=.5)#,dfig=(8,7))
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
                            cbfrac=.08, pad=.006, contours=True, scb_label='Flux', vmin=-.25, vmax=.25, gray=True)  # , vmin=-.5, vmax=.5)#,dfig=(8,7))
                    # hdu.data = fit*6.25/binsize#zw*1.25/(binsize*0.2) conversion to SB for 1" aperture (originally .2") for a zw z pix wide (zw*1.25 Angstrom)
                    # hdu.writeto(outname.replace('.fits', '.SB.fits'), clobber=True)

                if 0:#jpeg:
                    title = r'LAE%d oriented to LAE%d, z=%1.3f, shear $%1.3f\AA$.' % (id1, id2, red1, (red2-red1)*1.25)
                    astroim(outname.replace('.fits', '.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                            pad=.006, contours=True, scb_label='Flux', title=title, vmin=-3, vmax=3, x0=xw, y0=yw, gray=True)  # , vmin=-.5, vmax=.5)#,dfig=(8,7))
                    # makejpeg(outname.replace('.fits', '.IM.fits'), smooth, outname.replace('.fits', '.jpeg'))

        else:
            print 'No initial catalogue', cat

table.close()
