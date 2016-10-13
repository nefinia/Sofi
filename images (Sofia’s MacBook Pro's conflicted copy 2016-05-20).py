#!/usr/bin/env python
__author__ = 'nefinia'
from pyfits import getdata, PrimaryHDU
import h5py
import numpy as np
import os, sys
from math import sqrt
from tools_sofi import cic, cubex, makejpeg, astroim

fitcat = 'UDF'#'HDFS'
simon = True  # True: coordinates cubical False: cylindrical
norm = False
mask = True
jpeg = False
random = False
resample = True
smooth = True
cubesmooth = False
analysis = True
imtype = 'flux'  # 'mean'#
sb = True
sclip = 20
overwrite = True
theia = False


if cubesmooth:
    smooth = False
    ssmooth = ' -boxsm 3 '
else:
    ssmooth = ''

vartype = 'PROPVAR'
if simon:
    scalelims = '-0.2 0.2'
else:
    scalelims = '-0.1 0.1'

binsize = 2
zw0 = 1
if norm:
    xmin = -.2
    xmax = 1.2
    lenx = 100

ymax = 80
if simon:
    ymin = -ymax
else:
    ymin = 0

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
    sstype=''
    sstype2=''

if simon:
    sstype += 'cubical'
    sstype2 += 'cubical'
else:
    sstype += 'cylindrical'
    sstype2 += 'cylindrical'

if theia:
    folder = '/scratch/gallegos/MUSE/'
else:
    folder = '../../'


table = open('%s/%s/cats/LAE_pairs_info_%s_%s_%s.dat'%(folder, fitcat, ssnorm, ssmask, sstype), 'w')
table.write('#ID1 ID2 x1 x2 z1 z2 comdist mean std SB[erg/s/cm**2/arscec**2] 1sigma_SB combtype masked coord_type\n')


data = getdata('%s/%s/cats/lae_pairs.fits'%(folder, fitcat), 1)
id1 = data['id1'].astype(int)
id2 = data['id2'].astype(int)
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
xps = np.concatenate((xp1, [xp2[-1]]))
yps = np.concatenate((yp1, [yp2[-1]]))
ids = np.concatenate((id1, [id2[-1]]))


if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR.fit'
    filevar = '%s/%s/%s' % (folder, fitcat, cubename)
    filevar = '%s/%s/DATACUBE-HDFS-1.35-%s.fits' % (folder, fitcat, vartype)

if fitcat == 'UDF':
    cubename = 'DATACUBEFINALuser_20141021T055145_212058ad.fits'
    filevar = '%s/%s/%s' % (folder, fitcat, cubename)

datavar = getdata(filevar, 2)

good = np.array(
    [43, 71, 89, 92, 95, 112, 139, 144, 159, 181, 183, 186, 200, 202, 216, 218, 225, 232, 246, 271, 294, 308,
     311, 325, 334, 338, 422, 430, 437, 441, 449, 484, 489, 492, 499, 503, 514, 520, 543, 546, 547, 549, 551,
     553, 554, 555, 557, 560, 563, 568])

if 0:
    good2 = np.array(
        ["43-449", "43-513", "89-216", "89-308", "92-554", "139-181", "144-216", "144-294", "144-308", "186-334",
         "186-484", "216-294", "218-334", "218-484", "294-308", "334-338", "334-484", "499-549", "499-568"])
    close = np.where((["%d-%d" % (id1[i], id2[i]) in good2 for i in range(len(id1))]) & (dist <= 20) & (dist > 1))[0]
else:
    # close = np.where(([id1[i] in good for i in range(len(id1))]) & (dist <= 20) & (dist > .5))[0]
    #close = np.where((id1 == 144) & (id2 == 216))[0]
    close = np.where((dist <= 20) & (dist > .5) & (theta > 6))[0]
    #close = np.where((id1 >= 499))[0]
if mask:
    # Open mask file
    if fitcat == 'HDFS':
        fits = '%s/%s/DATACUBE-HDFS-1.35-STATVAR.IM.Objects_Id.fits' %(folder,fitcat)
    if fitcat == 'UDF':
        fits = '%s/%s/DATACUBEFINALuser_20141021T055145_212058ad.IM.Objects_Id.fits' %(folder,fitcat)
    data_cube, header_data_cube = getdata(fits, 0, header=True)
    data_cube = data_cube[0]
    # unmasked = np.where((data_cube == 0))

    if 1:
        # remove LAE from the mask catalogue
        laeid = {}
        for x, y, ii in zip(xps, yps, ids):
            if data_cube[y - 1, x - 1] > 0:
                laeid['%d'%ii] = data_cube[y - 1, x - 1]
            else:
                laeid['%d'%ii] = -1
                # print 'LAE found in mask! id', ii
        # laeids = np.where((data_cube == 0) | ([data_cube == i for i in laeid]))
        #unmasked = np.where([data_cube == i for i in laeid])
        #data_cube[unmasked[1], unmasked[2]] = 0

    # Open sky file
    sky = np.loadtxt('%s/%s/skyspec.txt' %(folder,fitcat)).T
    # z = sky[0].astype('int')
    # wav = sky[1].astype('float')
    skyflux = sky[2].astype('float')
    skymean = np.mean(skyflux)
    skystd = np.std(skyflux)
    zmask = np.where(abs(skyflux - skymean) < skystd)[0] + 1

for i1, i2, d, z1, z2 in zip(id1[close], id2[close], dist[close], zp1[close], zp2[close]):

    cat = '%s/%s/pairs/%s/cats/%s_pair%d-%d.h5' % (folder, fitcat, ssnorm, sstype, i1, i2)
    outname = '%s/%s/pairs/%s/%s/%s_pair%d-%d.fits' % (folder, fitcat, ssnorm, ssmask, sstype, i1, i2)

    if os.path.isfile(cat):
        if not os.path.isfile(outname) or overwrite:

            print 'Catalog name', cat
            print 'Output name', outname
            data = h5py.File(cat, 'r')

            if mask:
                print 'Defining good pixels'
                xpix = np.array(data['x'])
                ypix = np.array(data['y'])
                a = np.where((data_cube[ypix - 1, xpix - 1] == 0) | (data_cube[ypix - 1, xpix - 1] == laeid['%d'%i1]) | (data_cube[ypix - 1, xpix - 1] == laeid['%d'%i2]))[0]
                if len(a)==0:
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

            if not norm:
                xmin = min(x) * .99  # int(-.2*dpix)
                xmax = max(x) * 1.01  # int(1.2*dpix)
                # ymax = max(y)
                lenx = int((xmax - xmin) / binsize)

            leny = int((ymax - ymin + 1) / binsize)
            zmin = -20  # min(z)
            zmax = 20  # max(z)
            lenz = int(zmax - zmin) + 1
            zw = '%d:%d' % (lenz / 2 - zw0, lenz / 2 + zw0)

            px = (xmax - xmin) * np.arange(lenx + 1) / lenx + xmin
            py = (ymax - ymin) * np.arange(leny + 1) / leny + ymin

            # area = abs((px[1]-px[0])*(py[1]-py[0]))

            leny = len(py)
            ltot = 0

            if simon:
                inside = (z >= zmin) & (z <= zmax) & (y >= ymin) & (y <= ymax) & (f > -20)
                x = x[inside]
                y = y[inside]
                z = z[inside]
                f = f[inside]
                v = v[inside]
                fmean = np.nanmean(f)
                fstd = np.nanstd(f)
                good = abs(f - fmean) < sclip * fstd
                x = x[good]
                y = y[good]
                z = z[good]
                f = f[good]
                v = v[good]
            # check position of LAEs
            for i in range(lenx):
                if (0 >= px[i]) and (0 <= px[i + 1]):
                    posin = i
                if norm:
                    if (1 >= px[i]) and (1 <= px[i + 1]):
                        posout = i
                else:
                    if (dpix >= px[i]) and (dpix <= px[i + 1]):
                        posout = i

            v[v < 0] = float('NaN')
            if simon:#cubical coords
                pz = range(lenz)
                flux = np.zeros([lenx, leny, lenz]) - 999
                variance = np.zeros([lenx, leny, lenz]) - 999
                npix = np.zeros([lenx, leny, lenz])
                for k in pz:
                    print '%d-%d:' % (i1, i2), 'z', k + 1, 'of', lenz, 'mask', mask, 'norm', norm, 'resample %s'%resample
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
                        fk, ck = cic(fin[a], (xin[a] - xmin) * lenx / (xmax - xmin), lenx,
                                     (yin[a] - ymin) * leny / (ymax - ymin), leny, weighting=resample)
                        vk, ck = cic(vin[a], (xin[a] - xmin) * lenx / (xmax - xmin), lenx,
                                     (yin[a] - ymin) * leny / (ymax - ymin), leny, weighting=resample)
                        if imtype == 'mean':
                            flux[:, :, k] = fk[:, ::-1] * 2 / ck[:, ::-1].astype(float)
                            variance[:, :, k] = vk[:, ::-1] * 2 / ck[:, ::-1].astype(float) ** 2

                        if imtype == 'flux':
                            flux[:, :, k] = fk[:, ::-1] * 2
                            variance[:, :, k] = vk[:, ::-1] * 4 / ck[:, ::-1].astype(float)
                        npix[:, :, k] = ck[:, ::-1]

                    x = x[~zgood]
                    y = y[~zgood]
                    z = z[~zgood]
                    f = f[~zgood]
                    v = v[~zgood]

            else:#cylindrical coords
                fk, ck = cic(f, (x - xmin) / (xmax - xmin), lenx, (y - ymin) / (ymax - ymin), leny)
                flux[:, k] = fk[:, ::-1] * binsize ** 2 / ck[:, ::-1].astype(float)
                vk, ck = cic(v, (x - xmin) / (xmax - xmin), lenx, (y - ymin) / (ymax - ymin), leny)
                variance[:, k] = vk[:, ::-1] * binsize ** 2 / ck[:, ::-1].astype(float) ** 2
                npix[:, :] = ck[:, ::-1]

            flux[np.isinf(flux)] = float('NaN')
            flux[ck[:, ::-1] == 0] = float('NaN')

            variance[np.isinf(variance)] = float('NaN')
            variance[ck[:, ::-1] == 0] = float('NaN')

            if 0:
                if not norm: lenx2 = lenx

                if simon:
                    pz = range(lenz)
                    flux = np.zeros([lenx2, leny, lenz]) - 999
                    npix = np.zeros([lenx2, leny, lenz])
                    variance = np.zeros([lenx2, leny, lenz]) - 999

                else:
                    flux = np.zeros([lenx2, 2 * leny - 1]) - 999
                    npix = np.zeros([lenx2, 2 * leny - 1])
                    variance = np.zeros([lenx2, 2 * leny - 1]) - 999

                for i in range(lenx):
                    print 'Pair', i1, i2, 'x', i + 1, 'of', lenx, 'mask', mask, 'norm', norm, 'coords', sstype, 'zlen', lenz
                    for j in range(leny - 1):
                        inside = (x >= px[i]) & (x <= px[i + 1]) & (y >= py[j]) & (y <= py[j + 1])
                        lin = inside.sum()
                        ltot += lin

                        if simon:
                            if lin > 0:
                                zin = z[inside]
                                fin = f[inside]
                                vin = v[inside]

                                for k in pz:
                                    zgood = zin == (k + zmin)
                                    lzgood = zgood.sum()
                                    if lzgood > 0:
                                        if imtype == 'flux' or imtype == 'SB':
                                            fm = np.nansum(fin[zgood])
                                        else:
                                            fm = np.nanmean(fin[zgood])
                                        flux[i, leny - 1 - j, k] = fm
                                        npix[i, leny - 1 - j, k] = lzgood
                                        if imtype == 'flux' or imtype == 'SB':
                                            varmean = np.nanmean(vin[zgood])
                                        else:
                                            varmean = np.nanmean(vin[zgood]) / lzgood
                                        variance[i, leny - 1 - j, k] = varmean

                        else:
                            if lin > 0:
                                fm = np.mean(f[inside])
                                flux[i, leny - j - 1] = fm
                                flux[i, leny + j - 1] = fm
                                npix[i, leny - j - 1] = lin
                                npix[i, leny + j - 1] = lin
                                varmean = np.nanmean(v[inside]) / lin
                                variance[i, leny - j - 1] = varmean
                                variance[i, leny + j - 1] = varmean

                        x = x[~inside]
                        y = y[~inside]
                        if simon:
                            z = z[~inside]

                        f = f[~inside]
                        v = v[~inside]

                    print 'number inside', ltot, 'outside', len(x)
                    ltot = 0

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
            if imtype == 'mean':
                hdu.header['SB'] = (
                '10**(-20)*erg/s/cm**2/Angstrom/arcsec**2', 'Surface brightness for array values????')
            if imtype == 'flux':
                conv = 1.25e-20/(binsize*.2)#**2
                hdu.header['SB'] = (
                '%.2e*erg/s/cm**2/arcsec**2' % conv, 'Surface brightness for array values')

                xin = min(posin + 7, posout - 9)
                xout = max(posin + 7, posout - 9)
                if xin==xout:
                    xout+=1
                coolmean = np.nanmean(flux[xin: xout, leny/2 - 10:leny/2 + 10, lenz/2])
                coolstd = np.nanstd(flux[xin: xout, leny/2 - 10:leny/2 + 10, lenz/2])
                hdu.header['mean'] = ('%f'%coolmean, 'Flux mean in relevant region')
                hdu.header['std'] = ('%f'%coolstd, 'Flux std in relevant region')
                hdu.header['1sigma'] = ('%.2e erg/s/cm**2/arcsec**2' % (coolstd*conv), '1 sigma SB in relevant region')
            if not norm:
                hdu.header['x1'] = (posin, 'Horizontal position of LAE 1 in pixels')
                hdu.header['x2'] = (posout, 'Horizontal position of LAE 2 in pixels')
            hdu.header['pix2asec'] = (.2 * binsize, 'Conversion between pixels and arcseconds')

            if simon:
                hdu.header['TYPE'] = ('Projected+Shear', 'Type of coordinate transformation')
            else:
                hdu.header['TYPE'] = ('Cylindrical', 'Type of coordinate transformation')

            table.write('%d %d %d %d %d %d %f %f %f %.3e %.3e %s %s Projected+Shear\n'%(i1, i2, posin, posout, z1, z2, d, coolmean, coolstd, conv, coolstd*conv, imtype, mask))


            if overwrite:
                print "Warning: Overwriting output fits file"

            #flux fits file
            hdu.data = flux.T
            hdu.writeto(outname, clobber=overwrite)

            #variance cube
            hdu.data = variance.T
            hdu.header['COMMENT'] = 'Propagated variance in new coordinates'
            hdu.writeto(outname.replace('.fits', '.%s.fits' % vartype), clobber=overwrite)

            #Npixel cube
            hdu.data = npix.T
            hdu.header['COMMENT'] = 'Pixel map'
            hdu.writeto(outname.replace('.fits', '.NPIX.fits'), clobber=overwrite)

            #snr cube
            hdu.data = (flux / np.sqrt(variance)).T
            hdu.header['COMMENT'] = 'SNR with %s in new coordinates' % vartype
            hdu.writeto(outname.replace('.fits', '.SNR.fits'), clobber=overwrite)

            if simon:
                #flux image and snr image
                s = 'Cube2Im -cube %s[*,*,%s] -snrmap %s -varcube %s %s -imtype %s -writeNaN .true.' % (
                    outname, zw, outname.replace('.fits', '.SNR.IM.fits'),
                    outname.replace('.fits', '.%s.fits' % vartype), ssmooth, imtype)
            else:
                s = 'Cube2Im -cube %s %s' % (outname.replace('.fits', '.SNR.fits'), ssmooth)
            os.system(s)

            cont = False
            if jpeg:
                makejpeg(outname.replace('.fits', '.IM.fits'), smooth, outname.replace('.fits', '.jpeg'))
                makejpeg(outname.replace('.fits', '.SNR.IM.fits'), smooth, outname.replace('.fits', '.SNR.jpeg'), cont,  scalelims = '')

            if sb:
                hdu = PrimaryHDU()
                fit = getdata(outname.replace('.fits', '.IM.fits'))
                hdu.data = fit*(zw0*2+1)*3.125#*3*1.25/.4 conversion to SB for 1" aperture for a 3 z pix wide (3*1.25 Angstrom)
                hdu.writeto(outname.replace('.fits', '.SB.fits'), clobber=True)

        else:
            print 'Computation already done.'

            if jpeg:
                astroim(outname.replace('.fits', '.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08, pad=.006,dfig=(8,7))

                #makejpeg(outname.replace('.fits', '.IM.fits'), smooth, outname.replace('.fits', '.jpeg'))

            if 0:#sb:
                hdu = PrimaryHDU()
                fit = getdata(outname.replace('.fits', '.IM.fits'))
                hdu.data = fit*6.25/binsize#zw*1.25/(binsize*0.2) conversion to SB for 1" aperture (originally .2") for a zw z pix wide (zw*1.25 Angstrom)
                hdu.writeto(outname.replace('.fits', '.SB.fits'), clobber=True)

    else:
        print 'No initial catalogue', cat


table.close()