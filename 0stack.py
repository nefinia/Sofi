__author__ = 'gallegos'

from pyfits import getdata, PrimaryHDU, HDUList, getheader
import pyfits as pf
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import LogNorm
from tools_sofi import astroim, rdarg
from sys import argv
from scipy import stats
import random as rd

half = rdarg(argv, 'half', bool, True)
fitcat = rdarg(argv, key='fitcat', type=str, default='all')  # 'all'#'HDFS'
scalelims = rdarg(argv, 'scalelims', str, '-0.01 0.01')
nodist = rdarg(argv, 'nodist', bool, True)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)
if cubesmooth:
    smooth = False
    ssmooth = ' -boxsm 3 '
else:
    smooth = True
    ssmooth = ''
simon = rdarg(argv, 'cubical', bool, True)
random = rdarg(argv, 'random', bool, False)
nrand = rdarg(argv, 'nrand', int, 1)
tests = rdarg(argv, 'tests', bool, False)
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
flipx = rdarg(argv, 'flipx', bool, False)
flipy = rdarg(argv, 'flipy', bool, False)
rotate = rdarg(argv, 'rotate', bool, False)
weights = rdarg(argv, 'weights', bool, False)
mask = rdarg(argv, 'mask', bool, True)
jpeg = rdarg(argv, 'image', bool, True)
norm = rdarg(argv, 'norm', bool, False)
propvar = rdarg(argv, 'pvar', bool, False)
statvar = rdarg(argv, 'svar', bool, False)
imtype = rdarg(argv, 'imtype', str, 'mean')  # 'mean'
sclip = rdarg(argv, 'sclip', int, 5)
extraname = rdarg(argv, 'extraname', str, '')
if not extraname:
    extraname = ''
sb = rdarg(argv, 'sb', bool, False)
binsize = rdarg(argv, 'binsize', int, 2)
vmin = rdarg(argv, 'vmin', float, -1.2)
vmax = rdarg(argv, 'vmax', float, 1.2)
zw0 = rdarg(argv, 'zw0', int, 1)
distmin = rdarg(argv, 'dmin', list, [.5])
distmax = rdarg(argv, 'dmax', list, [20])
pdmin = rdarg(argv, 'pdmin', list, [20])
pdmax = rdarg(argv, 'pdmax', list, [80])
redmin = rdarg(argv, 'zmin', list, [2.9])
redmax = rdarg(argv, 'zmax', list, [4])
fluxmin = rdarg(argv, 'fmin', list, [100])
cubexmask = rdarg(argv, 'cubexmask', bool, False)
if cubexmask:
    extraname = '.Residuals'

if mask:
    sm = 'masked'
else:
    sm = 'unmasked'

if weights:
    sw = 'weighted'
else:
    sw = 'unweighted'

if norm:
    sn = 'normalized'
else:
    sn = 'not-normalized'

if random:
    st = 'random_'
else:
    st = ''

if simon:
    st += 'cubical'
else:
    st += 'cylindrical'

if norm:
    lenx = 40

if fitcat == 'all':
    # table = open('%s/%s/cats/LAE_pairs_info_%s_%s_%s.fits'%(folder, 'HDFS', sn, sm, st), 'r')

    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'HDFS')
    data = getdata(cat2, 1)
    ids1 = data['id1']
    ids2 = data['id2']

    laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'HDFS'), 1)
    flae = laedata['LYALPHA_FLUX']
    idlae = laedata['ID']
    flux_lae1 = []
    flux_lae2 = []
    for i, j in zip(ids1.astype(int), ids2.astype(int)):
        flux_lae1.append(np.float(flae[idlae == i]))
        flux_lae2.append(np.float(flae[idlae == j]))
    #sconf = data['sconf'].astype(int)
    sconf1 = data['sconf1'].astype(int)
    sconf2 = data['sconf2'].astype(int)
    pdists = data['x2'] - data['x1']
    zs1 = data['z1']
    redshift = (data['redshift1'] + data['redshift2']) / 2.
    dists = data['pi_Mpc']
    theta = data['theta']
    fitcats = ['HDFS'] * len(data)

    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'UDF')
    del data
    #HDUList.close()
    data = getdata(cat2, 1)
    ids1 = np.concatenate((ids1, data['id1']), 0)
    ids2 = np.concatenate((ids2, data['id2']), 0)

    laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'UDF'), 1)
    flae = laedata['LYALPHA_FLUX']
    idlae = laedata['ID']
    for i, j in zip(data['id1'].astype(int), data['id2'].astype(int)):
        flux_lae1.append(np.float(flae[idlae == i]))
        flux_lae2.append(np.float(flae[idlae == j]))
    #sconf = np.concatenate((sconf, data['sconf'].astype(int)), 0)
    sconf1 = np.concatenate((sconf1, data['sconf1'].astype(int)), 0)
    sconf2 = np.concatenate((sconf2, data['sconf2'].astype(int)), 0)
    pdists = np.concatenate((pdists, data['x2'] - data['x1']), 0)
    zs1 = np.concatenate((zs1, data['z1']), 0)
    redshift = np.concatenate((redshift, (data['redshift1'] + data['redshift2']) / 2.), 0)
    dists = np.concatenate((dists, data['pi_Mpc']), 0)
    theta = np.concatenate((theta, data['theta']), 0)
    nhdfs = len(fitcats)
    fitcats = np.concatenate((fitcats, ['UDF'] * len(data)), 0)
    flux_lae2 = np.array(flux_lae2)
    flux_lae1 = np.array(flux_lae1)
    #HDUList.close()
    del data

else:
    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, fitcat)
    data = getdata(cat2, 1)
    ids1 = data['id1']
    ids2 = data['id2']

    pdists = data['x2'] - data['x1']
    zs1 = data['z1']
    redshift = (data['redshift1'] + data['redshift2']) / 2.
    dists = data['pi_Mpc']
    theta = data['theta']
    fitcats = [fitcat] * len(data)
    del data


# UDF!!
ns = []

if simon:
    st = 'cubical'
else:
    st = 'cylindrical'
if random:
    st = 'random_' + st


def makejpeg(im, scalelims, smooth, imout=''):
    if imout == '':
        if im.find('.fits') == -1:
            print 'Standard .fits extension not found. Is it a fits file?'
            imout = im + '.jpeg'
        else:
            imout = im.replace('.fits', '.jpeg')
    # -scale limits %s
    scalelims = ''
    s = 'ds9 -zscale %s -cmap b -smooth %s %s -zoom to fit -saveimage jpeg %s 100 -exit' % (
        scalelims, smooth, im, imout)
    print 'Creating jpeg image: ', imout
    os.system(s)


h = 'theta '

hdu = PrimaryHDU()

for dmin, dmax, rmin, rmax, pmin, pmax, fmin in zip(distmin, distmax, redmin, redmax, pdmin, pdmax, fluxmin):
    # close = np.where((["%d-%d" % (ids1[i], ids2[i]) in good3 for i in range(len(ids1))]) & (redshift <= rmax) & (
    # redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax) & (theta > pmin))[0]
    # & (["%d-%d"%(ids1[i],ids2[i]) in good2 for i in range(len(ids1))]))[0]#
    close = np.where((redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax)
                     & (theta > pmin) & (flux_lae1 >= fmin) & (flux_lae2 >= fmin) & (sconf1 > 0) & (sconf2 > 0)
                     & (ids1 != 144) & (ids2 != 144))[0]  # 144 agn discarded & (sconf > 0)

    if len(close) > 1:
        ns.append([dmin, dmax, rmin, rmax, len(close)])
        print 'Number of pairs to be stacked', len(close)
        id1 = ids1[close]
        id2 = ids2[close]
        z1 = zs1[close]
        dist = dists[close]
        fcats = np.array(fitcats)[close]

        fits = []
        x1 = []
        x2 = []
        means = []
        stds = []
        randf = []
        if propvar:
            pvars = []
        pairs_lst = []

        for i in range(len(id1)):
            for j in range(nrand):
                #HDUList.close()
                #del data

                if j > 0 and random:
                    pair = '%s/%s/pairs/%s/%s/%s_pair%d-%d.%d.fits' % (folder, fcats[i], sn, sm, st, id1[i], id2[i], j)
                else:
                    pair = '%s/%s/pairs/%s/%s/%s_pair%d-%d.fits' % (folder, fcats[i], sn, sm, st, id1[i], id2[i])



                if os.path.isfile(pair):
                    pairs_lst.append(pair)

                    if cubexmask:
                        hdul = pf.open(pair.replace('.fits', '.Residuals.fits'))
                        fit = (hdul[0]).data#getdata(pair.replace('.fits', '.Residuals.fits'))
                        hdul.close()
                        fit[fit == -999] = np.nan
                        fits.append(fit)
                    else:
                        fit = getdata(pair)
                        fits.append(fit)

                    zw = "%d:%d" % (fit.shape[0] / 2 - zw0 + 1, fit.shape[0] / 2 + zw0 + 1)
                    fhead = getheader((pair))


                    if propvar:
                        hdul = pf.open(pair.replace('.fits', '.PROPVAR.fits'))
                        var = (hdul[0]).data#getdata(pair.replace('.fits', '.PROPVAR.fits'))
                        pvars.append(var)
                        hdul.close()

                    if flipx:
                        pairflip = pair.replace('.fits', '.flip.fits')
                        if cubexmask:
                            pairflip = pairflip.replace('.fits', '.Residuals.fits')
                        if not os.path.isfile(pairflip):
                            if simon:
                                hdu.data = fit[:, :, ::-1]
                                s = 'Cube2Im %s[*,*,%s] %s' % (pairflip, zw, ssmooth)
                                os.system(s)
                            else:
                                hdu.data = fit[:, ::-1]
                            hdu.writeto(pairflip, clobber=True)

                        if simon:
                            fits.append(fit[:, :, ::-1])  # - float(fhead['mean']))

                            if propvar: pvars.append(var[:, :, ::-1])
                        else:
                            fits.append(fit[:, ::-1])  # - float(fhead['mean']))

                            if propvar: pvars.append(var[:, ::-1])

                    if flipy and simon:
                        pairflip = pair.replace('.fits', '.flip.fits')
                        if not os.path.isfile(pairflip):
                            hdu.data = fit[:, ::-1, :]
                            hdu.writeto(pairflip, clobber=True)
                            s = 'Cube2Im %s[*,*,%s] %s' % (pairflip, zw, ssmooth)
                            os.system(s)

                        fits.append(fit[:, ::-1, :])  # - float(fhead['mean']))

                        if propvar: pvars.append(var[:, ::-1, :])

                    if rotate:
                        pairrot = pair.replace('.fits', '.rot.fits')
                        if cubexmask:
                            pairrot = pairrot.replace('.fits', '.Residuals.fits')
                        if not os.path.isfile(pairrot):
                            if simon:
                                hdu.data = fit[:, ::-1, ::-1]
                            else:
                                hdu.data = fit[::-1, ::-1]
                            hdu.writeto(pairrot, clobber=True)
                            s = 'Cube2Im %s[*,*,%s] %s' % (pairrot, zw, ssmooth)
                            os.system(s)

                        if simon:
                            fits.append(fit[:, ::-1, ::-1])
                            if propvar: pvars.append(var[:, ::-1, ::-1])

                        else:
                            fits.append(fit[::-1, ::-1])
                            if propvar: pvars.append(var[::-1, ::-1])

                    if nodist:
                        pair = '%s/%s/pairs/%s/%s/%s_pair%d-%d.fits' % (folder, fcats[i], sn, sm, st, id2[i], id1[i])

                        if os.path.isfile(pair):
                            pairs_lst.append(pair)
                            if cubexmask:
                                fit = getdata(pair.replace('.fits', '.Residuals.fits'))
                                fit[fit == -999] = np.nan
                                fits.append(fit)
                            else:
                                fit = getdata(pair)
                                fits.append(fit)
                    del fit
                    #hdu.close()
        if not propvar:
            pvars = fits

        if not norm and not nodist:
            if simon:
                fzlen, fylen, fxlen = np.array([f.shape for f in fits]).T
                xmin = -20#min(fxlen)
                xmax = max(fxlen)
                ymin = min(fylen)
                ymax = max(fylen)
                zmax = max(fzlen)
                lenz = zmax
                shift = 20
                shape = (zmax, ymax, xmax + shift)
                fits2 = []
                pvars2 = []
                for f, p, xl, yl, zl, i1, std in zip(fits, pvars, fxlen, fylen, fzlen, x1, stds):
                    fnew = np.zeros(shape, dtype=f.dtype) - 999
                    dxmin = max(0, shift - i1)
                    dxmin2 = max(0, i1 - shift)
                    dxmax = xl
                    fnew = np.zeros(shape, dtype=f.dtype) - 999

                    if weights:
                        fnew[:, :, dxmin:dxmax + dxmin] = f / std ** 2
                    else:
                        fnew[:, :, dxmin:dxmax + dxmin] = f
                    fits2.append(fnew)
                    if propvar:
                        pnew = np.zeros(shape, dtype=p.dtype) - 999
                        if weights:
                            pnew[:, :, dxmin:dxmax + dxmin] = p / std ** 2
                        else:
                            pnew[:, :, dxmin:dxmax + dxmin] = p
                        pvars2.append(pnew)

            else:
                fits2 = []
                pvars2 = []
                fylen, fxlen = np.array([f.shape for f in fits]).T
                xmin = min(fxlen)
                xmax = max(fxlen)
                ymin = min(fylen)
                ymax = max(fylen)
                shape = (ymax, xmax / 2 + 15)
                fits2 = []
                for f, p, xl, yl, i1, i2, std in zip(fits, pvars, fxlen, fylen, x1, x2, stds):
                    dxmin = max(0, shift - i1)
                    dxmin2 = max(0, i1 - shift)
                    dxmax = (i1 + i2) / 2
                    fnew = np.zeros(shape, dtype=f.dtype) - 999
                    if weights:
                        fnew[:yl, dxmin:dxmax + dxmin] = f[:yl, dxmin2:dxmax + dxmin2] / std ** 2
                    else:
                        fnew[:yl, dxmin:dxmax + dxmin] = f[:yl, dxmin2:dxmax + dxmin2]

                    fits2.append(fnew)
                    if propvar:
                        pnew = np.zeros(shape, dtype=p.dtype) - 999
                        if weights:
                            pnew[:yl, dxmin:dxmax + dxmin] = p[:yl, dxmin2:dxmax + dxmin2] / std ** 2
                        else:
                            pnew[:yl, dxmin:dxmax + dxmin] = p[:yl, dxmin2:dxmax + dxmin2]
                        pvars2.append(pnew)

        else:
            fits2 = np.copy(fits)
            if propvar: pvars2 = np.copy(pvars)

        fits2 = np.array(fits2)
        fits2[fits2 == -999] = np.nan
        nf2 = fits2.shape[0]
        if tests:
            sample1 = np.array(rd.sample(np.arange(nf2), nf2/2))
            sample2 = np.arange(fits2.shape[0])[~sample1]

        if weights:
            sweight = np.array(stds) ** -2
            a = np.array([np.isfinite(f) * sss for f, sss in zip(fits2, sweight)])
            wsum = np.nansum(a, axis=0)
            if tests:
                wsum_t1 = np.nansum(a[sample1, :], axis=0)
                wsum_t2 = np.nansum(a[sample2, :], axis=0)

        if propvar:
            pvars2 = np.array(pvars2)
            pvars2[pvars2 == -999] = float('NaN')
            if imtype == 'mean':
                if weights:
                    pvar = np.nanmean(pvars2, axis=0) / wsum
                else:
                    pvsum = np.nansum(np.isfinite(pvars2), axis=0)
                    pvar = np.nanmean(pvars2, axis=0) / pvsum
            else:
                pvsum = np.nansum(np.isfinite(pvars2), axis=0)
                pvar = np.nansum(pvars2, axis=0) / pvsum
        sigma = np.nanstd(fits2, 0)
        high_sigma = np.where(abs(fits2) > sclip * sigma)
        fits2[high_sigma] = np.nan

        if imtype == 'mean':
            if weights:
                f = np.nansum(fits2, axis=0) / wsum

                if tests:
                    f_t1 = np.nansum(fits2[sample1, :, :, :], axis=0) / wsum_t1  # np.nansum(a, axis=0)
                    f_t2 = np.nansum(fits2[sample2, :, :, :], axis=0) / wsum_t1  # np.nansum(a, axis=0)
            else:
                f = np.nanmean(fits2, axis=0)
                if tests:
                    f_t1 = np.nanmean(fits2[sample1, :, :, :], axis=0)
                    f_t2 = np.nanmean(fits2[sample2, :, :, :], axis=0)

        elif imtype == 'median':
            f = stats.nanmedian(fits2, axis=0)
            if tests:
                f_t1 = stats.nanmedian(fits2[sample1, :, :, :], axis=0)
                f_t2 = stats.nanmedian(fits2[sample2, :, :, :], axis=0)

        else:
            f = np.nansum(fits2, axis=0)
            if tests:
                f_t1 = np.nansum(fits2[sample1, :, :, :], axis=0)
                f_t2 = np.nansum(fits2[sample2, :, :, :], axis=0)

        nf = np.nansum(np.isfinite(fits2), axis=0)
        nstack = len(fits2)

        stackname = '%s/%s/stacks/%s/%s/%s_stack_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d%s.fits' % (
        folder, fitcat, sn, sm, st,
        dmin, dmax, pmin, pmax, rmin, rmax, fmin, extraname)

        hdu.header['N_images'] = (nstack, 'Number of LAEs stacked')
        hdu.header['zmin'] = (rmin, 'Minimum redshift')
        hdu.header['zmax'] = (rmax, 'Maximum redshift')
        hdu.header['Dmin'] = (dmin, 'Minimum pair distance in cMpc')
        hdu.header['Dmax'] = (dmax, 'Maximum pair distance in cMpc')
        hdu.header['PDmin'] = (pmin, 'Minimum projected distance in arcsec')
        hdu.header['PDmax'] = (pmax, 'Maximum projected distance in arcsec')

        if not nodist:
            if norm:
                # 14 and 85 are the position of the LAEs when using xlen=100 in normalized coords
                xin = 14 + 8
                xout = 85 / 2
                zl, yl, xl = f.shape
            else:
                xin = shift + 8
            xout = xmax / 2
        else:
            zl, yl, xl = f.shape
        stdbin = 4
        stdbins = stdbin * np.arange(xl / stdbin)
        stack_stds = []
        stack_flux = []
        stack_npix = []
        if tests:
            stack_stds_t1 = []
            stack_stds_t2 = []

        #Analysissssss
        for i in stdbins:
            stack_flux.append(np.nanmean(f[zl/2 - zw0: zl/2 + zw0, yl/2 - 3: yl/2 + 3, i:i + stdbin]))
            stack_npix.append(np.sum(nf[zl/2 - zw0: zl/2 + zw0, yl/2 - 3: yl/2 + 3, i:i + stdbin]))
            stack_stds.append(np.nanstd(np.concatenate((f[zl/2 - zw0: zl/2 + zw0, yl / 2 + 5:, i:i + stdbin], f[zl/2 - zw0: zl/2 + zw0, :yl / 2 - 4, i:i + stdbin]))))

            if tests:
                stack_stds_t1.append(np.nanstd(np.concatenate((f_t1[zl/2 - zw0: zl/2 + zw0, yl / 2 + 5:, i:i + stdbin], f_t1[zl/2 - zw0: zl/2 + zw0, :yl / 2 - 4, i:i + stdbin]))))
                stack_stds_t2.append(np.nanstd(np.concatenate((f_t2[zl/2 - zw0: zl/2 + zw0, yl / 2 + 5:, i:i + stdbin], f_t2[zl/2 - zw0: zl/2 + zw0, :yl / 2 - 4, i:i + stdbin]))))
        hdu.data = f
        hdu.writeto(stackname, clobber=True)
        s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (
            stackname, zw, stackname.replace('.fits', '.IM.fits'))
        print s
        os.system(s)
        title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $flux>%d$' % (
        nstack, dmin, dmax, pmin, pmax, rmin, rmax, fmin)

        f_snr = np.zeros(f.shape)
        for i, j in zip(stdbins, stack_stds):
            f_snr[:, :, i:i + stdbin] = f[:, :, i:i + stdbin] / j
        hdu.data = f_snr#np.nanmean(f_snr[zl / 2 - zw0:zl / 2 + zw0, :, :], 0)
        hdu.writeto(stackname.replace('.fits', '.SNR.fits'), clobber=True)
        s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (
            stackname.replace('.fits', '.SNR.fits'), zw, stackname.replace('.fits', '.SNR.IM.fits'))
        print s
        os.system(s)

        if tests:
            hdu.data = nf
            hdu.writeto(stackname.replace('.fits', '.NPIX.fits'), clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (
                stackname.replace('.fits', '.NPIX.fits'), zw, stackname.replace('.fits', '.NPIX.IM.fits'))
            print s
            os.system(s)

            f_snr_t1 = np.zeros(f.shape)
            f_snr_t2 = np.zeros(f.shape)
            for i, j, k in zip(stdbins, stack_stds_t1, stack_stds_t2):
                f_snr_t1[:, :, i:i + stdbin] = f_t1[:, :, i:i + stdbin] / j
                f_snr_t2[:, :, i:i + stdbin] = f_t2[:, :, i:i + stdbin] / k
            hdu.data = f_snr_t1#np.nanmean(f_snr_t1[zl / 2 - zw0:zl / 2 + zw0, :], 0)
            hdu.writeto(stackname.replace('.fits', '.SNR.halftest1.fits'), clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (
                stackname.replace('.fits', '.SNR.halftest1.fits'), zw,
                stackname.replace('.fits', '.SNR.halftest1.IM.fits'))
            print s
            os.system(s)
            hdu.data = f_snr_t2#np.nanmean(f_snr_t2[zl / 2 - zw0:zl / 2 + zw0, :], 0)
            hdu.writeto(stackname.replace('.fits', '.SNR.halftest2.fits'), clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (
                stackname.replace('.fits', '.SNR.halftest2.fits'), zw,
                stackname.replace('.fits', '.SNR.halftest2.IM.fits'))
            print s
            os.system(s)

        if propvar:
            hdu.data = pvar
            hdu.writeto(stackname.replace('.fits', '.PROPVAR.fits'), clobber=True)
            hdu.data = f / np.sqrt(pvar)
            hdu.writeto(stackname.replace('.fits', '.PROPVAR.SNR.fits'), clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -snrmap %s -varcube %s %s' % (
                stackname, zw, stackname.replace('.fits', '.PROPVAR.SNR.IM.fits'),
                stackname.replace('.fits', '.PROPVAR.fits'), ssmooth)
            print s
            os.system(s)
            if sb:
                fit = getdata(stackname.replace('.fits', '.IM.fits'))
                hdu.data = fit * 3.125
                hdu.writeto(stackname.replace('.fits', '.SB.IM.fits'), clobber=True)
                astroim(stackname.replace('.fits', '.SB.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                        pad=.006, dfig=(8, 7), title=title)
            if jpeg: astroim(stackname.replace('.fits', '.PROPVAR.SNR.IM.fits'), smooth=smooth, saveim=True, show=False,
                             cbfrac=.08,
                             pad=.006, dfig=(8, 7), contours=True, title=title)

        if statvar:
            var = np.nanvar(fits2, 0)
            hdu.data = var
            hdu.writeto(stackname.replace('.fits', '.STATVAR.fits'), clobber=True)
            hdu.data = f / np.sqrt(var)
            hdu.writeto(stackname.replace('.fits', '.STATVAR.SNR.fits'), clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -snrmap %s -varcube %s %s' % (
                stackname, zw, stackname.replace('.fits', '.STATVAR.SNR.IM.fits'),
                stackname.replace('.fits', '.STATVAR.fits'), ssmooth)
            print s
            os.system(s)
            if sb:
                fit = getdata(stackname.replace('.fits', '.IM.fits'))
                hdu.data = fit * 3.125
                hdu.writeto(stackname.replace('.fits', '.STATVAR.SB.fits'), clobber=True)
            if jpeg: astroim(stackname.replace('.fits', '.STATVAR.SNR.IM.fits'), smooth=smooth, saveim=True, show=False,
                             cbfrac=.08,
                             pad=.006, dfig=(8, 9), contours=True, scb_label='SNR', title=title)
        if jpeg:
            astroim(stackname.replace('.fits', '.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                    pad=.006, dfig=(8, 9), contours=True, scb_label='Flux', title=title, vmin=-.5, vmax=.5)
            astroim(stackname.replace('.fits', '.SNR.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                    pad=.006, dfig=(8, 7), contours=True, scb_label='SNR', title=title, nsigma=7, vmin=-3, vmax=3, std=0.3)

        if tests:

            if jpeg:
                astroim(stackname.replace('.fits', '.NPIX.IM.fits'), smooth=False, saveim=True, show=False, vmin=0,
                        vmax=np.amax(nf) * 3,
                        cbfrac=.08, pad=.006, dfig=(8, 7), contours=False, scb_label='Number of voxels', title=title)
                astroim(stackname.replace('.fits', '.SNR.halftest1.IM.fits'), smooth=smooth, saveim=True, show=False,
                        cbfrac=.08, pad=.006, dfig=(8, 7), contours=True, scb_label='Flux [1e-20 cgs]',
                        title=title)  # , std=1)
                astroim(stackname.replace('.fits', '.SNR.halftest2.IM.fits'), smooth=smooth, saveim=True, show=False,
                        cbfrac=.08, pad=.006, dfig=(8, 7), contours=True, scb_label='Flux [1e-20 cgs]',
                        title=title)  # , std=1)
            if os.path.isfile(stackname.replace(st, 'random_' + st)):

                frand = getdata(stackname.replace(st, 'random_' + st))#.replace('n%d.' % nstack, 'n%d.' % (nstack*nrand)))
                frand[np.isnan(frand)] = 0
                hdu.data = f - frand
                hdu.writeto(stackname.replace('.fits', '.randtest.fits'), clobber=True)
                randf = [np.nanmean(frand[zl / 2 - zw0: zl / 2 + zw0, yl / 2 - 3:yl / 2 + 3, i:i + stdbin]) for i in stdbins]
                s = 'Cube2Im %s[*,*,%s]' % (stackname.replace('.fits', '.randtest.fits'), zw)
                print s
                os.system(s)
                if jpeg: astroim(stackname.replace('.fits', '.randtest.IM.fits'), smooth=smooth, saveim=True,
                                 show=False, cbfrac=.08, pad=.006, dfig=(7, 7), contours=True, scb_label='SNR',
                                 title=title)  # , std=1)

            h = "theta flux std npix half1 half2 flux_rand"
            all = [stdbins*.4, stack_flux, stack_npix, stack_stds, stack_stds_t1, stack_stds_t2, randf]
        else:
            h = "theta flux npix std"
            all = [stdbins*.4, stack_flux, stack_npix, stack_stds]

        if not random:
            np.savetxt('%s/%s/cats/%s_%s_%s_%s_%s_%s_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d.dat' %
                    (folder, fitcat, extraname, sn, sm, sw, st, imtype, dmin, dmax, pmin, pmax, rmin, rmax, fmin),
                    np.matrix(all).T, header=h)
