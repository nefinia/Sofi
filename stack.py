__author__ = 'gallegos'

from pyfits import getdata, PrimaryHDU
import numpy as np
import os
from os.path import isfile, isdir
from tools_sofi import astroim, rdarg, stack, pdfim, analysis
from sys import argv
from random import randint

fshear = False
half = rdarg(argv, 'half', bool, True)
reject = rdarg(argv, 'reject', bool, False)
meancorr = rdarg(argv, 'meancorr', bool, True)
parallel = rdarg(argv, 'parallel', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
flipshear = rdarg(argv, 'flipshear', bool, False)
fitcat = rdarg(argv, key='fitcat', type=str, default='all')  # 'all'#'HDFS'
scalelims = rdarg(argv, 'scalelims', str, '-0.01 0.01')
nodist = rdarg(argv, 'nodist', bool, True)
makepdf = rdarg(argv, 'makepdf', bool, True)
cubesmooth = rdarg(argv, 'cubesmooth', bool, False)
if cubesmooth:
    smooth = False
    ssmooth = ' -boxsm 3 '
else:
    smooth = True
    ssmooth = ''
simon = rdarg(argv, 'cubical', bool, True)
random = rdarg(argv, 'random', bool, True)
nrand = rdarg(argv, 'nrand', int, 200)
ntest = rdarg(argv, 'ntest', int, 200)
jin = rdarg(argv, 'jin', int, 0)
cubexmask = rdarg(argv, 'cubexmask', bool, False)
tests = rdarg(argv, 'tests', bool, True)
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
flipx = rdarg(argv, 'flipx', bool, False)
flipy = rdarg(argv, 'flipy', bool, False)
randflipy = rdarg(argv, 'randflipy', bool, False)
rotate = rdarg(argv, 'rotate', bool, False)
weights = rdarg(argv, 'weights', bool, False)
mask = rdarg(argv, 'mask', bool, True)
makeim = rdarg(argv, 'makeim', bool, True)
norm = rdarg(argv, 'norm', bool, False)
propvar = rdarg(argv, 'pvar', bool, False)
statvar = rdarg(argv, 'svar', bool, False)
prename = rdarg(argv, 'prename', str, '')
extraname = rdarg(argv, 'extraname', str, '')
imtype = rdarg(argv, 'imtype', str, 'mean')  # 'mean'
if imtype == 'median': extraname = '.' + imtype
sclip = rdarg(argv, 'sclip', int, 3)
if not prename: prename = ''
if not extraname: extraname = ''
sb = rdarg(argv, 'sb', bool, False)
binsize = rdarg(argv, 'binsize', int, 2)
vmin = rdarg(argv, 'vmin', float, -.05)
vmax = rdarg(argv, 'vmax', float, .5)
xmin = rdarg(argv, 'xmin', float, -40)
xmax = rdarg(argv, 'xmax', float, 40)
ymin = rdarg(argv, 'ymin', float, -40)
ymax = rdarg(argv, 'ymax', float, 40)
zmin = rdarg(argv, 'zmin', float, -10)
zmax = rdarg(argv, 'zmax', float, 10)
zw0 = rdarg(argv, 'zw0', int, 1)
yw0 = rdarg(argv, 'yw0', int, 2)
distmin = rdarg(argv, 'dmin', list, [.5])
distmax = rdarg(argv, 'dmax', list, [20])
pdmin = rdarg(argv, 'pdmin', list, [16])
pdmax = rdarg(argv, 'pdmax', list, [80])
redmin = rdarg(argv, 'rmin', list, [2.9])
redmax = rdarg(argv, 'rmax', list, [4.0])
lummin = rdarg(argv, 'lmin', list, [0])
lummax = rdarg(argv, 'lmax', list, [99999])
npairsmin = rdarg(argv, 'nmin', list, [1])
npairsmax = rdarg(argv, 'nmax', list, [99999])

if cubexmask:
    prename = '.Residuals'
    vmin = -.5
    vmax = .5
if mask: sm = 'masked'
else: sm = 'unmasked'
if weights: sw = 'weighted'
else: sw = 'unweighted'
if norm: sn = 'normalized'
else: sn = 'not-normalized'
if simon: st = 'cubical'
else: st = 'cylindrical'
if norm: lenx = 40

if fitcat == 'all':
    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'HDFS')
    data = getdata(cat2, 1)
    ids1 = data['id1']
    ids2 = data['id2']

    laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'HDFS'), 1)
    flae = laedata['LYALPHA_FLUX']
    idlae = laedata['ID']
    lum_lae1 = []
    lum_lae2 = []
    for i, j in zip(ids1.astype(int), ids2.astype(int)):
        lum_lae1.append(np.float(flae[idlae == i]))
        lum_lae2.append(np.float(flae[idlae == j]))
    sconf1 = data['sconf1'].astype(int)
    sconf2 = data['sconf2'].astype(int)
    pdists = data['x2'] - data['x1']
    zs1 = data['z1']
    zs2 = data['z2']
    redshift = (data['redshift1'] + data['redshift2']) / 2.
    dists = data['pi_Mpc']
    theta = data['theta']
    fitcats = ['HDFS'] * len(data)
    dclose = np.where((dists > .5) & (dists <= 20))
    npairs1 = np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids1])
    npairs2 = np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids2])
    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'UDF')
    data = getdata(cat2, 1)
    ids1 = np.concatenate((ids1, data['id1']), 0)
    ids2 = np.concatenate((ids2, data['id2']), 0)
    laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'UDF'), 1)
    flae = laedata['LYALPHA_FLUX']
    idlae = laedata['ID']
    for i, j in zip(data['id1'].astype(int), data['id2'].astype(int)):
        lum_lae1.append(np.float(flae[idlae == i]))
        lum_lae2.append(np.float(flae[idlae == j]))
    sconf1 = np.concatenate((sconf1, data['sconf1'].astype(int)), 0)
    sconf2 = np.concatenate((sconf2, data['sconf2'].astype(int)), 0)
    pdists = np.concatenate((pdists, data['x2'] - data['x1']), 0)
    zs1 = np.concatenate((zs1, data['z1']), 0)
    zs2 = np.concatenate((zs2, data['z2']), 0)
    redshift = np.concatenate((redshift, (data['redshift1'] + data['redshift2']) / 2.), 0)
    dists = np.concatenate((dists, data['pi_Mpc']), 0)
    theta = np.concatenate((theta, data['theta']), 0)
    nhdfs = len(fitcats)
    fitcats = np.concatenate((fitcats, ['UDF'] * len(data)), 0)
    lum_lae2 = np.array(lum_lae2)
    lum_lae1 = np.array(lum_lae1)
    dclose = np.where((data['pi_Mpc'] > .5) & (data['pi_Mpc'] <= 20))
    temp = np.array([np.sum(data['id1'][dclose] == i) + np.sum(data['id2'][dclose] == i) for i in data['id1']])
    npairs1 = np.concatenate((npairs1, temp))
    temp = np.array([np.sum(data['id1'][dclose] == i) + np.sum(data['id2'][dclose] == i) for i in data['id2']])
    npairs2 = np.concatenate((npairs2, temp))

else:
    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, fitcat)
    data = getdata(cat2, 1)
    ids1 = data['id1']
    ids2 = data['id2']
    sconf1 = data['sconf1'].astype(int)
    sconf2 = data['sconf2'].astype(int)
    laedata = getdata('%s/%s/cats/laes.fits' % (folder, fitcat), 1)
    idlae = laedata['ID']
    flae = laedata['LYALPHA_FLUX']
    lum_lae1 = []
    lum_lae2 = []
    for i, j in zip(ids1.astype(int), ids2.astype(int)):
        lum_lae1.append(np.float(flae[idlae == i]))
        lum_lae2.append(np.float(flae[idlae == j]))
    pdists = data['x2'] - data['x1']
    zs1 = data['z1']
    zs2 = data['z2']
    redshift = (data['redshift1'] + data['redshift2']) / 2.
    dists = data['pi_Mpc']
    theta = data['theta']
    fitcats = [fitcat] * len(data)
    dclose = np.where((dists > .5) & (dists <= 20))
    npairs1 = np.array([np.sum(ids1[dclose] == i)+np.sum(ids2[dclose] == i) for i in ids1])
    npairs2 = np.array([np.sum(ids1[dclose] == i)+np.sum(ids2[dclose] == i) for i in ids2])

h = 'theta '
hdu = PrimaryHDU()

if nodist:
    temp = ids1
    ids1 = np.concatenate((temp, ids2), 0).astype(int)
    ids2 = np.concatenate((ids2, temp), 0).astype(int)
    temp = sconf1
    sconf1 = np.concatenate((temp, sconf2), 0).astype(int)
    sconf2 = np.concatenate((sconf2, temp), 0).astype(int)
    redshift = np.concatenate((redshift, redshift), 0)
    dists = np.concatenate((dists, dists), 0)
    theta = np.concatenate((theta, theta), 0)
    temp = lum_lae1
    lum_lae1 = np.concatenate((temp, lum_lae2), 0)
    lum_lae2 = np.concatenate((lum_lae2, temp), 0)
    fitcats = np.concatenate((fitcats, fitcats), 0)
    temp = zs1
    zs1 = np.concatenate((temp, zs2), 0)
    zs2 = np.concatenate((zs2, temp), 0)
    npairs = np.concatenate((npairs1, npairs2), 0)

for dmin, dmax, rmin, rmax, pmin, pmax, fmin, fmax, nmin, nmax in zip(distmin, distmax, redmin, redmax, pdmin, pdmax,
                                                                      lummin,
                                                                      lummax, npairsmin, npairsmax):

    if reject:
        rejected = ((ids1 == 97) & (ids2 == 138)) | ((ids1 == 137) & (ids2 == 170)) | ((ids1 == 238) & (ids2 == 563)) \
                   | ((ids1 == 170) & (ids2 == 137)) | ((ids1 == 391) & (ids2 == 149)) | ((ids1 == 293) & (ids2 == 163)) \
                   | ((ids1 == 293) & (ids2 == 180)) | ((ids1 == 780) & (ids2 == 230)) | ((ids1 == 780) & (ids2 == 252)) \
                   | ((ids1 == 391) & (ids2 == 506)) | ((ids1 == 391) & (ids2 == 753)) | ((ids1 == 494) & (ids2 == 605)) \
                   | ((ids1 == 780) & (ids2 == 753)) | ((ids1 == 391) & (ids2 == 629))

        close = np.where((redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax) \
                         & (theta > pmin) #& (lum_lae1 >= fmin) & (lum_lae1 < fmax) \
                         & (sconf1 > 0) & (sconf2 > 0) & (ids1 != 144) & (~rejected) \
                         & (npairs >= nmin) & (npairs < nmax))[0]  # 144 agn discarded & (sconf > 0)
    else:
        close = np.where((redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax) \
                         & (theta > pmin) #& (lum_lae1 >= fmin) & (lum_lae1 < fmax) \
                         & (sconf1 > 0) & (sconf2 > 0) & (ids1 != 144) \
                         & (npairs >= nmin) & (npairs < nmax))[0]  # 144 agn discarded & (sconf > 0)

    if len(close) > 1:
        print 'Number of subcubes to be stacked', len(close)
        id1 = ids1[close]
        id2 = ids2[close]
        z1 = zs1[close]
        z2 = zs2[close]
        fcats = np.array(fitcats)[close]
        lst = []

        random_lst = [[] for i in range(nrand)]

        for i in range(len(id1)):

            pair = '%s/%s/pairs/%s/%s/%s_pair%d-%d%s.fits' % (folder, fcats[i], sn, sm, st, id1[i], id2[i], prename)
            # print "pair", pair
            if os.path.isfile(pair): lst.append(pair)

            if random:
                for j in range(nrand):
                    pair = '%s/%s/pairs/%s/%s/random_%s_pair%d-%d%s.%d.fits' % (
                        folder, fcats[i], sn, sm, st, id1[i], id2[i], prename, j + jin)
                    if os.path.isfile(pair): random_lst[j].append(pair)

        nstack = len(lst)

        foldername = '%s/%s/stacks/%s_%s_%s_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d-%d%s%s' % (
            folder, fitcat, sn, sm, st, dmin, dmax,
            pmin, pmax, rmin, rmax, fmin, fmax, prename, extraname)

        if not isdir(foldername): os.system('mkdir %s' % foldername)

        lstname = '%s/stack.lst' % (foldername)
        flst = open(lstname, 'w')
        for l in lst: flst.write('%s\n' % l)
        flst.close()

        stackname = lstname.replace('.lst', '.fits')

        title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $%d<L<%d$' % (
            nstack, dmin, dmax, pmin, pmax, rmin, rmax, fmin, fmax)


        f, nf, var = stack(lst, stackname, imtype, 100, zw0, makeim, title, vmin, vmax, npix=True, var=True,
                           flipshear=fshear, corrmean=meancorr, flipy=flipy, overwrite=overwrite)

        if makepdf:
            imout = lstname.replace('.lst', '.pdf')
            if not os.path.isfile(imout) or overwrite:
                lstIM = [l.replace('.fits', '.IM.fits') for l in lst]
                pdfim(lstIM, imout=imout, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, contours=True)

        zl, yl, xl = f.shape
        zw = '%d:%d' % (zl / 2 - zw0 + 1, zl / 2 + zw0 + 1)
        stdbin = 1  # 6/binsize
        stdbins = stdbin * np.arange(xl / stdbin)

        if propvar:
            proplst = [l.relace('.fits', '.PROPVAR.fits') for l in lst]
            pvar = stack(lst, stackname, imtype='var', sclip=sclip, zw=zw0, makeim=makeim, title=title, vmin=vmin,
                         vmax=vmax, npix=False, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
                         overwrite=overwrite)

        if random:
            avr_right = []
            avr_left = []
            avr_top = []
            avr_bottom = []
            frs = []
            all = np.array([(stdbins - stdbins[-1] / 2.) * (binsize * .2)])
            h = 'theta'
            if not isdir(foldername+'/randoms'): os.system('mkdir %s/randoms' % foldername)
            for i in range(nrand):
                rlstname = '%s/randoms/random_stack.%d.lst' % (foldername, i)
                rstackname = rlstname.replace('.lst', '.fits')
                title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $lum>%d$' % (
                    len(random_lst[i]), dmin, dmax, pmin, pmax, rmin, rmax, fmin)
                if len(random_lst[i]) > 0:
                    fr, nfr = stack(random_lst[i], rstackname, imtype, sclip, zw0, makeim, title, vmin, vmax,
                                    npix=True, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
                                    overwrite=overwrite, std=0.0707064)
                    flst = open(rlstname, 'w')
                    for l in random_lst[i]: flst.write('%s\n' % l)
                    flst.close()
                    aa, _aa, hh = analysis(fr, nfr, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
                    h += ' SB%d' % i
                    all = np.concatenate((all, np.array(aa[1:])))
                    avr_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]) * 4)
                    avr_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]) * 4)
                    avr_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]) * 4)
                    avr_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]) * 4)
                    if len(random_lst[i]) == nstack: frs.append(rlstname.replace('.lst', '.fits'))
                    if makepdf:
                        pdfname = rlstname.replace('.lst', '.pdf')
                        if not os.path.isfile(pdfname) or overwrite:
                            lstIM = [l.replace('.fits', '.IM.fits') for l in random_lst[i]]
                            pdfim(lstIM, imout=pdfname, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, contours=True)

            if ntest > nrand:
                for i in range(nrand, ntest):
                    rlstname = '%s/randoms/random_stack.%d.lst' % (foldername, i)
                    rstackname = rlstname.replace('.lst', '.fits')
                    flst = open(rlstname, 'w')
                    full = np.where([len(random_lst[i]) == nstack for i in range(nrand)])[0]
                    lists = [l.replace('.fits','.%d.fits' % np.random.choice(full)) for l in lst]
                    for l in lists: flst.write('%s\n' % l)
                    flst.close()
                    title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $lum>%d$' % (
                        nrand, dmin, dmax, pmin, pmax, rmin, rmax, fmin)

                    fr, nfr = stack(lists, rstackname, imtype, sclip, zw0, makeim, title, vmin, vmax,
                                    npix=True, var=False, flipshear=fshear, corrmean=meancorr, flipy=flipy,
                                    overwrite=overwrite, std=0.0707064)

                    aa, _aa, hh = analysis(fr, nfr, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
                    h += ' SB%d' % (i + nrand)
                    all = np.concatenate((all, np.array(aa[1:])))
                    avr_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]) * 4)
                    avr_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]) * 4)
                    avr_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]) * 4)
                    avr_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]) * 4)
                    # frs.append(rlstname.replace('.lst', '.fits'))
                    if makepdf:
                        pdfname = rlstname.replace('.lst', '.pdf')
                        if not os.path.isfile(pdfname) or overwrite:
                            lstIM = [l.replace('.fits', '.IM.fits') for l in lists]
                            pdfim(lstIM, imout=pdfname, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, contours=True)

            superrandom = '%s/random_stack.fits' % (foldername)
            nfrs = [fff.replace('.fits', '.NPIX.fits') for fff in frs]
            print "Creating super-random file"
            frs = stack(frs, superrandom, imtype, sclip, zw0, makeim, title, vmin, vmax,
                        npix=False, var=False,
                        flipshear=fshear, flipy=flipy, overwrite=overwrite, std=0.0707064)  # sclip=100 -> no sigma clipping on the superrandom



            nfrs = stack(nfrs, superrandom.replace('.fits', '.NPIX.fits'), 'flux', sclip, zw0, makeim, title, vmin=None,
                         vmax=None, npix=False, var=False, corrmean=False, flipy=flipy, overwrite=overwrite)

            np.savetxt('%s/random_stack.dat' % foldername, np.matrix(all).T, header=h)

        if tests:
            avt_right = []
            avt_left = []
            avt_top = []
            avt_bottom = []
            fts = []
            all = np.array([(stdbins - stdbins[-1] / 2.) * (binsize * .2)])
            h = 'theta'
            if not isdir(foldername+'/tests'): os.system('mkdir %s/tests' % foldername)
            if randflipy:
                if not isdir(foldername+'/flipy'): os.system('mkdir %s/flipy' % foldername)

            for i in range(ntest):
                sample = np.random.choice(lst, nstack)
                print "Stacking bootstrapped sample of %d subcubes" % len(sample)
                title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $L>%d$' % (
                    len(sample), dmin, dmax, pmin, pmax, rmin, rmax, fmin)
                tlstname = '%s/tests/stack.t%d.lst' % (foldername, i)
                tstackname = tlstname.replace('.lst', '.fits')
                flst = open(tlstname, 'w')
                for l in sample: flst.write('%s\n' % l)
                flst.close()

                ftest, nftest = stack(sample, tstackname, imtype, sclip, zw0, makeim, title,
                                      vmin, vmax, npix=True, var=False, corrmean=meancorr, overwrite=overwrite)

                aa, _aa, hh = analysis(ftest, nftest, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
                h += ' SB%d' % i
                all = np.concatenate((all, np.array(aa[1:])))
                avt_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]) * 4)
                avt_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]) * 4)
                avt_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]) * 4)
                avt_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]) * 4)
                fts.append(tlstname.replace('.lst', '.fits'))
                if makepdf:
                    pdfname = tlstname.replace('.lst', '.pdf')
                    if not isfile(pdfname) or overwrite:
                        lstIM = [l.replace('.fits', '.IM.fits') for l in sample]
                        pdfim(lstIM, imout=pdfname, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, contours=True)

                if randflipy:
                    ylstname = '%s/flipy/stack.y%d.lst' % (foldername, i)
                    ystackname = ylstname.replace('.lst', '.fits')

                    fy, ny = stack(lst, ystackname, imtype, sclip, zw0, makeim, title, vmin, vmax, npix=True, var=False,
                                   corrmean=meancorr, randflipy=True, overwrite=overwrite)

            fts = stack(fts, '%s/stack.t.fits' % (foldername), imtype, sclip, zw0, makeim, title, vmin, vmax,
                        npix=False, var=False, flipy=flipy, overwrite=overwrite, std=0.0707064)
            h += '\n'
            np.savetxt('%s/stack.t.dat' % foldername, np.matrix(all).T, header=h)

            if random:
                s = 'CubeArit  %s - %s %s' % (
                    stackname, '%s/random_stack.fits' % (foldername), stackname.replace('.fits', '.randtest.fits'))
                print s
                os.system(s)
                s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (
                    stackname.replace('.fits', '.randtest.fits'), zw, stackname.replace('.fits', '.randtest.IM.fits'))
                print s
                os.system(s)
                astroim(stackname.replace('.fits', '.randtest.IM.fits'), smooth=smooth, saveim=makeim,
                         show=False,
                         cbfrac=.08, pad=.006, dfig=(8, 10), contours=True, scb_label=r'Flux [$10^{-20}\,\rm{erg/s/cm^2}$]',
                         title=title, vmin=vmin, vmax=vmax)


        if not random: frs = None
        all, _all, h = analysis(f, nf, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0, frs)
        avf_right = [np.nanmean(all[1][(all[0] <= 6) & (all[0] >= 2)]) * 4]
        avf_left = [np.nanmean(all[1][(all[0] <= -2) & (all[0] >= -6)]) * 4]
        avf_top = [np.nanmean(_all[1][(_all[0] <= 6) & (_all[0] >= 2)]) * 4]
        avf_bottom = [np.nanmean(_all[1][(_all[0] <= -2) & (_all[0] >= -6)]) * 4]
        np.savetxt('%s/stack_n%d.dat' % (foldername, nstack), np.matrix(all).T, header=h)
        np.savetxt('%s/stack_n%d_vertical.dat' % (foldername, nstack), np.matrix(_all).T, header=h)

        if tests and random:
            avp_right = []
            avp_left = []
            avp_top = []
            avp_bottom = []
            for l in lst:
                fp = getdata(l)
                if isfile(l.replace('.fits', '.NPIX.fits')): nfp = getdata(l.replace('.fits', '.NPIX.fits'))
                else: nfp = fp
                aa, _aa, hh = analysis(fp, nfp, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
                avp_right.append(np.nanmean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)]) * 4)
                avp_left.append(np.nanmean(aa[1][(aa[0] <= -2) & (aa[0] >= -6)]) * 4)
                avp_top.append(np.nanmean(_aa[1][(_aa[0] <= 6) & (_aa[0] >= 2)]) * 4)
                avp_bottom.append(np.nanmean(_aa[1][(_aa[0] <= -2) & (_aa[0] >= -6)]) * 4)
                del fp, nfp
            all = [avp_right, avp_top, avp_left, avp_bottom]
            np.savetxt('%s/subcubes.avs.dat' % foldername, np.matrix(all).T,
                       header='subcubes_r subcubes_t subcubes_l subcubes_b')

            all = [avf_right * ntest, avr_right, avt_right, avf_top * ntest, avr_top, avt_top,
                   avf_left * ntest, avr_left, avt_left, avf_bottom * ntest, avr_bottom, avt_bottom]
            np.savetxt('%s/stack.avs.dat' % foldername, np.matrix(all).T,
                       header='full_r random_r sets_r full_t random_t sets_t full_l random_l sets_l full_b random_b sets_b')

        if 1:
            f_snr = f * nf / np.sqrt(var)
            hdu.data = f_snr
            hdu.writeto(stackname.replace('.fits', '.SNR.fits'), clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -varcube %s -snrmap %s -imtype flux' % (stackname, zw,
                                                                     stackname.replace('.fits', '.VAR.fits'),
                                                                     stackname.replace('.fits', '.SNR.IM.fits'))

            print s
            os.system(s)

            astroim(stackname.replace('.fits', '.SNR.IM.fits'), smooth=smooth, saveim=makeim, show=False,
                         cbfrac=.08,
                         pad=.006, dfig=(8, 10), contours=True, scb_label='SNR', title=title, nsigma=10, vmin=-.1,
                         vmax=1,
                         std=0.3)

        pairs_data = [id1, id2, redshift[close], dists[close], theta[close],npairs[close]]
        np.savetxt('%s/stack_pairs.dat' % foldername, np.matrix(pairs_data).T,
                   header="#id1 id2 redshift comoving_distance theta n_neigbors")
