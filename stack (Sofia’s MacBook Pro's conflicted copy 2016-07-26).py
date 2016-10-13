__author__ = 'gallegos'

from pyfits import getdata, PrimaryHDU
import numpy as np
import os
from tools_sofi import astroim, rdarg, stack, pdfim, analysis
from sys import argv
from random import randint

half = rdarg(argv, 'half', bool, True)
overwrite = rdarg(argv, 'overwrite', bool, False)
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
nrand = rdarg(argv, 'nrand', int, 100)
ntest = rdarg(argv, 'ntest', int, 100)
jin = rdarg(argv, 'jin', int, 0)
cubexmask = rdarg(argv, 'cubexmask', bool, False)
tests = rdarg(argv, 'tests', bool, True)
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
prename = rdarg(argv, 'prename', str, '')
extraname = rdarg(argv, 'extraname', str, '')
imtype = rdarg(argv, 'imtype', str, 'mean')  # 'mean'
if imtype == 'median':
    extraname = '.' + imtype
sclip = rdarg(argv, 'sclip', int, 5)
if not prename:
    prename = ''
if not extraname:
    extraname = ''
sb = rdarg(argv, 'sb', bool, False)
binsize = rdarg(argv, 'binsize', int, 2)
vmin = rdarg(argv, 'vmin', float, -.5)
vmax = rdarg(argv, 'vmax', float, .5)
xmin = rdarg(argv, 'xmin', float, -40)
xmax = rdarg(argv, 'xmax', float, 40)
ymin = rdarg(argv, 'ymin', float, -40)
ymax = rdarg(argv, 'ymax', float, 40)
zmin = rdarg(argv, 'zmin', float, -10)
zmax = rdarg(argv, 'zmax', float, 10)
zw0 = rdarg(argv, 'zw0', int, 1)
yw0 = rdarg(argv, 'yw0', int, 3)
distmin = rdarg(argv, 'dmin', list, [.5])
distmax = rdarg(argv, 'dmax', list, [20])
pdmin = rdarg(argv, 'pdmin', list, [16])
pdmax = rdarg(argv, 'pdmax', list, [80])
redmin = rdarg(argv, 'redmin', list, [2.9])
redmax = rdarg(argv, 'redmax', list, [4.0])
fluxmin = rdarg(argv, 'fmin', list, [0])
if cubexmask:
    prename = '.Residuals'
    vmin = -.5
    vmax = .5

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

if simon:
    st = 'cubical'
else:
    st = 'cylindrical'

if norm:
    lenx = 40

if fitcat == 'all':
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
    sconf1 = data['sconf1'].astype(int)
    sconf2 = data['sconf2'].astype(int)
    pdists = data['x2'] - data['x1']
    zs1 = data['z1']
    redshift = (data['redshift1'] + data['redshift2']) / 2.
    dists = data['pi_Mpc']
    theta = data['theta']
    fitcats = ['HDFS'] * len(data)

    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'UDF')
    data = getdata(cat2, 1)
    ids1 = np.concatenate((ids1, data['id1']), 0)
    ids2 = np.concatenate((ids2, data['id2']), 0)

    laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'UDF'), 1)
    flae = laedata['LYALPHA_FLUX']
    idlae = laedata['ID']
    for i, j in zip(data['id1'].astype(int), data['id2'].astype(int)):
        flux_lae1.append(np.float(flae[idlae == i]))
        flux_lae2.append(np.float(flae[idlae == j]))
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
    flux_lae1 = []
    flux_lae2 = []
    for i, j in zip(ids1.astype(int), ids2.astype(int)):
        flux_lae1.append(np.float(flae[idlae == i]))
        flux_lae2.append(np.float(flae[idlae == j]))
    pdists = data['x2'] - data['x1']
    zs1 = data['z1']
    redshift = (data['redshift1'] + data['redshift2']) / 2.
    dists = data['pi_Mpc']
    theta = data['theta']
    fitcats = [fitcat] * len(data)

ns = []
h = 'theta '
hdu = PrimaryHDU()

for dmin, dmax, rmin, rmax, pmin, pmax, fmin in zip(distmin, distmax, redmin, redmax, pdmin, pdmax, fluxmin):

    close = np.where((redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax)
                     & (theta > pmin) & (flux_lae1 >= fmin) & (flux_lae2 >= fmin) & (sconf1 > 0)# & (sconf2 > 0)
                     & (ids1 != 144) & (ids2 != 144))[0]  # 144 agn discarded & (sconf > 0)

    if len(close) > 1:
        ns.append([dmin, dmax, rmin, rmax, len(close)])
        print 'Number of pairs to be stacked', len(close)
        id1 = ids1[close]
        id2 = ids2[close]
        z1 = zs1[close]
        fcats = np.array(fitcats)[close]
        lst = []

        random_lst = [[] for i in range(nrand)]

        for i in range(len(id1)):

            pair = '%s/%s/pairs/%s/%s/%s_pair%d-%d%s.fits' % (folder, fcats[i], sn, sm, st, id1[i], id2[i], prename)
            if os.path.isfile(pair): lst.append(pair)

            if nodist:
                pair = '%s/%s/pairs/%s/%s/%s_pair%d-%d%s.fits' % (folder, fcats[i], sn, sm, st, id2[i], id1[i], prename)
                if os.path.isfile(pair): lst.append(pair)
            if random:
                for j in range(nrand):
                    pair = '%s/%s/pairs/%s/%s/random_%s_pair%d-%d%s.%d.fits' % (
                        folder, fcats[i], sn, sm, st, id1[i], id2[i], prename, j+jin)
                    if os.path.isfile(pair): random_lst[j].append(pair)

                    if nodist:
                        pair = '%s/%s/pairs/%s/%s/random_%s_pair%d-%d%s.%d.fits' % (
                            folder, fcats[i], sn, sm, st, id2[i], id1[i], prename, j+jin)
                        if os.path.isfile(pair): random_lst[j].append(pair)

        nstack = len(lst)
        foldername = '%s/%s/stacks/%s_%s_%s_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d%s%s' % (folder, fitcat, sn, sm, st, dmin, dmax,
                                                                              pmin, pmax, rmin, rmax, fmin, prename, extraname)

        os.system('mkdir %s' % foldername)

        lstname = '%s/stack.lst' % (foldername)
        flst = open(lstname, 'w')
        for l in lst: flst.write('%s\n' % l)
        flst.close()

        stackname = lstname.replace('.lst', '.fits')

        title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $flux>%d$' % (
            nstack, dmin, dmax, pmin, pmax, rmin, rmax, fmin)

        if not os.path.isfile(stackname) or overwrite:
            f, nf, var = stack(lst, stackname, imtype, sclip, zw0, jpeg, title, vmin, vmax, npix=True, var=True)
            astroim(stackname.replace('.fits', '.NPIX.IM.fits'), smooth=False, vmin=np.amin(nf), vmax=np.amax(nf), show=False)
        else:
            print '%s already exists.' % stackname
            f = getdata(stackname)
            nf = getdata(stackname.replace('.fits', '.NPIX.fits'))

        if makepdf:
            imout = lstname.replace('.lst', '.pdf')
            if not os.path.isfile(imout) or overwrite:
                lstIM = [l.replace('.fits', '.IM.fits') for l in lst]
                pdfim(lstIM, imout=imout, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

        zl, yl, xl = f.shape
        zw = '%d:%d' % (zl / 2 - zw0 + 1, zl / 2 + zw0 + 1)
        stdbin = 1#6/binsize
        stdbins = stdbin * np.arange(xl / stdbin)

        if propvar:
            proplst = [l.relace('.fits', '.PROPVAR.fits') for l in lst]
            pvar = stack(lst, stackname, imtype='var', sclip=sclip, zw=zw0, jpeg=jpeg, title=title, vmin=vmin, vmax=vmax, npix=False, var=False)


        if random:
            avr = []
            frs = []
            all = np.array([(stdbins - stdbins[-1]/2.) * (binsize*.2)])
            h = 'theta'
            for i in range(nrand):
                rlstname = '%s/random_stack.%d.lst' % (foldername, i)
                rstackname = rlstname.replace('.lst', '.fits')
                flst = open(rlstname, 'w')
                for l in random_lst[i]: flst.write('%s\n' % l)
                flst.close()
                title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $flux>%d$' % (
                len(random_lst[i]), dmin, dmax, pmin, pmax, rmin, rmax, fmin)
                if not os.path.isfile(rstackname) or overwrite:
                    fr, nfr = stack(random_lst[i], rstackname, imtype, sclip, zw0, jpeg, title, vmin, vmax,
                                  npix=True, var=False)
                else:
                    print '%s already exists.' % rstackname
                    fr = getdata(rstackname)
                    nfr = getdata(rstackname.replace('.fits', '.NPIX.fits'))
                aa, _aa, hh = analysis(fr, nfr, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
                h += ' SB%d' % i
                all = np.concatenate((all, np.array(aa[1:])))
                avr.append(np.mean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)])*4)
                frs.append(rlstname.replace('.lst', '.fits'))
                if makepdf:
                    pdfname = rlstname.replace('.lst', '.pdf')
                    if not os.path.isfile(pdfname) or overwrite:
                        lstIM = [l.replace('.fits', '.IM.fits') for l in random_lst[i]]
                        pdfim(lstIM, imout=pdfname, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

            if ntest > nrand:
                diff = ntest - nrand
                for i in range(diff):
                    rlstname = '%s/random_stack.%d.lst' % (foldername, i+nrand)
                    rstackname = rlstname.replace('.lst', '.fits')
                    flst = open(rlstname, 'w')
                    lists = [random_lst[randint(0, nrand-1)][i] for i in range(nstack)]
                    for l in lists: flst.write('%s\n' % l)
                    flst.close()
                    title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $flux>%d$' % (
                    nrand, dmin, dmax, pmin, pmax, rmin, rmax, fmin)
                    if not os.path.isfile(rstackname) or overwrite:
                        fr, nfr = stack(lists, rstackname, imtype, sclip, zw0, jpeg, title, vmin, vmax,
                                      npix=True, var=False)
                    else:
                        print '%s already exists.' % rstackname
                        fr = getdata(rstackname)
                        nfr = getdata(rstackname.replace('.fits', '.NPIX.fits'))
                    aa, _aa, hh = analysis(fr, nfr, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
                    h += ' SB%d' % (i+nrand)
                    all = np.concatenate((all, np.array(aa[1:])))
                    avr.append(np.mean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)])*4)
                    frs.append(rlstname.replace('.lst', '.fits'))
                    if makepdf:
                        pdfname = rlstname.replace('.lst', '.pdf')
                        if not os.path.isfile(pdfname) or overwrite:
                            lstIM = [l.replace('.fits', '.IM.fits') for l in lists]
                            pdfim(lstIM, imout=pdfname, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

            superrandom = '%s/random_stack.fits' % (foldername)

            if not os.path.isfile(superrandom) or overwrite:
                frs = stack(frs, superrandom, imtype, 100, zw0, jpeg, title, vmin, vmax,
                                  npix=False, var=False)
            h += '\n'
            np.savetxt('%s/random_stack.dat' % foldername, np.matrix(all).T, header=h)

            superlst = np.concatenate(random_lst)
            rlstname = '%s/random_stack_b.lst' % foldername
            flst = open(rlstname, 'w')
            for l in superlst: flst.write('%s\n' % l)
            flst.close()
            superrandom = '%s/random_stack_b.fits' % foldername
            print "aaaaaaaaaaaaaaa"
            if not os.path.isfile(superrandom) or overwrite:
                frs = stack(superlst, superrandom, imtype, sclip, zw0, jpeg, title, vmin, vmax,
                                  npix=False, var=False)

        if tests:
            avt = []
            fts = []
            all = np.array([(stdbins - stdbins[-1]/2.) * (binsize*.2)])
            h = 'theta'
            for i in range(ntest):
                sample = np.random.choice(lst, nstack / 2)
                title = r'$N\,\,%d$, $%d<d<%d$, $%d<\theta<%d$, $%.1f<z<%.1f$, $flux>%d$' % (
                len(sample), dmin, dmax, pmin, pmax, rmin, rmax, fmin)
                tlstname = '%s/stack.t%d.lst' % (foldername, i)
                tstackname = tlstname.replace('.lst', '.fits')
                flst = open(tlstname, 'w')
                for l in sample: flst.write('%s\n' % l)
                flst.close()
                if not os.path.isfile(tstackname) or overwrite:
                    ftest, nftest = stack(sample, tstackname, imtype, sclip, zw0, jpeg, title,
                                 vmin, vmax, npix=True, var=False)
                else:
                    print '%s already exists.' % tstackname
                    ftest = getdata(tstackname)
                    nftest = getdata(tstackname.replace('.fits', '.NPIX.fits'))
                aa, _aa, hh = analysis(ftest, nftest, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0)
                h += ' SB%d' % i
                all = np.concatenate((all, np.array(aa[1:])))
                avt.append(np.mean(aa[1][(aa[0] <= 6) & (aa[0] >= 2)])*4)
                fts.append(tlstname.replace('.lst', '.fits'))
                if makepdf:
                    pdfname = tlstname.replace('.lst', '.pdf')
                    if not os.path.isfile(pdfname) or overwrite:
                        lstIM = [l.replace('.fits', '.IM.fits') for l in sample]
                        pdfim(lstIM, imout=pdfname, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
            fts = stack(fts, '%s/stack.t.fits' % (foldername), imtype, 100, zw0, jpeg, title, vmin, vmax,
                              npix=False, var=False)
            h += '\n'
            np.savetxt('%s/stack.t.dat' % foldername, np.matrix(all).T, header=h)

            if random:
                s = 'CubeArit  %s - %s %s' % (
                    stackname, '%s/random_stack.fits'  % (foldername), stackname.replace('.fits', '.randtest.fits'))
                print s
                os.system(s)
                s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (
                    stackname.replace('.fits', '.randtest.fits'), zw, stackname.replace('.fits', '.randtest.IM.fits'))
                print s
                os.system(s)
                if jpeg: astroim(stackname.replace('.fits', '.randtest.IM.fits'), smooth=smooth, saveim=True,
                                 show=False,
                                 cbfrac=.08, pad=.006, dfig=(8, 10), contours=True, scb_label='Flux [1e-20 cgs]',
                                 title=title, vmin=-.5, vmax=.5)

        all, _all, h = analysis(f, nf, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0, frs)
        avf = [np.mean(all[1][(all[0] <= 6) & (all[0] >= 2)])*4]
        np.savetxt('%s/stack_n%d.dat' % (foldername, nstack), np.matrix(all).T, header=h)
        np.savetxt('%s/stack_n%d_vertical.dat' % (foldername, nstack), np.matrix(_all).T, header=h)

        if tests and random:
            all = [avf*ntest, avr, avt]
            np.savetxt('%s/stack.avs.dat' % (foldername), np.matrix(all).T, header='full random sets')

        if 0:
            f_snr = f / np.sqrt(var)
            hdu.data = f_snr
            hdu.writeto(stackname.replace('.fits', '.SNR.fits'), clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -varcube %s -snrmap %s' % (stackname, zw,
                stackname.replace('.fits', '.VAR.fits'), stackname.replace('.fits', '.SNR.IM.fits'))

            print s
            os.system(s)

            if jpeg: astroim(stackname.replace('.fits', '.SNR.IM.fits'), smooth=smooth, saveim=True, show=False, cbfrac=.08,
                             pad=.006, dfig=(8, 10), contours=True, scb_label='SNR', title=title, nsigma=7, vmin=-3, vmax=3,
                            std=0.3)

        pairs_data = [id1, id2, redshift[close], dists[close], theta[close]]
        np.savetxt('%s/stack_pairs.dat' % (foldername), np.matrix(pairs_data).T, header="#id1 id2 redshift comoving_distance theta")

