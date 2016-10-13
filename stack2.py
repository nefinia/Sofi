__author__ = 'gallegos'

from pyfits import getdata, PrimaryHDU
import numpy as np
import os
from tools_sofi import astroim, rdarg, stack, pdfim, analysis
from sys import argv
from random import randint

if 0:
    fitcat = rdarg(argv, key='fitcat', type=str, default='all')  # 'all'#'HDFS'
    folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'

    if fitcat == 'all':
        cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'HDFS')
        data = getdata(cat2, 1)
        ids1 = data['id1']
        ids2 = data['id2']

        laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'HDFS'), 1)
        flae = laedata['Luminosity']
        idlae = laedata['ID']

        zs = laedata['Z']
        flux_lae1 = []
        flux_lae2 = []
        for i, j in zip(ids1.astype(int), ids2.astype(int)):
            flux_lae1.append(np.float(flae[idlae == i]))
            flux_lae2.append(np.float(flae[idlae == j]))
        sconf1 = data['sconf1'].astype(int)
        sconf2 = data['sconf2'].astype(int)
        pdists = data['x2'] - data['x1']
        zs1 = data['z1']
        zs2 = data['z2']
        redshift = (data['redshift1'] + data['redshift2']) / 2.
        dists = data['pi_Mpc']
        theta = data['theta']
        fitcats = ['HDFS'] * len(data)
        fitcat = ['HDFS'] * len(laedata)

        dclose = np.where((dists > .5) & (dists <= 20))
        npair = np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in idlae])
        npairs1 = np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids1])
        npairs2 = np.array([np.sum(ids1[dclose] == i) + np.sum(ids2[dclose] == i) for i in ids2])
        cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'UDF')
        data = getdata(cat2, 1)
        ids1 = np.concatenate((ids1, data['id1']), 0)
        ids2 = np.concatenate((ids2, data['id2']), 0)

        laedata = getdata('%s/%s/cats/laes.fits' % (folder, 'UDF'), 1)
        flaes = np.concatenate((flae, laedata['Luminosity']))
        idlaes = np.concatenate((idlae, laedata['ID']))
        flae = laedata['Luminosity']
        idlae = laedata['ID']
        zs = np.concatenate((zs, laedata['Z_MUSE']))

        for i, j in zip(data['id1'].astype(int), data['id2'].astype(int)):
            flux_lae1.append(np.float(flae[idlae == i]))
            flux_lae2.append(np.float(flae[idlae == j]))
        sconf1 = np.concatenate((sconf1, data['sconf1'].astype(int)), 0)
        sconf2 = np.concatenate((sconf2, data['sconf2'].astype(int)), 0)
        pdists = np.concatenate((pdists, data['x2'] - data['x1']), 0)
        zs1 = np.concatenate((zs1, data['z1']), 0)
        zs2 = np.concatenate((zs2, data['z2']), 0)
        redshift = np.concatenate((redshift, (data['redshift1'] + data['redshift2']) / 2.), 0)

        dclose = np.where((data['pi_Mpc'] > .5) & (data['pi_Mpc'] <= 20))
        temp = np.array([np.sum(data['id1'][dclose] == i) + np.sum(data['id2'][dclose] == i) for i in data['id1']])
        npairs1 = np.concatenate((npairs1, temp))
        temp = np.array([np.sum(data['id1'][dclose] == i) + np.sum(data['id2'][dclose] == i) for i in data['id2']])
        npairs2 = np.concatenate((npairs2, temp))
        temp = np.array([np.sum(data['id1'][dclose] == i) + np.sum(data['id2'][dclose] == i) for i in idlae])
        npair = np.concatenate((npair, temp))
        dists = np.concatenate((dists, data['pi_Mpc']), 0)
        theta = np.concatenate((theta, data['theta']), 0)
        nhdfs = len(fitcats)
        fitcats = np.concatenate((fitcats, ['UDF'] * len(data)), 0)
        fitcat = np.concatenate((fitcat, ['UDF'] * len(laedata)), 0)
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
        idlaes = idlae
        if fitcat == 'UDF':
            zs = laedata['Z_MUSE']
        else:
            zs = laedata['Z']
        flae = laedata['Luminosity']
        flux_lae1 = []
        flux_lae2 = []
        for i, j in zip(ids1.astype(int), ids2.astype(int)):
            if len(flae[idlae == i]) > 0:
                flux_lae1.append(np.float(flae[idlae == i]))
            else:
                flux_lae1.append(np.nan)

            if len(flae[idlae == j]) > 0:
                flux_lae2.append(np.float(flae[idlae == j]))
            else:
                flux_lae2.append(np.nan)
        pdists = data['x2'] - data['x1']
        zs1 = data['z1']
        zs2 = data['z2']
        redshift = (data['redshift1'] + data['redshift2']) / 2.
        dists = data['pi_Mpc']
        theta = data['theta']
        fitcats = [fitcat] * len(data)

if 0:
    ns = []
    h = 'theta '
    hdu = PrimaryHDU()
    nodist = True

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
        temp = flux_lae1
        flux_lae1 = np.concatenate((temp, flux_lae2), 0)
        flux_lae2 = np.concatenate((flux_lae2, temp), 0)
        fitcats = np.concatenate((fitcats, fitcats), 0)
        temp = zs1
        zs1 = np.concatenate((temp, zs2), 0)
        zs2 = np.concatenate((zs2, temp), 0)
        npairs = np.concatenate((npairs1, npairs2), 0)
        np.savetxt('%s/%s/cats/lae_npairs.dat' % (folder, fitcat), (npairs), fmt='%d')

    # specially discarded: 97-138 137-170 170-137 238-563
    dmin = .5
    dmax = 20
    pmin = 16
    pmax = 80
    rmin = 2.9
    rmax = 4.0
    fmin = 0
    fmax = 9999999

    rejected = ((ids1 == 97) & (ids2 == 138)) | ((ids1 == 137) & (ids2 == 170)) | ((ids1 == 238) & (ids2 == 563)) \
               | ((ids1 == 170) & (ids2 == 137)) | ((ids1 == 391) & (ids2 == 149)) | ((ids1 == 293) & (ids2 == 163)) \
               | ((ids1 == 293) & (ids2 == 180)) | ((ids1 == 780) & (ids2 == 230)) | ((ids1 == 780) & (ids2 == 252)) \
               | ((ids1 == 391) & (ids2 == 506)) | ((ids1 == 391) & (ids2 == 753)) | ((ids1 == 494) & (ids2 == 605)) \
               | ((ids1 == 780) & (ids2 == 753)) | ((ids1 == 391) & (ids2 == 629))

    close = np.where((redshift <= rmax) & (redshift > rmin) & (dists <= dmax) & (dists > dmin) & (theta <= pmax)
                     & (theta > pmin) & (sconf1 > 0) & (sconf2 > 0) & (ids1 != 144) & (ids2 != 144))[0]
    # & (~rejected) & (npairs >=6))[0]  # 144 agn discarded
    # & (flux_lae1 >= fmin) & (flux_lae2 >= fmin)
    # from scipy import stats

    mflux = np.nanmedian(flux_lae1[close])
    print mflux
    mred = np.nanmedian(redshift[close])
    print mred
    mtheta = np.nanmedian(theta[close])
    print mtheta
    mdists = np.nanmedian(dists[close])
    print mdists
    mnpairs = np.nanmedian(npairs[close])
    print mnpairs

if 1:
    "median will all number of pairs and no rejected subcubes"
    mflux = 181.82522583
    mred = 3.47335
    mtheta = 31.732
    mdists = 7.78095
    mnpairs = 6.0

if 0:
    "median will npairs>=6"
    mflux = 167.97416687
    mred = 3.7075
    mtheta = 34.165
    mdists = 7.7093
    mnpairs = 9.0

s = "python stack.py -tests True -random True -overwrite False -randflipy False -makepdf False -reject False -makeim True -ntest 300 "

extraname = "meancorr"

if 1:
    os.system(s + ' &')
    os.system(s + "-nmin %d -extraname _nmin%d.%s.half &" % (mnpairs, mnpairs, extraname))
    os.system(s + "-nmax %d -extraname _nmax%d.%s.half &" % (mnpairs, mnpairs, extraname))
    os.system(s + "-pdmin %f -extraname .%s.half &" % (mtheta, extraname))
    os.system(s + "-pdmax %f -extraname .%s.half &" % (mtheta, extraname))
    os.system(s + "-dmin %f -extraname .%s.half &" % (mdists, extraname))
    os.system(s + "-dmax %f -extraname .%s.half &" % (mdists, extraname))
    os.system(s + "-rmin %f -extraname .%s.half &" % (mred, extraname))
    os.system(s + "-rmax %f -extraname .%s.half &" % (mred, extraname))
    # os.system(s + "-lmin %f -extraname .half &" % mflux)
    # os.system(s + "-lmax %f -extraname .half &" % mflux)

if 1:
    for i in np.arange(9) + 1:
        os.system(s + "-nmin %d -nmax %d -extraname _nmin%d_nmax%d.%s &" % (i, i + 2, i, i + 1, extraname))

if 0:
    os.system(s + "-meancorr True -extraname .meancorr")

if 0:

    dl = []
    dc = []
    l = []
    for i in range(len(zs)):
        from tools_sofi import comdist

        dc.append(comdist(zs[i]))
        dl.append((1 + zs[i]) * dc[i])
        l.append(4. * np.pi * dl[i] ** 2 * flux_lae1[i] * 9.521e-14)  # In units of 1e40 ergs/s

    mat = np.array([idlaes, dl, l])

    np.savetxt('%s/%s/cats/lae_lums.dat' % (folder, fitcat), mat.T, fmt='%f', header='ID dl Luminosity')
