__author__ = 'gallegos'

from pyfits import getdata, PrimaryHDU
import numpy as np
import os
from tools_sofi import astroim, rdarg, stack, pdfim, analysis, angdist
from sys import argv
from random import randint
import matplotlib.mlab as mlab
from scipy.stats import norm
import matplotlib.pyplot as plt


def histplot(datos, drand, dmix, hname, title=None, gaussian=None):
    n, bins, patches = plt.hist(dmix, 16, normed=1, facecolor='blue', alpha=0.75)
    plt.xlabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=16)
    plt.savefig(hname.replace(".png", "_full.png"))
    plt.close()
    n, bins, patches = plt.hist(drand, 16, normed=1, facecolor='grey', alpha=0.75, label='Random Subsets')
    bins = np.append(bins, (np.arange(5) + 1) * abs(bins[1] - bins[0]) + bins[-1])
    print bins
    n, bins, patches = plt.hist(datos, bins, normed=1, facecolor='green', alpha=0.75, label='Oriented Subsets')
    if gaussian:
        (mu, sigma) = norm.fit(datos)
        (mu2, sigma2) = norm.fit(drand)
        # y = mlab.normpdf(bins, mu, sigma)
        y2 = mlab.normpdf(bins, mu2, sigma2)
        l = plt.plot(bins, y, 'r--', linewidth=2)
        l2 = plt.plot(bins, y2, 'r--', linewidth=2)
    if title:
        plt.title(title)
    plt.legend(prop={'size': 10})
    print hname
    plt.savefig(hname)
    plt.close()



halfplot = False
envplot = False
ksplot = False
radprof = False
zneig = True

fitcat = rdarg(argv, key='fitcat', type=str, default='all')  # 'all'#'HDFS'
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
variables = ['Flux', 'Environment', 'Redshift', "Com_Distance", "Proj_Distance"]
init = 'not-normalized_masked_cubical_'
n = len(variables)
medians = [199, 6, 3.5, 7, 31]
mins = [0, 1, 2.9, 0, 16]
maxs = [99999, 11, 4.0, 20, 80]
name = 'stack.avs.dat'

if halfplot:
    # max flux around 1700 (check units...1e-20?)
    halves = [[]] * 2

    f = [[[]] * n] * 2
    fout = open('%s/%s/cats/halfs.dat' % (folder, fitcat), 'w')
    fout.write('#id name half median SB_stack SB_mean SB_err\n')
    nmin = np.zeros((2, n))
    nmax = np.zeros((2, n))
    mu = np.zeros((2, n))
    mu2 = np.zeros((2, n))
    sigma = np.zeros((2, n))

    for i in np.arange(n - 1) + 1:
        for j in range(n):
            if j == i:
                nmin[0][j] = medians[j]
                nmax[0][j] = maxs[j]
                nmin[1][j] = mins[j]
                nmax[1][j] = medians[j]
                print nmin[0][j], medians[j], nmax[0][j], maxs[j], nmin[1][j], mins[j], nmax[0][j], medians[j], j, i
            else:
                nmin[1][j] = mins[j]
                nmax[1][j] = maxs[j]
                nmin[0][j] = mins[j]
                nmax[0][j] = maxs[j]
                print nmin[0][j], nmax[0][j], nmin[1][j], nmax[0][j], j, i

        foldernames = []
        s = ['min', 'max']
        s2 = ['>', '<']
        s3 = ['_ht_', '_lt_']

        nn = [medians[1]] * 2

        for k in range(2):

            if i == 1:
                ss = '%s/%s/stacks/not-normalized_masked_cubical_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d-%d_n%s%d.meancorr.half/%s' % (
                    folder, fitcat, nmin[k][3], nmax[k][3], nmin[k][4], nmax[k][4], nmin[k][2], nmax[k][2],
                    nmin[k][0], nmax[k][0], s[k], nn[k], name)
            else:
                ss = '%s/%s/stacks/not-normalized_masked_cubical_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d-%d.meancorr.half/%s' % (
                    folder, fitcat, nmin[k][3], nmax[k][3], nmin[k][4], nmax[k][4], nmin[k][2], nmax[k][2],
                    nmin[k][0], nmax[k][0], name)

            print ss
            ff = np.loadtxt(ss).transpose()
            fmix = []
            for ii in ff[2]:
                for jj in ff[1]:
                    fmix.append(ii-jj)
            sigma[k][i] = np.nanstd(fmix)
            mu[k][i] = np.nanmean(fmix)
            mu2[k][i] = ff[0][0]
            ss2 = '%s/%s/plots/%s' % (folder, fitcat, '%s%s%d.png' % (variables[i], s3[k], medians[i]))
            histplot(ff[2], ff[1], fmix, ss2, r'$%s%s%.1f$' % (variables[i], s2[k], medians[i]))
            sss = '%d %s %s %.1f %f %f %f\n' % (
            i, variables[i], s2[k], medians[i], mu2[k][i], mu[k][i], sigma[k][i])
            print sss
            fout.write(sss)

    ss = '%s/%s/stacks/not-normalized_masked_cubical_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d-%d/%s' % (
        folder, fitcat, mins[3], maxs[3], mins[4], maxs[4], mins[2], maxs[2],
        mins[0], maxs[0], name)
    print ss
    ff = np.loadtxt(ss).transpose()
    for ii in ff[2]:
        for jj in ff[1]:
            fmix.append(ii-jj)
    fsigma = np.nanstd(fmix)
    fmu = np.nanmean(fmix)
    fmu2 = ff[0][0]
    ss2 = '%s/%s/plots/%s' % (folder, fitcat, name.replace('.dat', '.png'))
    histplot(ff[2], ff[1], fmix, ss2)
    sss = '%d full - -1 %f %f %f\n' % (n, fmu2, fmu, fsigma)
    print sss
    fout.write(sss)

    fout.close()
    plt.errorbar([n+.5], [fmu], yerr=[fsigma], fmt='o', color='green')
    plt.scatter([n+.5], [fmu], color='green')
    labels = [r'$L$', r'$n_{\rm{neighbors}}$', r'$z$', r"$d_{\rm c}$", r"$d_{\rm{projected}}$", "Full"]
    x = np.arange(n+1) + .5
    plt.xticks(x, labels, rotation=0)
    x = np.arange(n) + .5
    y1 = mu[0]
    y2 = mu[1]
    yerr1 = sigma[0]
    yerr2 = sigma[1]
    plt.errorbar(x, y1, yerr=yerr1, fmt='o', color='red')
    plt.scatter(x, y1, label=r'$f<\~f}$', color='red')
    plt.errorbar(x, y2, yerr=yerr2, fmt='o', color='blue')
    plt.scatter(x, y2, label=r'$f>\~f}$', color='blue')
    plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=16)
    plt.legend()
    plt.savefig('%s/%s/plots/halfs.png' % (folder, fitcat))
    plt.close()


#####################
##Environment plot
####################
if envplot:
    n = 8
    x = np.arange(n) + 2
    randmean = np.zeros(n)
    sigma = np.zeros(n)
    mu = np.zeros(n)
    mu2 = np.zeros(n)

    for i in range(n):
        ss = '%s/%s/stacks/not-normalized_masked_cubical_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d-%d_nmin%d_nmax%d.meancorr/%s' % (
            folder, fitcat, mins[3], maxs[3], mins[4], maxs[4], mins[2], maxs[2],
            mins[0], maxs[0], i + 1, i + 2, name)
        ff = np.loadtxt(ss).transpose()
        fmix = []
        for ii in ff[2]:
            for jj in ff[1]:
                fmix.append(ii-jj)
        sigma[i] = np.nanstd(fmix)
        mu[i] = np.nanmean(fmix)
        mu2[i] = ff[0][0]
        ss2 = '%s/%s/plots/%s' % (folder, fitcat, '%s%d-%d.png' % (variables[1], i + 1, i + 2))
        histplot(ff[2], ff[1], fmix, ss2)#, r'$%d \leq \rm{N_{neighbors}} \leq%d$' % (i + 1, i + 2))
    y = mu - randmean
    yerr = sigma
    plt.errorbar(x, y, yerr=yerr, xerr=np.zeros(len(x)) + .5, fmt='o', color='red')
    plt.scatter(x, y, color='red')
    plt.xticks(np.arange(n + 2) + 1, rotation=0)
    plt.ylabel(r'$\rm{SB}\,\rm{[}10^{-20}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=16)
    plt.xlabel('# of neighbours', fontsize=16)
    plt.legend()
    plt.savefig('%s/%s/plots/environment.png' % (folder, fitcat))
    plt.close()
    # plt.show()


##### K-S test between random and oriented
#Cumulative histogram for both
if ksplot:
    ff = '%s/%s/stacks/not-normalized_masked_cubical_d%d-%d_pd%d-%d_z%.1f-%.1f_f%d-%d/' % (
        folder, fitcat, mins[3], maxs[3], mins[4], maxs[4], mins[2], maxs[2],
        mins[0], maxs[0])




if radprof:
    folder = "../../all/stacks/not-normalized_masked_cubical_d0-20_pd16-80_z2.9-4.0_f0-99999.nosclip/"
    # folder='not-normalized_masked_cubical_d0-20_pd16-80_z2.9-4.0_f0_nmin6.half/'

    sbdat = folder + "SBprof.dat"

    if os.path.isfile(sbdat):
        print sbdat, "already exists"
    else:
        fits = getdata(folder + "random_stack.fits")
        sets = []
        nsets = 100
        for i in range(nsets):
            sets.append(getdata(folder + 'random_stack.%d.fits' % i))
        sets = np.array(sets)

        zl, yl, xl = fits.shape
        s = 4
        s2 = float(s ** 2)
        fbig = np.zeros((zl, yl * s, xl * s))
        sbig = np.zeros((nsets, zl, yl * s, xl * s))

        for i in range(xl):
            print i + 1, "of", xl
            for j in range(yl):
                for k in range(zl):
                    fbig[k, j * s:s * (j + 1), i * s:s * (i + 1)] = fits[k, j, i] / s2
                    for n in range(nsets):
                        # print "set", n
                        sbig[n, k, j * s:s * (j + 1), i * s:s * (i + 1)] = sets[n, k, j, i] / s2

        y, x = np.ogrid[-yl * s / 2 + 1: yl * s / 2 + 1, -xl * s / 2 + 1: xl * s / 2 + 1]

        zw0 = 1
        zmean = 3.5

        d = angdist(zmean)
        conv0 = np.pi * d * 1000 / (180 * 3600.)
        n0 = 1
        zw0 = 1

        zcen = range(zl / 2 - zw0, zl / 2 + zw0 + 1)
        zblue = range(zl / 2 - 5)
        zred = range(zl / 2 - 5, zl)

        fcen = fbig[zcen, :]
        fblue = fbig[:zl / 2 - 5, :]
        fred = fbig[zl / 2 + 6:, :]

        fout = open(sbdat, "w")
        fout.write(
            "#rad_pix rad_ang rad_phys SB_rand std_rand SB_diff SB_sky std_sky SB_blue std_blue SB_red std_red\n")

        for rad in np.arange(30) * s:
            n = n0 * (2 * rad + 1)
            conv = 1 / (np.pi * ((rad / s + 1) ** 2 - (rad / s) ** 2) * .4)
            cool = (x ** 2 + y ** 2 >= rad ** 2) & (x ** 2 + y ** 2 < (rad + s) ** 2)
            sum = np.nansum(fcen[:, cool])
            N = np.sum(np.isfinite(fcen[:, cool]))
            std = np.nanstd(sbig[:, zl / 2 - zw0: zl / 2 + zw0 + 1, cool]) * N
            Nred = np.sum(np.isfinite(fred[:, cool]))
            mred = np.nansum(fred[:, cool])
            stdred = np.nanstd(sbig[:, zl / 2 + 6:, cool]) * N
            Nblue = np.sum(np.isfinite(fblue[:, cool]))
            mblue = np.nansum(fblue[:, cool])
            stdblue = np.nanstd(sbig[:, :zl / 2 - 5, cool]) * N
            msky = (mred + mblue) / 2.
            stdsky = np.sqrt(stdred ** 2 + stdblue ** 2)
            print rad / s + 1, sum, sum * conv, std
            fout.write("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n" % (
            rad / s + 1, (rad / s + 1) * .4, (rad / s + 1) * .4 * conv0,
            sum * conv, std * conv, (sum - msky) * conv, msky * conv, stdsky * conv,
            mblue * conv, stdblue * conv, mred * conv, stdred * conv))

        fout.close()

    zmean = 3.5

    d = angdist(zmean)
    conv0 = np.pi * d * 1000 / (180 * 3600.)



    a = np.loadtxt(sbdat).T
    x = a[1]

    y = a[3] * 1e-20
    sky = a[6] * 1e-20
    yerr = a[4] * 1e-20
    yerrb = yerr.copy()
    yerrb[(y - yerr) < 0] = y[(y - yerr) < 0] - 1e-22
    plt.plot(x, y)
    # plt.plot(x, x-x+1e-20)
    plt.xlim(0.5, 5)
    plt.errorbar(x, y, yerr=[yerrb, yerr], fmt='o', color='red')
    plt.scatter(x, y, label=r'$f<\~f}$', color='red')
    plt.ylabel(r'$\rm{SB}\,\rm{[}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$', fontsize=16)
    plt.xlabel(r'$\theta\,\rm{[arcsec]}$', fontsize=16)
    plt.loglog()
    plt.xticks([1, 2, 3, 4], ['1', '2', '3', '4'])
    plt.twiny()
    plt.xlabel(r'$\theta\,\rm{[kpc]}$', fontsize=16)
    plt.loglog()
    plt.xlim(0.5, 5)
    plt.xticks([1, 2, 3, 4], ['7.5', '15', '22.5', '30'])
    plt.savefig('../../all/plots/SB_prof.png')
    plt.show()
    plt.close()
    # plt.show()

    z_wish = 3

    zLutz, rhLutz, r19Lutz = np.loadtxt('lutz.dat', usecols=(1, 3, 6), unpack=1)
    sbCoeffLutz = ((1. + zLutz) / (1. + z_wish) / 4.) ** 4  # to scale / 4.
    cMpcCoeffLutz = (1. + zLutz) / (1. + z_wish)

    ####plot
    from math import exp


    def lyaW2016(r, rh, r19):
        return 1e19 * exp(-r19 / rh) * np.exp(-r / rh)


    imin = 20
    imax = 0
    imax2 = 3
    lineW = 3
    sbLutzmean = np.zeros((len(rhLutz), len(np.arange(0., 170., 1))))
    sbLutzmeanSB = np.zeros((len(rhLutz), len(np.arange(0., 170., 1))))

    for i in range(len(rhLutz)):
        b = np.arange(0., 170., 1)
        sb = lyaW2016(b, rhLutz[i], r19Lutz[i])
        sbLutzmean[i] = sb
        sbLutzmeanSB[i] = sb * sbCoeffLutz[i]
    plt.subplot()

    plt.plot(b, sbLutzmean.mean(axis=0), c='dodgerblue', lw=lineW, label='Wisotzki+16 (LAEs)')
    plt.fill_between(b, sbLutzmean[imin], sbLutzmean[[imax, imax2]].max(axis=0), facecolor='dodgerblue', alpha=0.3, lw=0,
                     edgecolor='none')


    plt.plot(b, sbLutzmeanSB.mean(axis=0), c='dodgerblue', lw=lineW )
    plt.fill_between(b,sbLutzmeanSB[imin], sbLutzmeanSB[[imax,imax2]].max(axis=0), facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')
    plt.loglog()
    plt.plot(b * cMpcCoeffLutz.mean(), sbLutzmeanSB.mean(axis=0), c='dodgerblue', lw=lineW) #, label='Wisotzki+16 (LAEs)')
    plt.fill_between(b * cMpcCoeffLutz.mean(),sbLutzmeanSB[imin], sbLutzmeanSB[[imax,imax2]].max(axis=0),
                     facecolor='dodgerblue', alpha=0.3, lw=0, edgecolor='none')
    plt.show()


if zneig:
    cat = '../../all/cats/laes.fits'
    data = getdata(cat, 1)
    npair = data['npairs']
    good = np.isfinite(npair)
    npair = npair[good]
    z = data['Z'][good]
    am = []
    mm = []
    t = []
    b = []
    bins = []
    bin = .3
    for i in np.arange(2.9, 4.1, bin):
        p = npair[(z > i) & (z <= i + bin)]
        mm.append(np.median(p))
        t.append(np.percentile(p, 75))
        b.append(np.percentile(p, 25))
        bins.append(i + bin / 2.)
        print mm, bin / 2.
    am = np.array([bins, mm, b, t]).T
    yerr = [np.array(mm) - np.array(b), np.array(t) - np.array(mm)]
    np.savetxt('../../all/cats/nvsz.dat', am, header='bins median n_25 n_75')
    plt.scatter(bins, mm)
    plt.errorbar(bins, mm, yerr=yerr, fmt='o', color='red')
    plt.axvline(x=3.5, linestyle='--')
    plt.xlabel("z")
    plt.ylabel(r"$n_{\rm{neighbors}}$", fontsize=18)
    plt.savefig("../../AAAPaper/zn.png")
    plt.show()