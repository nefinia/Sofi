#!/usr/bin/env python
from copy import copy
from pyfits import getdata, PrimaryHDU
from scipy import integrate

import numpy as np


def mean_sclip(flux, stds=None, nsigma=4, iterations=6):
    if stds is None: stds = np.nanstd(flux)
    good = np.abs(flux) < nsigma * stds
    flux = flux[good]
    for i in range(iterations-1):
        stds = np.nanstd(flux)
        good = np.abs(flux) < nsigma * stds
        flux = flux[good]
    return np.nanmean(flux), np.nanstd(flux)

def l2pix(l, l0 = 4750., pixsize = 1.25):
    """Convert wavelength to pixels"""
    # MUSE: 1 pixel = 1.25 Angstrom, 1st pixel = 4750 Angstrom
    return (l - l0) / pixsize + 1

def oneovere(z, omega_l=0.714):
    "Routine to compute 1/E(z) for angular diameter distance calculation"
    omega_m = 1. - omega_l
    return np.power(np.sqrt(omega_m * np.power(1 + z, 3.) + omega_l), -1)

def angdist(z, h0=69.6, omega_l=0.714):
    "Routine to compute angular diameter distance"
    c = 2.99792458e5
    dh = c / h0
    dm = integrate.quad(oneovere, 0, z, args=omega_l)

    return dh * dm[0] / (1 + z)

cubefolder = '/net/astrogate/export/astrodata/gallegos/'
fitcat = 'mosaic'
zmin = 1400
zw = 10
asec2radian = np.pi/6.48e5
asec2pix = 5 # ~5 pixels per asec in MUSE

do_datacube = False
if do_datacube:
    cubename = 'DATACUBE_UDF-MOSAIC.fits'
    cube = getdata('%s/%s/%s' % (cubefolder, fitcat, cubename))
    ttt, header_data_cube = getdata(cubename, 0, header=True)
    header = copy(header_data_cube)
    header_data_cube['NAXIS3'] -= zmin
    crval3 = header_data_cube['CRVAL3']
    header_data_cube['CRVAL3'] += zmin
    hdu = PrimaryHDU()
    hdu.header = header_data_cube
    hdu.data =  cube[1400:,:,:]
    hdu.writeto('%s/%s/%s' % (cubefolder, fitcat, cubename.replace('.fits', '.zmin%d.fits' % zmin)))

cubename = 'UDF-mosaic.zmin1400.csub.fits'
cube = getdata('%s/%s/%s' % (cubefolder, fitcat, cubename))
zl, yl, xl = cube.shape
cdfs = True
catname = 'gals.fits'

halpha_wav = 6564.60986328125

lines = ['CIII1907', 'CIII1909', 'HALPHA', 'HBETA', 'HDELTA', 'HEII1640', 'HEPSILON', 'HGAMMA', 'LYALPHA',
         'MGII2796', 'MGII2803', 'NEIII3869', 'NII6548', 'NII6584', 'OII3726', 'OII3729', 'OIII4959',
         'OIII5007', 'SII6717', 'SII6731']
_wlines = [1906.68, 1908.73, 6564.61, 4862.68, 4102.89, 1640.42, 3971.2, 4341.68, 1215.67, 2796.35,
          2803.53, 3870.16, 6549.85, 6585.28, 3727.09, 3729.88, 4960.3, 5008.24, 6718.29, 6732.67]
wlines = {}
for i in range(len(lines)): wlines[lines[i]] = _wlines[i]

if 0:#os.path.isfile(catname):
    galcat = getdata('../../%s/cats/%s' % (fitcat, catname), 1)
    xs, ys, zs = [galcat['x'], galcat['y'], galcat['z']]
    wavelength = galcat['HALPA_LBDA_OBS']
    redshift = galcat['Z_MUSE']
    ang_dists = angdist(redshift)

else:
    import astropy.io.fits as fits
    from astropy import wcs
    galcat = getdata('../../%s/cats/mosaic_dr1.team.20170829.fits' % fitcat, 1)
    cdfscat = getdata('../../%s/cats/cdfs.fits' % fitcat, 1)
    ra = np.concatenate((galcat['RA'], cdfscat['RAJ2000']))
    dec = np.concatenate((galcat['DEC'], cdfscat['DEJ2000']))
    filename = '%s/%s/%s' % (cubefolder, fitcat, 'DATACUBE_UDF-MOSAIC.fits')
    ttt, header_data_cube = getdata(filename, 0, header=True)
    cards = header_data_cube.cards
    bad = ['COMMENT' == b[0] for b in cards]
    for i in range(np.sum(bad)): header_data_cube.remove('COMMENT')
    hdulist = fits.open(filename)
    w = wcs.WCS(header_data_cube, hdulist)
    tt = np.zeros(len(ra))
    xs, ys, zs = w.all_world2pix(ra, dec, tt, 1)
    redshift = np.concatenate((galcat['Z_MUSE'], cdfscat['Z']))
    #ang_dists = angdist(redshift)
    ang_dists = []
    for z in redshift:
        ang_dists.append(angdist(z))
    kpc2pix = np.array(ang_dists) * asec2radian * 1000 * asec2pix
    wavelength = {}
    for line in lines:
        wavelength[line] = np.concatenate((galcat['%s_LBDA_OBS' % line], wlines[line]*(1+cdfscat['Z'])))

radii_kpc = np.array([0, 2.5, 7.5, 12.5, 23.5, 40.5, 75.5, 138.5, 251.5, 434.5, 769.5, 1102.5, 1517.5])
radticks = ['2', '5', '10', '18', '32', '58', '107', '195', '343', '602', '936', '1310']

flux_sdss = [np.nan, 27.5, 25.6, 5.1, 1.39, 1.04, .703, .326, .208, .155, .129, .107]
conv = 176.715 #(1.5)**2*3.1416/.04 conversion factor from muse pixel area to sdss fiber area
flux_sdss = np.array(flux_sdss)/conv
nrad = len(radii_kpc)

wavmin, wavmax = [4750.+zmin*1.25, 4750.+(zmin+zl-1)*1.25]

foutline = '../../%s/lines.dat' % fitcat

if 0:
    #isfiles = True
    #isfiles &= os.path.isfile(foutline)
    f_mean, f_median, f_std = [{},{},{}]
    fmean
    for line in lines:
        fgals_mean[line] = '../../%s/%s.mean.fits' % (fitcat, line)
        fgals_median[line] = '../../%s/%s.median.fits' % (fitcat, line)
        fgals_std[line] = '../../%s/%s.std.fits' % (fitcat, line)
        #isfiles &= os.path.isfile(fgals_mean[line]) & os.path.isfile(fgals_median[line]) & \
        #          os.path.isfile(fgals_std[line])
        #if os.path.isfile(fgals_mean[line])

_y, _x = np.ogrid[0: yl, 0: xl]

flux_mean = {}
flux_median = {}
flux_std = {}
nz = 2*zw+1
lines = ['HALPHA', 'HBETA', 'HGAMMA', 'HDELTA']
for line in lines:#[lines[2]]:
    print 'Computing radial profile for %s,' % line,
    wav = wavelength[line]
    goodgals = (wavmin<wav) & (wav<wavmax)
    print np.sum(goodgals), 'galaxies including',
    wav = wav[goodgals]
    x = xs[goodgals]
    y = ys[goodgals]
    z = np.round(l2pix(wav)).astype(int)-zmin-1
    N = len(x)
    print np.sum((x > 0) & (x < xl) & (y > 0) & (y < yl)), 'inside UDF-mosaic'

    fmean = np.zeros((nz, nrad-1))
    fmedian = np.zeros((nz, nrad-1))
    fstd = np.zeros((nz, nrad-1))

    for r in range(nrad - 1):
        fall = [[]]*nz
        for zz in range(nz):
            good = {}
            rgood = []
            print line,r,zz
            for j in range(N):
                #print j
                zzz = z[j]-zw+zz
                if (zzz>-1) and (zzz<zl):
                    radii = radii_kpc * kpc2pix[j]
                    rgood = ((x[j]-_x)**2+(y[j]-_y)**2 >= radii[r]) &\
                        ((x[j]-_x)**2+(y[j]-_y)**2 <= radii[r+1])
                    if len(rgood)>0:
                        if zzz in good.keys():
                            good[zzz] |= rgood
                        else: good[zzz] = rgood
                    else:
                        print 'nooo'
            for k in good.keys():
                fall[zz] = np.concatenate((fall[zz], cube[k, good[k]]))
        for zz in range(nz):
            fmean[zz,r], fstd[zz,r] = mean_sclip(fall[zz])
            fmedian[zz,r] = np.nanmedian(fall[zz])

    hdu = PrimaryHDU()
    hdu.data = fmean
    hdu.writeto('../../%s/%s.mean.fits' % (fitcat, line), clobber=True)
    hdu.data = fmedian
    hdu.writeto('../../%s/%s.median.fits' % (fitcat, line), clobber=True)
    hdu.data = fstd
    hdu.writeto('../../%s/%s.std.fits' % (fitcat, line), clobber=True)
    flux_mean[line] = fmean[nz/2,:]
    flux_median[line] = fmedian[nz/2,:]
    flux_std[line] = fstd[nz/2,:]

fout = open(foutline, 'w')
fout.write('# radius HALPHA_SDSS')
for line in lines: fout.write(' %s_mean %s_median %s_std' % (line, line, line))
fout.write('\n')
for i in range(nrad-1):
    fout.write('%s %f' % (radticks[i], flux_sdss[i]))
    for line in lines: fout.write(' %f %f %f' % (flux_mean[line][i], flux_median[line][i], flux_std[line][i]))
    fout.write('\n')
fout.close()

#get some spectra!
#convert to SB
#how much of this is from the continuum??