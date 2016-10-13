#!/usr/bin/env python
from scipy import integrate
import numpy as np
from math import *
from pyfits import getdata, PrimaryHDU
import h5py
import numpy as np
import os, sys
from math import sqrt
import astropy as ap
from astropy.io import fits
from pyfits import getdata, PrimaryHDU, getheader
# from pylab import *
import scipy.ndimage as ndimage


def deg2hms(ra, dec):
    ra_h = int(ra / 15)
    ra_m = int((ra / 15. - ra_h) * 60.)

    ra_s = ((ra / 15. - ra_h) * 60. - ra_m) * 60.
    dec_d = int(dec)
    dec_m = abs(int((dec - dec_d) * 60.))
    dec_s = abs((abs((dec - dec_d) * 60.) - dec_m) * 60.)

    ra_hms = [ra_h, ra_m, ra_s]
    dec_dms = [dec_d, dec_m, dec_s]

    print 'ra', ra, '->', ra_hms, 'dec', dec, '->', dec_dms
    return ra_hms, dec_dms


def lum_func(imag, mstar, limmag, alpha=-1.):
    "Routine to compute luminosity filter weight (see Eqn 10 in Rykoff et al. 2012)"

    ## integrate
    n = integrate.quad(schechter, 10, limmag, args=(alpha, mstar))

    wts = schechter(imag, alpha, mstar) / n[0]

    return wts


def schechter(m, alpha, mstar):
    "Schechter luminosity function"

    s1 = 10. ** (-0.4 * (m - mstar) * (alpha + 1))
    s2 = -10. ** (-0.4 * (m - mstar))

    return s1 * np.exp(s2)


def oneovere(z, omega_l=0.714):
    "Routine to compute 1/E(z) for angular diameter distance calculation"

    omega_m = 1. - omega_l

    return 1. / sqrt(omega_m * pow(1 + z, 3.) + omega_l)


def angdist(z, h0=69.6, omega_l=0.714):
    "Routine to compute angular diameter distance"

    c = 2.99792458e5
    dh = c / h0

    dm = integrate.quad(oneovere, 0, z, args=omega_l)

    return dh * dm[0] / (1 + z)


def comdist(z, h0=69.6, omega_l=0.714):
    "Routine to compute comoving distance"

    c = 2.99792458e5
    dh = c / h0

    dm = integrate.quad(oneovere, 0, z, args=omega_l)

    return dh * dm[0]


def distance(z1, ra1, dec1, z2, ra2, dec2):
    "Routine to compute several cosmological distances between 2 objects"

    # pi dist will be comoving!

    zdif = abs(z1 - z2)
    zmean = (z1 + z2) / 2.
    zmin = min(z1, z2)
    cd1 = comdist(z1)
    cd2 = comdist(z2)
    cdmean = (cd1 + cd2) / 2.
    cd = min(cd1, cd2)
    pi_v = zdif * 2.99792458e5
    pi_Mpc = abs(cd1 - cd2)  # / (1 + zmean)
    theta = sqrt(
        ((ra1 - ra2) * cos((dec1 + dec2) / 2. * pi / 180.)) ** 2 + (dec1 - dec2) ** 2)  # ang distance in degrees
    proj_dist = 1000. * theta * pi / 180. * cd / (1. + zmin)
    # phys_dist = sqrt(pi_Mpc**2+proj_dist**2) #sqrt((pi / 180. * theta * cd / (1 + minz)) ** 2 + cd ** 2)

    return theta * 3600., cd1, cd2, pi_Mpc, pi_v, proj_dist


def cic(value, x, nx, y=None, ny=1, weighting=True, wraparound=False):
    """ Interpolate an irregularly sampled field using Cloud in Cell
    method.

    This function interpolates an irregularly sampled field to a
    regular grid using Cloud In Cell (nearest grid point gets weight
    1-dngp, point on other side gets weight dngp, where dngp is the
    distance to the nearest grid point in units of the cell size).

    Inputs
    ------
    value: array, shape (N,)
        Sample weights (field values). For a temperature field this
        would be the temperature and the keyword average should be
        True. For a density field this could be either the particle
        mass (average should be False) or the density (average should
        be True).
    x: array, shape (N,)
        X coordinates of field samples, unit indices: [0,NX>.
    nx: int
        Number of grid points in X-direction.
    y: array, shape (N,), optional
        Y coordinates of field samples, unit indices: [0,NY>.
    ny: int, optional
        Number of grid points in Y-direction.
    wraparound: bool (False)
        If True, then values past the first or last grid point can
        wrap around and contribute to the grid point on the opposite
        side (see the Notes section below).

    Returns
    -------
    dens: ndarray, shape (nx, ny)
        The grid point values.

    Notes
    -----
    Example of default allocation of nearest grid points: nx = 4, * = gridpoint.

      0   1   2   3     Index of gridpoints
      *   *   *   *     Grid points
    |---|---|---|---|   Range allocated to gridpoints ([0.0,1.0> -> 0, etc.)
    0   1   2   3   4   posx

    Example of ngp allocation for wraparound=True: nx = 4, * = gridpoint.

      0   1   2   3        Index of gridpoints
      *   *   *   *        Grid points
    |---|---|---|---|--    Range allocated to gridpoints ([0.5,1.5> -> 1, etc.)
      0   1   2   3   4=0  posx


    References
    ----------
    R.W. Hockney and J.W. Eastwood, Computer Simulations Using Particles
        (New York: McGraw-Hill, 1981).

    Modification History
    --------------------
    IDL code written by Joop Schaye, Feb 1999.
    Avoid integer overflow for large dimensions P.Riley/W.Landsman Dec. 1999
    Translated to Python by Neil Crighton, July 2009.

    Examples
    --------
    #>>> nx = 20
    #>>> ny = 10
    #>>> posx = np.random.rand(size=1000)
    #>>> posy = np.random.rand(size=1000)
    #>>> value = posx**2 + posy**2
    #>>> field = cic(value, posx*nx, nx, posy*ny, ny)
    # plot surface
    """

    def findweights(pos, ngrid, ww=weighting):
        """ Calculate CIC weights.

        Coordinates of nearest grid point (ngp) to each value. """

        if wraparound:
            # grid points at integer values
            ngp = np.fix(pos + 0.5)
        else:
            # grid points are at half-integer values, starting at 0.5,
            # ending at len(grid) - 0.5
            ngp = np.fix(pos) + 0.5

        # Distance from sample to ngp.
        distngp = ngp - pos

        if ww:
            # weight for higher (right, w2) and lower (left, w1) ngp
            weight2 = np.abs(distngp)
            weight1 = 1.0 - weight2

        # indices of the nearest grid points
        if wraparound:
            ind1 = ngp
        else:
            ind1 = ngp - 0.5
        ind1 = ind1.astype(int)

        # print 'ind',ind1,'max min ind',max(ind1),min(ind1)

        ind2 = ind1 - 1

        # Correct points where ngp < pos (ngp to the left).
        ind2[distngp < 0] += 2

        # Note that ind2 can be both -1 and ngrid at this point,
        # regardless of wraparound. This is because distngp can be
        # exactly zero.
        bad = (ind2 == -1)
        ind2[bad] = ngrid - 1
        if not wraparound and ww:
            weight2[bad] = 0.
        bad = (ind2 == ngrid)
        ind2[bad] = 0
        bad1 = (ind1 == ngrid)
        ind1[bad1] = ngrid - 1

        if not wraparound and ww:
            weight2[bad] = 0.

        if wraparound:
            ind1[ind1 == ngrid] = 0

        if not ww:
            weight1 = ind1 - ind1
            weight2 = weight1

        return dict(weight=weight1, ind=ind1), dict(weight=weight2, ind=ind2)

    def update_field_vals(field, weight, count, value, ww=weighting, a=None, b=None, debug=False):
        """ This updates the field array (and the totweight array if
        average is True).

        The elements to update and their values are inferred from
        a,b,c and value.
        """
        # weight per coordinate
        if ww:
            weights = a['weight'] * b['weight']
            # Don't modify the input value array, just rebind the name.
            value = weights * value
            indices = []

        for i in range(len(value)):
            field[a['ind'][i]][b['ind'][i]] += value[i]
            weight[a['ind'][i]][b['ind'][i]] += weights[i]
            count[a['ind'][i]][b['ind'][i]] += 1

        if debug: print i, weights[i], value[i], field[a['ind'][i]][b['ind'][i]]

    nx, ny = (int(i) for i in (nx, ny))
    value = np.asarray(value)

    x1 = None
    x2 = None
    y1 = None
    y2 = None

    x1, x2 = findweights(np.asarray(x), nx)
    ind = []
    ind.append([x1, x2])
    if y != None:
        y1, y2 = findweights(np.asarray(y), ny)
        ind.append([y1, y2])

    # float32 to save memory for big arrays (e.g. 256**3)
    field = np.zeros(shape=(nx, ny))  # field = np.zeros(shape=(nx,ny,nz)).squeeze()
    weight = np.zeros(shape=(nx, ny))  # field = np.zeros(shape=(nx,ny,nz)).squeeze()
    count = np.zeros(shape=(nx, ny))

    update_field_vals(field, weight, count, value, weighting, x1, y1)
    update_field_vals(field, weight, count, value, weighting, x2, y1)
    update_field_vals(field, weight, count, value, weighting, x1, y2)
    update_field_vals(field, weight, count, value, weighting, x2, y2)

    return field, weight, count


def cubex(name, zw='*', imtype='flux', varname='', snrname='', imsnrname=''):
    if zw == '*':
        print 'Warning: selecting all wavelegth range...'

    if varname == '':
        varname = name.replace('.fits', '.VAR.fits')
    if snrname == '':
        snrname = name.replace('.fits', '.SNR.fits')
    if snrname == '':
        imsnrname = name.replace('.fits', '.SNR.IM.fits')

    s = 'CubEx-1.5 -cube %s -MultiExt .false. -ApplyFilter .true. -ApplyFilterVar .true. -FilterXYRad 1 -SN_Threshold 2 -MinNSpax 2' % (
        name, varname, snrname, imtype)
    print s
    os.system(s)
    s = 'Cube2Im -cube %s[*,*,%s] -varcube %s[*,*,zw] -snrmap %s -imtype %s' % (
        name, zw, varname, zw, imsnrname, imtype)
    print s
    os.system(s)


def makejpeg(im, smooth=True, imout='', cont=False, scalelims=''):
    if imout == '':
        if im.find('.fits') == -1:
            print 'Standard .fits extension not found. Is it a fits file?'
            imout = im + '.jpeg'
        else:
            imout = im.replace('.fits', '.jpeg')
    # -scale limits %s
    sscont = ''
    if cont: sscont = '-contour nlevels 10 -contour limits .6 6 -contour color black'

    s = 'ds9 -zscale %s -cmap b -nan red -smooth %s %s %s -zoom to fit -saveimage jpeg %s 100 -exit' % (
        scalelims, smooth, im, sscont, imout)
    print s
    os.system(s)


def astroim(im, smooth=True, xmin=None, xmax=None, ymin=None, ymax=None,
            vmin=None, vmax=None, xtmin=0, ytmin=None, contours=False, std=None,
            fsize=18, cbfrac=.047, pad=.03, dfig=None, cb_label=True,
            scb_label=r'Flux [$10^{-20}\,\rm{erg/s/cm^2}$]', saveim=False, clabelpad=None, ylabelpad=None,
            sbfix=False, units='theta', imout=None, show=True, title='', nsigma=6,
            x0=None, y0=None, res=600, gray=False, interpolation='gaussian', highq=False):
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, ListedColormap, BoundaryNorm

    hdu_list = fits.open(im)
    image_data = (hdu_list[0].data)[::-1, :]
    mask = np.isnan(image_data)
    image_data[mask] = 0  # -999
    image_data[image_data == -999] = 0  # -999

    if smooth:
        img = ndimage.gaussian_filter(image_data, sigma=1, order=0)
    else:
        img = image_data

    if sbfix:
        img /= 100.
    if xmin is None:
        xmin = 0
        xxmin=0
    else:
        xxmin=image_data.shape[1]/2+xmin
    if xmax is None:
        xmax = image_data.shape[1]
        xxmax=xmax
    else:
        xxmax = image_data.shape[1]/2+xmax
    if ymin is None:
        ymin = 0
        yymin = 0
    else:
        yymin = image_data.shape[0]/2+ymin
    if ymax is None:
        ymax = image_data.shape[0]
        yymax = ymax
    else:
        yymax = image_data.shape[0]/2+ymax

    fig, ax = plt.subplots(figsize=dfig)
    tn = 10
    x = xmax - tn * np.arange(xmax / tn + 1) + xmin - .5
    if xtmin != 0:
        dxt = xtmin - x[np.where([min(abs(x - xtmin)) == abs(i - xtmin) for i in x])[0]]
        x += dxt
    y = ymax - tn * np.arange(ymax / tn + 1) + ymin - .5

    xticks = ((xmax / 2 - x + xmin + .5) * .4).astype(int)[::-1]
    yticks = ((ymax / 2 - y + ymin + .5) * .4).astype(int)
    y = ymax - tn * np.arange(ymax / tn + 1) + ymin
    plt.xlim(xmin, xmax)
    plt.ylim(ymax, ymin)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.xaxis.label.set_fontsize(fsize)
    ax.yaxis.label.set_fontsize(fsize)
    plt.xticks(x, xticks, fontsize=fsize)
    plt.yticks(y, yticks, fontsize=fsize)
    plt.title(title, y=1.14)
    plt.xlabel(r'$\rm{\theta\,[arcsec]}$', fontsize=fsize)
    plt.ylabel(r'$\rm{\theta\,[arcsec]}$', fontsize=fsize, labelpad=ylabelpad)

    if x0 is not None and y0 is not None:
        plt.scatter([x0], [y0], marker='x')

    if not std:
        rad = 8
        y, x = np.ogrid[-ymax/2+1: ymax/2+1, -xmax/2+1: xmax/2+1]
        mask = x**2 + y**2 > rad**2
        std = np.nanstd(img[mask])
        print 'std with sigma clipping', std
    if vmin is None:
        vmin = -nsigma * std
    if vmax is None:
        vmax = nsigma * std
    if contours:
        levels = np.concatenate(((np.arange(nsigma)+1)*std, (np.arange(nsigma*2)**2+nsigma+1)*std))
        colors = ['white', 'lightgray', 'beige', 'yellow', 'orange', 'red', 'darkred']
        cset = plt.contour(img[yymin: yymax+1, xxmin: xxmax+1], levels, aspect='auto', colors='grey', extent=[xmin, xmax, ymin, ymax])
    print nsigma, std, vmin, vmax
    bounds = ((np.arange(nsigma*2)+1)*std-vmin)/(vmax-vmin)
    bounds[bounds>1] = 1
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0, 'white'),
                                                        (bounds[0], 'lightgrey'),
                                                        (bounds[1], 'beige'),
                                                        (bounds[2], 'yellow'),
                                                        (bounds[3], 'orange'),
                                                        (bounds[4], 'red'),
                                                        (bounds[5], 'darkred'),
                                                        (1, 'black')])
    if gray:
        cmap='gray_r'

    if 1: cmap = 'OrRd'
    plt.imshow(img[yymin: yymax+1, xxmin: xxmax+1], cmap=cmap, vmin=vmin, vmax=vmax, extent=[xmin, xmax, ymax, ymin], interpolation=interpolation)
    cbar = plt.colorbar(orientation='horizontal', pad=pad, fraction=cbfrac)#, boundaries=[levels[0], levels[-1]])
    tick_loc = ticker.MaxNLocator(nbins=7)
    cbar.locator = tick_loc
    cbar.update_ticks()
    cbar.ax.tick_params(labelsize=fsize)
    if cb_label and scb_label == '':
        if sbfix:
            aa = 18
        else:
            aa = 20
        scb_label = r'$\rm{SB}\,\rm{[}10^{-%s}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$' % aa
    cbar.set_label(scb_label, fontsize=(fsize - 1), labelpad=clabelpad)
    if saveim:
        if highq:
            if imout == None: imout = im.replace('.fits', '.eps')
            print 'Saving image', imout
            plt.savefig(imout, format='eps', dpi=1000)
        else:
            if imout == None: imout = im.replace('.fits', '.png')
            print 'Saving image', imout
            plt.savefig(imout)
    if show:
        plt.show()

    plt.close()


def pdfim(lst, smooth=True, xmin=0, xmax=41, ymin=0, ymax=41,
          vmin=None, vmax=None, xtmin=0, contours=False, std=None,
          fsize=18, cbfrac=.05, pad=.7, dfig=None, cb_label=True,
          scb_label='', clabelpad=None, ylabelpad=None,
          sbfix=False, units='theta', imout='test.pdf', title='', nsigma=5, parallel=True):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    cmap='OrRd'
    n = len(lst)
    ncol = 4
    nrow = 6
    jj = range(ncol)
    kk = range(nrow)
    i = 0
    nfig = 0

    with PdfPages(imout) as pdf:
        while i < n:

            f, ax = plt.subplots(ncol, nrow, sharex='col', sharey='row')
            nfig += 1
            for j in jj:
                for k in kk:
                    a = ax[j, k]
                    a.set_xticks([])
                    a.set_yticks([])
                    if i < n:
                        im = lst[i]
                        hdu_list = fits.open(im)
                        image_data = (hdu_list[0].data)[::-1, :]
                        mask = np.isnan(image_data)
                        image_data[mask] = 0
                        if smooth:
                            img = ndimage.gaussian_filter(image_data, sigma=1, order=0)
                        else:
                            img = image_data

                        if xmin is None:
                            xmin = 0
                            xmax = img.shape[1]
                            ymin = 0
                            ymax = img.shape[0]
                        a.set_xlim(xmin, xmax)
                        a.set_ylim(ymax, ymin)
                        a.xaxis.tick_top()
                        axtitle = im.split('_pair')[1].split('.')[0]
                        a.set_title(axtitle, fontsize=8)#, y=1.14)
                        if not std:
                            rad = 8
                            y, x = np.ogrid[-ymax/2+1: ymax/2+1, -xmax/2+1: xmax/2+1]
                            mask = x**2 + y**2 > rad**2
                            std = np.nanstd(img[mask])
                            print 'std with sigma clipping', std
                        if not vmin: vmin = -5 * std
                        if not vmax: vmax = 5 * std
                        a.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, extent=[xmin, xmax, ymax, ymin])
                        if contours:
                            levels = [(n + 1) * std for n in range(nsigma)]
                            a.contour(img, levels, aspect='auto', colors=['black'] * nsigma, extent=[xmin, xmax, ymin, ymax])
                        i += 1
            print "Writing PDF page", nfig
            pdf.savefig()
            plt.close()


def rdarg(argv, key, type=None, default=None):
    #print argv, key, type, default

    if len(argv) > 1:
        opt = np.where([a == '-%s' % key for a in argv])[0] + 1
        if len(opt) > 0:
            name = argv[int(opt)]
            if type is list:
                name = name.split()
                name = [float(i) for i in name]
            elif type is bool:
                print key, name
                name = eval(str(name))
            elif type is int:
                print key, name
                name = int(name)
            elif type is float:
                name = float(name)
            elif type is str:
                name = str(name)
            return name

    if default is not None:
        return default


def stack(lst, out, imtype='mean', sclip=5, zw0=1, makeim=True, title='', vmin=None, vmax=None, npix=True, var=True,
          output=True, smooth=True, flipshear=False, snr=True, corrmean=False, flipy=False, xmin=None, xmax=None,
          ymin=None, ymax=None, randflipy=False, overwrite=False, std=None, highq=False):
    from pyfits import getdata, PrimaryHDU

    def cubim(data, name, zw0, label, vmi=None, vma=None, contours=True, smooth=True):
        if not os.path.isfile(name) or overwrite:
            zl, xl, yl = data.shape
            zw = '%d:%d' % (zl / 2 - zw0 + 1, zl / 2 + zw0 + 1)
            hdu.data = data
            print 'cubim', name
            hdu.writeto(name, clobber=True)
            s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (name, zw, name.replace('.fits', '.IM.fits'))
            print s
            os.system(s)
        astroim(name.replace('.fits', '.IM.fits'), smooth=smooth, saveim=makeim, show=False, cbfrac=.08, pad=.006,
                    dfig=(8, 10), contours=contours, scb_label=label, title=title, vmin=vmi, vmax=vma, std=std, highq=highq)

    _output = []
    if not os.path.isfile(out) or overwrite:
        fits = []
        hdu = PrimaryHDU()
        if flipshear:
            r1=1
            r2=2
            for l in lst:
                ff = getdata(l)
                if r1 > r2:
                    ff = ff[::-1, :, :]
                fits.append(ff)
                if flipy:
                    fits.append(ff[:, ::-1, :])
                if randflipy:
                    if np.random.randint(2):
                        fits.append(ff[:, ::-1, :])
                del ff
        else:
            for l in lst:
                ff = getdata(l)
                fits.append(ff)
                if flipy:
                    fits.append(ff[:, ::-1, :])
                if randflipy:
                    if np.random.randint(2):
                        fits.append(ff[:, ::-1, :])
                del ff
        fits = np.array(fits).astype(float)
        bad = fits < -100
        if len(bad) > 0:
            fits[bad] = np.nan
        nl, zl, yl, xl = fits.shape
        #Estimate std in the region outise the central source
        rad = 8
        y, x = np.ogrid[-yl/2+1: yl/2+1, -xl/2+1: xl/2+1]
        mask = x**2 + y**2 > rad**2
        if corrmean:
            for i in range(nl):
                outmean = np.nanmean(fits[i, :, mask])
                fits[i] = fits[i] - outmean
                print "mean correction", outmean
        if std: sigma = std
        else: sigma = np.nanstd(fits[:, :, mask])
        print 'Stack sigma!!!!!!!!!!!!!!!!!!!!!!!!!!!', sigma
        high_sigma = np.abs(fits) > sclip * sigma
        fits[high_sigma & mask] = np.nan
        if imtype == 'mean':
            f = np.nanmean(fits, axis=0)
        if imtype == 'flux':
            f = np.nansum(fits, axis=0)
        if imtype == 'median':
            from scipy import stats
            f = stats.nanmedian(fits, axis=0)
        _output.append(f)
        if makeim: cubim(f, out, zw0, 'Flux [1e-20 cgs]', vmin, vmax)

        if npix:
            nf = np.nansum(np.isfinite(fits), axis=0)
            if makeim: cubim(nf, out.replace('.fits', '.NPIX.fits'), zw0, label='Number of voxels', vmi=0, contours=False, smooth=False)#, vmi=0, std=1000)#, vmi=np.amin(nf), vma=np.amax(nf))
            _output.append(nf)
        if var:
            fvar = np.nanvar(fits, axis=0)
            if makeim: cubim(fvar, out.replace('.fits', '.VAR.fits'), zw0, label='Variance')#, vmi=0, vma=30)# vmi=np.amin(var), vma=np.amax(var))
            _output.append(fvar)
        #hdu._close()
        del fits

    else:
        print out, 'already exists.'
        f = getdata(out)
        _output.append(f)
        if makeim: cubim(f, out, zw0, 'Flux [1e-20 cgs]', vmin, vmax)
        if npix:
            nf = getdata(out.replace('.fits', '.NPIX.fits'))
            if makeim: cubim(nf, out.replace('.fits', '.NPIX.fits'), zw0, label='Number of voxels', vmi=0, contours=False, smooth=False)#, vmi=np.amin(nf), vma=np.amax(nf))
            _output.append(nf)
        if var:
            fvar = getdata(out.replace('.fits', '.VAR.fits'))
            if makeim: cubim(fvar, out.replace('.fits', '.VAR.fits'), zw0, label='Variance')
            _output.append(fvar)

    if output:
        if len(_output) == 1:
            return _output[0]
        else:
            return _output


def analysis(f, nf, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0, frs=None):
    stack_stds = []
    stack_flux = []
    stack_npix = []
    randf = []
    _stack_stds = []
    _stack_flux = []
    _stack_npix = []
    _randf = []

    for i in stdbins:
        stack_flux.append(np.nanmean(f[zl / 2 - zw0: zl / 2 + zw0 + 1, yl / 2 - yw0: yl / 2 + yw0 + 1, i:i + stdbin]))
        _stack_flux.append(np.nanmean(f[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))

        stack_npix.append(np.sum(nf[zl / 2 - zw0: zl / 2 + zw0 + 1, yl / 2 - yw0: yl / 2 + yw0 + 1, i:i + stdbin]))
        _stack_npix.append(np.nanmean(nf[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))

        stack_stds.append(np.nanstd(
            np.concatenate((f[zl / 2 - zw0: zl / 2 + zw0 + 1, yl / 2 + yw0 + 2:, i:i + stdbin],
                            f[zl / 2 - zw0: zl / 2 + zw0 + 1, :yl / 2 - yw0 - 1, i:i + stdbin]))))
        _stack_stds.append(np.nanstd(
            np.concatenate((f[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 + yw0 + 2:],
                            f[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, :xl / 2 - yw0 - 1]))))
        if frs is not None:
            randf.append(np.nanmean(frs[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))
            _randf.append(np.nanmean(frs[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))

    conv1 = (zw0*2+1)*1.25*(yw0*2+1)*stdbin
    conv2 = (zw0*2+1)*1.25/(binsize*.2)
    stack_flux = np.array(stack_flux)
    stack_stds = np.array(stack_stds)
    _stack_flux = np.array(_stack_flux)
    _stack_stds = np.array(_stack_stds)
    if frs is not None:
        randf = np.array(randf)
        _randf = np.array(_randf)
        h = "theta SB npix 1sigma_SB SB_rand"
        all = [(stdbins - stdbins[-1]/2.) * (binsize*.2), stack_flux*conv1, stack_npix, stack_stds*conv2, randf*conv1]
        _all = [(stdbins - stdbins[-1]/2.) * (binsize*.2), _stack_flux*conv1, _stack_npix, _stack_stds*conv2, _randf*conv1]
    else:
        h = "theta SB npix 1sigma_SB"
        all = [(stdbins - stdbins[-1]/2.) * (binsize*.2), stack_flux*conv1, stack_npix, stack_stds*conv2]
        _all = [(stdbins - stdbins[-1]/2.) * (binsize*.2), _stack_flux*conv1, _stack_npix, _stack_stds*conv2]

    return all, _all, h