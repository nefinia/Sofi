#!/usr/bin/env python
__author__ = 'nefinia'

from pyfits import getdata
import h5py
import numpy as np
import os, subprocess, time
from math import *
from tools_sofi import rdarg
from random import randrange, uniform, randint, random, choice
from sys import argv


def b2s(s):
    if s: return '.true.'
    else: return '.false.'


# initial parameters
# random.seed()
csub = rdarg(argv, 'csub', bool, True)
daysmin = rdarg(argv, 'days', float, 10000)  # days = 1e-10#7
docat = rdarg(argv, 'docat', bool, True)
dorandom = rdarg(argv, 'random', bool, False)
dovar = rdarg(argv, 'dovar', bool, False)
extraname = rdarg(argv, 'extraname', str, '')
fitcat = rdarg(argv, key='fitcat', type=str, default='HDFS')  # 'UDF
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
folderout = rdarg(argv, 'folder', str, '/net/astrogate/export/astrodata/gallegos/')
i1 = rdarg(argv, 'id1', int, 0)
i2 = rdarg(argv, 'id2', int, 0)
idmin = rdarg(argv, 'idmin', int, 0)
idmax = rdarg(argv, 'idmax', int, 999999)
laecentered = rdarg(argv, 'laecentered', bool, True)
hdf5 = rdarg(argv, 'hdf5', bool, True)
norm = rdarg(argv, 'norm', bool, False)
noshear = rdarg(argv, 'noshear', bool, False)
overwrite = rdarg(argv, 'overwrite', bool, False)
nmin = rdarg(argv, 'nmin', int, 0)
nmax = rdarg(argv, 'nmax', int, 0)
randdir = rdarg(argv, 'rdir', bool, False)
redmin = rdarg(argv, 'rmin', list, [2.9])
redmax = rdarg(argv, 'rmax', list, [4.0])
simon = rdarg(argv, 'cubical', bool, True)
single = rdarg(argv, 'single', bool, False)
sizemin = rdarg(argv, 'size', int, 1)
r = rdarg(argv, 'rad', int, 50)
xw = rdarg(argv, 'xw', int, 80)
yw = rdarg(argv, 'yw', int, 80)
zw = rdarg(argv, 'zw', int, 10)

if fitcat == 'HDFS' or fitcat == 'mosaic': cut = 12
else: cut = 0
if not extraname: extraname = ''
if noshear: extraname += '.noshear'
if not laecentered: extraname += '.pair'

now = time.time()
byte2Mb = 9.5367431640625e-07


def outside(x, y):
    xmax = 326
    ymax = 331
    if x < 1 or x > xmax or y < 1 or y > ymax: return True
    else: return False


if fitcat == 'HDFS':
    cat = '%s/%s/cats/laes.fits' % (folder, fitcat)
    data = getdata(cat, 1)
    ids = data['ID']
    zs = data['redshift']
    ra = data['RaHMS']
    dec = data['DecDMS']
    sconf = data['Sofi_confidence']

if fitcat == 'UDF':
    cat = '%s/%s/cats/laes.fits' % (folder, fitcat)
    data = getdata(cat, 1)
    ids = data['ID']
    zs = data['Z_MUSE']
    sconf = data['Sofi_confidence']

if fitcat == 'mosaic':
    cat = '%s/%s/cats/laes.fits' % (folder, fitcat)
    data = getdata(cat, 1)
    ids = data['ID']
    zs = data['Z_MUSE']
    sconf = np.zeros(len(ids))+2#data['Sofi_confidence']

if fitcat == 'mosaic10':
    cat = '%s/%s/cats/laes.fits' % (folder, 'UDF')
    data = getdata(cat, 1)
    ids = data['ID']
    zs = data['Z_MUSE']
    sconf = data['Sofi_confidence']

if fitcat == 'EAGLE':
    cat = '%s/%s/cats/laes.fits' % (folder, 'UDF')
    data = getdata(cat, 1)
    ids = data['ID']
    zs = data['Z_MUSE']
    sconf = np.zeros(len(ids))+2

if fitcat == 'mosaic10':
    cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, 'UDF')
    fitcat = 'mosaic'
else: cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, fitcat)
data = getdata(cat2, 1)

if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR.fits'
    filevar = '%s/%s/%s' % (folder, fitcat, cubename)

if fitcat == 'UDF':
    cubename = 'DATACUBEFINALuser_20141021T055145_72f74684.fits'
    filevar = '%s/%s/%s' % (folder, fitcat, 'DATACUBEFINALuser_20141021T055145_212058ad.fits')
    # cubename = 'DATACUBEFINALuser_20141021T055145_212058ad.fits'

if fitcat == 'mosaic' or fitcat == 'mosaic10':
    cubename = 'mosaic.z1300.fits'#####
    folder = folderout

if csub: filename = '%s/%s/%s' % (folder, fitcat, cubename.replace('.fits', '.csub.fits'))
else: filename = '%s/%s/%s' % (folder, fitcat, cubename)

if laecentered:
    id1 = np.concatenate((data['id1'], data['id2']), 0)
    id2 = np.concatenate((data['id2'], data['id1']), 0)
    sconf1 = np.concatenate((data['sconf1'], data['sconf2']), 0)
    sconf2 = np.concatenate((data['sconf2'], data['sconf1']), 0)
    x1 = np.concatenate((data['x1'], data['x2']), 0)
    y1 = np.concatenate((data['y1'], data['y2']), 0)
    z1 = np.concatenate((data['z1'], data['z2']), 0)
    x2 = np.concatenate((data['x2'], data['x1']), 0)
    y2 = np.concatenate((data['y2'], data['y1']), 0)
    z2 = np.concatenate((data['z2'], data['z1']), 0)
    r1 = np.concatenate((data['redshift1'], data['redshift2']), 0)
    r2 = np.concatenate((data['redshift2'], data['redshift1']), 0)
    dist = np.concatenate((data['pi_Mpc'], data['pi_Mpc']), 0)
    theta = np.concatenate((data['theta'], data['theta']), 0)
else:
    id1 = data['id1']
    id2 = data['id2']
    x1 = data['x1']
    y1 = data['y1']
    z1 = data['z1']
    x2 = data['x2']
    y2 = data['y2']
    z2 = data['z2']
    r1 = data['redshift1']
    r2 = data['redshift2']
    dist = data['pi_Mpc']
    theta = data['theta']
    sconf1 = data['sconf1']
    sconf2 = data['sconf2']


if single: idpairs = np.where((id1 == i1) & (id2 == i2))[0]
else:
    idpairs = np.where((dist <= 20) & (dist >= .5) & (theta > 6) & (sconf1 > 0) & (sconf2 > 0)
                   & (r1 < redmax) & (r2 < redmax) & (r1 > redmin) & (r2 > redmin)
                   & (id1 > idmin) & (id1 < idmax))[0]  # & (r1 < 4) & (r2 < 4)

if nmax > nmin: idpairs = idpairs[nmin:nmax]

npairs = len(idpairs)
shear = (z2 - z1)[idpairs]
# idpairs = idpairs[::-1]
#idpairs = np.where(((id1 == 221) & (id2 == 287)))[0]# | ((id1 == 216) & (id2 == 144)))[0]############

if simon: sstype = 'cubical'
else: sstype = 'cylindrical'
if dorandom: sstype = 'random_' + sstype
if norm: ssnorm = 'normalized'
else: ssnorm = 'not-normalized'

angs = []
shears = []

dxs = []
dys = []
dds = []
ddds = []
mat = []
fout2 = open('%s/%s/pairs/%s/%s%s.dat' % (folderout, fitcat, ssnorm, sstype, extraname), 'w')
fout2.write('#id1 id2 angle shear\n')

ss1 = []
ss2 = []

posnans = np.isnan(x1) | np.isnan(y1) | np.isnan(z1) | np.isnan(x2) | np.isnan(y2) | np.isnan(z2)

for i in idpairs:
    print 'Pair %d %d' % (id1[i], id2[i]), 'hdf5', hdf5, 'overwrite', overwrite

    if posnans[i]: print 'Position nan! i1', x1[i], y1[i], z1[i], 'i2', x2[i], y2[i], z2[i]
    else:
        output = '%s/%s/cats/%s_pair%d-%d%s' % (folderout, fitcat, sstype, id1[i], id2[i], extraname)
        print output

        if hdf5:
            output += '.h5'

        if dorandom:
            if laecentered:
                ang = random() * 2 * pi  # uniform(0, 2*pi)#(pi * .2, 2 * pi * .8)
                if randdir:
                    b = (-1) ** randint(0, 1)
                    xx2 = x1[i] + (b * x2[i] - x1[i]) * cos(ang) - (-y2[i] - y1[i]) * sin(ang)
                    yy2 = y1[i] + (b * x2[i] - x1[i]) * sin(ang) + (-y2[i] - y1[i]) * cos(ang)
                else:
                    j = choice(shear)
                    xx2 = x1[i] + (x2[i] - x1[i]) * cos(ang) - (y2[i] - y1[i]) * sin(ang)
                    yy2 = y1[i] + (x2[i] - x1[i]) * sin(ang) + (y2[i] - y1[i]) * cos(ang)
                    print "z2 old", z2[i]
                    z2[i] = z1[i] + j
                    print "z2 new", z2[i]
            else:
                out = True
                while out:
                    ang = random() * 2 * pi
                    xx2 = x1[i] + (x2[i] - x1[i]) * cos(ang) - (y2[i] - y1[i]) * sin(ang)
                    yy2 = y1[i] + (x2[i] - x1[i]) * sin(ang) + (y2[i] - y1[i]) * cos(ang)

                    out = outside(x1[i], y1[i]) or outside(xx2, yy2)
                    if (randrange(1, 100) % 100 == 1):
                        print 'aaaaa', ang, xx2, yy2

            fout2.write('%d %d %.3f %d\n' % (id1[i], id2[i], ang * 180. / pi, j))

            if noshear:
                z2[i] = z1[i]
            s1 = '"%d %d %d"' % (np.round(x1[i]), np.round(y1[i]), np.round(z1[i]))
            s2 = '"%d %d %d"' % (np.round(xx2), np.round(yy2), np.round(z2[i]))
            x2[i] = np.round(xx2)
            y2[i] = np.round(yy2)
        else:
            angs.append(atan2(y2[i] - y1[i], x2[i] - x1[i]))
            if noshear:
                z2[i] = z1[i]
            shears.append(z1[i] - z2[i])
            if fitcat == 'mosaic':
                s1 = '"%d %d %d"' % (x1[i], y1[i], z1[i])
                s2 = '"%d %d %d"' % (x2[i], y2[i], z2[i])
            else:
                s1 = '"%d %d %d"' % (x1[i], y1[i], z1[i])
                s2 = '"%d %d %d"' % (x2[i], y2[i], z2[i])

        if simon:
            coords = 'cub'
        else:
            coords = 'cyl'

        sc = 'Coord-conv -cube %s -PosIn %s -PosEnd %s -print .true. -cut %d -single %s -dovar %s -norm %s -coords %s' % \
                 (filename, s1, s2, cut, b2s(laecentered), b2s(dovar), b2s(norm), coords)
        if simon:
            sc += ' -xw %d -yw %d -zw %d' % (xw, yw, zw)
        else:
            sc += ' -rad %d' % r

        if dovar: sc += ' -varFile %s' % filevar
        print sc
        print 'output:', output
        dorun = False

        if not os.path.isfile(output):
            print "Output file does not exist."
            dorun = True

        else:
            ftime = os.stat(output).st_mtime
            fsize = os.path.getsize(output)

            if now - ftime > daysmin * 86400:
                print "File is older than %d days" % daysmin
                dorun = True
                os.system('rm %s' % output)
            else:
                print "File is younger than  %d days!!!!!!!!!!!" % daysmin

            if fsize * byte2Mb < sizemin:
                print "File is smaller than %d Mb" % sizemin
                dorun = True
                os.system('rm %s' % output)
            else:
                print "File is bigger than %d Mb!!!!!!!!!!" % sizemin

            #((now - ftime > days * 86400) or (fsize * byte2Mb < 20.)) and

            if overwrite:
                print "Overwriting output file"
                dorun = True
                os.system('rm %s' % output)


        if dorun:
            temp = subprocess.check_output(sc, shell=True)
            nh = 6
            if simon: nh += 1
            if dovar: nh += 1

            if hdf5:
                with h5py.File(output, 'w') as h:
                    dpix = sqrt((x1[i] - x2[i]) ** 2 + (y1[i] - y2[i]) ** 2)
                    data = np.array(temp.split())

                    data = data.reshape(len(data) / nh, nh).T
                    h.create_dataset('x', data=data[0].astype('int'))
                    h.create_dataset('y', data=data[1].astype('int'))
                    h.create_dataset('z', data=data[2].astype('int'))
                    h.create_dataset('px', data=data[3].astype('float'))
                    h.create_dataset('py', data=data[4].astype('float'))
                    if simon:
                        h.create_dataset('pz', data=data[5].astype('int'))
                        h.create_dataset('flux', data=data[6].astype('float'))
                    else:
                        h.create_dataset('flux', data=data[5].astype('float'))
                    if dovar: h.create_dataset('variance', data=data[-1].astype('float'))
                    h.create_dataset('dpix', data=dpix)
            else:
                fout = open(output, 'w')
                fout.write(' '.join(h) + "\n")
                fout.write(temp)
                fout.close
        else:
            print "Ouput file will not be produced\n\n"

fout2.close()

if not docat:
    x2 = x2[idpairs]
    x1 = x1[idpairs]
    y2 = y2[idpairs]
    y1 = y1[idpairs]
    import matplotlib.pyplot as plt

    plt.scatter(x2 - x1, y2 - y1, s=2)
    plt.show()
    # plt.hist(np.array(angs)*180/np.pi)
    # plt.savefig('%s/%s/pairs/%s/cats/%s_ang_dist%s.png' % (folder, fitcat, ssnorm, sstype, extraname), format='png')
    # plt.close()
