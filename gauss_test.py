#!/usr/bin/env python
__author__ = 'nefinia'

from pyfits import getdata
import h5py
import numpy as np
import os, subprocess
from math import *
from tools_sofi import rdarg
from random import randrange, uniform, randint, random, choice
from sys import argv

#random.seed()
fitcat = rdarg(argv, key='fitcat', type=str, default='Test')  # 'UDF
folder = rdarg(argv, 'folder', str, '../../')  # '/scratch/gallegos/MUSE/'
folderout = rdarg(argv, 'folder', str, '/net/astrogate/export/astrodata/gallegos/cats/')
dorandom = rdarg(argv, 'random', bool, False)
docat = rdarg(argv, 'docat', bool, True)
hdf5 = True
norm = rdarg(argv, 'norm', bool, False)
nodist = rdarg(argv, 'nodist', bool, True)
noshear = rdarg(argv, 'noshear', bool, False)
randdir = rdarg(argv, 'rdir', bool, False)

if nodist:
    snodist = '.true.'
else:
    snodist = '.false.'

simon = rdarg(argv, 'cubical', bool, True)
half = rdarg(argv, 'half', bool, True)
if half:
    shalf = '.true.'
else:
    shalf = '.false.'
overwrite = rdarg(argv, 'overwrite', bool, False)
csub = True
r = rdarg(argv, 'rad', int, 50)
xw = rdarg(argv, 'xw', int, 80)
yw = rdarg(argv, 'yw', int, 80)
zw = rdarg(argv, 'zw', int, 10)

if fitcat == 'HDFS':
    cut = 12
else:
    cut = 0
extraname = rdarg(argv, 'extraname', str, '')
if not extraname:
    extraname = ''

if noshear:
    extraname = '.noshear' + extraname

def outside(x, y):
    xmax = 326
    ymax = 331
    if x < 1 or x > xmax or y < 1 or y > ymax:
        return True
    else:
        return False



def gauss(x,y,z):


np.matrix()

cat2 = '%s/%s/cats/lae_pairs.fits' % (folder, fitcat)
data = getdata(cat2, 1)

if fitcat == 'HDFS':
    cubename = 'DATACUBE-HDFS-1.35-PROPVAR.fits'

if fitcat == 'UDF':
    cubename = 'DATACUBEFINALuser_20141021T055145_212058ad.fits'

filevar = '%s/%s/%s' % (folder, fitcat, cubename)

if csub:
    filename = '%s/%s/%s' % (folder, fitcat, cubename.replace('.fits', '.csub.fits'))
else:
    filename = '%s/%s/%s' % (folder, fitcat, cubename)


if half:
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

idpairs = np.where((dist <= 20) & (dist >= .5) & (theta > 6))[0]# & (r1 < 4) & (r2 < 4))[0]
#idpairs = np.where(((id1 == 144) & (id2 == 216)))[0]# | ((id1 == 216) & (id2 == 144)))[0]
npairs = len(idpairs)
shear = (z2 - z1)[idpairs]
#idpairs = np.where(((id1 == 275) & (id2 == 170)))[0]# | ((id1 == 216) & (id2 == 144)))[0]
#idpairs = np.where(((id1 == 144) & (id2 == 216)))[0]# | ((id1 == 216) & (id2 == 144)))[0]
#idpairs = idpairs[::-1]

if simon:
    sstype = 'cubical'
else:
    sstype = 'cylindrical'
if dorandom:
    sstype = 'random_' + sstype
if norm:
    ssnorm = 'normalized'
else:
    ssnorm = 'not-normalized'

angs = []
shears = []

dxs = []
dys = []
dds = []
ddds = []
mat = []
fout2 = open('%s/%s/pairs/%s/cats/%s%s.dat' % (folder, fitcat, ssnorm, sstype, extraname), 'w')
fout2.write('#id1 id2 angle shear\n')
for i in idpairs:
    print 'Pair %d %d' % (id1[i], id2[i]), 'hdf5', hdf5, 'overwrite', overwrite

    output = '%s/%s/%s_pair%d-%d%s' % (folderout, fitcat, sstype, id1[i], id2[i], extraname)

    if hdf5:
        output += '.h5'

    if os.path.isfile(output) and not overwrite:
        print 'File %s already exists' % output

    else:
        if dorandom:
            if nodist:
                ang = random()*2*pi#uniform(0, 2*pi)#(pi * .2, 2 * pi * .8)
                if randdir:
                    b = (-1)**randint(0,1)
                    xx2 = x1[i] + (b*x2[i] - x1[i]) * cos(ang) - (-y2[i] - y1[i]) * sin(ang)
                    yy2 = y1[i] + (b*x2[i] - x1[i]) * sin(ang) + (-y2[i] - y1[i]) * cos(ang)
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
                    ang = random()*2*pi
                    xx2 = x1[i] + (x2[i] - x1[i]) * cos(ang) - (y2[i] - y1[i]) * sin(ang)
                    yy2 = y1[i] + (x2[i] - x1[i]) * sin(ang) + (y2[i] - y1[i]) * cos(ang)

                    if half:
                        out = outside(x1[i], y1[i]) or outside(xx2/2., yy2/2.)
                    else:
                        out = outside(x1[i], y1[i]) or outside(xx2, yy2)
                    if(randrange(1,100)%100 == 1):
                        print 'aaaaa',ang, xx2, yy2

            fout2.write('%d %d %.3f %d\n' % (id1[i], id2[i], ang*180./pi, j))

            if noshear:
                z2[i] = z1[i]
            s1 = '"%d %d %d"' % (np.round(x1[i]), np.round(y1[i]), np.round(z1[i]))
            s2 = '"%d %d %d"' % (np.round(xx2), np.round(yy2), np.round(z2[i]))
            x2[i] = np.round(xx2)
            y2[i] = np.round(yy2)
        else:
            angs.append(atan2(y2[i]-y1[i], x2[i]-x1[i]))
            if noshear:
                z2[i] = z1[i]
            shears.append(z1[i]-z2[i])
            s1 = '"%d %d %d"' % (x1[i], y1[i], z1[i])
            s2 = '"%d %d %d"' % (x2[i], y2[i], z2[i])

        if docat:
            if norm:
                snorm = '.true.'
            else:
                snorm = '.false.'

            if simon:
                coords = 'cub'
            else:
                coords = 'cyl'

            if half:
                if simon:
                    sc = 'Coord-conv -cube %s -PosIn %s -PosEnd %s -xw %d -yw %d -zw %d -print .true. -norm %s -coords %s -varFile %s -cut %d -half .false. -nodist %s' % (
                        filename, s1, s2, xw, yw, zw, snorm, coords, filevar, cut, snodist)
                else:
                    sc = 'Coord-conv -cube %s -PosIn %s -PosEnd %s -rad %d -print .true. -norm %s -coords %s -varFile %s -cut %d' % (
                        filename, s1, s2, r, snorm, coords, filevar, cut)
            else:
                if simon:
                    sc = 'Coord-conv -cube %s -PosIn %s -PosEnd %s -xw %d -yw %d -zw %d -print .true. -norm %s -coords %s -varFile %s -cut %d -half .false. -nodist %s' % (
                        filename, s1, s2, xw, yw, zw, snorm, coords, filevar, cut, snodist)
                else:
                    sc = 'Coord-conv -cube %s -PosIn %s -PosEnd %s -rad %d -print .true. -norm %s -coords %s -varFile %s -cut %d' % (
                        filename, s1, s2, r, snorm, coords, filevar, cut)
            print sc
            print 'output:', output

            if os.path.isfile(output) and overwrite:
                # print "Overwriting output file"
                os.system('rm %s' % output)

            temp = subprocess.check_output(sc, shell=True)

            if simon:
                nh = 8
            else:
                nh = 7

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
                        h.create_dataset('variance', data=data[7].astype('float'))
                    else:
                        h.create_dataset('flux', data=data[5].astype('float'))
                        h.create_dataset('variance', data=data[6].astype('float'))

                    h.create_dataset('dpix', data=dpix)

            else:
                fout = open(output, 'w')
                fout.write(' '.join(h) + "\n")
                fout.write(temp)
                fout.close


#mat = np.array([id1[idpairs], id2[idpairs], angs, np.array(shears)*180/np.pi])
#np.savetxt('%s/%s/pairs/%s/cats/%s%s.dat' % (folder, fitcat, ssnorm, sstype, extraname), nparmat.T,
           #header='id1 id2 angle shear')

fout2.close()

if not docat:

    x2=x2[idpairs]
    x1=x1[idpairs]
    y2=y2[idpairs]
    y1=y1[idpairs]
    import matplotlib.pyplot as plt
    plt.scatter(x2-x1, y2-y1, s=2)
    plt.show()
    #plt.hist(np.array(angs)*180/np.pi)
    #plt.savefig('%s/%s/pairs/%s/cats/%s_ang_dist%s.png' % (folder, fitcat, ssnorm, sstype, extraname), format='png')
    #plt.close()