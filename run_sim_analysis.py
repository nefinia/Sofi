#!/usr/bin/env python
__author__ = 'gallegos'
import h5py
import numpy as np

def run(snap, coord, i, xmin, xmax, yin, ymax, zmin, zmax):
    print 'Pair %d-%d, snap%d %s' % (id[0][cool[i]], id[1][cool[i]], snap, coord)
    cat = '%s%d-%d.h5' % (pairsname, id[0][cool[i]], id[1][cool[i]])
    data = h5py.File(cat, 'r')
    pz = np.array(data['pz'])
    zgood = (pz >= zmin) & (pz <= zmax)
    f = np.array(data['flux'])[zgood]
    px = np.array(data['px'])[zgood]
    py = np.array(data['py'])[zgood]
    right = (px > xmin) & (px < xmax[i]) & (py > ymin) & (py < ymax)
    left = (px < -xmin) & (px > -xmax[i]) & (py > ymin) & (py < ymax)
    down = (px > ymin) & (px < ymax) & (py < -xmin) & (py > -xmax[i])
    up = (px > ymin) & (px < ymax) & (py > xmin) & (py < xmax[i])
    rand = left | down | up
    f_ = np.nanmean(f[right])
    f_rand = np.nanmean(f[rand])
    fout.write('%d %d %f %f\n' % (id[0][cool[i]], id[1][cool[i]], f_, f_rand))
    for _z in zr:
        # remember f[zgood]..... so the sample is smaller already
        rightz = right & (pz == _z)
        randz = rand & (pz == _z)
        f_ = np.nanmean(f[rightz])
        f_rand = np.nanmean(f[randz])
        foutz[str(_z)].write(
            '%d %d %f %f\n' % (id[0][cool[i]], id[1][cool[i]], f_, f_rand))