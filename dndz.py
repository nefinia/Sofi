import numpy as np
a=np.loadtxt('../../fLLS/prochaska14.dat', dtype=str).T
z=np.array([float(i) for i in a[1]])
q=np.array([i for i in a[0]])
zps = [[1.7, 2.5],[2.5, 3], [3, 3.3], [3.3, 3.7], [3.2, 4.3], [3.8, 6]]

if 0:
	for zp in zps:
		cool = (z>zp[0]) & (z<zp[1])
		zmin = np.amin(z[cool])
		zmax = np.amax(z[cool])
		dz = zmax-zmin
		zmean = np.mean(z[cool])
		nlls = np.sum(cool)
		nq = len(np.unique(q[cool]))
		print zp, 'zmean', zmean, 'dndz', nlls/float(nq)/dz, 'nlls', nlls, 'nq', nq

redshifts = {'7': 5.487, '8': 5.037, '9': 4.485, '10': 3.984, '11': 3.528, '12': 3.017, '13': 2.478, '14': 2.237, '15': 2.012}


def l(z, ls, zs, a):
	return ls*np.power((1+z)/(1.+zs), a)*1.1

ls = 1.46
_ls = .11
a = 1.7
_a = .22
zs = 3
for snap in range(7, 16):
	dndz = []
	s = str(snap)
	z = redshifts[s]
	print 'snap', s, 'z', z
	dndz.append(l(z, ls, zs, a))
	dndz.append(l(z, ls, zs, a+_a))
	dndz.append(l(z, ls, zs, a-_a))
	dndz.append(l(z, ls+_ls, zs, a))
	dndz.append(l(z, ls-_ls, zs, a))
	dndz.append(l(z, ls-_ls, zs, a+_a))
	dndz.append(l(z, ls-_ls, zs, a-_a))
	dndz.append(l(z, ls+_ls, zs, a+_a))
	dndz.append(l(z, ls+_ls, zs, a-_a))
	lmax = np.amax(dndz)
	lmin = np.amin(dndz)
	print dndz[0], lmax, lmin, (lmax+lmin)/2.
