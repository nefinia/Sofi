import numpy as np
import h5py
import params
from tools_sofi import rdarg  # , hubblefrom tools_sofi import cic, cubex, makejpeg, astroim,
from sys import argv

snaps = rdarg(argv, 'snap', list, [8, 9, 10, 11, 12, 13, 14], int)
_ssthr = rdarg(argv, 'ssthr', float, 6.73e-3)#None) or 1e10
model = rdarg(argv, 'model', str, 'HM01')#'HM12')#

def read_dataset(att, snap, itype=0, nfiles=16):
	""" Read a selected dataset, itype is the PartType and att is the attribute name. """

	# Output array.
	print 'read dataset snap', snap
	data = []
	sname = params.snames[str(snap)]
	# Loop over each file and extract the data.
	for i in range(nfiles):
		f = h5py.File('/scratch/gallegos/EAGLE/R/snapshot_%s/snap_%s.%i.hdf5'
					  % (sname, sname, i), 'r')
		tmp = f['PartType%i/%s' % (itype, att)][...]
		data.append(tmp)

		# Get conversion factors.
		cgs = f['PartType%i/%s' % (itype, att)].attrs.get('CGSConversionFactor')
		aexp = f['PartType%i/%s' % (itype, att)].attrs.get('aexp-scale-exponent')
		hexp = f['PartType%i/%s' % (itype, att)].attrs.get('h-scale-exponent')

		# Get expansion factor and Hubble parameter from the header.
		a = f['Header'].attrs.get('Time')
		h = f['Header'].attrs.get('HubbleParam')

		f.close()

	# Combine to a single array.
	if len(tmp.shape) > 1:
		data = np.vstack(data)
	else:
		data = np.concatenate(data)

	# Convert to physical.
	if data.dtype != np.int32 and data.dtype != np.int64:
		data = np.multiply(data, cgs * a ** aexp * h ** hexp, dtype='f8')


	return data


att = 'Temperature', 'Density'

for snap in snaps:
	s = str(snap)
	gamma_bkg = params.gamma_bkg[model][s][0]
	sname = params.snames[str(snap)]
	nfiles = 16
	for i in range(nfiles):
		Temp, Dens = read_dataset('Temperature', snap), read_dataset('Density', snap)
		Temp[Dens > .1] = 1.e4
		Temp[Temp < 100] = 100.
		UVBSSThr_phys = _ssthr * ((gamma_bkg * 1.e12) ** (2. / 3.))
		ThisGamma = gamma_bkg*(0.98*(1.+(Dens/UVBSSThr_phys)**1.64)**(-2.28)+0.02*(1.+Dens/UVBSSThr_phys)**(-0.84))
		aA = 1.269e-13 * (315614. / Temp) ** 1.503 / (1 + (315614. / Temp / 0.522) ** 0.47) ** 1.923
		LT = 1.17e-10 * np.sqrt(Temp) * np.exp(-157809. / Temp) / (1 + np.sqrt(Temp * 1.e-5))
		A = aA + LT
		B = 2. * aA + ThisGamma / Dens + LT
		C = aA
		nHI = .5*(B - np.sqrt(B * B - 4. * A * C)) / A
