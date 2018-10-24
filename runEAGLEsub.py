import numpy as np
import os
from pyfits import getdata
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

fitcat = 'EAGLE'
coords = ['x', 'y', 'z']
snaps = [10, 11, 12]
foldercat = '../../'
superscript = open('subEAGLE.sh', 'w')
overwrite = False
imov = True


comm.send(idata, dest=1)

if 0:

	ncores = 36


	def b2s(b):
		if b: return 'True'
		else: return 'False'

	name = 'subEAGLE.sh'
	print name
	script = open(name, 'w')
	n = 0
	for snap in snaps:
		for coord in coords:
			cat = '%s/%s/cats/gals_snap%d.fits' % (foldercat, fitcat, snap)
			data = getdata(cat, 1)
			ids = data['ID']
			for i in ids:
				script.write('python subcubes.py -overwrite %s -imov %s -fitcat %s -snap %d -coord %s -single True -id1 %d &\n'
							 % (b2s(overwrite), b2s(imov), fitcat, snap, coord, i))
				n += 1
			#script.write('wait')
			script.close()
			superscript.write('bash %s\nwait\n' % name)

	superscript.close()