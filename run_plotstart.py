import os
overwrite = False
fin = '/net/astrogate/export/astrodata/saeed/EAGLE/RefL0025N0752/'

names = ['010_z003p984', '011_z003p528', '012_z003p017', '013_z002p478', '014_z002p237', '015_z002p012']
layers = [34, 29, 25, 20, 18, 16]
types1 = ['NHI', 'NHII']
types2 = ['stars', 'SFR']
script = open('run_plotstart.sh', 'w')
for name, ls in zip(names, layers)[:3]:
	input = '%s/snapshot_%s/init-evol-512-4096-8/SO.snap_%s_512_4096_8_init_evol.ts0000' % (fin, name, name)
	input2 = '%s/snapshot_%s/chombo-512-4096-8/snap_%s_512_4096_8.chombo.hdf5' % (fin, name, name)
	fout = '%s/snapshot_%s/layers/' % (fin, name)
	fbin = '%s/snapshot_%s/layers/PlotStArt/' % (fin, name)
	script.write('fbin=%s\ninp=%s\ninp2=%s\nfout=%s\n\n' % (fbin, input, input2, fout))
	fl = 1 / float(ls)
	for l in range(ls):
		f = (l + 1) * fl
		i = f - fl

		for t in types1:
			outs = ['%s/snap%s_512_4096_8_x_%d.%s.fits' % (fout, name, l, t),
					'%s/snap%s_512_4096_8_y_%d.%s.fits' % (fout, name, l, t),
					'%s/snap%s_512_4096_8_z_%d.%s.fits' % (fout, name, l, t)]
			if not os.path.isfile(outs[0]):
				script.write(
					'$fbin./PlotStART -inp "$inp" -out "%s" -var %s -ftype fits -proj 1 '
					'-box_le "%0.3f,0,0" -box_re "%0.3f,1,1" > "$fout/log_x_%d.%s.txt" &\n' % (
					outs[0].replace(fout, '$fout'), t, i, f, l, t))
			if not os.path.isfile(outs[1]):
				script.write(
					'$fbin./PlotStART -inp "$inp" -out "%s" -var %s -ftype fits -proj 2 '
					'-box_le "0,%0.3f,0" -box_re "1,%0.3f,1" > "$fout/log_y_%d.%s.txt" &\n' % (
					outs[1].replace(fout, '$fout'), t, i, f, l, t))
			if not os.path.isfile(outs[2]):
				script.write(
					'$fbin./PlotStART -inp "$inp" -out "%s" -var %s -ftype fits -proj 3 '
					'-box_le "0,0,%0.3f" -box_re "1,1,%0.3f" > "$fout/log_z_%d.%s.txt" &\n' % (
					outs[2].replace(fout, '$fout'), t, i, f, l, t))
		if 0:#for type in types2:
			script.write(
				'$fbin./PlotStART -inp "$inp2" -out "$fout/snap%s_x_%d.%s.fits" -var %s -ftype fits -proj 1 '
				'-box_le "%0.3f,0,0" -box_re "%0.3f,1,1" > "$fout/log_x_%d.%s.txt" &\n' % (
				name, l, type, type, i, f, l, type))
			script.write(
				'$fbin./PlotStART -inp "$inp2" -out "$fout/snap%s_y_%d.%s.fits" -var %s -ftype fits -proj 2 '
				'-box_le "0,%0.3f,0" -box_re "1,%0.3f,1" > "$fout/log_y_%d.%s.txt" &\n' % (
				name, l, type, type, i, f, l, type))
			script.write(
				'$fbin./PlotStART -inp "$inp2" -out "$fout/snap%s_z_%d.%s.fits" -var %s -ftype fits -proj 3 '
				'-box_le "0,0,%0.3f" -box_re "1,1,%0.3f" > "$fout/log_z_%d.%s.txt" &\n' % (
				name, l, type, type, i, f, l, type))
		script.write('wait\n\n')
script.close()
