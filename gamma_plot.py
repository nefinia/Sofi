import numpy as np
import matplotlib.pyplot as plt
fin = '../../UVB/'

if 1:
	g, gup, glo, z = {}, {}, {}, {}
	aut = 'Faucher-Giguere 2019', r'Khaire & Srianand $2019$', 'Haardt & Madau 2012', \
		  'Fumagalli et al. 2017', 'Gaikward et al. 2017a', 'Kollmeier et al. 2014', \
		  'Viel et al. 2017', 'Khaire & Srianand 2019', 'Becker & Bolton 2013', \
		  "D'Aloisio et al. 2018", 'Wyithe et al. 2011', 'Calverley et al. 2011', \
		  'This work'

	# Data
	#units 1e-12 1/s
	i=0
	a = np.loadtxt(fin+'fg19.dat').T
	g[aut[i]], z[aut[i]], gup[aut[i]], glo[aut[i]] = a[1], a[0], np.zeros(len(a[0])), np.zeros(len(a[0]))
	i+=1
	del a
	a = np.loadtxt(fin+'ks19.dat').T
	g[aut[i]], z[aut[i]], gup[aut[i]], glo[aut[i]] = a[1], a[0], np.zeros(len(a[0])), np.zeros(len(a[0]))
	del a
	i+=1
	a = np.loadtxt(fin+'hm12.dat').T
	g[aut[i]], z[aut[i]], gup[aut[i]], glo[aut[i]] = a[1], a[0], np.zeros(len(a[0])), np.zeros(len(a[0]))
	del a
	i+=1
	#Fumagali17
	g[aut[i]] = [.07]
	gup[aut[i]] = [.01]
	glo[aut[i]] = [.01]
	z[aut[i]] = [.0]
	i+=1
	#Gaikward17
	g[aut[i]] = [.066, .1, .145, .21]
	gup[aut[i]] = [.015, .021, .037, .052]
	glo[aut[i]] = [.015, .021, .037, .052]
	z[aut[i]] = [.1125, .2, .3, .4]
	i+=1
	#Kollmeier14
	g[aut[i]] = [.178]
	gup[aut[i]] = [.0]
	glo[aut[i]] = [.0]
	z[aut[i]] = [.1]
	i+=1
	#Viel17
	g[aut[i]] = [.071]
	gup[aut[i]] = [.101]
	glo[aut[i]] = [.049]
	z[aut[i]] = [.1]
	i+=1
	#Khaire19
	g[aut[i]] = [.0585, .0756, .1135, .1479, .1418]
	gup[aut[i]] = [.017, .017, .032, .044, .053]
	glo[aut[i]] = [.018, .018, .032, .052, .057]
	z[aut[i]] = [.03, .1, .2, .3, .41]
	i+=1
	#Becker and Bolton13
	g[aut[i]] = [1.03514217, 0.85901352, 0.78886012, 0.79983426, 0.84722741, 0.95719407, 0.93540567]
	gup[aut[i]] = [0.368, 0.297, 0.275, 0.282, 0.309, 0.364, 0.404]
	glo[aut[i]] = [0.2955, 0.2236, 0.1918, 0.19, 0.2, 0.2344, 0.2686]
	z[aut[i]] = [2.4, 2.8, 3.2, 3.6, 4., 4.4, 4.75]
	i+=1

	#D'Aloisio18
	g[aut[i]] = [.6, .55, .51, .5, .48, .29]
	gup[aut[i]] = [.08, .1, .11, .13, .15, .11]
	glo[aut[i]] = [.2, .2, .19, .19, .12, .12]
	z[aut[i]] = [4.8, 5., 5.2, 5.4, 5.6, 5.8]
	i+=1

	#Whyithe11
	g[aut[i]] = [.47, .18]
	gup[aut[i]] = [.3, .18]
	glo[aut[i]] = [.2, .09]
	z[aut[i]] = [5., 6.]
	i+=1

	#Caverley+11
	g[aut[i]] = [.7079, .1445]
	gup[aut[i]] = [.315, .074]
	glo[aut[i]] = [.21, .049]
	z[aut[i]] = [5., 6.]
	i+=1

	########################################
	##Gallego+19 work... check this!!!!!!!!!
	########################################

	#protocluster detections

	g[aut[i]] = [3.53, 4.57, 3.19]
	gup[aut[i]] = [1.5, 1.1, 1.26]
	glo[aut[i]] = gup[aut[i]]
	z[aut[i]] = [3.47, 3.71, 4.51]
	# 3.467<z<3.476, mean 3.47 detection r6-12 zw2 zoff-4 rho 6
	# 3.68<z<3.74, mean z 3.71 detection with r6-14 zw2 zoff-2 rho 4
	#4.496<z<4.53 mean 4.51 detection with r6-8 zw2 zoff1 rho ???


	#detections
	g[aut[i]] = [2.67, 2.02]
	gup[aut[i]] = [.69, .62]
	glo[aut[i]] = gup[aut[i]]
	z[aut[i]] = [4.01, 4.47]
	# 2.8<z<3.2 mean 3.07 detection 6<r<12 zw2 zoff-1 gamma .87+-.32
	# 3.8<z<5, mean z 4.47 detection with r6-10 zw2 zoff-2
	# 3.68<z<3.74, mean z 3.71 detection with r6-14 zw2 zoff-2
	# 3.5<z<5, mean z 4.28 detection with r6-10 zw2 zoff-2 gamma 1.26 +- .39
	# 3.8<z<4.8, mean z 4.39 detection with r6-13 zw2 zoff-2
	# 3.8<z<5.5, mean z 4.63 detection with r6-8 zw2 zoff-2 gamma 1.06 +- .39
	# 4.496<z<4.53 mean 4.51 detection with r6-8 zw2 zoff1


	#upper limits
	g[aut[i]] = [1.31, 1.28, 1.34]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.23, 3.47, 3.64]
	# all at 6<r<18

	#upper limits
	g[aut[i]] = [1.18, 1.33, 1.4]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.23, 3.47, 3.64]
	# z 3.23 at 6<r<12 2.8<z<3.5 gamma 1.18
	# z 3.47 at 6<r<18 3.2<z<3.8 gamma 1.33
	# z 3.64 at 6<r<18 3.2<z<4.2 gamma 1.4


	#upper limits UDF
	g[aut[i]] = [1.17, 1.55, 1.55, 1.54, 1.88]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.23, 3.46, 3.59, 4.43, 4.99]
	# z 3.23 at 6<r<20 2.8<z<3.5 gamma 1.17
	# z 3.46 at 6<r<18 3.2<z<3.8 gamma 1.55
	# z 3.59 at 6<r<18 3.2<z<4.2 gamma 1.55
	# z 4.43 at 6<r<8 3.8<z<5 gamma 1.54
	# z 4.99 at 5<r<20 4.5<z<5.5 gamma 1.88

	#upper limits UDF+mosaic
	g[aut[i]] = [1.2, .97, 1.29, 1.29, 1.29, 1.02, 1.25]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.07, 3.23, 3.47, 3.64, 4, 4.45, 4.96]
	# z 3.22 at 6<r<12 2.8<z<3.5 gamma 1.2
	#check the rest
	# z 3.46 at 6<r<18 3.2<z<3.8 gamma 1.55
	# z 3.59 at 6<r<18 3.2<z<4.2 gamma 1.55
	# z 4.43 at 6<r<8 3.8<z<5 gamma 1.54
	# z 4.99 at 5<r<20 4.5<z<5.5 gamma 1.88


	#upper limits UDF+mosaic
	g[aut[i]] = [1.2, .97, 1.29, 1.29, 1.02, 1.25]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.07, 3.23, 3.56, 4, 4.45, 4.96]
	#4.96 	4.5<z<5.5
	#4.45 	3.8<z<5
	#4 		3.5<z<4.5
	#3.64 	3.2<z<4.2
	#3.47	3.2<z<3.8
	#3.23	2.8<z<3.2



	#upper limits UDF+mosaic
	g[aut[i]] = [1.2, 1.29, 1.02, 1.25]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.07, 3.64, 4.45, 4.96]
	#4.96 	4.5<z<5.5
	#4.45 	3.8<z<5
	#4 		3.5<z<4.5
	#3.64 	3.2<z<4.2
	#3.07	2.8<z<3.2


	#upper limits UDF
	g[aut[i]] = [.75, .63, 1.39]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.23, 3.59, 3.97]
	# z 3.23 at 6<r<20 2.8<z<3.5 gamma .76 .56
	# z 3.46 at 5<r<20 3.2<z<3.8 gamma .91
	# z 3.59 at 10<r<20 3.2<z<4.2 gamma .84 .42
	# z 3.97 at 6<r<20 3.2<z<4.2 gamma

	#upper limits UDF
	g[aut[i]] = [.85, .68, 1.09, .8]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.26, 3.52, 4.17, 4.89]
	#zprob = [[2.8, 3.5], [3.2, 3.8], [3.8, 5], [4.5, 5.5]]#[2.8, 3.2]

	#upper limits UDF
	g[aut[i]] = [.56, 1.09, .8]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.37, 4.17, 4.89]
	#zprob = [[2.8, 4], [3.5, 5], [4.5, 5.5]]#[2.8, 3.2]

	#upper limits UDF
	g[aut[i]] = [.68, 1.09, .8]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.25, 4.17, 4.89]
	#zprob = [[2.8, 3.5], [3.5, 5], [4.5, 5.5]]#[2.8, 3.2]


	#upper limits UDF+mosaic
	g[aut[i]] = [.96, .75, 1.28, .73]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.07, 3.47, 4.18, 4.96]
	#zprob = [[2.9, 3.2], [3.2, 3.8], [3.8, 4.5], [4.5, 5.5]]
	#3.07 6<r<20 .96
	#3.47 8<r<14 .67, 6<r<20 .75
	#4.18 6<r<13 1.14, 6<r<20 1.28
	#4.88 6<r<10 .61, 6<r<20 .77
	#4.96 6<r<15 .64, 6<r<20 .73


	#upper limits UDF+mosaic
	g[aut[i]] = [.99, .74, .79]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.07, 3.47, 5.01]


	#detections
	g[aut[i]] = [2.67, 2.02]
	gup[aut[i]] = [.69, .62]
	glo[aut[i]] = gup[aut[i]]
	z[aut[i]] = [4.14, 4.5]


	#upper limits UDF+mosaic
	g[aut[i]] = [.76, .77, .74]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.16, 3.9, 4.88]




	#new upper limits UDF+mosaic with 2020 catalogs and simulations
	
	
	
	# for 6<r<12
	#zprob = [[2.9, 3.3], [3.3, 4.3], [4.3, 5.3], [5.3, 6.6]]
	g[aut[i]] = [.93, .75, .71]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.09, 3.71, 4.76]

	# for 6<r<12
	#zprob = [[2.9, 3.4], [3.4, 4.4], [4.4, 5.4], [5.4, 6.6]]
	g[aut[i]] = [.79, .91, .66]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.15, 3.79, 4.82]

	#zprob = [[2.9, 3.5], [3.5, 4.5], [4.5, 5.5], [5.5, 6.6]]
	g[aut[i]] = [.67, 1.28, 1.08]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.23, 3.93, 4.91]


	#choosing best
	
	#[[2.9, 3.5], [3, 4], [3.3, 4.3], [3.5, 4.5], [3.7, 4.7], [4, 5], [4.4, 5.4], [4.5, 5.5]
	
	g[aut[i]] = [.67, .51, .75, 1.28, 1.30, 1.2, .66, 1.08]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.23, 3.46, 3.71, 3.93, 4.14, 4.5, 4.82, 4.91]
	
	# choosing best
	
	# [[2.9, 3.4], [3.4, 3.8], [3.8, 4.6], [4.6, 5.6]
	
	g[aut[i]] = [.79, .77, 1.67, .98]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]]) * .3)
	z[aut[i]] = [3.15, 3.59, 4.26, 5.02]

	#best with no galaxy masking
	g[aut[i]] = [.71, .80, .77]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]]) * .3)
	z[aut[i]] = [3.15, 3.91, 5.02]
	zrange = [[2.9, 3.4], [3.4, 4.6], [4.6, 5.6]]
	zerr = np.array([[abs(zz-zr[0]), abs(zz-zr[1])] for zz, zr in zip(z[aut[i]], zrange)]).T



	#best with galaxy masking
	g[aut[i]] = [.98, 1., .87]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]]) * .3)
	z[aut[i]] = [3.15, 3.91, 5.02]
	zrange = [[2.9, 3.4], [3.4, 4.6], [4.6, 5.6]]
	zerr = np.array([[abs(zz-zr[0]), abs(zz-zr[1])] for zz, zr in zip(z[aut[i]], zrange)]).T


	# for 6<r<20
	#zprob = [[2.9, 3.4], [3.4, 4.4], [4.4, 5.4], [5.4, 6.6]]
	g[aut[i]] = [.7, .97, .76]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.15, 3.79, 4.82]
	zrange = [[2.9, 3.4], [3.4, 4.4], [4.4, 5.4]]
	zerr = np.array([[abs(zz-zr[0]), abs(zz-zr[1])] for zz, zr in zip(z[aut[i]], zrange)]).T


	# for 6<r<20
	g[aut[i]] = [.71, .87, .62]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.15, 3.84, 4.91]
	zrange = [[2.9, 3.4], [3.4, 4.5], [4.5, 5.5]]
	zerr = np.array([[abs(zz-zr[0]), abs(zz-zr[1])] for zz, zr in zip(z[aut[i]], zrange)]).T

	#Desired redshift 3.1 mean redshift(s) {'UDF': 3.1023664820934016, 'mosaic': 3.0864228705083074} ngal {'UDF': 26, 'mosaic': 109}
# gamma 0.71

	# for 6<r<20 with mask2d
	g[aut[i]] = [.69, .91, .54]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.15, 3.84, 4.91]
	zrange = [[2.9, 3.4], [3.4, 4.5], [4.5, 5.5]]
	zerr = np.array([[abs(zz-zr[0]), abs(zz-zr[1])] for zz, zr in zip(z[aut[i]], zrange)]).T

	# for 6<r<20 with mask2d, removing protos
	g[aut[i]] = [.69, 1.1, .54]
	nprob = len(g[aut[i]])
	gup[aut[i]] = np.zeros(nprob)
	glo[aut[i]] = (np.array(g[aut[i]])*.3)
	z[aut[i]] = [3.15, 3.84, 4.91]
	zrange = [[2.9, 3.4], [3.4, 4.5], [4.5, 5.5]]
	zerr = np.array([[abs(zz-zr[0]), abs(zz-zr[1])] for zz, zr in zip(z[aut[i]], zrange)]).T


	#Do plot

	lss = '--', '-.', ':'
	colors = 'chocolate', 'orange', 'lime', 'blue', 'yellow', 'purple', 'teal', 'cyan', 'brown'
	m = "H", "d", "s", "^", "p", "o", "h", "D", "d"

	fig, ax = plt.subplots(1)
	naut = len(aut)
	p = [[]]*naut
	for i in 0,1,2:
		print aut[i]
		p[i] = ax.plot(z[aut[i]], g[aut[i]], label=aut[i], ls=lss[i], color='grey', zorder=i)
	i = 3
	print aut[i]
	p[i] = ax.errorbar(z[aut[i]], g[aut[i]], label=aut[i], yerr=[glo[aut[i]], gup[aut[i]]], color='orchid', fmt="none", zorder=i, lw=4)

	imin, imax = 4, 12
	for i in range(imin, imax):
		zord = i
		print aut[i], zord
		plt.errorbar(z[aut[i]], g[aut[i]], yerr=[glo[aut[i]], gup[aut[i]]], color='black', xerr=0, fmt="none", capsize=4)
		p[i] = ax.scatter(z[aut[i]], g[aut[i]], label=aut[i], color=colors[i-imin], zorder=zord+1, marker=m[i-imin], edgecolor='black', s=50)
	
	
	
	#This work
	i = naut-1
	print aut[i]
	p[i] = ax.scatter(z[aut[i]], g[aut[i]], color='red', zorder=14, marker='v', s=100, label=aut[i])
	plt.errorbar(z[aut[i]], g[aut[i]], yerr=glo[aut[i]], xerr=zerr, uplims=np.zeros(nprob)+1, color='red', ls="none", zorder=14)

	if 0:#detections
		
		p[i] = ax.scatter([3.84], [1.2], color='red', zorder=14, marker='*', s=200)
		plt.errorbar([3.84], [1.2], yerr=[.45], color='red', ls="none", zorder=14, xerr=0, fmt="none", capsize=4)

	plt.xlim(-.05, 6.1)
	plt.ylim(.02, 2)

	h,l = ax.get_legend_handles_labels()
	#h = [h[13], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11], h[12]]
	#l = [l[13], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12]]
	h0 = h[:3]
	l0 = l[:3]
	h1, l1 = [], []
	h1.append(h[naut-1])
	l1.append(l[naut - 1])
	for i in range(3, naut-1):
		h1.append(h[i])
		l1.append(l[i])
	#h = [h[12], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11]]
	#l = [l[12], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11]]

	plt.legend(h, l)
	leg1 = ax.legend(h0, l0, loc=[.435, .02], title='UVB Models')
	leg2 = ax.legend(h1, l1, loc=[.2, 0.02])
	ax.add_artist(leg1)
	plt.semilogy()
	plt.yticks([.1, 1], ['0.1', '1.0'])
	plt.xlabel('z')
	plt.ylabel(r'$\Gamma_{\rm{HI}}\,[10^{-12}\rm{s^{-1}}]$')
	#plt.show()
	ext='pdf'
	plt.savefig(fin+'UVB.'+ext, ext=ext, pdi=100)
	plt.legend(h1[5:], l1[5:])
	plt.xlim(2, 6.1)
	#plt.yscale('linear')
	a = [.1, .2, .3, .5, .7, 1., 1.4, 2.]
	b = ['%.1f' % aa for aa in a]
	plt.yticks(a, b)
	plt.ylim(.087, 2)
	plt.savefig(fin+'UVB2.'+ext, ext=ext, pdi=100)
	plt.close()

if 0:
	from tools_sofi import UVBcalc
	import scipy.interpolate as inter

	z = np.arange(3,6,.1)
	div = []
	for zz in z:
		a=UVBcalc(zz)
		div = a[2] / a[0]/3.84e-22
		print zz, div

if 0:
	z2div = inter.interp1d(z, div)

	#T = np.arange(1, 5, .1)
	def e(x): return .686-.106*np.log10(x)-.009*x**-.44


if 0:
	z = [3.2, 3.9, 4.9]
	flls = [.05, .1, .25]
	gf = [4.1, 6.1, 11]
	plt.figure()
	plt.scatter(z, flls)
	plt.xlabel('redshift')
	plt.ylabel(r'$\rm{f_{LLS}}$ upper limit')
	plt.savefig('../../Figures/flls_uplim.pdf', fmt='pdf')
	plt.close()

	plt.figure(num=None, figsize=(5,3.5))
	plt.scatter(z, gf, marker='v', color='red')
	plt.yticks([3,6,9,12])
	plt.xticks([3,3.5,4,4.5,5])
	plt.xlabel('redshift')
	plt.ylabel(r'$\rm{\Gamma_{NHI}}\times\rm{f_{LLS}}$ upper limit [$10^{-12}$ s$^{-1}$]')
	plt.tight_layout()
	plt.savefig('../../Figures/gammaflls.pdf', fmt='pdf')
	plt.close()

