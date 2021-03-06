import numpy as np

import params

fin = '../../UVB/'

l, h, w, z = {}, {}, {}, {}
aut = 'Crighton+19 fit', 'Prochaska+10 fit', 'Ribaudo+11', 'Fumagalli+13', 'Crighton+20', 'Rahmati+15', 'This work'


def dndz(z, ls=1.46, zs=3., a=1.7):
	return ls*((1+z)/(1+zs))**a

# Data


########################################
##Crighton+19 fit... check this!!!!!!!!!
########################################
i = 0
z[aut[i]] = np.arange(0, 7, .01)
l[aut[i]] = dndz(z[aut[i]], 1.46, 3., 1.7)

########################################
##Prochaska+10 fit... check this!!!!!!!!!
########################################
i += 1
z[aut[i]] = np.arange(3, 6, .01)
l[aut[i]] = dndz(z[aut[i]], 1.9, 3.7, 5.1)


# HST Ribaudo+11
i += 1
l[aut[i]] = [0.37,
			 0.52,
			 0.64,
			 0.86,
			 1.80,
			 1.99,
			 2.57]
h[aut[i]] = [[0.1, 0.1],
			 [0.14, 0.14],
			 [0.21, 0.21],
			 [0.3, 0.3],
			 [0.26, 0.26],
			 [0.29, 0.29],
			 [0.44, 0.44]]
w[aut[i]] = [[0.242, 1.078],
			 [1.078, 1.544],
			 [1.544, 1.947],
			 [1.947, 2.594],
			 [3.500, 3.663],
			 [3.663, 3.925],
			 [3.925, 4.907]]
z[aut[i]] = [0.895,
			 1.157,
			 1.751,
			 2.257,
			 3.593,
			 3.799,
			 4.125]

# MagE Fumagalli+13
i += 1
l[aut[i]] = [1.21]

h[aut[i]] = [.28]

w[aut[i]] = [[2.6, 3]]
z[aut[i]] = [2.8]

# GGG
i += 1

l[aut[i]] = [2.21,
			 2.2,
			 2.95]

h[aut[i]] = [0.4,
			 0.47,
			 0.58]

w[aut[i]] = [[3.75, 4.4],
			 [4.4, 4.7],
			 [4.7, 5.4]]
z[aut[i]] = [4.075,
			 4.55,
			 5.05]




########################################
##Ramati+15
########################################
i +=1
z[aut[i]] = [1, 2, 3, 4, 5]
l[aut[i]] = [1.02, 1.22, 1.86, 3.26, 5.49]
# rescaled to latest cosmology
l[aut[i]] = [1.08, 1.27, 1.93, 3.39, 5.70]



########################################
##Gallego+19 work... check this!!!!!!!!!
########################################
i +=1
z[aut[i]] = [params.redshifts[s] for s in range(6, 15)]

#No self shielding
l[aut[i]] = [12.38, 9.1, 5.02, 3.67, 2.77, 1.52, 1.22]

#Self shielding
l[aut[i]] = [31.6, 19.35, 13.4, 7.5, 4.99, 3.08, 2.03]

#Self shielding Gamma x 10
l[aut[i]] = [10.7, 5.9, 7.6, 3.2, 2.9, 1.9, 1.6]


#Self shielding
l[aut[i]] = [20.1, 13.9, 10.1, 7.3, 5.1, 3.7]

#Case A, Gamma x 2
l[aut[i]] = [11.7, 8.3, 6.2, 4.3, 2.9, 2.]


#Case A, Gamma x 1, new values
l[aut[i]] = [11.7, 2.9, 1.7, 1.2, 0.8, 2.]

# variable ssthr
l[aut[i]] = [6.39, 4.51, 2.81, 1.95, 1.56, 1.14, 0.89]

#fixed ssthr 6.73e-3
l[aut[i]] = [5.52, 4.18, 3.34, 2.52, 2.1, 1.64, 1.42, 0.94, 0.86]



# Do plot3
import matplotlib.pyplot as plt

lss = '--', '-.', ':'
colors = 'chocolate', 'lime', 'blue', 'orange', 'yellow', 'purple', 'teal', 'cyan', 'brown'
m = "H", "d", "s", "^", "p", "o", "h", "D", "d"

fig, ax = plt.subplots(1)
naut = len(aut)
dotsize = 60

for i in [0, 1]:
	ax.plot(z[aut[i]], l[aut[i]], label=aut[i], ls=lss[i], color='grey', zorder=i)

for i in range(2, naut-2):
	zord = i
	print aut[i], zord
	yerr = np.array(h[aut[i]]).T
	xerr = np.abs(np.array(w[aut[i]]).T-z[aut[i]])
	plt.errorbar(z[aut[i]], l[aut[i]], yerr=yerr, xerr=xerr, color='black', fmt="none", capsize=4)
	ax.scatter(z[aut[i]], l[aut[i]], label=aut[i], color=colors[i-2], zorder=zord + 1, marker=m[i],
					  edgecolor='black', s=dotsize)


i=naut-2
plt.errorbar(z[aut[i]], l[aut[i]], color='black', fmt="none", capsize=4)
ax.scatter(z[aut[i]], l[aut[i]], label=aut[i], color=colors[i], zorder=zord + 1, marker=m[i],
				  edgecolor='black', s=dotsize)



i=naut-1
plt.errorbar(z[aut[i]], l[aut[i]], color='black', fmt="none", capsize=4)
ax.scatter(z[aut[i]], l[aut[i]], label=aut[i], color='yellow', zorder=zord + 1, marker=m[i],
				  edgecolor='black', s=dotsize)


fsize=14
plt.semilogy()
plt.xlim(0, 6.2)
plt.ylim(.05, 7)
plt.yticks([.1, .3, 1, 3], ['0.1', '0.3', '1.0', '3.0'], fontsize=fsize)
plt.xticks([0,1,2,3,4,5,6], fontsize=fsize)

plt.legend(loc=4, fontsize=fsize)
plt.xlabel('z', fontsize=fsize)
plt.ylabel(r'$\ell\rm{(z)[\tau>2]}$', fontsize=fsize)
# plt.show()
ext = 'pdf'
plt.savefig('../../Figures/dndz.' + ext, ext=ext, pdi=100)
plt.close()
