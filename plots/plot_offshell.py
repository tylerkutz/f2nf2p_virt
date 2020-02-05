
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

if len(sys.argv) < 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Min-He3+H3] [Min-H3] [Min-He3]"
	exit(-1)


PARs = []
for fi in sys.argv[1:]:
	if '.py' in fi: continue
	if '.txt' in fi: continue

	of_a = 0;
	#of_b = 0;
	with open(fi,"rb") as f:
		for line in f:
			if 'of_a' in line and "+/-" in line: of_a = float(line.strip().split("=")[-1].split("+/-")[0])
			#if 'of_b' in line and "+/-" in line: of_b = float(line.strip().split("=")[-1].split("+/-")[0])
	#PARs.append( [ of_a, of_b ] )
	PARs.append( [ of_a ] )

print PARs
dx, dnu = 0.05, 0.05
# generate 2 2d grids for the x & y bounds
nu, x = np.mgrid[slice(-3, 0 + dnu, dnu),
                slice(0, 1 + dx, dx)]


fig, axs = plt.subplots(3,1, figsize=(19, 18))
fig.subplots_adjust(hspace=.8)
fig.subplots_adjust(wspace=.8)

labs=['He-3 + H-3','H-3','He-3']
j=0
for i in range(len(PARs)):
	axs[i].tick_params(axis='both', which='major', labelsize=10)

	#off = 1. + nu*nu*(PARs[i][0] + PARs[i][1]*x)
	off = 1. + nu*nu*(PARs[i][0])
	im = axs[i].contourf(x + dx/2. , nu + dnu/2. , off ,cmap="RdBu")

	divider = make_axes_locatable(axs[i])
	cax = divider.append_axes('right', size='5%', pad=0.05)


	fig.colorbar(im,cax=cax)
	axs[i].set_xlabel('x',fontsize=16)
	axs[i].set_ylabel('nu',fontsize=16)

	#plt.xticks(fontsize=16)
	#plt.yticks(fontsize=16)
	axs[i].set_title('Offshell function for '+labs[i],fontsize=17)

plt.savefig('offshell-linearOff-LC.pdf')
plt.show()
