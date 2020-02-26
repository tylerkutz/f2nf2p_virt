import numpy as np
import matplotlib.pyplot as plt
from matplotlib  import cm
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from scipy.stats import norm
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys


if len(sys.argv) != 5:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Plot opt] [Min-He3+H3] [Min-H3] [Min-He3]"
	print "\t [0 = const offshell]"
	print "\t [1 = linear offshell]"
	print "\t [2 = isodep offshell]"
	exit(-1)

DIM = -1
OPT = -1
if int(sys.argv[1]) == 0:
	DIM = 6
elif int(sys.argv[1]) == 1:
	DIM = 7
	OPT = 0
elif int(sys.argv[1]) == 2:
	DIM = 7
	OPT = 1


PARs = []
for fi in sys.argv[4:5]:
	if '.py' in fi: continue
	if '.txt' in fi: continue
	np_a = 0;
	np_b = 0;
	np_c = 0;
	of_a = 0;
	of_b = 0;
	N_he3 = 0;
	N_h3 = 0;
	chi2 = 0;
	read_chi2 = False
	with open(fi,"rb") as f:
		for line in f:
			#if len(PARs) == 500: break
			if read_chi2:
				chi2 = float(line.strip())
				if DIM == 6 and chi2 < 50:
					PARs.append( [ np_a, np_b, np_c, of_a ,  N_he3 , N_h3 , chi2 ] )
				elif DIM == 7 and chi2 < 50:
					PARs.append( [ np_a, np_b, np_c, of_a , of_b , N_he3 , N_h3 , chi2 ] )
				read_chi2 = False
				continue
			if '+/-' in line or 'fixed' in line: continue
			if 'Current chi2' in line:
				read_chi2 = True
				continue
			if 'np_a' in line: np_a = float(line.strip().split(":")[-1])
			if 'np_b' in line: np_b = float(line.strip().split(":")[-1])
			if 'np_c' in line: np_c = float(line.strip().split(":")[-1])
			if DIM == 6:
				if 'of_a' in line: of_a = float(line.strip().split(":")[-1])
			elif DIM == 7 and OPT == 0:
				if 'of_a' in line: of_a = float(line.strip().split(":")[-1])
				if 'of_b' in line: of_b = float(line.strip().split(":")[-1])
			elif DIM == 7 and OPT == 1:
				if 'of_p' in line: of_a = float(line.strip().split(":")[-1])
				if 'of_n' in line: of_b = float(line.strip().split(":")[-1])
			if 'N_he3'in line: N_he3= float(line.strip().split(":")[-1])
			if 'N_h3' in line: N_h3 = float(line.strip().split(":")[-1])

PARs = np.asarray(PARs)


fig, axs = plt.subplots(DIM, DIM, figsize=(19, 18))
fig.subplots_adjust(hspace=.5)
fig.subplots_adjust(wspace=.5)


widths = [[],[],[]]
if DIM == 6:
	labels = ["np_a","np_b","np_c","of_a","N_he3","N_h3"]
elif DIM == 7 and OPT == 0:
	labels = ["np_a","np_b","np_c","of_a","of_b","N_he3","N_h3"]
elif DIM == 7 and OPT == 1:
	labels = ["np_a","np_b","np_c","of_p","of_n","N_he3","N_h3"]

'''
for i in range(DIM):
	for j in range(DIM):
		
		if j < i : 
			axs[j][i].set_visible(False)
			continue
		if j == i: 
			axs[j][i].set_visible(False)
			continue


		x = PARs[:,i]
		y = PARs[:,j]
		z = PARs[:,-1]
		xi = np.linspace(	min(x),	max(x), 200 )
		yi = np.linspace(	min(y),	max(y), 200 )
		#xi = np.linspace(	np.percentile(x,34),	np.percentile(x,84), 100 )
		#yi = np.linspace(	np.percentile(y,34),	np.percentile(y,84), 100 )
		zi = griddata( (x,y), z, (xi[None,:], yi[:,None])	 ,method='linear')
		

		CS = axs[j][i].contourf(	xi,yi,zi,15,cmap=plt.cm.jet )
		cbar = fig.colorbar(CS, ax=axs[j][i])
		#cbar.ax.set_yticklabels(['<0.5', '1', '>1.5'])

		#axs[j][i].scatter(PARs[:,j],PARs[:,i],s=20,c=PARs[:,-1], marker = 'o', cmap = cm.jet );
		axs[j][i].tick_params(axis='both', which='major', labelsize=10)
		start, end = axs[j][i].get_xlim()
		axs[j][i].xaxis.set_ticks(np.round(np.linspace(start, end, 2),4))
		start, end = axs[j][i].get_ylim()
		axs[j][i].yaxis.set_ticks(np.round(np.linspace(start, end, 2),4))
		axs[j][i].set_xlabel(labels[i],fontsize=16)
		axs[j][i].set_ylabel(labels[j],fontsize=16)
		axs[j][i].grid(True)

'''
#plt.savefig('test.pdf')

plt.figure(2)
x = PARs[:,2]
y = PARs[:,3] - PARs[:,4]
z = PARs[:,-1]
xi = np.linspace(       min(x), max(x), 200 )
yi = np.linspace(       min(y), max(y), 200 )
zi = griddata( (x,y), z, (xi[None,:], yi[:,None])        ,method='linear')
plt.contourf( xi,yi,zi,15,cmap=plt.cm.jet )
plt.colorbar()
plt.show()



#plt.show()
