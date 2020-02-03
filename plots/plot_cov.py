import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from scipy.stats import norm
import sys

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
def confidence_ellipse(cov,mean_x,mean_y, ax, n_std=3.0, facecolor='none', **kwargs):
	pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
	# Using a special case to obtain the eigenvalues of this
	# two-dimensionl dataset.
	ell_radius_x = np.sqrt(1 + pearson)
	ell_radius_y = np.sqrt(1 - pearson)
	ellipse = Ellipse((0, 0),
		width=ell_radius_x * 2,
		height=ell_radius_y * 2,
		facecolor=facecolor,
		**kwargs)

	# Calculating the stdandard deviation of x from
	# the squareroot of the cov[i][i] and multiplying
	# with the given number of standard deviations.
	scale_x = np.sqrt(cov[0, 0]) * n_std

	# calculating the stdandard deviation of y ...
	scale_y = np.sqrt(cov[1, 1]) * n_std

	transf = transforms.Affine2D() \
		 .rotate_deg(45) \
		 .scale(scale_x, scale_y) \
		 .translate(mean_x, mean_y)

	ellipse.set_transform(transf + ax.transData)

	ax.add_patch(ellipse)

	#return ax.add_patch(ellipse)
	#return ell_radius_x,ell_radius_y
	return scale_x,scale_y,pearson

if len(sys.argv) != 7:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Cov-He3+H3] [Cov-H3] [Cov-He3] [Min-He3+H3] [Min-H3] [Min-He3]"
	exit(-1)

COVs = []
for fi in sys.argv[1:4]:
	if '.py' in fi: continue
	cov = []
	ctr1 = 0
	ctr2 = 0
	with open(fi,"rb") as f:
		arr = []
		for i in range(2): next(f)
		for line in f:
			if '*' in line: continue
			val = float(line.strip())
			arr.append(val)
			ctr2+=1
			if ctr2 == 7:
				ctr2 = 0
				ctr1 += 1
				cov.append(arr)
				arr = []

	COVs.append(cov)

PARs = []
for fi in sys.argv[4:]:
	if '.py' in fi: continue
	np_a = 0;
	np_b = 0;
	np_c = 0;
	#np_d = 0;
	of_a = 0;
	of_b = 0;
	N_he3 = 0;
	N_h3 = 0;
	with open(fi,"rb") as f:
		for line in f:
			if 'np_a' in line and "+/-" in line: np_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_b' in line and "+/-" in line: np_b = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_c' in line and "+/-" in line: np_c = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_a' in line and "+/-" in line: of_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_b' in line and "+/-" in line: of_b = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'N_he3'in line and "+/-" in line: N_he3= float(line.strip().split("=")[-1].split("+/-")[0])
			if 'N_h3' in line and "+/-" in line: N_h3 = float(line.strip().split("=")[-1].split("+/-")[0])
			
	PARs.append( [ np_a, np_b, np_c, of_a , of_b , N_he3 , N_h3 ] )


fig, axs = plt.subplots(7, 7, figsize=(19, 11))
fig.subplots_adjust(hspace=.8)
fig.subplots_adjust(wspace=.8)


labels = ["np_a","np_b","np_c","of_a","of_b","N_he3","N_h3"]
ii = 0
for cov in COVs:
	cols = ['orange','red','blue']
	labs = ['He-3 and H-3','H-3','He-3']
	for i in range(len(cov)):
		for j in range(len(cov[i])):
			print labs[ii],i,j
			if j < i : 
				axs[j][i].set_visible(False)
				continue
			if j == i: 
				if cov[i][i] == 0 : continue
				mean = PARs[ii][i];
				#mean = min_pars_he3_h3[i];
				std = np.sqrt(cov[i][i])
				x_axis = np.linspace(mean-10*std,mean+10*std,100)
				axs[j][i].plot(x_axis, gaussian(x_axis,mean,std),color=cols[ii])

				axs[i][j].set_title(labels[j])
				axs[j][i].set_xlim(mean-4*std,mean+4*std)
				axs[j][i].set_ylim(0,1)
				
				start, end = axs[j][i].get_xlim()
				axs[j][i].xaxis.set_ticks(np.round(np.linspace(start, end, 3),2))
				start, end = axs[j][i].get_ylim()
				axs[j][i].yaxis.set_ticks(np.round(np.linspace(start, end, 3),2))

				continue
			if cov[i][j] == 0 or  cov[j][j] == 0: continue
			subcov = np.asarray([	[ cov[i][i] , cov[i][j] ], [ cov[j][i] , cov[j][j] ]	])
			X = PARs[ii][i]	
			Y = PARs[ii][j]	

			[x,y,pear] = confidence_ellipse(subcov, X , Y , axs[j][i] , edgecolor=cols[ii])
			axs[j][i].set_xlim(X-1.5*x,X+1.5*x)
			axs[j][i].set_ylim(Y-1.5*y,Y+1.5*y)
			axs[j][i].tick_params(axis='both', which='major', labelsize=10)
			start, end = axs[j][i].get_xlim()
			axs[j][i].xaxis.set_ticks(np.round(np.linspace(start, end, 3),2))
			start, end = axs[j][i].get_ylim()
			axs[j][i].yaxis.set_ticks(np.round(np.linspace(start, end, 3),2))
			axs[j][i].set_xlabel(labels[i],fontsize=16)
			axs[j][i].set_ylabel(labels[j],fontsize=16)
			axs[j][i].grid()
			axs[j][i].set_title("R: %.2f" % np.round(pear,2) ,fontsize=12)

			#if x_min[j][i] > (PARs[ii][2]) - x		: x_min[j][i] = (PARs[ii][2]) - x*1.5
			#if y_min[j][i] > (PARs[ii][3]) - y		: y_min[j][i] = (PARs[ii][3]) - y*1.5
			#if x_max[j][i] < (PARs[ii][2]) + x		: x_max[j][i] = (PARs[ii][2]) + x*1.5
			#if y_max[j][i] < (PARs[ii][3]) + y		: y_max[j][i] = (PARs[ii][3]) + y*1.5
			#axs[j][i].set_ylim(y_min[j][i],y_max[j][i])
			#axs[j][i].set_xlim(x_min[j][i],x_max[j][i])
	ii+=1
plt.savefig("correlation.pdf")
plt.show()



