import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from scipy.stats import norm
import sys

#def gaussian(x, mu, sig):
#    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
def gaussian(x, mu, sig):
	return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
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

if len(sys.argv) != 8:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [dim of cov] [Cov-He3+H3] [Cov-H3] [Cov-He3] [Min-He3+H3] [Min-H3] [Min-He3]"
	exit(-1)

COVs = []
for fi in sys.argv[2:]:
	if '.py' in fi: continue
	if '.out' in fi: continue
	cov = []
	ctr1 = 0
	ctr2 = 0
	with open(fi,"rb") as f:
		arr = []
		for i in range(2): next(f)
		for line in f:
			if '*' in line: continue
			val = float(line.strip().split(" ")[0])
			arr.append(val)
			ctr2+=1
			if ctr2 == int(sys.argv[1]):
				ctr2 = 0
				ctr1 += 1
				cov.append(arr)
				arr = []

	COVs.append(cov)

PARs = []
for fi in sys.argv[2:]:
	if '.py' in fi: continue
	if '.txt' in fi: continue
	np_a = 0;
	np_b = 0;
	np_c = 0;
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
	#PARs.append( [ np_a, np_b, np_c, of_a ,  N_he3 , N_h3 ] )


fig, axs = plt.subplots(int(sys.argv[1]), int(sys.argv[1]), figsize=(19, 18))
fig.subplots_adjust(hspace=.8)
fig.subplots_adjust(wspace=.8)


widths = [[],[],[]]
labels = ["np_a","np_b","np_c","of_a","of_b","N_he3","N_h3"]
#labels = ["np_a","np_b","np_c","of_a","N_he3","N_h3"]
ii = 0

x_min = []
x_max = []
for i in range(len(COVs[0])):
	this_min = []
	this_max = []
	for j in range(len(COVs[0][i])):
		this_min.append(1e4)
		this_max.append(-1e4)
	x_min.append(this_min)
	x_max.append(this_max)
y_min = []
y_max = []
for i in range(len(COVs[0])):
	this_min = []
	this_max = []
	for j in range(len(COVs[0][i])):
		this_min.append(1e4)
		this_max.append(-1e4)
	y_min.append(this_min)
	y_max.append(this_max)

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
				std = np.sqrt(cov[i][i])
				widths[ii].append(std)
				x_axis = np.linspace(mean-10*std,mean+10*std,100)
				axs[j][i].plot(x_axis, gaussian(x_axis,mean,std),color=cols[ii])
				
				axs[i][j].set_title(labels[j])
				axs[j][i].set_xlim(mean-4*std,mean+4*std)
				#axs[j][i].set_ylim(0,1)
				
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
			print "\t",X-1.5*x,X+1.5*x
			print "\t",Y-1.5*y,Y+1.5*y

			if y_min[j][i] > (Y-1.5*y): y_min[j][i] = Y-1.5*y
			if y_max[j][i] < (Y+1.5*y): y_max[j][i] = Y+1.5*y
			if x_min[j][i] > (X-1.5*x): x_min[j][i] = X-1.5*x
			if x_max[j][i] < (X+1.5*x): x_max[j][i] = X+1.5*x
			if ii == 0:
				axs[j][i].set_title("R: %.2f" % np.round(pear,2) ,fontsize=12)

	ii+=1
for i in range(len(COVs[0])):
	for j in range(len(COVs[0][i])):
			if j < i : continue
			if j == i: continue
			axs[j][i].set_xlim(x_min[j][i],x_max[j][i])
			axs[j][i].set_ylim(y_min[j][i],y_max[j][i])
			axs[j][i].tick_params(axis='both', which='major', labelsize=10)
			start, end = axs[j][i].get_xlim()
			axs[j][i].xaxis.set_ticks(np.round(np.linspace(start, end, 3),2))
			start, end = axs[j][i].get_ylim()
			axs[j][i].yaxis.set_ticks(np.round(np.linspace(start, end, 3),2))
			axs[j][i].set_xlabel(labels[i],fontsize=16)
			axs[j][i].set_ylabel(labels[j],fontsize=16)
			axs[j][i].grid()



plt.savefig("correlation-linearOff.pdf")


# Do 1D plot of just the c parameter:
plt.figure(3)
xs = np.arange(0,1,0.001)
on = np.ones(len(xs))
for i in range(3):
	std = widths[i][2]
	mean = PARs[i][2]
	x_axis = np.linspace(mean-10*std,mean+10*std,100)
	plt.plot(x_axis, gaussian(x_axis,mean,std),color=cols[i],label=labs[i],linewidth=3)
	#ys = gaussian( xs, PARs[i][2], widths[i][0])
	#plt.plot(xs,ys,color=cols[i],linewidth=3,label=labs[i])
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.xlabel('F2n/F2p(x=1)',fontsize=16)
plt.ylabel('PDF',fontsize=16)
plt.xlim([0,1])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('f2nf2p-linearOff-1D.pdf',bbox_inches="tight")


plt.show()



