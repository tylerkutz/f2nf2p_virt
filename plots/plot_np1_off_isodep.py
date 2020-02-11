import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from scipy.stats import norm
import sys

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
for fi in sys.argv[1:]:
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
			if ctr2 == 7:
				ctr2 = 0
				ctr1 += 1
				cov.append(arr)
				arr = []

	COVs.append(cov)

PARs = []
for fi in sys.argv[1:]:
	if '.py' in fi: continue
	if '.txt' in fi: continue
	np_a = 0;
	np_b = 0;
	np_c = 0;
	of_p = 0;
	of_n = 0;
	#of_a0 = 0;
	#of_a1 = 0;
	N_he3 = 0;
	N_h3 = 0;
	with open(fi,"rb") as f:
		for line in f:
			if 'np_a' in line and "+/-" in line: np_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_b' in line and "+/-" in line: np_b = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_c' in line and "+/-" in line: np_c = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_p' in line and "+/-" in line: of_p = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_n' in line and "+/-" in line: of_n = float(line.strip().split("=")[-1].split("+/-")[0])
			#if 'of_a0' in line and "+/-" in line: of_a0 = float(line.strip().split("=")[-1].split("+/-")[0])
			#if 'of_a1' in line and "+/-" in line: of_a1 = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'N_he3'in line and "+/-" in line: N_he3= float(line.strip().split("=")[-1].split("+/-")[0])
			if 'N_h3' in line and "+/-" in line: N_h3 = float(line.strip().split("=")[-1].split("+/-")[0])
			
	#PARs.append( [ np_a, np_b, np_c, np_d, of_a ] )
	PARs.append( [ np_a, np_b, np_c, of_p , of_n , N_he3 , N_h3 ] )
	#PARs.append( [ np_a, np_b, np_c, of_a0 , of_a1 , N_he3 , N_h3 ] )


widths = [[],[],[]]

# Figure 1
fig, ax = plt.subplots(1, 1, figsize=(19, 11))
# Create covariance matrix just for F2n/F2p(x=1) vs offshell parameter
ii = 0
cols = ['orange','red','blue']
labs = ['He-3 and H-3','H-3','He-3']
y_min = 1e4; x_min = 1e4;
y_max = -1e4; x_max = -1e4;
for cov in COVs:
	# cov[0] = np_a
	# cov[1] = np_b
	# cov[2] = np_c
	# cov[3] = of_p
	# cov[4] = of_n
	cov_c_c = cov[2][2]
	cov_c_o = cov[2][3] - cov[2][4]
	cov_o_c = cov[3][2] - cov[4][2]
	cov_o_o = cov[3][3] - cov[3][4] - cov[4][3] + cov[4][4]
	#cov_c_o = cov[2][3] 
	#cov_o_c = cov[3][2] 
	#cov_o_o = cov[3][3] 
	new_cov = np.asarray([ [cov_c_c,cov_c_o],[cov_o_c,cov_o_o] ])	
	x = PARs[ii][2]
	y = PARs[ii][3] - PARs[ii][4]
	#y = PARs[ii][3] 

	[tmx,tmy,pear] = confidence_ellipse(new_cov, x , y , ax ,n_std=1, edgecolor=cols[ii],linewidth=3)
	widths[ii].append(tmx)
	[tmx,tmy,pear] = confidence_ellipse(new_cov, x , y , ax ,n_std=2, edgecolor=cols[ii],linewidth=3)
	widths[ii].append(tmx)
	[X,Y,pear] = confidence_ellipse(new_cov, x , y , ax ,n_std=3, edgecolor=cols[ii],linewidth=3,label=labs[ii])
	widths[ii].append(X)


	#if x_min > (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) - x	: x_min = (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) - x*1.5
	#if y_min > (PARs[ii][4]) - y				: y_min = (PARs[ii][4]) - y*1.5
	#if x_max < (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) + x 	: x_max = (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) + x*1.5
	#if y_max < (PARs[ii][4]) + y				: y_max = (PARs[ii][4]) + y*1.5
	if x_min > (x) - X	: x_min = (x) - X*1.5
	if y_min > (y) - Y	: y_min = (y) - Y*1.5
	if x_max < (x) + X 	: x_max = (x) + X*1.5
	if y_max < (y) + Y	: y_max = (y) + Y*1.5
	ii+=1

plt.scatter([0.436479],[0],color='black',s=50,zorder=99,label='Isoscalar Result')
ax.set_ylim(y_min,y_max)
ax.set_xlim(x_min,x_max)
ax.set_xlabel('F2n/F2p at x=1',fontsize=16)
ax.set_ylabel('Proton offshell - Neutron offshell',fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=16)
#ax.set_ylabel('Offshell a0')
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.tight_layout()
plt.savefig("npx1_off_diff.pdf",bbox_inches="tight")

# Figure 2
fig2, ax2 = plt.subplots(1, 1, figsize=(19, 11))
# Create covariance matrix just for F2n/F2p(x=1) vs offshell parameter
ii = 0
cols = ['orange','red','blue']
labs = ['He-3 and H-3','H-3','He-3']
y_min = 1e4; x_min = 1e4;
y_max = -1e4; x_max = -1e4;
for cov in COVs:
	# cov[0] = np_a
	# cov[1] = np_b
	# cov[2] = np_c
	# cov[3] = of_p
	# cov[4] = of_n
	cov_c_c = cov[2][2]
	cov_c_o = cov[2][3] + cov[2][4]
	cov_o_c = cov[3][2] + cov[4][2]
	cov_o_o = cov[3][3] + cov[3][4] + cov[4][3] + cov[4][4]
	#cov_c_o = cov[2][4]
	#cov_o_c = cov[4][2]
	#cov_o_o = cov[4][4]
	new_cov = np.asarray([ [cov_c_c,cov_c_o],[cov_o_c,cov_o_o] ])	
	x = PARs[ii][2]
	y = PARs[ii][3] + PARs[ii][4]
	#y = PARs[ii][4]

	[tm,tm,pear] = confidence_ellipse(new_cov, x , y , ax2 ,n_std=1, edgecolor=cols[ii],linewidth=3)
	[tm,tm,pear] = confidence_ellipse(new_cov, x , y , ax2 ,n_std=2, edgecolor=cols[ii],linewidth=3)
	[X,Y,pear] = confidence_ellipse(new_cov, x , y , ax2 ,n_std=3, edgecolor=cols[ii],linewidth=3,label=labs[ii])


	#if x_min > (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) - x	: x_min = (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) - x*1.5
	#if y_min > (PARs[ii][4]) - y				: y_min = (PARs[ii][4]) - y*1.5
	#if x_max < (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) + x 	: x_max = (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) + x*1.5
	#if y_max < (PARs[ii][4]) + y				: y_max = (PARs[ii][4]) + y*1.5
	if x_min > (x) - X	: x_min = (x) - X*1.5
	if y_min > (y) - Y	: y_min = (y) - Y*1.5
	if x_max < (x) + X 	: x_max = (x) + X*1.5
	if y_max < (y) + Y	: y_max = (y) + Y*1.5
	ii+=1

plt.scatter([0.436479],[2.*(-1.25864)],color='black',s=50,zorder=99,label='Isoscalar Result')
ax2.set_ylim(y_min,y_max)
ax2.set_xlim(x_min,x_max)
ax2.set_xlabel('F2n/F2p at x=1',fontsize=16)
ax2.set_ylabel('Proton offshell + Neutron offshell',fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
#ax2.set_ylabel('Offshell a1')
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.tight_layout()
plt.savefig("npx1_off_sum.pdf",bbox_inches="tight")
#plt.savefig("npx1_off_a1.pdf",bbox_inches="tight")

'''
# Figure 3
def gaussian(x, mu, sig):
	return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
plt.figure(3)
xs = np.arange(0,1,0.001)
on = np.ones(len(xs))
for i in range(3):
	ys = gaussian( xs, PARs[i][2], widths[i][0])
	plt.plot(xs,ys,color=cols[i],linewidth=3,label=labs[i])
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.savefig("npx1_offLin_1d.pdf",bbox_inches="tight")
'''

#plt.savefig("npx1_off.pdf",bbox_inches="tight")
plt.show()
