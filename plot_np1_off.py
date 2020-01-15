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
			if ctr2 == 5:
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
	np_d = 0;
	of_a = 0;
	with open(fi,"rb") as f:
		for line in f:
			if 'np_a' in line and "+/-" in line: np_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_b' in line and "+/-" in line: np_b = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_c' in line and "+/-" in line: np_c = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_d' in line and "+/-" in line: np_d = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_a' in line and "+/-" in line: of_a = float(line.strip().split("=")[-1].split("+/-")[0])
	PARs.append( [ np_a, np_b, np_c, np_d, of_a ] )


print PARs
print len(PARs)
fig, ax = plt.subplots(1, 1, figsize=(19, 11))
# Create covariance matrix just for F2n/F2p(x=1) vs offshell parameter
ii = 0
cols = ['orange','red','blue']
labs = ['He-3 and H-3','H-3','He-3']
y_min = 1e4; x_min = 1e4;
y_max = -1e4; x_max = -1e4;
for cov in COVs:
	cov_abc_o = cov[0][4] + cov[1][4] + cov[2][4]
	cov_o_abc = cov[4][0] + cov[4][1] + cov[4][2]
	cov_abc_abc = 0
	for i in range(3):
		cov_abc_abc += (cov[i][0] +cov[i][1] + cov[i][2] )
	cov_o_o = cov[4][4]
	new_cov = np.asarray([	[cov_abc_abc,cov_abc_o],[cov_o_abc,cov_o_o] ])
	
	x = PARs[ii][0]+PARs[ii][1]+PARs[ii][2]
	y = PARs[ii][4]
	[tm,tm,pear] = confidence_ellipse(new_cov, x , y , ax ,n_std=1, edgecolor=cols[ii],linewidth=3)
	[tm,tm,pear] = confidence_ellipse(new_cov, x , y , ax ,n_std=2, edgecolor=cols[ii],linewidth=3)
	[x,y,pear] = confidence_ellipse(new_cov, x , y , ax ,n_std=3, edgecolor=cols[ii],linewidth=3,label=labs[ii])

	if x_min > (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) - x	: x_min = (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) - x*1.5
	if y_min > (PARs[ii][4]) - y				: y_min = (PARs[ii][4]) - y*1.5
	if x_max < (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) + x 	: x_max = (PARs[ii][0]+PARs[ii][1]+PARs[ii][2]) + x*1.5
	if y_max < (PARs[ii][4]) + y				: y_max = (PARs[ii][4]) + y*1.5
	ii+=1

ax.set_ylim(y_min,y_max)
ax.set_xlim(x_min,x_max)
ax.set_xlabel('F2n/F2p at x=1')
ax.set_ylabel('Offshell Par a')
plt.legend(numpoints=1,loc='best')
plt.savefig("npx1_off.pdf",bbox_inches="tight")
plt.show()
