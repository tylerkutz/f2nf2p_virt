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



cov = []
ctr1 = 0
ctr2 = 0
DIM = 7
OPT = 1
with open('../output/he-3-h-3-newframework-isodep.txt',"rb") as f:
	arr = []
	for i in range(2): next(f)
	for line in f:
		if '*' in line: continue
		val = float(line.strip().split(" ")[0])
		arr.append(val)
		ctr2+=1
		if ctr2 == DIM:
			ctr2 = 0
			ctr1 += 1
			cov.append(arr)
			arr = []
np_a = 0;
np_b = 0;
np_c = 0;
of_a = 0;
of_b = 0;
N_he3 = 0;
N_h3 = 0;
with open('../output/he-3-h-3-newframework-isodep-minimization.out',"rb") as f:
	for line in f:
		if 'np_a' in line and "+/-" in line: np_a = float(line.strip().split("=")[-1].split("+/-")[0])
		if 'np_b' in line and "+/-" in line: np_b = float(line.strip().split("=")[-1].split("+/-")[0])
		if 'np_c' in line and "+/-" in line: np_c = float(line.strip().split("=")[-1].split("+/-")[0])
		if DIM == 6:
			if 'of_a' in line and "+/-" in line: of_a = float(line.strip().split("=")[-1].split("+/-")[0])
		elif DIM == 7 and OPT == 0:
			if 'of_a' in line and "+/-" in line: of_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_b' in line and "+/-" in line: of_b = float(line.strip().split("=")[-1].split("+/-")[0])
		elif DIM == 7 and OPT == 1:
			if 'of_p' in line and "+/-" in line: of_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_n' in line and "+/-" in line: of_b = float(line.strip().split("=")[-1].split("+/-")[0])
		if 'N_he3'in line and "+/-" in line: N_he3= float(line.strip().split("=")[-1].split("+/-")[0])
		if 'N_h3' in line and "+/-" in line: N_h3 = float(line.strip().split("=")[-1].split("+/-")[0])
PARs = [ np_a, np_b, np_c, of_a , of_b , N_he3 , N_h3 ]

fig, axs = plt.subplots(1, 1)
for i in range(len(cov)):
	for j in range(len(cov[i])):
		if j != 4: continue
		if i != 3: continue
		if j < i : 
			axs[j][i].set_visible(False)
			continue
		if j == i: 
			axs[j][i].set_visible(False)
			continue
		if cov[i][j] == 0 or  cov[j][j] == 0: continue
		print i,j
		X = PARs[i]
		Y = PARs[j]
		subcov = np.asarray([	[ cov[i][i] , cov[i][j] ], [ cov[j][i] , cov[j][j] ]	])
		print subcov
		[x,y,pear] = confidence_ellipse(subcov, X , Y , axs , edgecolor='Red')
		axs.set_xlim(X-1.5*x,X+1.5*x)
		axs.set_ylim(Y-1.5*y,Y+1.5*y)

A = []; B = []
with open("../build/error/contour-test-Up14-isodep.txt","rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		A.append( float(parse[0]) )
		B.append( float(parse[1]) )
axs.plot( A, B ,color='blue',linewidth=3)
axs.scatter( [-1.25864],[-1.25864],color='black',s=40)

plt.savefig('testing-contour.pdf')
plt.show()
