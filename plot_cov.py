import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from scipy.stats import norm
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

fig, axs = plt.subplots(5, 5, figsize=(19, 11))
fig.subplots_adjust(hspace=.8)
fig.subplots_adjust(wspace=.8)

cov_he3_h3 = []
cov_h3 = []
arr = []
ctr1 = 0
ctr2 = 0
with open("build/test-he3-h3.txt","rb") as f:
	for line in f:
		if '#' in line: continue
		val = float(line.strip())
		arr.append(val)
		ctr2+=1
		if ctr2 == 6:
			ctr2 = 0
			ctr1 += 1
			cov_he3_h3.append(arr)
			arr = []
with open("build/test-no-h3.txt","rb") as f:
	for line in f:
		if '#' in line: continue
		val = float(line.strip())
		arr.append(val)
		ctr2+=1
		if ctr2 == 6:
			ctr2 = 0
			ctr1 += 1
			cov_h3.append(arr)
			arr = []

# He3 and H3 data:
min_pars_he3_h3 = [-1.49198,1.35858,0.694808,1.29278,-2.56522]

# He3 data only:
min_pars_h3 = [-1.28977,0.791046,0.777979,0.905046,1.50513]

init_pars = [-1.21721713,0.8622478,0.82047886,0.96399233]

labels = ["np_a","np_b","np_c","np_d","o_a"]
for i in range(len(cov_he3_h3)):
	for j in range(len(cov_he3_h3[i])):
		if j < i : 
			axs[j][i].set_visible(False)
			continue
		if j == i: 
			mean = min_pars_he3_h3[i];
			std = np.sqrt(cov_he3_h3[i][i])
			x_axis = np.linspace(mean-10*std,mean+10*std,100)
			axs[j][i].plot(x_axis, gaussian(x_axis,mean,std))

			axs[i][j].set_title(labels[j])
			axs[j][i].set_xlim(mean-4*std,mean+4*std)
			axs[j][i].set_ylim(0,1)
			continue
		
		subcov_he3_h3 = np.asarray([	[ cov_he3_h3[i][i] , cov_he3_h3[i][j] ], [ cov_he3_h3[j][i] , cov_he3_h3[j][j] ]	])

		[x,y,pear] = confidence_ellipse(subcov_he3_h3, min_pars_he3_h3[i] ,min_pars_he3_h3[j] , axs[j][i] , edgecolor='red')
		axs[j][i].set_xlim(min_pars_he3_h3[i]-1.5*x,min_pars_he3_h3[i]+1.5*x)
		axs[j][i].set_ylim(min_pars_he3_h3[j]-1.5*y,min_pars_he3_h3[j]+1.5*y)
		axs[j][i].set_xlabel(labels[i])
		axs[j][i].set_ylabel(labels[j])
		axs[j][i].set_title("Pearson R: %.2f" % np.round(pear,2) )
		#axs[i][j].set_xlim(-3*x,3*x)
		#axs[i][j].set_ylim(-3*y,3*y)
		#axs[i][j].set_xlim(-10,10)
plt.savefig("correlation.pdf")



fig2, ax2 = plt.subplots(1, 4, figsize=(19, 11))

cov_he3_h3_abc_o = cov_he3_h3[0][4] + cov_he3_h3[1][4] + cov_he3_h3[2][4]
cov_he3_h3_o_abc = cov_he3_h3[4][0] + cov_he3_h3[4][1] + cov_he3_h3[4][2]
cov_he3_h3_abc_abc = 0
for i in range(3):
	cov_he3_h3_abc_abc += (cov_he3_h3[i][0] +cov_he3_h3[i][1] + cov_he3_h3[i][2] )
cov_he3_h3_o_o = cov_he3_h3[4][4]
new_cov_he3_h3 = np.asarray([	[cov_he3_h3_abc_abc,cov_he3_h3_abc_o],[cov_he3_h3_o_abc,cov_he3_h3_o_o] ])

cov_h3_abc_o = cov_h3[0][4] + cov_h3[1][4] + cov_h3[2][4]
cov_h3_o_abc = cov_h3[4][0] + cov_h3[4][1] + cov_h3[4][2]
cov_h3_abc_abc = 0
for i in range(3):
	cov_h3_abc_abc += (cov_h3[i][0] +cov_h3[i][1] + cov_h3[i][2] )
cov_h3_o_o = cov_h3[4][4]
new_cov_h3 = np.asarray([	[cov_h3_abc_abc,cov_h3_abc_o],[cov_h3_o_abc,cov_h3_o_o] ])


[x,y,pear] = confidence_ellipse(new_cov_he3_h3, min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2] , min_pars_he3_h3[-1] , ax2[0] , edgecolor='red')
ax2[0].set_xlim((min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2])-1.5*x,(min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2])+1.5*x)
ax2[0].set_ylim(min_pars_he3_h3[-1]-1.5*y,min_pars_he3_h3[-1]+1.5*y)
ax2[0].set_title("He3 and H3 - Pearson R: %.2f" % np.round(pear,2) )

[x,y,pear] = confidence_ellipse(new_cov_h3, min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2] , min_pars_h3[-1] , ax2[1] , edgecolor='blue')
ax2[1].set_xlim((min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2])-1.5*x,(min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2])+1.5*x)
ax2[1].set_ylim(min_pars_h3[-1]-1.5*y,min_pars_h3[-1]+1.5*y)
ax2[1].set_title("He3 - Pearson R: %.2f" % np.round(pear,2) )

[x1,y1,pear1] = confidence_ellipse(new_cov_he3_h3, 0,0 , ax2[2] , edgecolor='red')
[x2,y2,pear2] = confidence_ellipse(new_cov_h3, 0,0 , ax2[2] , edgecolor='blue')
ax2[2].set_xlim(-max(x1,x2)*1.5,max(x1,x2)*1.5)
ax2[2].set_ylim(-max(y1,y2)*1.5,max(y1,y2)*1.5)

[x1,y1,pear1] = confidence_ellipse(new_cov_he3_h3, min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2] , min_pars_he3_h3[-1] , ax2[3] , edgecolor='red')
[x2,y2,pear2] = confidence_ellipse(new_cov_h3, min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2] , min_pars_h3[-1] , ax2[3] , edgecolor='blue')
ax2[3].set_xlim(0.1,0.65)
ax2[3].set_ylim(-3.5,4)
ax2[3].set_xlabel("F2n/F2p at x=1")
ax2[3].set_ylabel("Offshell par a")


plt.savefig("before-after.pdf")


plt.figure(4)
x_axis = np.linspace(0,1,1000)
plt.plot(x_axis, min_pars_he3_h3[0] + min_pars_he3_h3[1]*x_axis + min_pars_he3_h3[2]*np.exp( min_pars_he3_h3[3]*(np.ones(len(x_axis))-x_axis)) ,color='red',linewidth=3,label='He3 and H3')
plt.plot(x_axis, min_pars_h3[0] + min_pars_h3[1]*x_axis + min_pars_h3[2]*np.exp( min_pars_h3[3]*(np.ones(len(x_axis))-x_axis)) ,color='blue',linewidth=3,label='He3 Only')
plt.plot(x_axis, init_pars[0] + init_pars[1]*x_axis + init_pars[2]*np.exp( init_pars[3]*(np.ones(len(x_axis))-x_axis)) ,color='black',linestyle='--',linewidth=3,label='Initial Par')
plt.plot([0.9,1],[min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2]+2.*np.sqrt(cov_h3_abc_abc),min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2]+2.*np.sqrt(cov_h3_abc_abc)],color='blue')
plt.plot([0.9,1],[min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2]-2.*np.sqrt(cov_h3_abc_abc),min_pars_h3[0]+min_pars_h3[1]+min_pars_h3[2]-2.*np.sqrt(cov_h3_abc_abc)],color='blue')
plt.plot([0.9,1],[min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2]+2.*np.sqrt(cov_he3_h3_abc_abc),min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2]+2.*np.sqrt(cov_he3_h3_abc_abc)],color='red')
plt.plot([0.9,1],[min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2]-2.*np.sqrt(cov_he3_h3_abc_abc),min_pars_he3_h3[0]+min_pars_he3_h3[1]+min_pars_he3_h3[2]-2.*np.sqrt(cov_he3_h3_abc_abc)],color='red')
plt.ylabel('F2n/F2p')
plt.legend(numpoints=1,loc='best')
plt.xlabel('x')
plt.ylim([0.1,1])
plt.savefig("f2nf2p.pdf")
plt.show()
'''
for ax in axs:
	for subax in ax:
		# Make a corner plot from this!
		sub_cov = 
		#confidence_ellipse(cov, 0,0,ax, edgecolor='red')
		subax.set_xlim(-10,10)
		subax.set_ylim(-10,10)

'''
