import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from scipy.stats import norm
import sys
from scipy.interpolate import griddata

OF_p = []
OF_n = []
Chi2 = []
with open("test-contour.txt","rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		of_p = float(parse[0])
		of_n = float(parse[1])
		chi2 = float(parse[2])

		if chi2 == -1: continue
		if chi2 > 80: continue
		
		print chi2

		OF_p.append(of_p)
		OF_n.append(of_n)
		Chi2.append(chi2)

Chi2 = np.asarray(Chi2) - np.min(Chi2)

xi = np.linspace(-5,0,20)
yi = np.linspace(-5,15,20)
zi = griddata( (OF_p,OF_n), Chi2, (xi[None,:], yi[:,None]))

plt.contourf(xi,yi,zi, cmap=plt.cm.jet)
#plt.scatter(OF_p,OF_n,c=Chi2)
plt.xlabel('of_p',fontsize=15)
plt.ylabel('of_n',fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
cbar = plt.colorbar()
cbar.ax.set_title('Delta Chi2', fontsize=15 )
plt.scatter([-1.25864],[-1.25864],s=30,color='black')
plt.tight_layout()
plt.savefig('newcontour.pdf')
plt.show()
